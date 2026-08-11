"""Async volume processing functions for histology data."""

from __future__ import annotations

import asyncio
import json
import logging
from pathlib import Path
from typing import Any

import ants
import SimpleITK as sitk
from aind_registration_utils.annotations import expand_compacted_image
from ants.core import ANTsImage
from ants.utils import to_sitk

from aind_ibl_ephys_alignment_preprocessing._async.concurrency import (
    Limits,
    io_to_thread_on,
    to_thread_logged,
)
from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION
from aind_ibl_ephys_alignment_preprocessing._timing import timed
from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    OutputDirs,
    ReferencePaths,
    ReferenceVolumes,
)

try:
    import numpy as np
except ImportError:  # pragma: no cover
    pass

logger = logging.getLogger(__name__)


async def compress_reorient_nrrd_file_async(
    input_path: Path,
    output_path: Path,
    limits: Limits,
    force_orientation: str | None = None,
) -> None:
    """Async version of compress/reorient NRRD file."""
    logger.info("Reading %s for compression and reorientation", input_path)
    img = await io_to_thread_on(limits, str(input_path), sitk.ReadImage, str(input_path))
    orientation_code = sitk.DICOMOrientImageFilter.GetOrientationFromDirectionCosines(img.GetDirection())
    if force_orientation is not None and orientation_code != force_orientation:
        logger.info("Reorienting %s from %s to %s", input_path, orientation_code, force_orientation)
        out_img = await to_thread_logged(sitk.DICOMOrient, img, force_orientation)
    else:
        out_img = img
    temp_output_path = output_path.with_suffix(".temp.nrrd")
    logger.info("Writing %s for compression and reorientation", input_path)
    await io_to_thread_on(
        limits, str(temp_output_path), sitk.WriteImage, out_img, str(temp_output_path), useCompression=True
    )
    logger.info("Replacing %s for compression and reorientation", input_path)
    await io_to_thread_on(limits, str(output_path), temp_output_path.replace, output_path)


#: Identifies the sidecar's shape, so a reader can tell what it is holding.
_GEOMETRY_SCHEMA = "anatomical-header/1"

#: Physical frame of ``origin`` and the direction columns. Fixed by ITK, stated
#: because a consumer reaching for nibabel would otherwise assume RAS.
_GEOMETRY_SPACE = "left-posterior-superior"

#: ITK is unit-agnostic, so the numbers alone cannot say. These volumes inherit
#: micrometres from the OME-Zarr scales, which matters because
#: ``AnatomicalHeader`` documents millimetres: same class, different unit, and
#: nothing raises if they are mixed.
_GEOMETRY_UNITS = "micrometer"


def image_geometry(img: sitk.Image) -> dict[str, Any]:
    """Serialize a grid's geometry: everything index -> physical needs, nothing more.

    Built from :class:`aind_anatomical_utils.anatomical_volume.AnatomicalHeader`
    so the payload is exactly that class's constructor arguments, and a consumer
    reconstructs with

    .. code-block:: python

        header = AnatomicalHeader(
            origin=tuple(h["origin"]),
            spacing=tuple(h["spacing"]),
            direction=np.array(h["direction"]).reshape(3, 3),
            size_ijk=tuple(h["size_ijk"]),
        )

    rather than re-deriving a convention from loose numbers.

    The field names carry the convention because getting it wrong is silent:

    - ``size_ijk`` and ``spacing`` are per **index axis** ``(i, j, k)``, which is
      SimpleITK's ``GetSize()`` order and therefore the **reverse** of the numpy
      ``array.shape``. Read as a shape, it transposes the volume.
    - ``origin`` and the direction's **columns** are physical LPS.
    - ``direction`` is the 3x3 flattened **row-major**, ready for
      ``SetDirection``; its columns are the index axes' LPS unit vectors.

    Together: ``physical = origin + direction @ (spacing * index)``.

    Parameters
    ----------
    img : sitk.Image
        Image whose geometry to capture. Must be the image *as written*: the
        orientation conversion changes origin and direction, so geometry taken
        beforehand describes a grid that does not exist on disk.

    Returns
    -------
    dict
        Self-describing payload; the constructor arguments live under ``header``.
    """
    from aind_anatomical_utils.anatomical_volume import AnatomicalHeader

    header = AnatomicalHeader.from_sitk(img)
    return {
        "schema": _GEOMETRY_SCHEMA,
        "space": _GEOMETRY_SPACE,
        "units": _GEOMETRY_UNITS,
        "header": {
            "origin": [float(v) for v in header.origin],
            "spacing": [float(v) for v in header.spacing],
            "direction": [float(v) for v in header.direction_tuple()],
            "size_ijk": [int(v) for v in header.size_ijk],
        },
    }


async def convert_img_to_direction_and_write_async(
    img: sitk.Image,
    dst_path: Path | str,
    limits: Limits,
    direction: str = _BLESSED_DIRECTION,
) -> dict[str, Any]:
    """Async convert image orientation and write to disk.

    Returns
    -------
    dict
        The written image's geometry, per :func:`image_geometry`.
    """
    name = Path(dst_path).name
    logger.info("[Histology] Converting image for %s to %s orientation", dst_path, direction)
    with timed("histology.reorient", volume=name):
        img_oriented = await to_thread_logged(sitk.DICOMOrient, img, direction)
    logger.info("[Histology] Writing image for %s to disk", dst_path)
    with timed("histology.write", volume=name) as record:
        await io_to_thread_on(limits, str(dst_path), sitk.WriteImage, img_oriented, str(dst_path), useCompression=True)
        record["mvox"] = int(np.prod(img_oriented.GetSize()) // 10**6)
    logger.info("[Histology] Done writing image for %s to disk", dst_path)
    return image_geometry(img_oriented)


async def copy_registration_channel_ccf_reorient_async(
    asset_info: AssetInfo,
    outputs: OutputDirs,
    limits: Limits,
) -> None:
    """Async copy precomputed CCF registration to results."""
    logger.info("[CCF Copy] Copying precomputed CCF registration to results")
    if not asset_info.registration_in_ccf_precomputed.exists():
        raise FileNotFoundError(
            f"Precomputed registration in CCF not found: {asset_info.registration_in_ccf_precomputed}"
        )
    ccf_img = await io_to_thread_on(
        limits,
        str(asset_info.registration_in_ccf_precomputed),
        sitk.ReadImage,
        str(asset_info.registration_in_ccf_precomputed),
    )
    logger.info("[CCF Copy] Read precomputed CCF registration image")
    ccf_img_dest = str(outputs.histology_ccf / "histology_registration.nrrd")
    await convert_img_to_direction_and_write_async(ccf_img, ccf_img_dest, limits)
    logger.info("[CCF Copy] Completed: histology_registration.nrrd in CCF space")


async def write_registration_channel_images_async(
    asset_info: AssetInfo,
    outputs: OutputDirs,
    limits: Limits,
    *,
    level: int = 3,
    output_voxel_size_um: float,
    opened_zarr: tuple[Any, dict[str, Any]] | None = None,
) -> tuple[Path, Path]:
    """Async write registration-channel outputs to CCF and image space."""
    from aind_zarr_utils.pipeline_transformed import base_and_pipeline_zarr_to_sitk
    from aind_zarr_utils.zarr import _open_zarr

    reg_zarr = asset_info.zarr_volumes.registration
    zarr_name = Path(reg_zarr).stem
    logger.info("[Histology] Reading registration channel from zarr: %s at level %d", zarr_name, level)
    if opened_zarr is None:
        with timed("histology.zarr_open", volume=zarr_name):
            zarr_node, zarr_metadata = await to_thread_logged(_open_zarr, reg_zarr)
    else:
        zarr_node, zarr_metadata = opened_zarr

    metadata = asset_info.zarr_volumes.metadata
    processing = asset_info.zarr_volumes.processing
    with timed("histology.zarr_read", volume=zarr_name, level=level) as record:
        raw_img, pipeline_raw_img = await to_thread_logged(
            base_and_pipeline_zarr_to_sitk,
            reg_zarr,
            metadata,
            processing,
            level=level,
            opened_zarr=(zarr_node, zarr_metadata),
        )
        record["mvox"] = int(np.prod(raw_img.GetSize()) // 10**6)
    logger.info("[Histology] Registration channel loaded: raw + pipeline-transformed images")
    with timed("histology.resample", volume=zarr_name):
        raw_img = resample_to_isotropic(raw_img, output_voxel_size_um, "registration")
        pipeline_raw_img = resample_to_isotropic(pipeline_raw_img, output_voxel_size_um, "registration-pipeline")
    raw_img_dst = outputs.histology_img / "histology_registration.nrrd"
    bugged_img_dst = outputs.histology_img / "histology_registration_pipeline.nrrd"
    logger.info("[Histology] Registration channel conversion to %s + write started", _BLESSED_DIRECTION)
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            convert_img_to_direction_and_write_async(raw_img, raw_img_dst, limits),
            name="write-registration-raw",
        )
        pipeline_task = tg.create_task(
            convert_img_to_direction_and_write_async(pipeline_raw_img, bugged_img_dst, limits),
            name="write-registration-pipeline",
        )
    # The pipeline volume's pixels are never read -- it is voxel-identical to the
    # raw one and differs only in origin, and its only consumer maps indices to
    # physical points. Emit that geometry as a sidecar so the volume itself can
    # be retired once consumers read this instead.
    geometry_dst = bugged_img_dst.with_suffix(".json")
    geometry_dst.write_text(json.dumps(pipeline_task.result(), indent=2))
    logger.info("[Histology] Wrote pipeline-image geometry sidecar: %s", geometry_dst.name)
    return raw_img_dst, bugged_img_dst


async def apply_ccf_inverse_tx_then_fix_domain_async(
    ccf_space_img_moving: ANTsImage,
    pipeline_space_fixed_img: ANTsImage,
    correct_hist_domain_img: ANTsImage,
    asset_info: AssetInfo,
    limits: Limits,
    **kwargs: Any,
) -> ANTsImage:
    """Async version of CCF inverse transform with domain repair."""
    pt_tx_str = asset_info.pipeline_registration_chains.pt_tx_str
    pt_tx_inverted = asset_info.pipeline_registration_chains.pt_tx_inverted
    interpolator = str(kwargs.get("interpolator", "linear"))
    # ``timed`` is a *sync* context manager, so it cannot share the ``async
    # with``. Nesting it inside the semaphore is also what we want: it then
    # measures the warp rather than the wait for a registration slot.
    async with limits.registration:
        with timed("histology.warp", interpolator=interpolator):
            ccf_img_in_hist_space: ANTsImage = await to_thread_logged(
                ants.apply_transforms,
                fixed=pipeline_space_fixed_img,
                moving=ccf_space_img_moving,
                transformlist=pt_tx_str,
                whichtoinvert=pt_tx_inverted,
                **kwargs,
            )
    ccf_img_in_hist_space.set_spacing(correct_hist_domain_img.spacing)
    ccf_img_in_hist_space.set_origin(correct_hist_domain_img.origin)
    ccf_img_in_hist_space.set_direction(correct_hist_domain_img.direction)
    return ccf_img_in_hist_space


def resample_to_isotropic(img: sitk.Image, target_um: float, label: str) -> sitk.Image:
    """Resample an intensity volume onto an isotropic *target_um* grid.

    Reading a multiscale level does not give a predictable grid: whether a volume
    has a level near the target depends on when it was stitched, so without this
    each mouse lands on whatever its pyramid happened to offer and cross-animal
    comparisons are differently sampled. Resampling puts every mouse on one grid.

    Two properties this must preserve, both of which are easy to lose:

    - **Anti-aliasing.** ``sitk.Resample`` point-samples; it applies no low-pass.
      Going 14.4 -> 30 µm without smoothing aliases the fluorescence. A Gaussian
      with physical sigma of half the target spacing is applied first, per axis
      and only where that axis is actually being coarsened.
    - **Dtype.** Smoothing needs floats, so the image round-trips through float32
      and is cast back. The channels use their full uint16 depth and must not
      silently become float on disk.

    A volume already at or coarser than the target on every axis is returned
    unchanged -- the target never manufactures resolution. Reaching a *common*
    isotropic grid can still interpolate up on one axis when that axis alone is
    coarser than the target; that is deliberate, and small for the volumes in
    hand.

    Labels are deliberately not routed through here. They reach image space via
    the ANTs warp with ``genericLabel``, whose fixed image is already on this
    grid, so they are resampled once, by a label-aware interpolator, rather than
    twice.

    Parameters
    ----------
    img : sitk.Image
        Intensity volume to resample.
    target_um : float
        Isotropic target spacing, micrometres.
    label : str
        Volume name, for the log line.

    Returns
    -------
    sitk.Image
        The resampled volume, or *img* unchanged when it is already coarser.
    """
    spacing = img.GetSpacing()
    if min(spacing) >= target_um:
        logger.info(
            "[%s] already coarser than %.1f um (spacing %s); not resampling",
            label,
            target_um,
            tuple(round(s, 2) for s in spacing),
        )
        return img

    size = img.GetSize()
    target = (target_um,) * img.GetDimension()
    new_size = [max(1, int(round(n * s / t))) for n, s, t in zip(size, spacing, target, strict=True)]

    # Low-pass only the axes being coarsened; sigma = half the new spacing is the
    # usual cutoff at the new Nyquist frequency. Applied one axis at a time
    # because a mixed grid (some axes coarsening, some not) is normal here, and
    # SmoothingRecursiveGaussian rejects a zero sigma on any axis rather than
    # treating it as "leave this one alone".
    pixel_id = img.GetPixelID()
    work = sitk.Cast(img, sitk.sitkFloat32)
    for axis, (source, want) in enumerate(zip(spacing, target, strict=True)):
        if want > source:
            work = sitk.RecursiveGaussian(work, sigma=target_um / 2.0, direction=axis)

    resampled = sitk.Resample(
        work,
        new_size,
        sitk.Transform(),
        sitk.sitkLinear,
        img.GetOrigin(),
        target,
        img.GetDirection(),
        0.0,
        sitk.sitkFloat32,
    )
    logger.info(
        "[%s] resampled %s @ %s -> %s @ %.1f um isotropic",
        label,
        tuple(size),
        tuple(round(s, 2) for s in spacing),
        tuple(new_size),
        target_um,
    )
    out: sitk.Image = sitk.Cast(resampled, pixel_id)
    return out


#: Bounds of the signed/unsigned integer pixel types we narrow to, so a cast is
#: only made when it provably cannot alter a voxel.
_DTYPE_RANGE: dict[int, tuple[int, int]] = {
    sitk.sitkUInt8: (0, 2**8 - 1),
    sitk.sitkInt16: (-(2**15), 2**15 - 1),
    sitk.sitkUInt16: (0, 2**16 - 1),
    sitk.sitkInt32: (-(2**31), 2**31 - 1),
    sitk.sitkUInt32: (0, 2**32 - 1),
}


def _narrow_dtype(img: sitk.Image, pixel_type: int, label: str) -> sitk.Image:
    """Cast *img* to *pixel_type* when its values fit, else leave it alone.

    Two of the written volumes carry a dtype inherited from how they were built
    rather than from what they hold: the expanded label image takes int64 from
    the LUT array it is indexed out of, and the CCF template arrives as uint32
    holding a 10-bit intensity scale. Each is a straight doubling of the bytes
    pushed through the warp, the reorientation and the write.

    The cast is guarded rather than assumed, because a silent overflow here would
    corrupt structure ids into other valid structure ids -- wrong in a way that
    still looks like a brain. If the range does not fit, the image is returned
    unchanged and the reason logged.

    Note this is deliberately not applied to the fluorescence channels: they use
    their full uint16 depth, and the GUI windows contrast from the source values
    rather than from a pre-quantized image.

    Parameters
    ----------
    img : sitk.Image
        Image to narrow.
    pixel_type : int
        Target SimpleITK pixel type.
    label : str
        Volume name, for the log line.

    Returns
    -------
    sitk.Image
        The narrowed image, or *img* unchanged when the values do not fit.
    """
    low, high = _DTYPE_RANGE[pixel_type]
    stats = sitk.MinimumMaximumImageFilter()
    stats.Execute(img)
    actual_low, actual_high = stats.GetMinimum(), stats.GetMaximum()
    if actual_low < low or actual_high > high:
        logger.warning(
            "[%s] leaving dtype as %s: values [%s, %s] do not fit %s",
            label,
            img.GetPixelIDTypeAsString(),
            actual_low,
            actual_high,
            sitk.GetPixelIDValueAsString(pixel_type),
        )
        return img
    logger.info(
        "[%s] narrowing %s -> %s (values [%s, %s])",
        label,
        img.GetPixelIDTypeAsString(),
        sitk.GetPixelIDValueAsString(pixel_type),
        actual_low,
        actual_high,
    )
    narrowed: sitk.Image = sitk.Cast(img, pixel_type)
    return narrowed


async def transform_ccf_to_image_space_async(
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
    limits: Limits,
) -> None:
    """Async transform CCF template into native image space."""
    logger.info("[CCF Transform] Starting CCF template -> image space transform")
    ccf_in_hist_img = await apply_ccf_inverse_tx_then_fix_domain_async(
        refs.ccf_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
        limits=limits,
    )
    ccf_in_hist_img_path = outputs.histology_img / "ccf_in_mouse.nrrd"
    ccf_in_hist_sitk = await to_thread_logged(to_sitk, ccf_in_hist_img)
    del ccf_in_hist_img
    ccf_in_hist_sitk = _narrow_dtype(ccf_in_hist_sitk, sitk.sitkUInt16, "ccf_in_mouse")
    await convert_img_to_direction_and_write_async(ccf_in_hist_sitk, ccf_in_hist_img_path, limits)
    logger.info("[CCF Transform] Completed: %s", ccf_in_hist_img_path.name)


async def transform_ccf_labels_to_image_space_async(
    asset_info: AssetInfo,
    ref_paths: ReferencePaths,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
    limits: Limits,
) -> None:
    """Async transform lateralized CCF labels into native image space."""
    logger.info("[CCF Labels] Starting CCF labels -> image space transform")
    ccf_labels_lateralized_25 = await to_thread_logged(
        ants.image_read,
        str(ref_paths.ccf_labels_lateralized_25),
        pixeltype=None,
    )
    unq_vals = np.load(str(ref_paths.ccf_labels_lateralized_25_unq_vals))["unique_labels"]
    ccf_labels_in_hist_img = await apply_ccf_inverse_tx_then_fix_domain_async(
        ccf_labels_lateralized_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
        limits=limits,
        interpolator="genericLabel",
    )
    del ccf_labels_lateralized_25
    ccf_labels_sitk = await to_thread_logged(to_sitk, ccf_labels_in_hist_img)
    del ccf_labels_in_hist_img
    ccf_labels_expanded = expand_compacted_image(ccf_labels_sitk, unq_vals)
    ccf_labels_expanded = _narrow_dtype(ccf_labels_expanded, sitk.sitkInt32, "labels_in_mouse")
    ccf_labels_in_hist_img_path = outputs.histology_img / "labels_in_mouse.nrrd"
    await convert_img_to_direction_and_write_async(ccf_labels_expanded, ccf_labels_in_hist_img_path, limits)
    logger.info("[CCF Labels] Completed: %s", ccf_labels_in_hist_img_path.name)


async def process_additional_channel_pipeline_async(
    zarr_path: str,
    pipeline_histology_space_img: ANTsImage,
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    outputs: OutputDirs,
    limits: Limits,
    scratch_root: Path,
    level: int = 3,
    *,
    output_voxel_size_um: float,
    emit_qc: bool = False,
) -> None:
    """Async process a single additional OME-Zarr channel.

    Writes the image-space channel NRRD always; the CCF-space warp (GUI-unused
    QC) only when *emit_qc* is True.
    """
    from aind_zarr_utils.zarr import zarr_to_sitk

    ch_str = Path(zarr_path).stem
    logger.info("[Channel %s] Starting processing", ch_str)
    with timed("histology.zarr_read", volume=ch_str, level=level) as record:
        img_raw = await to_thread_logged(zarr_to_sitk, zarr_path, asset_info.zarr_volumes.metadata, level=level)
        record["mvox"] = int(np.prod(img_raw.GetSize()) // 10**6)
    logger.info("[Channel %s] read from zarr complete", ch_str)
    with timed("histology.resample", volume=ch_str):
        img_raw = resample_to_isotropic(img_raw, output_voxel_size_um, ch_str)
    channel_dst = outputs.histology_img / f"{ch_str}.nrrd"
    await convert_img_to_direction_and_write_async(img_raw, channel_dst, limits)
    logger.info("[Channel %s] converted to %s and written to disk", ch_str, _BLESSED_DIRECTION)

    if not emit_qc:
        logger.info("[Channel %s] Completed (image-space only; QC off)", ch_str)
        return

    ants_hist_img = await io_to_thread_on(limits, str(channel_dst), ants.image_read, str(channel_dst), pixeltype=None)
    logger.info("[Channel %s] read into ANTs complete", ch_str)
    ants.copy_image_info(pipeline_histology_space_img, ants_hist_img)

    logger.debug("[Channel %s] Applying ANTs transform to CCF", ch_str)
    async with limits.registration:
        ch_in_ccf = await to_thread_logged(
            ants.apply_transforms,
            refs.ccf_25,
            ants_hist_img,
            asset_info.pipeline_registration_chains.img_tx_str,
            whichtoinvert=asset_info.pipeline_registration_chains.img_tx_inverted,
        )
    ch_in_ccf_dst = outputs.histology_ccf / f"histology_{ch_str}.nrrd"
    ch_in_ccf_tmp_dst = scratch_root / f"histology-{ch_str}-ccf.nrrd"
    logger.info("[Registered channel %s] writing to disk", ch_str)
    await io_to_thread_on(limits, str(ch_in_ccf_tmp_dst), ants.image_write, ch_in_ccf, str(ch_in_ccf_tmp_dst))
    del ch_in_ccf
    try:
        logger.info("[Registered channel %s] compressing and reorienting", ch_str)
        await compress_reorient_nrrd_file_async(
            ch_in_ccf_tmp_dst, ch_in_ccf_dst, limits, force_orientation=_BLESSED_DIRECTION
        )
    finally:
        ch_in_ccf_tmp_dst.unlink(missing_ok=True)
    logger.info("[Channel %s] Completed: %s + histology_%s.nrrd", ch_str, channel_dst.name, ch_str)


async def process_additional_channels_pipeline_async(
    pipeline_histology_space_img: ANTsImage,
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    outputs: OutputDirs,
    limits: Limits,
    scratch_root: Path,
    level: int = 3,
    *,
    output_voxel_size_um: float,
    emit_qc: bool = False,
) -> None:
    """Async dispatch all additional channels in parallel."""
    async with asyncio.TaskGroup() as tg:
        for zarr_path in asset_info.zarr_volumes.additional:
            ch_name = Path(zarr_path).stem
            tg.create_task(
                process_additional_channel_pipeline_async(
                    zarr_path,
                    pipeline_histology_space_img,
                    asset_info,
                    refs,
                    outputs,
                    limits,
                    scratch_root=scratch_root,
                    level=level,
                    output_voxel_size_um=output_voxel_size_um,
                    emit_qc=emit_qc,
                ),
                name=f"channel-{ch_name}",
            )
