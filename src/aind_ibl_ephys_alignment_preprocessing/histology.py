"""Volume processing functions for histology data (synchronous).

Handles CCF reorientation, registration channel export, additional channel
transforms, and CCF-to-image-space inverse transforms.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import ants
import numpy as np
import SimpleITK as sitk
from aind_anatomical_utils.anatomical_volume import AnatomicalHeader
from aind_registration_utils.annotations import expand_compacted_image
from ants.core import ANTsImage
from ants.utils import to_sitk
from numpy.typing import DTypeLike

from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION
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


#: Identifies the sidecar's shape, so a reader can tell what it is holding.
_GEOMETRY_SCHEMA = "anatomical-header/1"
#: Physical frame of ``origin`` and the direction columns. Fixed by ITK, stated
#: because a consumer reaching for nibabel would otherwise assume RAS.
_GEOMETRY_SPACE = "left-posterior-superior"
#: ITK is unit-agnostic, so the numbers alone cannot say. ``aind_zarr_utils``
#: reads OME-Zarr with ``scale_unit="millimeter"``, so the pyramid's micrometre
#: scales arrive converted and every world coordinate here -- including the
#: stored ANTs warps -- is in millimetres.
_GEOMETRY_UNITS = "millimeter"


def image_geometry(img: sitk.Image | AnatomicalHeader) -> dict[str, Any]:
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
    header = img if isinstance(img, AnatomicalHeader) else AnatomicalHeader.from_sitk(img)
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


#: Image spacing throughout this pipeline is **millimetres**: ``aind_zarr_utils``
#: reads OME-Zarr with ``scale_unit="millimeter"``, so the micrometre scales in
#: the pyramid metadata arrive already converted, and the stored ANTs warps are
#: in the same frame. Targets are expressed in micrometres because that is the
#: natural unit for talking about these volumes, so they convert here -- once.
_UM_PER_MM = 1000.0
#: Smallest sane extent, in voxels, along any axis after resampling. A target
#: mistaken for the wrong unit is off by 1000x and silently collapses a volume to
#: a single voxel rather than failing, so it is checked instead of trusted.
_MIN_RESAMPLED_EXTENT = 8


def _mirror_pipeline_grid(pipeline: sitk.Image, base_before: sitk.Image, base_after: sitk.Image) -> AnatomicalHeader:
    """Carry the pipeline image through the base image's resample, in index space.

    ``base`` and ``pipeline`` are the *same voxels* under two headers: the base
    header comes from the Zarr scales, the pipeline header is the domain
    ``aind_zarr_utils`` reconstructs so the stored ANTs transforms apply in the
    frame they were trained in. Two things depend on that pairing staying exact:
    the CCF template and labels are warped onto the *pipeline* grid and then
    stamped with the *base* header, which is only meaningful when both grids have
    the same size; and the GUI overlays the results on the base volume voxel for
    voxel.

    Resampling the two independently breaks it. Their spacings need not be
    identical -- the pipeline's is the corrected level-0 spacing scaled by
    ``2**level``, the base's comes from the Zarr metadata -- so each computes a
    different output size from the same target, and nothing downstream notices.

    So the resample is done once, on the base, and the pipeline image inherits
    the result: the same voxels, with its header rescaled by the same per-axis
    factor. Origin is kept, which is the convention ``apply_pipeline_overlays``
    itself uses across levels (it scales spacing by ``2**level`` and passes the
    level-0 origin through unchanged), so index ``i`` of the resampled pipeline
    image lands where index ``i * factor`` of the original did.

    Parameters
    ----------
    pipeline : sitk.Image
        Pipeline-domain image, before resampling.
    base_before, base_after : sitk.Image
        The base image either side of the resample; their spacing ratio is the
        index-space factor to apply.

    Returns
    -------
    sitk.Image
        The resampled voxels carrying the rescaled pipeline header.
    """
    if base_after is base_before:
        return AnatomicalHeader.from_sitk(pipeline)

    factor = [after / before for before, after in zip(base_before.GetSpacing(), base_after.GetSpacing(), strict=True)]
    out = AnatomicalHeader(
        origin=pipeline.GetOrigin(),
        spacing=tuple(s * f for s, f in zip(pipeline.GetSpacing(), factor, strict=True)),
        direction=np.array(pipeline.GetDirection()).reshape(3, 3),
        size_ijk=base_after.GetSize(),
    )
    logger.info(
        "[registration-pipeline] grid mirrored from the base resample: %s @ %s mm (factor %s)",
        out.size_ijk,
        tuple(round(s, 6) for s in out.spacing),
        tuple(round(f, 4) for f in factor),
    )
    return out


def _blessed_header(header: AnatomicalHeader, label: str) -> AnatomicalHeader:
    """Return *header* in :data:`_BLESSED_DIRECTION`, without touching any pixels.

    The ANTs images this feeds were historically obtained by writing a volume and
    reading it back, and the write reorients -- so their headers were the
    *reoriented* ones. ``regrid_to`` reproduces that exactly (same permuted size,
    same origin as ``sitk.DICOMOrient``) from headers alone.

    It requires an axis-aligned header and raises otherwise, so a volume with an
    oblique direction falls back to orienting a stub image, which costs one
    voxel rather than the volume.
    """
    try:
        return header.regrid_to(_BLESSED_DIRECTION)
    except ValueError:
        logger.warning(
            "[%s] header is not axis-aligned; reorienting via a stub instead of regrid_to",
            label,
        )
        oriented = sitk.DICOMOrient(header.as_sitk_stub(), _BLESSED_DIRECTION)
        return AnatomicalHeader(
            origin=oriented.GetOrigin(),
            spacing=oriented.GetSpacing(),
            direction=np.array(oriented.GetDirection()).reshape(3, 3),
            size_ijk=tuple(int(n) for n in oriented.GetSize()),
        )


def ants_warp_domain(header: AnatomicalHeader, label: str, np_eltype: DTypeLike) -> ANTsImage:
    """Build the full-size, zero-filled ANTs image a warp uses as its ``fixed``.

    ``ants.apply_transforms`` reads ``fixed`` only for its *content*-free
    properties -- a zero-filled domain gives byte-identical output to the real
    volume -- but two of its properties are load-bearing and neither is
    negotiable:

    - **size**, which becomes the output size, so this must be full extent;
    - **pixel type**, because the documentation is explicit that "the output
      will have the same pixel type as this image". There is no separate
      output-type argument. A domain narrower than the moving data clips it
      silently: a ``uint8`` domain truncates the compacted label indices
      (~1341 distinct) and the CCF template (0-516) to 255, destroying both
      warped volumes without an error.

    So *np_eltype* is required rather than defaulted -- pass the dtype of the
    volume this domain replaces, which keeps the warp's output type exactly what
    the previous write-and-read-back produced.
    """
    return _blessed_header(header, label).as_ants(np_eltype)


def ants_domain_stub(header: AnatomicalHeader, label: str) -> ANTsImage:
    """Build the 1x1x1 ANTs image used only for its spacing/origin/direction.

    Pixel type is irrelevant here: this one is never a warp's ``fixed``, only the
    source of a spacing/origin/direction stamp.
    """
    blessed = _blessed_header(header, label)
    stub = AnatomicalHeader(
        origin=blessed.origin,
        spacing=blessed.spacing,
        direction=blessed.direction,
        size_ijk=(1, 1, 1),
    )
    return stub.as_ants()


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

    Notes
    -----
    ``ants.resample_image`` covers the grid arithmetic, but the volumes are
    SimpleITK at this point and a round trip through ANTs invites the axis-order
    and interop hazards that conversion is known for; it also does no low-pass,
    which is the part that actually matters here.

    The origin is kept, so output voxel ``(0, 0, 0)`` keeps the input's first
    voxel centre. The alternative is anchoring the outer corner, for which
    ``aind_anatomical_utils.anatomical_volume.fix_corner_compute_origin`` exists.
    The two differ by half the spacing change -- sub-micrometre here -- and
    either is physically consistent, because ``sitk.Resample`` places content by
    physical point regardless. It matters not at all for the CCF template and
    labels, which are warped onto *this* grid and so inherit whichever choice is
    made.

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
    target_mm = target_um / _UM_PER_MM
    if min(spacing) >= target_mm:
        logger.info(
            "[%s] already coarser than %.1f um (spacing %s mm); not resampling",
            label,
            target_um,
            tuple(round(s, 5) for s in spacing),
        )
        return img

    size = img.GetSize()
    target = (target_mm,) * img.GetDimension()
    new_size = [max(1, int(round(n * s / t))) for n, s, t in zip(size, spacing, target, strict=True)]
    if min(new_size) < _MIN_RESAMPLED_EXTENT:
        raise ValueError(
            f"[{label}] resampling to {target_um} um would give {tuple(new_size)} voxels from "
            f"{tuple(size)} at spacing {tuple(spacing)}. Spacing in this pipeline is in "
            f"millimetres, so a target that collapses the volume almost certainly means the "
            f"target and the image are in different units."
        )

    # Low-pass only the axes being coarsened; sigma = half the new spacing is the
    # usual cutoff at the new Nyquist frequency. Applied one axis at a time
    # because a mixed grid (some axes coarsening, some not) is normal here, and
    # SmoothingRecursiveGaussian rejects a zero sigma on any axis rather than
    # treating it as "leave this one alone".
    pixel_id = img.GetPixelID()
    work = sitk.Cast(img, sitk.sitkFloat32)
    for axis, (source, want) in enumerate(zip(spacing, target, strict=True)):
        if want > source:
            work = sitk.RecursiveGaussian(work, sigma=target_mm / 2.0, direction=axis)

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
        "[%s] resampled %s @ %s mm -> %s @ %.1f um (%.5f mm) isotropic",
        label,
        tuple(size),
        tuple(round(s, 5) for s in spacing),
        tuple(new_size),
        target_um,
        target_mm,
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


def convert_img_direction_and_write(
    img: sitk.Image,
    output_path: Path,
    direction: str = _BLESSED_DIRECTION,
) -> dict[str, Any]:
    """Convert an image to the specified DICOM orientation and write as compressed NRRD.

    Parameters
    ----------
    img : sitk.Image
        Input SimpleITK image.
    output_path : Path
        Destination NRRD file path.
    direction : str
        Target DICOM orientation code (default ``"IRP"``).
    """
    img_oriented = sitk.DICOMOrient(img, direction)
    sitk.WriteImage(img_oriented, str(output_path), useCompression=True)
    return image_geometry(img_oriented)


def copy_registration_channel_ccf_reorient(
    asset_info: AssetInfo,
    outputs: OutputDirs,
) -> None:
    """Copy precomputed CCF registration to results as reoriented NRRD.

    Parameters
    ----------
    asset_info : AssetInfo
        Asset metadata with precomputed registration path.
    outputs : OutputDirs
        Output directory tree.
    """
    if not asset_info.registration_in_ccf_precomputed.exists():
        raise FileNotFoundError(
            f"Precomputed registration in CCF not found: {asset_info.registration_in_ccf_precomputed}"
        )
    ccf_img = sitk.ReadImage(str(asset_info.registration_in_ccf_precomputed))
    img_in_ccf_dst = outputs.histology_ccf / "histology_registration.nrrd"
    convert_img_direction_and_write(ccf_img, img_in_ccf_dst)


def write_registration_channel_images(
    asset_info: AssetInfo,
    outputs: OutputDirs,
    *,
    level: int = 3,
    output_voxel_size_um: float,
    opened_zarr: tuple[Any, dict[str, Any]] | None = None,
) -> tuple[Path, AnatomicalHeader, AnatomicalHeader, DTypeLike]:
    """Write the registration channel and return the grids the warps need.

    Returns the written volume's path, the *base* and *pipeline* headers on the
    resampled grid, and the volume's dtype. The pipeline volume is not written:
    nothing reads its pixels, and the warps need only its geometry, which the
    sidecar carries.

    The dtype comes back because a warp's ``fixed`` image dictates the output
    pixel type. Returning the volume's own keeps the warps producing exactly
    what the previous write-and-read-back produced; the binding constraint is
    that it must hold the *moving* data's range, which uint16 does for both the
    compacted label indices and the CCF template.

    Parameters
    ----------
    asset_info : AssetInfo
        Asset metadata.
    outputs : OutputDirs
        Output directory tree.
    level : int
        Multiscale zarr level.
    opened_zarr : tuple | None
        Pre-opened zarr node and metadata.

    Returns
    -------
    tuple[Path, Path]
        ``(raw_img_path, pipeline_img_path)`` of written NRRD files.
    """
    from aind_zarr_utils.pipeline_transformed import base_and_pipeline_zarr_to_sitk
    from aind_zarr_utils.zarr import _open_zarr

    reg_zarr = asset_info.zarr_volumes.registration
    if opened_zarr is None:
        zarr_node, zarr_metadata = _open_zarr(reg_zarr)
    else:
        zarr_node, zarr_metadata = opened_zarr

    metadata = asset_info.zarr_volumes.metadata
    processing = asset_info.zarr_volumes.processing
    raw_img, pipeline_raw_img = base_and_pipeline_zarr_to_sitk(
        reg_zarr,
        metadata,
        processing,
        level=level,
        opened_zarr=(zarr_node, zarr_metadata),
    )
    resampled = resample_to_isotropic(raw_img, output_voxel_size_um, "registration")
    pipeline_header = _mirror_pipeline_grid(pipeline_raw_img, raw_img, resampled)
    del pipeline_raw_img
    base_header = AnatomicalHeader.from_sitk(resampled)
    warp_dtype = sitk.GetArrayViewFromImage(resampled).dtype
    raw_img_dst = outputs.histology_img / "histology_registration.nrrd"
    convert_img_direction_and_write(resampled, raw_img_dst)
    del resampled, raw_img
    sidecar = outputs.histology_img / "histology_registration_pipeline.json"
    sidecar.write_text(json.dumps(image_geometry(_blessed_header(pipeline_header, "registration-pipeline")), indent=2))
    return raw_img_dst, base_header, pipeline_header, warp_dtype


def process_additional_channels_pipeline(
    pipeline_histology_space_img: ANTsImage,
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    outputs: OutputDirs,
    scratch_root: Path,
    level: int = 3,
    *,
    output_voxel_size_um: float,
    emit_qc: bool = False,
) -> None:
    """Process non-alignment OME-Zarr channels into image space (+ CCF for QC).

    The image-space channel NRRD (consumed by the GUI) is always written. The
    CCF-space channel volume — a full-volume ANTs warp the GUI never reads — is
    produced only when *emit_qc* is True.

    Parameters
    ----------
    pipeline_histology_space_img : ANTsImage
        Pipeline-space histology image for spatial metadata copying.
    asset_info : AssetInfo
        Asset metadata.
    refs : ReferenceVolumes
        Reference volumes (CCF template).
    outputs : OutputDirs
        Output directory tree.
    scratch_root : Path
        Temporary directory for intermediate files.
    level : int
        Multiscale zarr level.
    emit_qc : bool
        Also warp each channel into CCF space (GUI-unused QC output).
    """
    from aind_zarr_utils.zarr import zarr_to_sitk

    for zarr_path in asset_info.zarr_volumes.additional:
        ch_str = Path(zarr_path).stem
        img_raw = zarr_to_sitk(zarr_path, asset_info.zarr_volumes.metadata, level=level)
        img_raw = resample_to_isotropic(img_raw, output_voxel_size_um, Path(zarr_path).stem)
        channel_dst = outputs.histology_img / f"{ch_str}.nrrd"
        convert_img_direction_and_write(img_raw, channel_dst)
        del img_raw

        if not emit_qc:
            continue  # skip the GUI-unused CCF-space volume warp

        ants_hist_img = ants.image_read(str(channel_dst), pixeltype=None)
        ants.copy_image_info(pipeline_histology_space_img, ants_hist_img)

        ch_in_ccf = ants.apply_transforms(
            refs.ccf_25,
            ants_hist_img,
            asset_info.pipeline_registration_chains.img_tx_str,
            whichtoinvert=asset_info.pipeline_registration_chains.img_tx_inverted,
        )
        ch_in_ccf_dst = outputs.histology_ccf / f"histology_{ch_str}.nrrd"
        ch_in_ccf_tmp_dst = scratch_root / f"histology-{ch_str}-ccf.nrrd"
        ants.image_write(ch_in_ccf, str(ch_in_ccf_tmp_dst))
        del ch_in_ccf
        try:
            compress_reorient_nrrd_file(
                ch_in_ccf_tmp_dst,
                ch_in_ccf_dst,
                force_orientation=_BLESSED_DIRECTION,
            )
        finally:
            ch_in_ccf_tmp_dst.unlink(missing_ok=True)


def apply_ccf_inverse_tx_then_fix_domain(
    ccf_space_img_moving: ANTsImage,
    pipeline_space_fixed_img: ANTsImage,
    correct_hist_domain_img: ANTsImage,
    asset_info: AssetInfo,
    **kwargs: Any,
) -> ANTsImage:
    """Apply inverse pipeline (CCF -> histology) transform then repair image domain.

    Parameters
    ----------
    ccf_space_img_moving : ANTsImage
        Image in CCF/template space to move into histology space.
    pipeline_space_fixed_img : ANTsImage
        Image in the pipeline's (buggy) histology space.
    correct_hist_domain_img : ANTsImage
        Reference histology image with correct spacing/origin/direction.
    asset_info : AssetInfo
        Pipeline registration chain paths.
    **kwargs
        Forwarded to ``ants.apply_transforms``.

    Returns
    -------
    ANTsImage
        Transformed image with corrected spatial domain.
    """
    pt_tx_str = asset_info.pipeline_registration_chains.pt_tx_str
    pt_tx_inverted = asset_info.pipeline_registration_chains.pt_tx_inverted
    ccf_space_img_in_hist_space: ANTsImage = ants.apply_transforms(
        fixed=pipeline_space_fixed_img,
        moving=ccf_space_img_moving,
        transformlist=pt_tx_str,
        whichtoinvert=pt_tx_inverted,
        **kwargs,
    )
    ccf_space_img_in_hist_space.set_spacing(correct_hist_domain_img.spacing)
    ccf_space_img_in_hist_space.set_origin(correct_hist_domain_img.origin)
    ccf_space_img_in_hist_space.set_direction(correct_hist_domain_img.direction)
    return ccf_space_img_in_hist_space


def compress_reorient_nrrd_file(
    input_path: Path,
    output_path: Path,
    force_orientation: str | None = None,
) -> None:
    """Re-compress and optionally reorient an NRRD file.

    Parameters
    ----------
    input_path : Path
        Source NRRD file.
    output_path : Path
        Destination NRRD file.
    force_orientation : str | None
        If set, reorient to this DICOM code before writing.
    """
    img = sitk.ReadImage(str(input_path))
    orientation_code = sitk.DICOMOrientImageFilter.GetOrientationFromDirectionCosines(img.GetDirection())
    if force_orientation is not None and orientation_code != force_orientation:
        logger.info("Reorienting %s from %s to %s", input_path, orientation_code, force_orientation)
        out_img = sitk.DICOMOrient(img, force_orientation)
    else:
        out_img = img
    temp_output_path = output_path.with_suffix(".temp.nrrd")
    sitk.WriteImage(out_img, str(temp_output_path), useCompression=True)
    temp_output_path.replace(output_path)


def transform_ccf_to_image_space(
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
) -> None:
    """Transform CCF template into native image space.

    Parameters
    ----------
    asset_info : AssetInfo
        Pipeline registration chain paths.
    refs : ReferenceVolumes
        Reference volumes.
    raw_hist_img : ANTsImage
        Histology image with correct spatial domain.
    pipeline_hist_domain_img : ANTsImage
        Histology image in pipeline (buggy) domain.
    outputs : OutputDirs
        Output directory tree.
    """
    ccf_in_hist_img = apply_ccf_inverse_tx_then_fix_domain(
        refs.ccf_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
    )
    ccf_in_hist_img_path = outputs.histology_img / "ccf_in_mouse.nrrd"
    ccf_in_hist_sitk = to_sitk(ccf_in_hist_img)
    del ccf_in_hist_img
    ccf_in_hist_sitk = _narrow_dtype(ccf_in_hist_sitk, sitk.sitkUInt16, "ccf_in_mouse")
    convert_img_direction_and_write(ccf_in_hist_sitk, ccf_in_hist_img_path)


def transform_ccf_labels_to_image_space(
    asset_info: AssetInfo,
    ref_paths: ReferencePaths,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
) -> None:
    """Transform lateralized CCF labels into native image space.

    Parameters
    ----------
    asset_info : AssetInfo
        Pipeline registration chain paths.
    ref_paths : ReferencePaths
        Reference data paths.
    raw_hist_img : ANTsImage
        Histology image with correct spatial domain.
    pipeline_hist_domain_img : ANTsImage
        Histology image in pipeline (buggy) domain.
    outputs : OutputDirs
        Output directory tree.
    """
    ccf_labels_lateralized_25 = ants.image_read(
        str(ref_paths.ccf_labels_lateralized_25),
        pixeltype=None,
    )
    unq_vals = np.load(str(ref_paths.ccf_labels_lateralized_25_unq_vals))["unique_labels"]
    ccf_labels_in_hist_img = apply_ccf_inverse_tx_then_fix_domain(
        ccf_labels_lateralized_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
        interpolator="genericLabel",
    )
    del ccf_labels_lateralized_25
    ccf_labels_sitk = to_sitk(ccf_labels_in_hist_img)
    del ccf_labels_in_hist_img
    ccf_labels_expanded = expand_compacted_image(ccf_labels_sitk, unq_vals)
    ccf_labels_expanded = _narrow_dtype(ccf_labels_expanded, sitk.sitkInt32, "labels_in_mouse")
    ccf_labels_in_hist_img_path = outputs.histology_img / "labels_in_mouse.nrrd"
    convert_img_direction_and_write(ccf_labels_expanded, ccf_labels_in_hist_img_path)
