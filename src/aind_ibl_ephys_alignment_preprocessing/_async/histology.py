"""Async volume processing functions for histology data."""

from __future__ import annotations

import asyncio
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


async def convert_img_to_direction_and_write_async(
    img: sitk.Image,
    dst_path: Path | str,
    limits: Limits,
    direction: str = _BLESSED_DIRECTION,
) -> None:
    """Async convert image orientation and write to disk."""
    logger.info("[Histology] Converting image for %s to %s orientation", dst_path, direction)
    img_oriented = await to_thread_logged(sitk.DICOMOrient, img, direction)
    logger.info("[Histology] Writing image for %s to disk", dst_path)
    await io_to_thread_on(limits, str(dst_path), sitk.WriteImage, img_oriented, str(dst_path), useCompression=True)
    logger.info("[Histology] Done writing image for %s to disk", dst_path)


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
    opened_zarr: tuple[Any, dict[str, Any]] | None = None,
) -> tuple[Path, Path]:
    """Async write registration-channel outputs to CCF and image space."""
    from aind_zarr_utils.pipeline_transformed import base_and_pipeline_zarr_to_sitk
    from aind_zarr_utils.zarr import _open_zarr

    reg_zarr = asset_info.zarr_volumes.registration
    zarr_name = Path(reg_zarr).stem
    logger.info("[Histology] Reading registration channel from zarr: %s at level %d", zarr_name, level)
    if opened_zarr is None:
        zarr_node, zarr_metadata = await to_thread_logged(_open_zarr, reg_zarr)
    else:
        zarr_node, zarr_metadata = opened_zarr

    metadata = asset_info.zarr_volumes.metadata
    processing = asset_info.zarr_volumes.processing
    raw_img, pipeline_raw_img = await to_thread_logged(
        base_and_pipeline_zarr_to_sitk,
        reg_zarr,
        metadata,
        processing,
        level=level,
        opened_zarr=(zarr_node, zarr_metadata),
    )
    logger.info("[Histology] Registration channel loaded: raw + pipeline-transformed images")
    raw_img_dst = outputs.histology_img / "histology_registration.nrrd"
    bugged_img_dst = outputs.histology_img / "histology_registration_pipeline.nrrd"
    logger.info("[Histology] Registration channel conversion to %s + write started", _BLESSED_DIRECTION)
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            convert_img_to_direction_and_write_async(raw_img, raw_img_dst, limits),
            name="write-registration-raw",
        )
        tg.create_task(
            convert_img_to_direction_and_write_async(pipeline_raw_img, bugged_img_dst, limits),
            name="write-registration-pipeline",
        )
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
    async with limits.registration:
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
) -> None:
    """Async process a single additional OME-Zarr channel."""
    from aind_zarr_utils.zarr import zarr_to_sitk

    ch_str = Path(zarr_path).stem
    logger.info("[Channel %s] Starting processing", ch_str)
    img_raw = await to_thread_logged(zarr_to_sitk, zarr_path, asset_info.zarr_volumes.metadata, level=level)
    logger.info("[Channel %s] read from zarr complete", ch_str)
    channel_dst = outputs.histology_img / f"{ch_str}.nrrd"
    await convert_img_to_direction_and_write_async(img_raw, channel_dst, limits)
    logger.info("[Channel %s] converted to %s and written to disk", ch_str, _BLESSED_DIRECTION)
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
                ),
                name=f"channel-{ch_name}",
            )
