"""Synchronous pipeline orchestrator."""

from __future__ import annotations

import logging
import shutil
import tempfile
from pathlib import Path

import ants
import pandas as pd
from aind_zarr_utils.pipeline_transformed import base_and_pipeline_anatomical_stub
from aind_zarr_utils.zarr import _open_zarr
from iblatlas.atlas import AllenAtlas

from aind_ibl_ephys_alignment_preprocessing.discovery import (
    determine_desired_level,
    find_asset_info,
    prepare_result_dirs,
)
from aind_ibl_ephys_alignment_preprocessing.ephys import run_ephys_for_recording
from aind_ibl_ephys_alignment_preprocessing.histology import (
    copy_registration_channel_ccf_reorient,
    process_additional_channels_pipeline,
    transform_ccf_labels_to_image_space,
    transform_ccf_to_image_space,
    write_registration_channel_images,
)
from aind_ibl_ephys_alignment_preprocessing.manifest import build_datapackage, write_datapackage
from aind_ibl_ephys_alignment_preprocessing.probes import process_manifest_row
from aind_ibl_ephys_alignment_preprocessing.types import (
    ManifestRow,
    PipelineConfig,
    ProcessResult,
    ReferencePaths,
    ReferenceVolumes,
)

logger = logging.getLogger(__name__)


def run_pipeline(config: PipelineConfig) -> list[ProcessResult]:
    """Run the full preprocessing pipeline synchronously.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.

    Returns
    -------
    list[ProcessResult]
        Per-probe processing results.
    """
    ref_paths = ReferencePaths.from_config(config)
    ref_imgs = ReferenceVolumes.from_paths(ref_paths)

    # Keep a manifest snapshot for reproducibility
    shutil.copy(config.manifest_csv, config.results_root / "manifest.csv")

    manifest_df = pd.read_csv(config.manifest_csv)
    mouse_id: str = str(manifest_df["mouseid"].astype("string").iat[0])

    out = prepare_result_dirs(mouse_id, config.results_root)
    asset_info = find_asset_info(config)

    node, zarr_metadata = _open_zarr(asset_info.zarr_volumes.registration)
    level = determine_desired_level(zarr_metadata, desired_voxel_size_um=config.desired_voxel_size_um)

    copy_registration_channel_ccf_reorient(asset_info, out)
    raw_img_path, pipeline_img_path = write_registration_channel_images(
        asset_info, out, level=level, opened_zarr=(node, zarr_metadata)
    )
    pipeline_img_ants = ants.image_read(str(pipeline_img_path), pixeltype=None)
    raw_img_ants = ants.image_read(str(raw_img_path), pixeltype=None)

    scratch_root = Path(config.scratch_root) if config.scratch_root is not None else Path(tempfile.mkdtemp())
    scratch_root.mkdir(parents=True, exist_ok=True)

    process_additional_channels_pipeline(
        pipeline_img_ants,
        asset_info,
        ref_imgs,
        out,
        scratch_root=scratch_root,
        level=level,
    )
    transform_ccf_to_image_space(asset_info, ref_imgs, raw_img_ants, pipeline_img_ants, out)
    transform_ccf_labels_to_image_space(asset_info, ref_paths, raw_img_ants, pipeline_img_ants, out)

    raw_img_stub, raw_img_stub_buggy, _ = base_and_pipeline_anatomical_stub(
        asset_info.zarr_volumes.registration,
        asset_info.zarr_volumes.metadata,
        asset_info.zarr_volumes.processing,
        opened_zarr=(node, zarr_metadata),
    )

    processed_recordings: set[str] = set()
    processed_results: list[ProcessResult] = []
    ibl_atlas = AllenAtlas(25, hist_path=ref_paths.ibl_atlas_histology_path)

    if config.skip_ephys:
        logger.info("Ephys processing disabled via skip_ephys flag")

    for _, row in manifest_df.iterrows():
        mr = ManifestRow.from_series(row)
        result = process_manifest_row(
            mr, asset_info, raw_img_stub, raw_img_stub_buggy, ibl_atlas, out, config.data_root
        )
        processed_results.append(result)
        if not result.wrote_files:
            logger.warning("Did not write files for %s: %s", mr.sorted_recording, result.skipped_reason)
            continue
        if not config.skip_ephys:
            run_ephys_for_recording(
                mr, out, config.data_root, processed_recordings, num_parallel_jobs=config.num_parallel_jobs
            )

    manifest_rows = [ManifestRow.from_series(row) for _, row in manifest_df.iterrows()]
    dp = build_datapackage(mouse_id, manifest_rows, processed_results, asset_info, out, config)
    dp_path = write_datapackage(dp, config.results_root / mouse_id)
    logger.info("Wrote datapackage manifest to %s", dp_path)

    return processed_results
