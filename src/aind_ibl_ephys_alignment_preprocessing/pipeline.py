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
from aind_ibl_ephys_alignment_preprocessing.manifest import (
    build_datapackage,
    infer_process_results_from_outputs,
    producer_asset_overrides,
    write_datapackage,
)
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
    dp_path = write_datapackage(
        dp,
        config.results_root / mouse_id,
        asset_roots=[config.data_root],
        asset_overrides=producer_asset_overrides(asset_info, config),
    )
    logger.info("Wrote datapackage manifest to %s", dp_path)

    return processed_results


def regenerate_datapackage(
    config: PipelineConfig,
    *,
    validate: bool = True,
    source_results: Path | None = None,
) -> Path:
    """Rewrite ``datapackage.json`` from existing preprocessing outputs only.

    If *source_results* is provided, copy the existing immutable output tree
    into ``config.results_root`` first. This supports Code Ocean workflows where
    a previous results asset is mounted read-only under ``/data`` and the new
    run must write corrected metadata under ``/results``.
    """
    manifest_df = pd.read_csv(config.manifest_csv)
    mouse_id: str = str(manifest_df["mouseid"].astype("string").iat[0])
    config.results_root.mkdir(parents=True, exist_ok=True)
    if source_results is not None:
        copy_existing_results(source_results, config.results_root, mouse_id)
    shutil.copy(config.manifest_csv, config.results_root / "manifest.csv")
    out = prepare_result_dirs(mouse_id, config.results_root)
    asset_info = find_asset_info(config)
    manifest_rows = [ManifestRow.from_series(row) for _, row in manifest_df.iterrows()]
    processed_results = infer_process_results_from_outputs(manifest_rows, out)
    dp = build_datapackage(mouse_id, manifest_rows, processed_results, asset_info, out, config)
    dp_path = write_datapackage(
        dp,
        config.results_root / mouse_id,
        validate=validate,
        asset_roots=[config.data_root],
        asset_overrides=producer_asset_overrides(asset_info, config),
    )
    logger.info("Regenerated datapackage manifest at %s", dp_path)
    return dp_path


def copy_existing_results(source_results: Path, results_root: Path, mouse_id: str) -> Path:
    """Copy a prior mouse output tree into *results_root* and return its destination.

    *source_results* may be either the prior results asset root, containing a
    ``<mouse_id>/`` child, or the mouse output directory itself.
    """
    source_mouse_root = _resolve_existing_mouse_root(Path(source_results), mouse_id)
    destination = Path(results_root) / mouse_id
    if source_mouse_root.resolve() == destination.resolve():
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Copying existing mouse results from %s to %s", source_mouse_root, destination)
    shutil.copytree(source_mouse_root, destination, dirs_exist_ok=True, copy_function=shutil.copy2)
    return destination


def _resolve_existing_mouse_root(source_results: Path, mouse_id: str) -> Path:
    """Resolve source results path to the directory containing datapackage.json."""
    source_results = Path(source_results)
    direct_datapackage = source_results / "datapackage.json"
    if direct_datapackage.is_file():
        return source_results
    nested = source_results / mouse_id
    if (nested / "datapackage.json").is_file() or nested.is_dir():
        return nested
    raise FileNotFoundError(
        "Could not find existing mouse results for "
        f"{mouse_id!r} under {source_results}. Expected either "
        f"{source_results}/datapackage.json or {nested}/datapackage.json."
    )
