"""Async pipeline orchestrator."""

from __future__ import annotations

import asyncio
import logging
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any

import ants
import pandas as pd
from aind_zarr_utils.zarr import _open_zarr

from aind_ibl_ephys_alignment_preprocessing._async.concurrency import Limits, io_to_thread_on, to_thread_logged
from aind_ibl_ephys_alignment_preprocessing._async.ephys import (
    _asyncio_exception_handler,
    run_manifest_subprocess_sync,
)
from aind_ibl_ephys_alignment_preprocessing._async.histology import (
    copy_registration_channel_ccf_reorient_async,
    process_additional_channels_pipeline_async,
    transform_ccf_labels_to_image_space_async,
    transform_ccf_to_image_space_async,
    write_registration_channel_images_async,
)
from aind_ibl_ephys_alignment_preprocessing.discovery import (
    determine_desired_level,
    find_asset_info,
    prepare_result_dirs,
)
from aind_ibl_ephys_alignment_preprocessing.types import (
    PipelineConfig,
    ProcessResult,
    ReferencePaths,
    ReferenceVolumes,
)

logger = logging.getLogger(__name__)

IO_THREADS = 40


async def _create_volumes_async(
    asset_info: Any,
    ref_imgs: ReferenceVolumes,
    ref_paths: ReferencePaths,
    out: Any,
    node: Any,
    zarr_metadata: dict[str, Any],
    limits: Limits,
    scratch_root: Path,
    desired_voxel_size_um: float = 25.0,
) -> None:
    """Create all volume outputs (registration channel + additional + CCF transforms)."""
    logger.info("[Histology] Starting volume processing")
    level = determine_desired_level(zarr_metadata, desired_voxel_size_um=desired_voxel_size_um)
    num_additional = len(asset_info.zarr_volumes.additional)
    logger.info(
        "[Histology] Processing registration channel (level %d) + %d additional channel(s)", level, num_additional
    )
    raw_img_path, pipeline_img_path = await write_registration_channel_images_async(
        asset_info, out, limits, level=level, opened_zarr=(node, zarr_metadata)
    )
    logger.info(
        "[Histology] Registration channel export complete: raw=%s, pipeline=%s",
        raw_img_path.name,
        pipeline_img_path.name,
    )
    async with asyncio.TaskGroup() as tg:
        pipeline_img_ants_task = tg.create_task(
            io_to_thread_on(limits, str(pipeline_img_path), ants.image_read, str(pipeline_img_path), pixeltype=None),
            name="load-ants-pipeline-img",
        )
        raw_img_ants_task = tg.create_task(
            io_to_thread_on(limits, str(raw_img_path), ants.image_read, str(raw_img_path), pixeltype=None),
            name="load-ants-raw-img",
        )
    pipeline_img_ants = pipeline_img_ants_task.result()
    raw_img_ants = raw_img_ants_task.result()
    logger.info(
        "[Histology] Starting parallel processing: %d additional channel(s), CCF template + labels transforms",
        num_additional,
    )
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            process_additional_channels_pipeline_async(
                pipeline_img_ants, asset_info, ref_imgs, out, limits, scratch_root=scratch_root, level=level
            ),
            name="process-additional-channels",
        )
        tg.create_task(
            transform_ccf_to_image_space_async(asset_info, ref_imgs, raw_img_ants, pipeline_img_ants, out, limits),
            name="transform-ccf-template-to-image",
        )
        tg.create_task(
            transform_ccf_labels_to_image_space_async(
                asset_info, ref_paths, raw_img_ants, pipeline_img_ants, out, limits
            ),
            name="transform-ccf-labels-to-image",
        )
    logger.info("[Histology] All volume processing complete")


async def run_pipeline_async(config: PipelineConfig, max_workers: int = 40) -> list[ProcessResult]:
    """Run the full preprocessing pipeline asynchronously.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.
    max_workers : int
        Number of thread pool workers.

    Returns
    -------
    list[ProcessResult]
        Per-probe processing results.
    """
    loop = asyncio.get_running_loop()
    loop.set_debug(True)
    loop.set_exception_handler(_asyncio_exception_handler)
    loop.set_default_executor(ThreadPoolExecutor(max_workers=max_workers))

    ref_paths = ReferencePaths.from_config(config)

    scratch_root = Path(config.scratch_root) if config.scratch_root is not None else Path(tempfile.mkdtemp())
    scratch_root.mkdir(parents=True, exist_ok=True)
    scratch_root_str = str(scratch_root)
    results_root_str = str(config.results_root)
    data_root_str = str(config.data_root)

    limits = Limits(
        scratch_root=scratch_root_str,
        results_root=results_root_str,
        data_root=data_root_str,
    )

    # Keep a manifest snapshot for reproducibility
    shutil.copy(config.manifest_csv, config.results_root / "manifest.csv")

    manifest_df = pd.read_csv(config.manifest_csv)
    mouse_id: str = str(manifest_df["mouseid"].astype("string").iat[0])
    num_probes = len(manifest_df)
    logger.info("[Orchestrator] Starting histology and ephys processing for mouse: %s", mouse_id)
    logger.info("[Orchestrator] Manifest contains %d probe(s)", num_probes)

    async with asyncio.TaskGroup() as tg:
        ref_imgs_task = tg.create_task(ReferenceVolumes.from_paths_async(ref_paths), name="load-ref-volumes")
        asset_info_task = tg.create_task(to_thread_logged(find_asset_info, config), name="find-asset-info")
    ref_imgs = ref_imgs_task.result()
    asset_info = asset_info_task.result()

    out = prepare_result_dirs(mouse_id, config.results_root)

    manifest_pool = ProcessPoolExecutor(max_workers=1)
    node, zarr_metadata = _open_zarr(asset_info.zarr_volumes.registration)

    skip_ephys_msg = " (ephys disabled)" if config.skip_ephys else ""
    logger.info(
        "[Orchestrator] Launching 3 parallel task groups: volumes, manifest (%d probes) in subprocess, CCF copy%s",
        num_probes,
        skip_ephys_msg,
    )
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            _create_volumes_async(
                asset_info,
                ref_imgs,
                ref_paths,
                out,
                node,
                zarr_metadata,
                limits,
                scratch_root=scratch_root,
                desired_voxel_size_um=config.desired_voxel_size_um,
            ),
            name=f"create-volumes-{mouse_id}",
        )

        async def _run_manifest_in_subprocess() -> list[ProcessResult]:
            return await loop.run_in_executor(
                manifest_pool,
                run_manifest_subprocess_sync,
                manifest_df,
                asset_info,
                ref_paths,
                out,
                config,
                limits.max_ephys,
                limits.max_manifest_rows,
                limits.max_scratch,
                limits.max_results,
                limits.max_data,
            )

        manifest_task = tg.create_task(
            _run_manifest_in_subprocess(),
            name=f"process-manifest-subprocess-{mouse_id}",
        )

        tg.create_task(
            copy_registration_channel_ccf_reorient_async(asset_info, out, limits),
            name=f"copy-ccf-registration-{mouse_id}",
        )
    logger.info("[Orchestrator] All parallel tasks completed")
    manifest_pool.shutdown(wait=True)
    processed_results = manifest_task.result()
    num_succeeded = sum(1 for r in processed_results if r.wrote_files)
    num_failed = len(processed_results) - num_succeeded
    logger.info("[Orchestrator] Pipeline complete: %d succeeded, %d failed", num_succeeded, num_failed)
    return processed_results
