"""Async ephys coordination (process pool with single-flight dedup)."""

from __future__ import annotations

import asyncio
import logging
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Any

import pandas as pd
from aind_zarr_utils.pipeline_transformed import base_and_pipeline_anatomical_stub
from aind_zarr_utils.zarr import _open_zarr
from iblatlas.atlas import AllenAtlas

from aind_ibl_ephys_alignment_preprocessing._async.concurrency import (
    EphysCoordinator,
    Limits,
    _run_ephys_sync,
    io_to_thread_on,
    to_thread_logged,
)
from aind_ibl_ephys_alignment_preprocessing._async.probes import process_manifest_row_safe_async
from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    ManifestRow,
    OutputDirs,
    PipelineConfig,
    ProcessResult,
    ReferencePaths,
)

logger = logging.getLogger(__name__)

EPROCS = int(os.environ.get("EPROCS", "8"))
IO_THREADS = 40


async def process_manifest_async(
    manifest_df: pd.DataFrame,
    asset_info: AssetInfo,
    ibl_atlas: AllenAtlas,
    out: OutputDirs,
    node: Any,
    zarr_metadata: dict[str, Any],
    config: PipelineConfig,
    ephys: EphysCoordinator,
    limits: Limits,
) -> list[ProcessResult]:
    """Async manifest processing: all probes + ephys in parallel.

    Parameters
    ----------
    manifest_df : pd.DataFrame
        Manifest data.
    asset_info : AssetInfo
        Asset metadata.
    ibl_atlas : AllenAtlas
        Allen atlas instance.
    out : OutputDirs
        Output directory tree.
    node : Any
        OME-Zarr node.
    zarr_metadata : dict
        Zarr metadata.
    config : PipelineConfig
        Pipeline configuration.
    ephys : EphysCoordinator
        Ephys coordinator for single-flight dedup.
    limits : Limits
        Concurrency limits.

    Returns
    -------
    list[ProcessResult]
        Per-probe processing results.
    """
    num_probes = len(manifest_df)
    logger.info("[Manifest] Starting manifest processing: %d probe(s)", num_probes)
    raw_img_stub, raw_img_stub_buggy, _ = await to_thread_logged(
        base_and_pipeline_anatomical_stub,
        asset_info.zarr_volumes.registration,
        asset_info.zarr_volumes.metadata,
        asset_info.zarr_volumes.processing,
        opened_zarr=(node, zarr_metadata),
    )

    if config.skip_ephys:
        logger.info("[Manifest] Ephys processing disabled via skip_ephys flag")

    row_tasks: list[tuple[ManifestRow, asyncio.Task[ProcessResult]]] = []
    logger.info("[Manifest] Creating parallel tasks for all %d probe(s)", num_probes)
    async with asyncio.TaskGroup() as tg:
        for _, row in manifest_df.iterrows():
            mr = ManifestRow.from_series(row)
            t = tg.create_task(
                process_manifest_row_safe_async(
                    mr, asset_info, raw_img_stub, raw_img_stub_buggy, ibl_atlas, out, limits, config.data_root
                ),
                name=f"probe-{mr.probe_id}-{mr.recording_id}",
            )
            row_tasks.append((mr, t))
            if not config.skip_ephys:
                tg.create_task(
                    ephys.ensure(
                        key=mr.sorted_recording,
                        run_sync=_run_ephys_sync,
                        mr=mr,
                        out=out,
                        data_root=config.data_root,
                    ),
                    name=f"ephys-ensure-{mr.sorted_recording}",
                )

    processed_results: list[ProcessResult] = []
    for mr, rt in row_tasks:
        result = rt.result()
        processed_results.append(result)
        if not result.wrote_files:
            logger.warning("Did not write files for %s: %s", mr.sorted_recording, result.skipped_reason)
    num_succeeded = sum(1 for r in processed_results if r.wrote_files)
    num_failed = len(processed_results) - num_succeeded
    logger.info("[Manifest] Completed: %d succeeded, %d failed", num_succeeded, num_failed)
    return processed_results


def _asyncio_exception_handler(loop: asyncio.AbstractEventLoop, context: dict[str, Any]) -> None:
    msg = context.get("message")
    exc = context.get("exception")
    logging.error("Asyncio exception: %s", msg or exc, exc_info=exc)
    if exc:
        loop.call_soon_threadsafe(lambda: (_ for _ in ()).throw(exc))


def run_manifest_subprocess_sync(
    manifest_df: pd.DataFrame,
    asset_info: AssetInfo,
    ref_paths: ReferencePaths,
    out: OutputDirs,
    config: PipelineConfig,
    max_ephys: int | None,
    max_manifest_rows: int | None,
    max_scratch: int | None,
    max_results: int | None,
    max_data: int | None,
) -> list[ProcessResult]:
    """Subprocess entry point for manifest processing.

    Recreates heavy objects (AllenAtlas, EphysCoordinator, Limits) in the
    subprocess, then runs the async manifest processing in a new event loop.
    """

    async def _async_main() -> list[ProcessResult]:
        loop = asyncio.get_running_loop()
        loop.set_debug(True)
        loop.set_exception_handler(_asyncio_exception_handler)
        loop.set_default_executor(ThreadPoolExecutor(max_workers=IO_THREADS))

        scratch_root = str(config.scratch_root) if config.scratch_root else "/scratch"
        results_root = str(config.results_root)
        data_root = str(config.data_root)

        limits = Limits(
            max_ephys=max_ephys,
            max_registration=None,
            max_manifest_rows=max_manifest_rows,
            max_scratch=max_scratch,
            max_results=max_results,
            max_data=max_data,
            scratch_root=scratch_root,
            results_root=results_root,
            data_root=data_root,
        )

        logger.info("[Subprocess] Loading AllenAtlas in subprocess")
        ibl_atlas = await io_to_thread_on(
            limits,
            str(ref_paths.ibl_atlas_histology_path),
            AllenAtlas,
            25,
            hist_path=ref_paths.ibl_atlas_histology_path,
        )

        logger.info("[Subprocess] Opening zarr in subprocess")
        node, zarr_metadata = await to_thread_logged(_open_zarr, asset_info.zarr_volumes.registration)

        ephys_pool = ProcessPoolExecutor(max_workers=EPROCS)
        ephys = EphysCoordinator(pool=ephys_pool, max_inflight=2)

        try:
            logger.info("[Subprocess] Starting manifest processing")
            result = await process_manifest_async(
                manifest_df, asset_info, ibl_atlas, out, node, zarr_metadata, config, ephys, limits
            )
            return result
        finally:
            ephys_pool.shutdown(wait=True)

    return asyncio.run(_async_main())
