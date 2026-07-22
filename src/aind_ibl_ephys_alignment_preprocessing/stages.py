"""Per-stage entry points for the Code Ocean fan-out pipeline.

The monolithic orchestrators (:func:`~aind_ibl_ephys_alignment_preprocessing.pipeline.run_pipeline`
and its async twin) run every step in one process. A Code Ocean Nextflow
pipeline instead runs each step as an independent node, so the work must be
sliced into functions a DAG can dispatch separately. This module provides that
slicing as thin, importable wrappers over the existing per-step functions:

- :func:`stage_discover` (1x) -- the **viability gate**: for each probe decide
  up front whether both its ephys (sorting output) and its histology track
  annotation exist, then emit only the viable work -- a filtered ``manifest.csv``
  for ``histology``/``pack`` and one fan-out config per viable
  ``(recording, ephys_collection)`` (via ``role_dispatch``). The pipeline's
  Flatten edge stages each config to its own ``ephys`` worker. Because the
  decision is made here, the ``histology`` and ``ephys`` stages do not re-check.
- :func:`stage_histology` (1x/mouse) -- registration volumes, additional
  channels, CCF-to-image transforms, and the per-probe coordinate conversion
  (``xyz_picks``). Coordinate conversion depends only on the header-only
  ``sitk`` stubs + the ANTs transforms in the SmartSPIM asset -- **not** on the
  histology volume outputs -- so it rides inside this per-mouse node and no stub
  ever crosses a process boundary.
- :func:`stage_ephys` (1x/``(rec, collection)``) -- find this worker's config
  and extract the one stream to IBL ALF.
- :func:`stage_pack` (1x/mouse) -- assemble ``datapackage.json`` from the
  outputs the upstream nodes wrote.

The Nextflow DAG provides parallelism *across nodes* -- so cross-stage and
cross-recording concurrency (which the monolith did by hand in ``_async``) is
subsumed, and ``discover``/``ephys``/``pack`` are plain synchronous wrappers.
But the DAG does not parallelize *within* a node, so :func:`stage_histology`
keeps the monolith's intra-node async: the heavy ANTs full-volume warps overlap
each other and the per-probe coordinate transforms. It reuses the ``_async``
building blocks and is driven via :func:`asyncio.run` so the dispatcher can call
every stage the same synchronous way. ``ephys``'s async was purely the
cross-recording axis -- now one Nextflow node per ``(rec, collection)`` -- while
its within-extraction threading lives inside ``extract_continuous`` /
``extract_spikes``; ``run_ephys_for_stream`` additionally overlaps those two
independent extractions so a fan-out node stays busy.

Metadata (``processing.json``) emission is a capsule-boundary concern (it needs
the optional ``aind-data-schema`` dependency) and is layered on in the capsule's
entry point via ``aind_code_ocean_pipeline_utils.step``; these library
functions stay free of it.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
from aind_code_ocean_pipeline_utils.role_dispatch import (
    find_stream_config,
    write_stream_configs,
)

from aind_ibl_ephys_alignment_preprocessing.discovery import prepare_result_dirs
from aind_ibl_ephys_alignment_preprocessing.ephys import has_sorting_output, run_ephys_for_stream
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, PipelineConfig, ProcessResult

logger = logging.getLogger(__name__)

# Marker injected into every ephys fan-out config so a worker can recognize its
# own config among CO's unpredictably-nested staged inputs (see role_dispatch).
EPHYS_STREAM_MARKER = "_ibl_ephys_stream_config"


def _ephys_unit_name(recording_id: str, ephys_collection: str | None) -> str:
    """Build a stable fan-out unit name for one ``(recording, collection)``."""
    collection = ephys_collection if ephys_collection is not None else "all"
    return f"{recording_id}__{collection}"


def _track_annotation_present(data_root: Path, mr: ManifestRow) -> tuple[bool, str | None]:
    """Whether the probe's Neuroglancer track annotation exists and has points.

    A light check (no coordinate transforms): the annotation JSON must be found
    under ``data_root`` and its layer named by ``histology_track_id`` must carry
    annotation points. Mirrors the two track-related skip paths in
    :func:`~aind_ibl_ephys_alignment_preprocessing.probes.process_manifest_row`
    (``Annotation not found`` / ``Probe points not found``) so ``discover`` can
    make the same call up front.
    """
    if mr.annotation_format != "json":
        return False, "only JSON annotations supported"
    ann_path = next(data_root.glob(f"*/{mr.probe_file}.json"), None)
    if ann_path is None:
        return False, f"annotation file not found: */{mr.probe_file}.json"

    from aind_s3_cache.json_utils import get_json
    from aind_zarr_utils.neuroglancer import neuroglancer_annotations_to_indices

    probe_id = str(mr.histology_track_id)
    try:
        points_by_layer, _ = neuroglancer_annotations_to_indices(get_json(str(ann_path)), layer_names=[probe_id])
    except Exception as exc:  # malformed JSON / missing or unreadable layer
        return False, f"track layer {probe_id!r} unreadable: {exc}"
    points = points_by_layer.get(probe_id)
    if points is None or len(points) == 0:
        return False, f"no track points for layer {probe_id!r}"
    return True, None


def _probe_viability(config: PipelineConfig, mr: ManifestRow) -> tuple[bool, str | None]:
    """Decide, up front, whether a probe is worth preprocessing.

    A probe is viable only when **both** its ephys and its histology track are
    usable: the upstream spike sorting produced output for the collection
    (``has_sorting_output``) **and** its Neuroglancer track annotation exists
    with points. Consolidating this at the launcher means histology and ephys
    workers never have to re-decide it -- they process only what ``discover``
    deemed viable. Mirrors the monolith's per-row skip, made once and up front.

    Returns
    -------
    tuple[bool, str or None]
        ``(viable, reason)``; ``reason`` is the skip cause when not viable.
    """
    if not config.skip_ephys and mr.ephys_collection is not None:
        if not has_sorting_output(config.data_root / mr.sorted_recording, str(mr.ephys_collection)):
            return False, "no spike-sorting output (bad sorting)"
    return _track_annotation_present(config.data_root, mr)


def stage_discover(config: PipelineConfig) -> list[Path]:
    """Gate probe viability up front, then emit the viable fan-out work.

    The launcher is the single decision point for *what to process*: for each
    manifest row it checks :func:`_probe_viability` (ephys output + track
    annotation both present) and

    - writes a **filtered ``manifest.csv``** (viable rows only) into
      ``config.results_root`` for ``histology`` and ``pack`` to consume, and
    - writes one schema-tagged ephys ``config.json`` per unique viable
      ``(recording, collection)`` (via ``role_dispatch.write_stream_configs``),
      which the pipeline's Flatten edge fans out to :func:`stage_ephys` workers.

    Skipped probes are logged with a reason. Because the decision is made here,
    ``stage_histology`` and ``stage_ephys`` do not re-check viability -- they
    trust the filtered manifest / configs. No ephys configs are emitted when
    ``skip_ephys`` is set.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration. ``data_root`` must have the
        sorted assets (for ``has_sorting_output``) and the annotation JSON (for
        the track check) mounted.

    Returns
    -------
    list[pathlib.Path]
        Paths of the written ephys ``config.json`` files, one per viable unit.
    """
    manifest_df = pd.read_csv(config.manifest_csv)
    rows = [ManifestRow.from_series(row) for _, row in manifest_df.iterrows()]

    viable: list[bool] = []
    for mr in rows:
        ok, reason = _probe_viability(config, mr)
        if not ok:
            logger.warning(
                "[discover] skipping probe %s (%s/%s): %s",
                mr.probe_id,
                mr.recording_id,
                mr.ephys_collection,
                reason,
            )
        viable.append(ok)

    config.results_root.mkdir(parents=True, exist_ok=True)
    filtered_manifest = config.results_root / "manifest.csv"
    manifest_df[pd.Series(viable, index=manifest_df.index)].to_csv(filtered_manifest, index=False)
    logger.info("[discover] %d/%d probes viable -> %s", sum(viable), len(rows), filtered_manifest)

    if config.skip_ephys:
        logger.info("[discover] ephys disabled (skip_ephys); no fan-out configs emitted")
        return []

    seen: set[tuple[str, str | None]] = set()
    items: list[dict[str, object]] = []
    for mr, ok in zip(rows, viable):
        if not ok or mr.ephys_collection is None:
            continue
        key = (mr.recording_id, mr.ephys_collection)
        if key in seen:
            continue
        seen.add(key)
        items.append(
            {
                "name": _ephys_unit_name(mr.recording_id, mr.ephys_collection),
                "mouseid": str(mr.mouseid),
                "sorted_recording": str(mr.sorted_recording),
                "recording_id": mr.recording_id,
                "ephys_collection": mr.ephys_collection,
                "surface_finding": str(mr.surface_finding) if mr.surface_finding is not None else None,
            }
        )

    written = write_stream_configs(
        items,
        results_dir=config.results_root,
        schema_marker=EPHYS_STREAM_MARKER,
        name_key="name",
    )
    logger.info("[discover] wrote %d ephys fan-out config(s)", len(written))
    return written


def stage_histology(config: PipelineConfig) -> list[ProcessResult]:
    """Produce histology volumes and per-probe coordinate outputs for one mouse.

    Runs the registration-channel export, additional-channel processing, the
    CCF-to-image-space transforms, and the per-probe coordinate conversion
    (``xyz_picks`` etc.) -- everything the monolith does except ephys extraction
    and datapackage assembly. Probes whose upstream spike sorting failed are
    skipped (no usable ephys means the track is dropped anyway), avoiding the
    expensive per-probe transforms.

    Unlike ``discover``/``ephys``/``pack``, this stage keeps the monolith's
    **intra-node async parallelism**: the heavy ANTs full-volume warps
    (additional channels, CCF template, CCF labels) overlap each other and the
    per-probe coordinate transforms. That concurrency is *within one node*, so
    the Nextflow DAG -- which parallelizes across nodes -- does not provide it;
    running this stage serially would be a real wall-time regression. The stage
    reuses the ``_async`` building blocks and is driven synchronously via
    :func:`asyncio.run` so the ``--stage`` dispatcher stays uniform.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.

    Returns
    -------
    list[ProcessResult]
        Per-probe coordinate-conversion results.
    """
    import asyncio
    import tempfile
    from concurrent.futures import ThreadPoolExecutor

    from aind_zarr_utils.pipeline_transformed import base_and_pipeline_anatomical_stub
    from aind_zarr_utils.zarr import _open_zarr
    from iblatlas.atlas import AllenAtlas

    from aind_ibl_ephys_alignment_preprocessing._async.concurrency import Limits, to_thread_logged
    from aind_ibl_ephys_alignment_preprocessing._async.ephys import _asyncio_exception_handler
    from aind_ibl_ephys_alignment_preprocessing._async.histology import (
        copy_registration_channel_ccf_reorient_async,
    )
    from aind_ibl_ephys_alignment_preprocessing._async.pipeline import _create_volumes_async
    from aind_ibl_ephys_alignment_preprocessing._async.probes import process_manifest_row_safe_async
    from aind_ibl_ephys_alignment_preprocessing.discovery import find_asset_info
    from aind_ibl_ephys_alignment_preprocessing.types import ReferencePaths, ReferenceVolumes

    async def _coords(
        manifest_df: pd.DataFrame,
        asset_info: Any,
        stub: Any,
        stub_buggy: Any,
        ibl_atlas: Any,
        out: Any,
        limits: Limits,
    ) -> list[ProcessResult]:
        """Convert every probe row's coordinates concurrently (no ephys).

        Viability (ephys output + track present) was already decided by
        ``stage_discover``, so every row in the (filtered) manifest is processed;
        ``process_manifest_row_safe_async`` still returns a graceful skip for any
        residual per-row issue.
        """
        row_tasks: list[tuple[ManifestRow, asyncio.Task[ProcessResult]]] = []
        async with asyncio.TaskGroup() as tg:
            for _, row in manifest_df.iterrows():
                mr = ManifestRow.from_series(row)
                task = tg.create_task(
                    process_manifest_row_safe_async(
                        mr,
                        asset_info,
                        stub,
                        stub_buggy,
                        ibl_atlas,
                        out,
                        limits,
                        config.data_root,
                        emit_qc=config.emit_qc,
                    )
                )
                row_tasks.append((mr, task))

        results: list[ProcessResult] = []
        for mr, row_task in row_tasks:
            result = row_task.result()
            results.append(result)
            if not result.wrote_files:
                logger.warning(
                    "[histology] Did not write coords for %s: %s", mr.sorted_recording, result.skipped_reason
                )
        return results

    async def _run() -> list[ProcessResult]:
        loop = asyncio.get_running_loop()
        loop.set_exception_handler(_asyncio_exception_handler)
        loop.set_default_executor(ThreadPoolExecutor(max_workers=40))

        ref_paths = ReferencePaths.from_config(config)
        scratch_root = Path(config.scratch_root) if config.scratch_root is not None else Path(tempfile.mkdtemp())
        scratch_root.mkdir(parents=True, exist_ok=True)
        limits = Limits(
            scratch_root=str(scratch_root),
            results_root=str(config.results_root),
            data_root=str(config.data_root),
        )

        manifest_df = pd.read_csv(config.manifest_csv)
        mouse_id = str(manifest_df["mouseid"].astype("string").iat[0])

        async with asyncio.TaskGroup() as tg:
            ref_imgs_task = tg.create_task(ReferenceVolumes.from_paths_async(ref_paths))
            asset_info_task = tg.create_task(to_thread_logged(find_asset_info, config))
        ref_imgs = ref_imgs_task.result()
        asset_info = asset_info_task.result()

        out = prepare_result_dirs(mouse_id, config.results_root)
        node, zarr_metadata = _open_zarr(asset_info.zarr_volumes.registration)

        # Stubs (header-only) and the atlas up front. The stubs read the opened
        # zarr node, so compute them before the volume tasks touch it
        # concurrently; coords needs only the stubs + atlas, never the node.
        async with asyncio.TaskGroup() as tg:
            stub_task = tg.create_task(
                to_thread_logged(
                    base_and_pipeline_anatomical_stub,
                    asset_info.zarr_volumes.registration,
                    asset_info.zarr_volumes.metadata,
                    asset_info.zarr_volumes.processing,
                    opened_zarr=(node, zarr_metadata),
                )
            )
            atlas_task = tg.create_task(to_thread_logged(AllenAtlas, 25, hist_path=ref_paths.ibl_atlas_histology_path))
        raw_img_stub, raw_img_stub_buggy, _ = stub_task.result()
        ibl_atlas = atlas_task.result()

        # Volumes (use the zarr node) run concurrently with per-probe coords
        # (use the stubs + atlas) and the emit_qc-gated CCF-space copy.
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
                    emit_qc=config.emit_qc,
                ),
                name=f"histology-volumes-{mouse_id}",
            )
            coords_task = tg.create_task(
                _coords(manifest_df, asset_info, raw_img_stub, raw_img_stub_buggy, ibl_atlas, out, limits),
                name=f"histology-coords-{mouse_id}",
            )
            if config.emit_qc:  # CCF-space registration volume is GUI-unused QC
                tg.create_task(
                    copy_registration_channel_ccf_reorient_async(asset_info, out, limits),
                    name=f"histology-ccf-copy-{mouse_id}",
                )

        logger.info("[histology] Completed volumes + coords for mouse %s", mouse_id)
        return coords_task.result()

    return asyncio.run(_run())


def find_ephys_stream_config(data_root: Path) -> dict[str, Any]:
    """Locate this ephys worker's staged fan-out config under ``data_root``.

    Thin wrapper over ``role_dispatch.find_stream_config`` with the marker
    :data:`EPHYS_STREAM_MARKER`. Exposed so a caller (e.g. the capsule's
    ``processing.json`` node-naming) can read the config once and hand it to
    :func:`stage_ephys` via ``stream_config`` rather than reading it twice.

    Parameters
    ----------
    data_root : Path
        Where the pipeline staged this worker's ``config.json`` (``/data``).

    Returns
    -------
    dict
        The parsed fan-out config for this ``(recording, collection)`` unit.
    """
    _, cfg = find_stream_config(data_root, schema_marker=EPHYS_STREAM_MARKER)
    return cfg


def stage_ephys(config: PipelineConfig, *, stream_config: dict[str, Any] | None = None) -> None:
    """Extract one ``(recording, collection)`` slice to IBL ALF.

    Locates this worker's fan-out config under ``config.data_root`` (the
    schema-tagged ``config.json`` written by :func:`stage_discover` and staged
    by the pipeline), then runs stream-filtered ephys extraction. ``discover``
    only emits configs for viable units (sorting output present), so this stage
    does not re-check ``has_sorting_output`` -- it trusts the config it is given.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration. ``data_root`` must be where the
        staged config and the sorted asset are mounted.
    stream_config : dict or None
        The already-read fan-out config for this unit (from
        :func:`find_ephys_stream_config`). When ``None`` it is read here. Pass
        it to avoid a second filesystem walk when the caller already read it.
    """
    cfg = stream_config if stream_config is not None else find_ephys_stream_config(config.data_root)
    mouse_id = str(cfg["mouseid"])
    sorted_recording = str(cfg["sorted_recording"])
    recording_id = str(cfg["recording_id"])
    ephys_collection = cfg.get("ephys_collection")
    surface_raw = cfg.get("surface_finding")
    surface_finding = Path(str(surface_raw)) if surface_raw else None

    out = prepare_result_dirs(mouse_id, config.results_root)
    run_ephys_for_stream(
        sorted_recording,
        recording_id,
        str(ephys_collection) if ephys_collection is not None else None,
        surface_finding,
        out,
        config.data_root,
        num_parallel_jobs=config.num_parallel_jobs,
    )
    logger.info("[ephys] Completed %s/%s", recording_id, ephys_collection)


def stage_pack(config: PipelineConfig, *, source_results: Path | None = None) -> Path:
    """Assemble ``datapackage.json`` from the upstream nodes' outputs.

    Delegates to
    :func:`~aind_ibl_ephys_alignment_preprocessing.pipeline.regenerate_datapackage`,
    which infers per-probe success from the assembled output tree (it does not
    rerun histology or ephys). In a pipeline the upstream results are mounted
    read-only under ``/data``; pass that path as ``source_results`` so they are
    copied into ``config.results_root`` before the datapackage is written.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.
    source_results : Path or None
        Prior results asset/root to stage into ``config.results_root`` first
        (typically the mounted upstream outputs). ``None`` assumes the outputs
        already live under ``config.results_root``.

    Returns
    -------
    pathlib.Path
        Path to the written ``datapackage.json``.
    """
    from aind_ibl_ephys_alignment_preprocessing.pipeline import regenerate_datapackage

    return regenerate_datapackage(config, source_results=source_results)
