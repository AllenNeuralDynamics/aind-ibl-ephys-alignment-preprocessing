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

from aind_ibl_ephys_alignment_preprocessing.coverage import (
    CoverageError,
    RunCoverage,
    build_coverage,
)
from aind_ibl_ephys_alignment_preprocessing.discovery import prepare_result_dirs
from aind_ibl_ephys_alignment_preprocessing.ephys import (
    find_session_dir,
    find_sorted_session_dir,
    find_sorted_session_dirs,
    has_sorting_output,
    read_sorted_input_recording,
    run_ephys_for_stream,
)
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, PipelineConfig, ProcessResult

logger = logging.getLogger(__name__)

# Marker injected into every ephys fan-out config so a worker can recognize its
# own config among CO's unpredictably-nested staged inputs (see role_dispatch).
EPHYS_STREAM_MARKER = "_ibl_ephys_stream_config"


def _ephys_unit_name(recording_id: str, ephys_collection: str | None) -> str:
    """Build a stable fan-out unit name for one ``(recording, collection)``."""
    collection = ephys_collection if ephys_collection is not None else "all"
    return f"{recording_id}__{collection}"


def _ephys_config_item(mr: ManifestRow) -> dict[str, object]:
    """Build one ephys fan-out config payload from a manifest row.

    Shared by :func:`stage_discover` (whole-mouse) and
    :func:`stage_ephys_launch` (per-sort) so the config schema stays identical.
    """
    return {
        "name": _ephys_unit_name(mr.recording_id, mr.ephys_collection),
        "mouseid": str(mr.mouseid),
        "sorted_recording": str(mr.sorted_recording),
        "recording_id": mr.recording_id,
        "ephys_collection": mr.ephys_collection,
        "surface_finding": str(mr.surface_finding) if mr.surface_finding is not None else None,
    }


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


def _probe_viability(
    config: PipelineConfig,
    mr: ManifestRow,
    sorted_dirs: dict[str, Path] | None = None,
) -> tuple[bool, str | None, list[str]]:
    """Decide, up front, whether a probe is worth preprocessing.

    A probe is viable only when **both** its ephys and its histology track are
    usable: the upstream spike sorting produced output for the collection
    (``has_sorting_output``) **and** its Neuroglancer track annotation exists
    with points. Consolidating this at the launcher means histology and ephys
    workers never have to re-decide it -- they process only what ``discover``
    deemed viable. Mirrors the monolith's per-row skip, made once and up front.

    Parameters
    ----------
    config : PipelineConfig
        Resolved pipeline configuration.
    mr : ManifestRow
        The row to judge.
    sorted_dirs : dict[str, Path], optional
        ``{recording_id: sorted session dir}`` from
        :func:`~aind_ibl_ephys_alignment_preprocessing.ephys.find_sorted_session_dirs`,
        built once by the caller because it walks ``data_root``. When omitted it is
        built here, which is correct but costs one walk per row.

    Returns
    -------
    tuple[bool, str or None, list[str]]
        ``(viable, reason, notes)``. *reason* is the skip cause when not viable.
        *notes* records non-fatal oddities for the coverage record -- a row can
        be viable and still worth flagging.
    """
    notes: list[str] = []
    if not config.skip_ephys and mr.ephys_collection is not None:
        if sorted_dirs is None:
            sorted_dirs = find_sorted_session_dirs(config.data_root)
        # Match the sorting by the recording it was derived from, not by the pinned
        # name. A sorted asset mounts under its *asset* name, which is free to differ
        # from the manifest's ``sorted_recording`` -- published captures have shipped
        # with the recording time missing from the name. Resolving by name then finds
        # nothing and the whole session reads as unsorted; that is how 791094 lost all
        # 10 probes of 2025-10-09_14-26-31 while the run still exited 0.
        sorted_dir = sorted_dirs.get(mr.recording_id)
        if sorted_dir is None:
            # No sorting identified itself as this recording's. Fall back to the old
            # name lookup so an asset with an unreadable data_description still works.
            sorted_dir = find_session_dir(config.data_root, str(mr.sorted_recording))
            if sorted_dir is not None:
                logger.warning(
                    "[discover] %s: no mounted sorting names this recording in its "
                    "data_description; matched %s by directory name instead",
                    mr.recording_id,
                    sorted_dir,
                )
                notes.append(f"sorting matched by directory name ({sorted_dir.name}), not by input recording")
        if not has_sorting_output(sorted_dir, str(mr.ephys_collection)):
            return False, "no spike-sorting output (bad sorting)", notes
    ok, reason = _track_annotation_present(config.data_root, mr)
    return ok, reason, notes


def stage_discover(
    config: PipelineConfig,
    *,
    producer_record: dict[str, Any] | None = None,
) -> list[Path]:
    """Gate probe viability up front, then emit the viable fan-out work.

    The launcher is the single decision point for *what to process*: for each
    manifest row it checks :func:`_probe_viability` (ephys output + track
    annotation both present) and

    - writes a **filtered ``manifest.csv``** (viable rows only) into
      ``config.results_root`` for ``histology`` and ``pack`` to consume,
    - writes one schema-tagged ephys ``config.json`` per unique viable
      ``(recording, collection)`` (via ``role_dispatch.write_stream_configs``),
      which the pipeline's Flatten edge fans out to :func:`stage_ephys` workers,
      and
    - writes a **coverage record** (``coverage.json``) naming every row it read,
      the unit each surviving row went into, and why each dropped row was
      dropped -- then **refuses** unless ``config.allow_partial``.

    The refusal is the point. Rows discarded here leave no trace in the filtered
    manifest or the fan-out set, so every downstream check compares two sets that
    have already shrunk together and agrees with itself; that is how a whole
    session was lost with exit 0. See
    :mod:`aind_ibl_ephys_alignment_preprocessing.coverage`.

    Parameters
    ----------
    config : PipelineConfig
        Resolved pipeline configuration.
    producer_record : dict[str, Any], optional
        This launcher's provenance shard (see
        :func:`aind_code_ocean_pipeline_utils.records.make_record`). When given, it
        is side-written into each fan-out unit's ``provenance/`` so an ephys worker
        can infer ``discover`` as its parent by frontier. A Flatten fan-out
        distributes per-item slices but does not broadcast the launcher's shared
        ``provenance/``, so this per-unit copy is how the edge is carried across.

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

    Raises
    ------
    aind_ibl_ephys_alignment_preprocessing.coverage.CoverageError
        If any manifest row was dropped and ``config.allow_partial`` is unset,
        or if no row survived at all (fatal either way).
    """
    manifest_df = pd.read_csv(config.manifest_csv)
    rows = [ManifestRow.from_series(row) for _, row in manifest_df.iterrows()]

    # Built once: it walks data_root, and every row consults it.
    sorted_dirs = {} if config.skip_ephys else find_sorted_session_dirs(config.data_root)

    # Keep the verdict *with* its row. Dropping the reason on the floor -- it used
    # to go to a log line and nowhere else -- is what made a lost session require
    # run-log archaeology instead of a file read.
    judged: list[tuple[ManifestRow, bool, str | None, list[str]]] = []
    for mr in rows:
        ok, reason, notes = _probe_viability(config, mr, sorted_dirs)
        if not ok:
            logger.warning(
                "[discover] skipping probe %s (%s/%s): %s",
                mr.probe_id,
                mr.recording_id,
                mr.ephys_collection,
                reason,
            )
        judged.append((mr, ok, reason, notes))
    viable = [ok for _, ok, _, _ in judged]

    config.results_root.mkdir(parents=True, exist_ok=True)
    filtered_manifest = config.results_root / "manifest.csv"
    manifest_df[pd.Series(viable, index=manifest_df.index)].to_csv(filtered_manifest, index=False)
    logger.info("[discover] %d/%d probes viable -> %s", sum(viable), len(rows), filtered_manifest)

    # Assign units before writing anything, so the coverage record can name the
    # unit each surviving row was folded into -- and so it is written on the
    # ``skip_ephys`` path too, which is precisely where no downstream fan-out
    # exists to notice a gap.
    seen: set[tuple[str, str | None]] = set()
    items: list[dict[str, object]] = []
    unit_by_row: dict[int, str] = {}
    if not config.skip_ephys:
        for index, (mr, ok, _, _) in enumerate(judged):
            if not ok or mr.ephys_collection is None:
                continue
            unit_by_row[index] = _ephys_unit_name(mr.recording_id, mr.ephys_collection)
            key = (mr.recording_id, mr.ephys_collection)
            if key in seen:
                continue
            seen.add(key)
            items.append(_ephys_config_item(mr))

    coverage = build_coverage(
        mouse_id=str(rows[0].mouseid) if rows else "",
        manifest_csv=config.manifest_csv,
        judged=judged,
        unit_by_row=unit_by_row,
        skip_ephys=config.skip_ephys,
    )
    coverage_path = coverage.write(config.results_root)
    logger.info(
        "[discover] coverage: %d/%d rows viable, %d unit(s) -> %s",
        coverage.rows_viable,
        coverage.rows_requested,
        len(coverage.units_emitted),
        coverage_path,
    )

    written: list[Path] = []
    if config.skip_ephys:
        logger.info("[discover] ephys disabled (skip_ephys); no fan-out configs emitted")
    else:
        written = write_stream_configs(
            items,
            results_dir=config.results_root,
            schema_marker=EPHYS_STREAM_MARKER,
            name_key="name",
            producer_record=producer_record,
        )
        logger.info("[discover] wrote %d ephys fan-out config(s)", len(written))

    # Refuse *after* writing, so the record that explains the refusal survives it.
    _assert_coverage(coverage, allow_partial=config.allow_partial)
    return written


def _assert_coverage(coverage: RunCoverage, *, allow_partial: bool) -> None:
    """Stop the run when it would cover less than the manifest asked for.

    Mirrors the trigger's ``_require_full_coverage`` so the two layers enforce one
    policy, and is deliberately the *default*: a dropped row does not crash
    anything on its own, it just makes the run finish green over less work.

    Zero viable rows stays fatal regardless of *allow_partial* -- ``RUNBOOK.md``
    already calls it the known failure mode, and "process nothing successfully"
    is never what a partial run means.

    An **empty manifest** is fatal on its own terms, checked before anything else.
    It has no dropped rows to report, so every gap-based check passes it: an
    earlier version guarded the zero-viable branch with ``rows_requested`` to
    avoid a confusing message, and bought a hole in exactly the failure class
    this function exists to close -- a header-only or mis-generated manifest ran
    green over nothing at all. A manifest that lists no work is a generation
    failure upstream, never a request for an empty run.
    """
    if not coverage.rows_requested:
        raise CoverageError(
            f"{coverage.manifest} lists no manifest rows, so there is nothing to process. "
            "An empty manifest is an upstream generation failure; allow_partial does not apply "
            "(there is no part to run)."
        )
    if not coverage.rows_viable:
        raise CoverageError(
            f"no viable probes in {coverage.manifest} ({coverage.rows_requested} row(s) requested); "
            "nothing to process. Reasons: " + "; ".join(coverage.gaps())
        )
    gaps = coverage.gaps()
    if not gaps:
        return
    listed = "; ".join(gaps)
    if allow_partial:
        logger.warning("[discover] proceeding with partial coverage (allow_partial): %s", listed)
        return
    raise CoverageError(
        f"discover covers {coverage.rows_viable}/{coverage.rows_requested} manifest row(s): {listed}. "
        "Fix the manifest or the assets, or set allow_partial to accept a partial run."
    )


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
                    output_voxel_size_um=config.output_voxel_size_um,
                    emit_qc=config.emit_qc,
                ),
                name=f"histology-volumes-{mouse_id}",
            )
            coords_task = tg.create_task(
                _coords(manifest_df, asset_info, raw_img_stub, raw_img_stub_buggy, ibl_atlas, out, limits),
                name=f"histology-coords-{mouse_id}",
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

    Output is written to ``config.results_root/<unit_name>/<mouse_id>/...``
    where ``unit_name`` is this fan-out unit's ``<recording>__<collection>``.
    The per-unit namespace keeps parallel ephys instances' top-level ``/results``
    names disjoint so the downstream Collect (``pack``) node does not abort on an
    "input file name collision"; ``pack``'s merge unions the nested trees back.

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

    # Namespace this worker's output by its fan-out unit so parallel ephys
    # instances never share a top-level ``/results`` name. Every worker would
    # otherwise write ``/results/<mouse_id>/``, which makes the downstream
    # Collect (``pack``) node abort at staging with an "input file name
    # collision" -- after every worker has already run. Writing under
    # ``<unit_name>/`` keeps the names disjoint; ``pack``'s layout-agnostic
    # merge (``merge_pipeline_outputs``) finds the nested ``<mouse_id>/`` tree
    # and unions them back together.
    unit_name = _ephys_unit_name(recording_id, str(ephys_collection) if ephys_collection is not None else None)
    out = prepare_result_dirs(mouse_id, config.results_root / unit_name)
    run_ephys_for_stream(
        sorted_recording,
        recording_id,
        str(ephys_collection) if ephys_collection is not None else None,
        surface_finding,
        out,
        config.data_root,
        num_parallel_jobs=config.num_parallel_jobs,
    )
    logger.info("[ephys] Completed %s/%s -> %s", recording_id, ephys_collection, unit_name)


def stage_ephys_launch(config: PipelineConfig) -> list[Path]:
    """Emit fan-out configs for the one sort mounted here (per-sort ephys launcher).

    The split architecture runs one ephys pipeline per sort. This launcher reads
    ``discover``'s filtered ``manifest.csv`` and writes one ephys ``config.json``
    per unique ``(recording, collection)`` -- but only for recordings whose
    **sorted asset is actually mounted** under ``config.data_root``. A per-sort
    ephys pipeline mounts exactly one sorted asset, so that scopes the fan-out to
    this sort's probes; the pipeline's Flatten edge then dispatches one
    :func:`stage_ephys` worker per config. No viability re-check (``discover``
    already filtered) and no manifest is written.

    Parameters
    ----------
    config : PipelineConfig
        Resolved configuration. ``manifest_csv`` must point at ``discover``'s
        filtered manifest; ``data_root`` must have this sort's sorted asset mounted.

    Returns
    -------
    list[pathlib.Path]
        Paths of the written ephys ``config.json`` files (this sort's units).
    """
    manifest_df = pd.read_csv(config.manifest_csv)
    rows = [ManifestRow.from_series(row) for _, row in manifest_df.iterrows()]

    # Scope to the sort mounted here. A per-sort ephys pipeline runs as pipeline
    # *nodes* (fixed slots: /data/sorted), so the asset name -- and thus the
    # manifest's ``sorted_recording`` -- is not in the mount path and cannot be
    # matched by directory basename. Identify the mounted sort by content instead:
    # its ``data_description.json`` ``input_data_name`` equals ``recording_id``.
    sorted_dir = find_sorted_session_dir(config.data_root)
    if sorted_dir is None:
        raise FileNotFoundError(f"[ephys-launch] no spike-sorted asset (spikesorted/) mounted under {config.data_root}")
    mounted_recording_id = read_sorted_input_recording(sorted_dir)
    if mounted_recording_id is None:
        raise ValueError(f"[ephys-launch] could not read input_data_name from {sorted_dir}/data_description.json")

    scoped = [mr for mr in rows if mr.ephys_collection is not None and mr.recording_id == mounted_recording_id]
    distinct_sorts = {str(mr.sorted_recording) for mr in scoped}
    if len(distinct_sorts) > 1:
        raise ValueError(
            f"[ephys-launch] recording_id {mounted_recording_id!r} maps to multiple sorts "
            f"{sorted(distinct_sorts)}; cannot tell which is mounted from asset content alone"
        )

    seen: set[tuple[str, str | None]] = set()
    items: list[dict[str, object]] = []
    for mr in scoped:
        key = (mr.recording_id, mr.ephys_collection)
        if key in seen:
            continue
        seen.add(key)
        items.append(_ephys_config_item(mr))

    config.results_root.mkdir(parents=True, exist_ok=True)
    written = write_stream_configs(
        items,
        results_dir=config.results_root,
        schema_marker=EPHYS_STREAM_MARKER,
        name_key="name",
    )
    if not written:
        raise RuntimeError(
            f"[ephys-launch] wrote 0 fan-out configs for mounted sort {mounted_recording_id!r}; "
            f"manifest recording_ids present: {sorted({mr.recording_id for mr in rows})}"
        )
    logger.info(
        "[ephys-launch] wrote %d fan-out config(s) for mounted sort %s",
        len(written),
        mounted_recording_id,
    )
    return written


def stage_ephys_collect(config: PipelineConfig, *, merge_from: Path | None = None) -> None:
    """Union one sort's fanned-out ephys ALF into ``/results`` (per-sort collector).

    The ephys pipeline's Collect node. Each :func:`stage_ephys` worker wrote a
    disjoint ``<unit>/<mouse_id>/`` subtree; this union-merges every ``<mouse_id>/``
    tree found under ``merge_from`` (default ``config.data_root``) into
    ``config.results_root/<mouse_id>/`` for capture as this sort's ephys
    intermediate asset. No datapackage is built -- that is the mouse-level
    :func:`stage_pack`.

    Parameters
    ----------
    config : PipelineConfig
        Resolved configuration. ``manifest_csv`` supplies the mouse id.
    merge_from : Path or None
        Mount holding the fanned-out worker outputs (default ``config.data_root``).
    """
    from aind_ibl_ephys_alignment_preprocessing.datapackage import merge_pipeline_outputs

    source = merge_from if merge_from is not None else config.data_root
    manifest_df = pd.read_csv(config.manifest_csv)
    mouse_id = str(manifest_df["mouseid"].astype("string").iat[0])
    merge_pipeline_outputs(source, config.results_root, mouse_id)
    logger.info("[ephys-collect] merged ephys outputs for mouse %s", mouse_id)


def stage_pack(
    config: PipelineConfig,
    *,
    source_results: Path | None = None,
    merge_from: Path | None = None,
    asset_resolution: str | None = None,
) -> Path:
    """Assemble ``datapackage.json`` from the upstream nodes' outputs.

    This is the pipeline's **fan-in** node. The histology node and each ephys
    fan-out node write disjoint slices of one ``<mouse_id>/`` tree; ``pack``
    unions them back together and builds the datapackage over the complete
    output (inferring per-probe success from the assembled tree -- it does not
    rerun histology or ephys).

    Two staging modes, mutually exclusive:

    - ``merge_from`` -- the pipeline fan-in path. All upstream captured results
      are mounted read-only under this root (``/data``), each an independent
      ``<mouse_id>/`` subtree (possibly one "indexed folders" level deep).
      :func:`~aind_ibl_ephys_alignment_preprocessing.datapackage.merge_pipeline_outputs`
      union-merges every one it finds into ``config.results_root/<mouse_id>/``
      before the datapackage is written.
    - ``source_results`` -- the single-tree regenerate path (rewrite metadata
      for one prior results asset). Delegated to ``regenerate_datapackage``.

    With neither set, the outputs are assumed to already live under
    ``config.results_root``. ``config.manifest_csv`` must point at ``discover``'s
    filtered manifest so the datapackage covers only viable probes.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.
    source_results : Path or None
        Single prior results asset/root to stage in first (regenerate path).
    merge_from : Path or None
        Mount holding all upstream node outputs to union-merge (pipeline path).
    asset_resolution : str or None
        The launcher's resolved-asset JSON, merged into the run record (see
        :func:`_carry_coverage_forward`). Asset ids are resolved before any
        capsule is launched and are invisible from inside one, so this is the
        only way they can reach the output.

    Returns
    -------
    pathlib.Path
        Path to the written ``datapackage.json``.
    """
    from aind_ibl_ephys_alignment_preprocessing.datapackage import merge_pipeline_outputs
    from aind_ibl_ephys_alignment_preprocessing.pipeline import regenerate_datapackage

    if merge_from is not None:
        if source_results is not None:
            raise ValueError("stage_pack: pass either merge_from or source_results, not both")
        manifest_df = pd.read_csv(config.manifest_csv)
        mouse_id = str(manifest_df["mouseid"].astype("string").iat[0])
        merge_pipeline_outputs(merge_from, config.results_root, mouse_id)
        dp_path = regenerate_datapackage(config, source_results=None)
    else:
        dp_path = regenerate_datapackage(config, source_results=source_results)

    _carry_coverage_forward(config, asset_resolution)
    return dp_path


def _carry_coverage_forward(config: PipelineConfig, asset_resolution: str | None) -> None:
    """Copy ``discover``'s run record into this run's output, with the asset ids added.

    Two reasons this is not merely a nicety:

    - **Durability.** ``discover``'s ``coverage.json`` lives in an intermediate
      asset the launcher *deletes on success*, so without this step the record of
      what a successful run covered is destroyed by that run succeeding. What
      survives is a computation log nobody thinks to fetch -- which is precisely
      why diagnosing the 791094 loss cost what it did.
    - **Completion.** The asset ids only exist out at the launcher; discover
      never sees one.

    Deliberately *beside* ``datapackage.json`` rather than inside it.
    ``datapackage.json`` is the alignment GUI's consumption contract -- what to
    read and where -- and nothing here is consumed by it. Folding this in would
    move a schema version, and therefore every consumer's version gate, for a
    question no consumer asks.

    Best-effort: a run that has already done the work must not fail over its own
    bookkeeping. Absence is logged, not raised -- the monolith has no discover
    stage, and a rerun against an older ``discover`` asset finds no record.
    """
    try:
        # discover's mount: the directory the filtered manifest resolved from.
        discovered = RunCoverage.find(config.manifest_csv.parent)
        if discovered is None:
            logger.info("[pack] no discover coverage record found; run record not carried forward")
            return
        coverage = discovered.with_assets(asset_resolution)
        path = coverage.write(config.results_root)
        logger.info(
            "[pack] run record -> %s (%d/%d rows viable, %d unit(s) with resolved assets)",
            path,
            coverage.rows_viable,
            coverage.rows_requested,
            len(coverage.unit_assets),
        )
    except Exception as exc:  # bookkeeping must never fail a completed run
        logger.warning("[pack] could not write the run record: %s", exc)
