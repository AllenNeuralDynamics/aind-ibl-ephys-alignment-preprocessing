"""Ephys extraction wrappers (synchronous)."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, OutputDirs

logger = logging.getLogger(__name__)


def find_session_dir(data_root: Path, name: str, *, max_depth: int = 4) -> Path | None:
    """Locate a session/recording directory named ``name`` under ``data_root``.

    Walks ``data_root`` (following symlinks -- Code Ocean stages pipeline inputs
    via ``/tmp`` symlink chains) for a directory whose basename is exactly
    ``name``, so callers never assume a flat ``data_root/name`` layout. This
    survives combined data assets (members nested under a wrapper), role-split
    assets (raw and sorted mounted under unrelated parents), and Code Ocean's
    extra input nesting -- replacing the old ``data_root / name`` direct join and
    the sibling-path assumption in ``extract_continuous``. A matched directory is
    not descended into (session trees hold millions of files) and the walk is
    depth-bounded for the same reason.

    Parameters
    ----------
    data_root : Path
        Root under which inputs are mounted (e.g. ``/data``).
    name : str
        Exact directory basename to find (a recording/session folder name).
    max_depth : int, default 4
        Maximum depth below ``data_root`` to descend before pruning a branch.

    Returns
    -------
    Path or None
        The matching directory, or ``None`` if none is found.

    Raises
    ------
    ValueError
        If ``name`` matches two or more distinct real directories (ambiguous).
    """
    base = Path(data_root)
    matches: dict[str, Path] = {}
    for dirpath, dirnames, _filenames in os.walk(base, followlinks=True):
        try:
            depth = len(Path(dirpath).relative_to(base).parts)
        except ValueError:
            depth = 0
        if os.path.basename(dirpath) == name:
            matches.setdefault(os.path.realpath(dirpath), Path(dirpath))
            dirnames[:] = []  # never descend into a matched (huge) session tree
            continue
        if depth >= max_depth:
            dirnames[:] = []
    if not matches:
        return None
    if len(matches) > 1:
        found = sorted(str(p) for p in matches.values())
        raise ValueError(f"ambiguous session dir {name!r} under {base}: {found}")
    return next(iter(matches.values()))


#: Directory suffixes whose contents are an opaque container, never a session.
#:
#: A zarr-backed NWB file is a *directory*, and its internals mimic the structure
#: this search looks for: an NWB processing module is literally named ``ecephys``,
#: so ``<x>.nwb/processing`` satisfies the raw-session marker. A sorted asset
#: carrying three NWB files therefore contributed three phantom "raw sessions"
#: alongside the real one and made the lookup ambiguous (771432, 2026-08-11).
#: Structural search has to skip containers whose insides are not a filesystem
#: layout we own.
_OPAQUE_DIR_SUFFIXES = (".nwb",)


def _find_marked_dirs(data_root: Path, markers: tuple[str, ...], *, max_depth: int = 4) -> dict[str, Path]:
    """Return session dirs identified by a child marker directory (realpath-deduped).

    Code Ocean *pipelines* mount inputs under fixed slot names (``/data/sorted``,
    ``/data/raw``) rather than the asset's own name, so a session cannot be found
    by matching its asset name against a directory basename (see
    :func:`find_session_dir`, which the capsule stages still use). Instead, locate
    it structurally: a directory holding one of ``markers`` as a child (e.g.
    ``spikesorted`` for a sorted asset, ``ecephys_compressed`` for a raw one). A
    matched directory is not descended into and the walk is depth-bounded, matching
    :func:`find_session_dir`'s traversal discipline.

    Directories named ``*.nwb`` are not descended into either -- see
    :data:`_OPAQUE_DIR_SUFFIXES`.
    """
    base = Path(data_root)
    matches: dict[str, Path] = {}
    for dirpath, dirnames, _filenames in os.walk(base, followlinks=True):
        here = Path(dirpath)
        try:
            depth = len(here.relative_to(base).parts)
        except ValueError:
            depth = 0
        # Prune opaque containers before anything else, so neither the marker test
        # nor the depth budget is ever spent inside one.
        dirnames[:] = [d for d in dirnames if not d.endswith(_OPAQUE_DIR_SUFFIXES)]
        if any((here / m).is_dir() for m in markers):
            matches.setdefault(os.path.realpath(dirpath), here)
            dirnames[:] = []  # never descend into the (huge) session tree
            continue
        if depth >= max_depth:
            dirnames[:] = []
    return matches


def _read_data_description_name(session_dir: Path) -> str | None:
    """Return the AIND ``data_description.json`` ``name`` for a session dir, or None."""
    dd = Path(session_dir) / "data_description.json"
    try:
        payload = json.loads(dd.read_text())
    except (OSError, ValueError):
        return None
    name = payload.get("name") if isinstance(payload, dict) else None
    return str(name) if name else None


def find_sorted_session_dir(data_root: Path, *, max_depth: int = 4) -> Path | None:
    """Locate the one spike-sorted session dir mounted under ``data_root`` by content.

    Identifies the sorted asset by its ``spikesorted`` child rather than its asset
    name, so it resolves under a fixed pipeline slot (``/data/sorted``). A per-sort
    ephys pipeline mounts exactly one sorted asset; ``discover`` and raw mounts have
    no ``spikesorted`` child and are ignored.

    Returns
    -------
    Path or None
        The sorted session directory, or ``None`` if none is mounted.

    Raises
    ------
    ValueError
        If two distinct sorted assets are mounted (not a per-sort fan-out layout).
    """
    matches = _find_marked_dirs(data_root, ("spikesorted",), max_depth=max_depth)
    if len(matches) > 1:
        found = sorted(str(p) for p in matches.values())
        raise ValueError(f"ambiguous spike-sorted session under {data_root}: {found}")
    return next(iter(matches.values()), None)


def find_raw_session_dir(data_root: Path, *, recording_id: str | None = None, max_depth: int = 4) -> Path | None:
    """Locate the raw ecephys session dir mounted under ``data_root`` by content.

    Identifies a raw asset by an ``ecephys_compressed`` / ``ecephys`` child, so it
    resolves under a fixed pipeline slot (``/data/raw``). When ``recording_id`` is
    given and more than one raw asset is mounted (e.g. a surface-finding recording
    alongside the main one), the raw whose ``data_description.json`` ``name`` equals
    ``recording_id`` is selected; this is how the main raw is told apart from a
    surface raw. ``None`` is returned when no raw asset is mounted, letting the
    converter fall back to its legacy sibling-of-sorted lookup for the monolith/RR
    path.

    Raises
    ------
    ValueError
        If several raw assets are mounted and none (or more than one) matches
        ``recording_id`` -- the main raw cannot be disambiguated.
    """
    matches = _find_marked_dirs(data_root, ("ecephys_compressed", "ecephys"), max_depth=max_depth)
    if len(matches) <= 1:
        return next(iter(matches.values()), None)
    if recording_id is not None:
        named = [p for p in matches.values() if _read_data_description_name(p) == recording_id]
        if len(named) == 1:
            return named[0]
    found = sorted(str(p) for p in matches.values())
    raise ValueError(f"ambiguous raw ecephys session under {data_root} (recording_id={recording_id!r}): {found}")


def read_sorted_input_recording(sorted_dir: Path) -> str | None:
    """Return a sorted asset's raw input recording name (``recording_id``).

    Reads ``input_data_name`` from the sorted asset's ``data_description.json`` --
    the raw recording the sorting was derived from, which equals ``ManifestRow``'s
    ``recording_id``. Under a fixed pipeline slot the asset name (and thus the
    manifest's ``sorted_recording``) is not recoverable from the mount path, so this
    is how the launcher tells *which* sort is mounted. Returns ``None`` if the field
    is missing or the file is unreadable.
    """
    dd = Path(sorted_dir) / "data_description.json"
    try:
        payload = json.loads(dd.read_text())
    except (OSError, ValueError):
        return None
    name = payload.get("input_data_name") if isinstance(payload, dict) else None
    return str(name) if name else None


def resolve_surface_finding(data_root: Path, surface_finding: Path | str) -> Path:
    """Resolve a surface-finding session folder, tolerant of input nesting.

    Surface-finding data currently lives *inside* another mounted asset (it is not
    registered as a separate asset), so combining the ecephys assets relocates it
    from ``data_root/<name>`` to somewhere nested. Prefer the literal
    ``data_root / surface_finding`` join (unchanged behaviour for a flat layout);
    if that path does not exist, fall back to a name-walk for the folder's
    basename (:func:`find_session_dir`), which finds it under the combined-asset
    nesting. Returns the direct join unchanged when neither resolves, so the
    caller/converter surfaces the missing input rather than this helper.

    Parameters
    ----------
    data_root : Path
        Root under which inputs are mounted.
    surface_finding : Path or str
        The manifest's ``surface_finding`` value (a session folder path/name).

    Returns
    -------
    Path
        The resolved surface-finding session folder.
    """
    direct = data_root / surface_finding
    if direct.exists():
        return direct
    found = find_session_dir(data_root, Path(surface_finding).name)
    return found if found is not None else direct


def has_sorting_output(recording_folder: Path | None, ephys_collection: str) -> bool:
    """Return whether a sorted asset holds postprocessed output for a collection.

    aind-ephys-ibl-gui-conversion builds each collection's ALF table from a
    ``postprocessed/experiment*_<stream>_recording*`` analyzer under the sorted
    asset, where ``<stream>`` embeds the ephys-collection token (e.g.
    ``...ProbeA-AP``). When no such analyzer exists -- the usual sign of failed
    upstream spike sorting -- the converter writes ``sorting_error.txt`` and
    produces no usable ephys, so the probe is not worth processing.

    Parameters
    ----------
    recording_folder : Path or None
        The sorted-recording asset directory (resolved via
        :func:`find_session_dir`), or ``None`` when it could not be located --
        treated the same as missing output (returns *False*).
    ephys_collection : str
        The ephys-collection token (the ALF output folder name / probe stream).

    Returns
    -------
    bool
        *True* if at least one non-LFP postprocessed analyzer matches the
        collection.
    """
    if recording_folder is None:
        return False
    postprocessed = recording_folder / "postprocessed"
    if not postprocessed.is_dir():
        return False
    return any(p.is_dir() and "-LFP" not in p.name for p in postprocessed.glob(f"*{ephys_collection}*"))


def run_ephys_for_recording(
    row: ManifestRow,
    outputs: OutputDirs,
    data_root: Path,
    processed: set[str],
    num_parallel_jobs: int = 4,
) -> None:
    """Run ephys extraction once per unique ``sorted_recording``.

    Creates the results folder under ``results_root/<mouseid>/<recording_id>/``
    and invokes ``extract_continuous`` + ``extract_spikes``.

    Parameters
    ----------
    row : ManifestRow
        Manifest row with ``sorted_recording`` and optional ``surface_finding``.
    outputs : OutputDirs
        Output directory tree.
    data_root : Path
        Root directory containing input data.
    processed : set[str]
        Set of already-processed ``sorted_recording`` strings for idempotency.
    num_parallel_jobs : int
        Number of parallel workers for ``compute_rms`` in ``extract_continuous``.
    """
    # Imported lazily so lightweight consumers (e.g. pre-flight validation)
    # can use ``has_sorting_output`` without pulling in spikeinterface.
    from aind_ephys_ibl_gui_conversion.ephys import extract_continuous, extract_spikes

    sorted_rec = str(row.sorted_recording)
    if sorted_rec in processed:
        return
    processed.add(sorted_rec)

    recording_id = row.recording_id
    mouse_root = outputs.tracks_root.parent
    results_folder = mouse_root / recording_id
    results_folder.mkdir(parents=True, exist_ok=True)

    recording_folder = find_session_dir(data_root, sorted_rec)
    if recording_folder is None:
        raise FileNotFoundError(f"sorted recording {sorted_rec!r} not found under {data_root}")
    # Resolve the raw session folder by name (not by stripping the sorted path),
    # so raw and sorted need not be siblings. None -> converter falls back to its
    # legacy sibling-split, preserving the monolith/RR behaviour.
    session_folder = find_session_dir(data_root, recording_id)

    if row.surface_finding is not None:
        extract_continuous(
            recording_folder,
            results_folder,
            probe_surface_finding=resolve_surface_finding(data_root, row.surface_finding),
            num_parallel_jobs=num_parallel_jobs,
            session_folder=session_folder,
        )
    else:
        extract_continuous(
            recording_folder,
            results_folder,
            num_parallel_jobs=num_parallel_jobs,
            session_folder=session_folder,
        )

    extract_spikes(recording_folder, results_folder, session_folder=session_folder)


def run_ephys_for_stream(
    sorted_recording: str,
    recording_id: str,
    ephys_collection: str | None,
    surface_finding: Path | None,
    outputs: OutputDirs,
    data_root: Path,
    num_parallel_jobs: int = 4,
) -> None:
    """Run ephys extraction for a single ``(recording, ephys_collection)`` slice.

    Unlike :func:`run_ephys_for_recording`, which processes every stream of a
    recording in one call, this restricts extraction to the one stream named by
    ``ephys_collection`` -- the unit of fan-out for the pipeline's ``ephys``
    stage. The collection token is passed straight through as ``stream_to_use``;
    the converter matches it against either the full neo stream name or the
    derived probe/collection token (see
    ``aind_ephys_ibl_gui_conversion.recording_utils._stream_matches``). When
    ``ephys_collection`` is ``None`` every stream is processed, matching
    :func:`run_ephys_for_recording`.

    Output lands under ``results_root/<mouseid>/<recording_id>/<probe_name>/``.
    Because the converter namespaces by probe name, two slices of the same
    recording write disjoint subtrees and never collide -- important for the
    pipeline's Collect fan-in.

    ``extract_continuous`` and ``extract_spikes`` are run **concurrently**: they
    read disjoint inputs (raw ``ecephys_compressed`` vs ``postprocessed``) and
    write disjoint files, so overlapping them keeps a fan-out node's cores busy
    where the monolith ran them serially per recording.

    Parameters
    ----------
    sorted_recording : str
        Name of the spike-sorting folder (relative to ``data_root``).
    recording_id : str
        Recording identifier; names the per-recording results folder.
    ephys_collection : str or None
        ALF/ephys collection token identifying the single stream to process, or
        ``None`` to process every stream in the recording.
    surface_finding : Path or None
        Optional surface-finding file path fragment (relative to ``data_root``).
    outputs : OutputDirs
        Output directory tree.
    data_root : Path
        Root directory containing input data.
    num_parallel_jobs : int
        Number of parallel workers for ``compute_rms`` in ``extract_continuous``.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from aind_ephys_ibl_gui_conversion.ephys import extract_continuous, extract_spikes

    mouse_root = outputs.tracks_root.parent
    results_folder = mouse_root / recording_id
    results_folder.mkdir(parents=True, exist_ok=True)

    # This is the pipeline fan-out worker: exactly one sorted + one raw asset are
    # mounted under fixed slots (/data/sorted, /data/raw), so resolve by content
    # (spikesorted/ marker) -- the asset name, and thus ``sorted_recording``, is not
    # in the mount path. Fall back to name-matching for the monolith/RR path, where
    # assets mount under their own names.
    recording_folder = find_sorted_session_dir(data_root)
    if recording_folder is None:
        recording_folder = find_session_dir(data_root, sorted_recording)
    if recording_folder is None:
        raise FileNotFoundError(
            f"no spike-sorted session (spikesorted/) or {sorted_recording!r} found under {data_root}"
        )
    # Raw session resolved by content too (ecephys_compressed/ marker), told apart
    # from a surface-finding raw by recording_id. None -> converter sibling-split.
    session_folder = find_raw_session_dir(data_root, recording_id=recording_id)
    if session_folder is None:
        session_folder = find_session_dir(data_root, recording_id)

    def _continuous() -> None:
        if surface_finding is not None:
            extract_continuous(
                recording_folder,
                results_folder,
                stream_to_use=ephys_collection,
                probe_surface_finding=resolve_surface_finding(data_root, surface_finding),
                num_parallel_jobs=num_parallel_jobs,
                session_folder=session_folder,
            )
        else:
            extract_continuous(
                recording_folder,
                results_folder,
                stream_to_use=ephys_collection,
                num_parallel_jobs=num_parallel_jobs,
                session_folder=session_folder,
            )

    def _spikes() -> None:
        extract_spikes(
            recording_folder,
            results_folder,
            stream_to_use=ephys_collection,
            session_folder=session_folder,
        )

    # Overlap the two independent extractions. Threads suffice -- the heavy
    # numeric work in each releases the GIL, and each already parallelizes
    # internally (block threads / num_parallel_jobs) -- and they avoid
    # process-pool pickling. Switch to a process pool if spikeinterface
    # global-state contention ever appears.
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(_continuous), executor.submit(_spikes)]
        for future in as_completed(futures):
            future.result()  # re-raise the first failure as soon as it occurs
