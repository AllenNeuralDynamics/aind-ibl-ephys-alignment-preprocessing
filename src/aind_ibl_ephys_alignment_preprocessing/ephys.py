"""Ephys extraction wrappers (synchronous)."""

from __future__ import annotations

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

    extract_spikes(recording_folder, results_folder)


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

    recording_folder = find_session_dir(data_root, sorted_recording)
    if recording_folder is None:
        raise FileNotFoundError(f"sorted recording {sorted_recording!r} not found under {data_root}")
    # Raw session resolved by name, independent of the sorted folder's parent, so
    # combined/role-split asset layouts work. None -> converter sibling-split.
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
        extract_spikes(recording_folder, results_folder, stream_to_use=ephys_collection)

    # Overlap the two independent extractions. Threads suffice -- the heavy
    # numeric work in each releases the GIL, and each already parallelizes
    # internally (block threads / num_parallel_jobs) -- and they avoid
    # process-pool pickling. Switch to a process pool if spikeinterface
    # global-state contention ever appears.
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(_continuous), executor.submit(_spikes)]
        for future in as_completed(futures):
            future.result()  # re-raise the first failure as soon as it occurs
