"""Ephys extraction wrappers (synchronous)."""

from __future__ import annotations

import logging
from pathlib import Path

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, OutputDirs

logger = logging.getLogger(__name__)


def has_sorting_output(recording_folder: Path, ephys_collection: str) -> bool:
    """Return whether a sorted asset holds postprocessed output for a collection.

    aind-ephys-ibl-gui-conversion builds each collection's ALF table from a
    ``postprocessed/experiment*_<stream>_recording*`` analyzer under the sorted
    asset, where ``<stream>`` embeds the ephys-collection token (e.g.
    ``...ProbeA-AP``). When no such analyzer exists -- the usual sign of failed
    upstream spike sorting -- the converter writes ``sorting_error.txt`` and
    produces no usable ephys, so the probe is not worth processing.

    Parameters
    ----------
    recording_folder : Path
        The sorted-recording asset directory (``data_root / sorted_recording``).
    ephys_collection : str
        The ephys-collection token (the ALF output folder name / probe stream).

    Returns
    -------
    bool
        *True* if at least one non-LFP postprocessed analyzer matches the
        collection.
    """
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

    recording_folder = data_root / sorted_rec

    if row.surface_finding is not None:
        extract_continuous(
            recording_folder,
            results_folder,
            probe_surface_finding=data_root / str(row.surface_finding),
            num_parallel_jobs=num_parallel_jobs,
        )
    else:
        extract_continuous(recording_folder, results_folder, num_parallel_jobs=num_parallel_jobs)

    extract_spikes(recording_folder, results_folder)
