"""Ephys extraction wrappers (synchronous)."""

from __future__ import annotations

import logging
from pathlib import Path

from aind_ephys_ibl_gui_conversion.ephys import extract_continuous, extract_spikes

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, OutputDirs

logger = logging.getLogger(__name__)


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
