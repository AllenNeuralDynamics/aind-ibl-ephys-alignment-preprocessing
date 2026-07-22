"""Tests for the fan-out pipeline stage entry points.

Only :func:`stage_discover` and the unit-naming helper are exercised here --
they depend on nothing heavier than pandas + ``role_dispatch``. The
histology/ephys/pack stages need ants/iblatlas/SimpleITK and real assets, so
they are covered by integration runs rather than these smoke tests.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from aind_code_ocean_pipeline_utils.role_dispatch import find_stream_config

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.stages import (
    EPHYS_STREAM_MARKER,
    _ephys_unit_name,
    stage_discover,
)

_MANIFEST_HEADER = "mouseid,sorted_recording,probe_file,probe_id,probe_name,surface_finding"


def _write_manifest(path: Path, rows: list[str]) -> None:
    """Write a minimal manifest CSV with the given data rows."""
    path.write_text("\n".join([_MANIFEST_HEADER, *rows]) + "\n")


def test_public_stage_functions_exist() -> None:
    """The four pipeline stages are importable from the module."""
    for name in ("stage_discover", "stage_histology", "stage_ephys", "stage_pack"):
        assert callable(getattr(stages, name))


def test_ephys_unit_name_uses_collection_and_falls_back() -> None:
    """The unit name joins recording + collection, with an 'all' fallback."""
    assert _ephys_unit_name("rec_2025-10-08", "ProbeA") == "rec_2025-10-08__ProbeA"
    assert _ephys_unit_name("rec_2025-10-08", None) == "rec_2025-10-08__all"


def test_stage_discover_writes_one_config_per_unique_unit(tmp_path: Path) -> None:
    """One ephys config per unique ``(recording, collection)``; duplicates collapse."""
    manifest = tmp_path / "manifest.csv"
    results = tmp_path / "results"
    # Two collections of one recording + one collection of another + a
    # duplicate of the first row -> three unique fan-out units.
    _write_manifest(
        manifest,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track2,ProbeB,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track3,ProbeA,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
        ],
    )
    config = SimpleNamespace(manifest_csv=manifest, results_root=results)

    written = stage_discover(config)

    assert len(written) == 3
    assert all(p.name == "config.json" for p in written)


def test_stage_discover_config_payload_is_complete(tmp_path: Path) -> None:
    """Each written config carries the marker and every field a worker needs."""
    manifest = tmp_path / "manifest.csv"
    results = tmp_path / "results"
    _write_manifest(
        manifest,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,surface_791094",
        ],
    )
    config = SimpleNamespace(manifest_csv=manifest, results_root=results)

    written = stage_discover(config)
    payload = json.loads(written[0].read_text())

    assert payload[EPHYS_STREAM_MARKER]
    assert payload["mouseid"] == "791094"
    assert payload["sorted_recording"] == "ecephys_791094_2025-10-08_sorted_x"
    assert payload["recording_id"] == "ecephys_791094_2025-10-08"
    assert payload["ephys_collection"] == "ProbeA"
    assert payload["surface_finding"] == "surface_791094"


def test_stage_discover_output_is_worker_discoverable(tmp_path: Path) -> None:
    """A worker's find_stream_config locates exactly one config it wrote."""
    manifest = tmp_path / "manifest.csv"
    results = tmp_path / "results"
    _write_manifest(
        manifest,
        [
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track3,ProbeA,",
        ],
    )
    config = SimpleNamespace(manifest_csv=manifest, results_root=results)

    stage_discover(config)

    # The pipeline stages each config subdir under /data; here results IS the
    # stage-in root. A worker must find exactly its one config by the marker.
    _, cfg = find_stream_config(results, schema_marker=EPHYS_STREAM_MARKER)
    assert cfg["recording_id"] == "ecephys_791094_2025-10-09"
    assert cfg["ephys_collection"] == "ProbeA"


def test_stage_discover_handles_missing_surface_finding(tmp_path: Path) -> None:
    """A blank surface_finding column serializes as null, not the string 'nan'."""
    manifest = tmp_path / "manifest.csv"
    results = tmp_path / "results"
    _write_manifest(
        manifest,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
        ],
    )
    config = SimpleNamespace(manifest_csv=manifest, results_root=results)

    written = stage_discover(config)
    payload = json.loads(written[0].read_text())

    assert payload["surface_finding"] is None
