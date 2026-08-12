"""Tests for the run record: manifest coverage plus the assets a run was built on.

Deliberately separate from ``test_datapackage.py``. The datapackage is the
alignment GUI's consumption contract; this is a record of the run that produced
it, and the two are versioned and read by different parties.
"""

from __future__ import annotations

import json
from pathlib import Path

from aind_ibl_ephys_alignment_preprocessing.coverage import (
    COVERAGE_FILENAME,
    RowCoverage,
    RunCoverage,
)

UNIT = "ecephys_791094_2025-10-08__ProbeA"


def _record(**overrides) -> RunCoverage:
    """A two-row record: one viable and emitted, one dropped."""
    fields = dict(
        mouse_id="791094",
        manifest="/data/annotations/manifest.csv",
        rows_requested=2,
        rows_viable=1,
        units_emitted=[UNIT],
        rows=[
            RowCoverage(
                recording_id="ecephys_791094_2025-10-08",
                ephys_collection="ProbeA",
                viable=True,
                unit=UNIT,
            ),
            RowCoverage(
                recording_id="ecephys_791094_2025-10-09",
                ephys_collection="ProbeB",
                viable=False,
                reason="no spike-sorting output (bad sorting)",
            ),
        ],
    )
    fields.update(overrides)
    return RunCoverage(**fields)


def _resolution(**overrides) -> str:
    """Exactly the shape the trigger's ``_asset_resolution_arg`` emits."""
    payload = {
        "units": {
            UNIT: {
                "raw": {"id": "raw-id", "name": "ecephys_791094_2025-10-08"},
                "sorted": {"id": "sorted-id", "name": "ecephys_791094_2025-10-08_sorted_x"},
            }
        },
        "smartspim": {"id": "spim-id", "name": "SmartSPIM_791094_stitched"},
    }
    payload.update(overrides)
    return json.dumps(payload)


def test_round_trips_through_the_file(tmp_path: Path):
    """The record survives the write/read that carries it between stages."""
    path = _record().write(tmp_path)

    assert path.name == COVERAGE_FILENAME
    assert RunCoverage.read(path) == _record()


def test_find_locates_a_mounted_discover_asset(tmp_path: Path):
    """pack looks one level down too -- that is where a captured asset lands."""
    _record().write(tmp_path / "discover")

    assert RunCoverage.find(tmp_path / "discover") is not None
    assert RunCoverage.find(tmp_path) is not None  # via the */coverage.json glob


def test_find_returns_none_when_there_is_no_record(tmp_path: Path):
    """The monolith has no discover stage; absence is not an error."""
    assert RunCoverage.find(tmp_path) is None


def test_find_ignores_an_unreadable_record(tmp_path: Path):
    """A corrupt record must not fail a run that already did the work."""
    (tmp_path / COVERAGE_FILENAME).write_text("{not json")

    assert RunCoverage.find(tmp_path) is None


def test_with_assets_merges_the_resolved_ids():
    """The launcher's half joins the producer's half on the unit name."""
    coverage = _record().with_assets(_resolution())

    assert set(coverage.unit_assets) == set(coverage.units_emitted)
    assert coverage.unit_assets[UNIT].sorted.id == "sorted-id"
    assert coverage.unit_assets[UNIT].raw.id == "raw-id"
    assert coverage.smartspim.id == "spim-id"
    # ...and the producer's half is untouched by the merge.
    assert coverage.rows_requested == 2
    assert coverage.is_complete is False


def test_with_assets_ignores_a_malformed_resolution():
    """Provenance is worth recording, never worth failing a completed run over."""
    for bad in ("", "   ", "not json", "[1, 2, 3]", None):
        coverage = _record().with_assets(bad)
        assert coverage.unit_assets == {}
        assert coverage.rows_requested == 2


def test_with_assets_warns_when_the_join_key_does_not_match(caplog):
    """A resolution keyed by the ``stream_`` directory would answer nothing.

    Three layers agree on one key -- discover's ``units_emitted``, the trigger's
    resolved-asset map, and this record. A silent mismatch produces a file that
    validates and joins to nothing, so it is called out.
    """
    with caplog.at_level("WARNING"):
        coverage = _record().with_assets(json.dumps({"units": {f"stream_{UNIT}": {"raw": {"id": "raw-id"}}}}))

    assert f"stream_{UNIT}" in caplog.text
    assert "never emitted" in caplog.text
    assert coverage.unit_assets  # still recorded, so the mismatch is inspectable


def test_complete_run_reports_complete():
    """The happy path: every requested row survived."""
    assert _record(rows_requested=1, rows_viable=1, rows=[]).is_complete is True


def test_row_key_is_shank_aware():
    """Rows are keyed per shank, because units are not."""
    row = RowCoverage(recording_id="rec", ephys_collection="ProbeA", ephys_shank=1, viable=False)
    assert row.key == "rec/ProbeA/shank1"

    unshanked = RowCoverage(recording_id="rec", ephys_collection="ProbeA", viable=False)
    assert unshanked.key == "rec/ProbeA/shank-"


def test_gap_message_names_the_track_and_the_reason():
    """A gap line has to be actionable on its own, in a log or an exception."""
    row = RowCoverage(
        recording_id="rec",
        ephys_collection="ProbeB",
        histology_track_id="Track2",
        viable=False,
        reason="no spike-sorting output (bad sorting)",
    )

    message = row.gap_message()
    assert "Track2" in message
    assert "no spike-sorting output" in message
