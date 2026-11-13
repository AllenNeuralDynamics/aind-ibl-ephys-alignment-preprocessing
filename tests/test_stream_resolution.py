"""Probe-and-shank manifest rows resolve onto the per-shank streams a quad-base probe records."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.coverage import COVERAGE_FILENAME, CoverageError, RunCoverage
from aind_ibl_ephys_alignment_preprocessing.stream_resolution import (
    SortedStream,
    StreamResolutionError,
    default_locator,
    read_resolved_manifest,
    resolve_manifest_streams,
    resolve_split_stream,
    sorted_streams,
)

RECORDING = "ecephys_863973_2026-08-26_14-55-23"
SORTED = f"{RECORDING}_sorted-opto_2026-10-10_04-21-47"
SERIAL_A, SERIAL_B = "23328424092", "23328408402"
HEADER = "mouseid,sorted_recording,probe_file,probe_name,probe_shank,probe_id,surface_finding,annotation_format"


def _analyzer(postprocessed: Path, collection: str, probes: list[dict]) -> None:
    info = (
        postprocessed / f"experiment1_Record Node 103#Neuropix-PXI-111.{collection}_recording1.zarr" / "recording_info"
    )
    info.mkdir(parents=True)
    (info / ".zattrs").write_text(json.dumps({"probegroup": {"probe_ids": [], "probes": probes}}))


def _probe(serial: str, shank: int | list[int], *, contact_ids: bool = True) -> dict:
    shanks = [shank] if isinstance(shank, int) else shank
    per_contact = [s for s in shanks for _ in range(4)]
    return {
        "annotations": {"serial_number": serial, "part_number": "NP2021"},
        "contact_ids": [f"s{s}e{i}" if contact_ids else str(i) for i, s in enumerate(per_contact)],
        "shank_ids": [str(s) for s in per_contact],
    }


def _mount_quad_sort(data_root: Path, probes: dict[str, str]) -> Path:
    """Mount a sorted asset recording each named probe as four per-shank streams."""
    root = data_root / SORTED
    (root / "spikesorted").mkdir(parents=True)
    (root / "data_description.json").write_text(json.dumps({"schema_version": "2.3.0", "source_data": [RECORDING]}))
    for probe, serial in probes.items():
        for position in range(1, 5):
            _analyzer(root / "postprocessed", f"{probe}-{position}", [_probe(serial, position - 1)])
    return root


def _manifest(rows: list[tuple[str, int | str]]) -> pd.DataFrame:
    lines = [HEADER] + [
        f"863973,{SORTED},ng,{probe},{shank},20260826_{probe}_Shank{shank if shank == '' else int(shank) + 1},,json"
        for probe, shank in rows
    ]
    from io import StringIO

    return pd.read_csv(StringIO("\n".join(lines)))


def _streams(*specs: tuple[str, str | None, set[int]]) -> list[SortedStream]:
    return [SortedStream(c, serial, frozenset(shanks)) for c, serial, shanks in specs]


QUAD_A = _streams(*[(f"ProbeA-{n}", SERIAL_A, {n - 1}) for n in range(1, 5)])


def test_sorted_streams_read_the_recorded_shank_of_each_stream(tmp_path: Path) -> None:
    root = _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    assert sorted_streams(root) == QUAD_A


def test_shank_ids_stand_in_when_contact_ids_carry_no_shank(tmp_path: Path) -> None:
    _analyzer(tmp_path / "postprocessed", "ProbeA-3", [_probe(SERIAL_A, 2, contact_ids=False)])
    assert sorted_streams(tmp_path) == _streams(("ProbeA-3", SERIAL_A, {2}))


@pytest.mark.parametrize("shank", [0, 1, 2, 3])
def test_probe_and_shank_resolve_to_the_stream_recording_that_shank(shank: int) -> None:
    resolution = resolve_split_stream(QUAD_A, "ProbeA", shank)
    assert resolution is not None
    assert (resolution.probe, resolution.shank, resolution.collection, resolution.serial) == (
        "ProbeA",
        shank,
        f"ProbeA-{shank + 1}",
        SERIAL_A,
    )


@pytest.mark.parametrize("collection", ["ProbeA-2", "ProbeA-2-AP", "ProbeZ"])
def test_rows_naming_a_stream_or_an_unknown_probe_are_left_alone(collection: str) -> None:
    assert resolve_split_stream(QUAD_A, collection, 1) is None


@pytest.mark.parametrize(
    ("streams", "shank", "message"),
    [
        (QUAD_A, None, "needs probe_shank"),
        (QUAD_A[:2] + _streams(("ProbeA-3", SERIAL_B, {2})), 0, "different probes"),
        (_streams(("ProbeA-1", SERIAL_A, {0, 1}), ("ProbeA-2", SERIAL_A, {1})), 0, "expected one shank"),
        (_streams(("ProbeA-1", SERIAL_A, set()), ("ProbeA-2", SERIAL_A, {1})), 1, "no readable shank"),
        (_streams(("ProbeA-1", SERIAL_A, {1}), ("ProbeA-2", SERIAL_A, {0})), 0, "contradicting its name"),
        (QUAD_A[:2], 3, "no stream of ProbeA records shank 3"),
    ],
    ids=["no-shank", "mixed-serials", "multi-shank-stream", "unreadable-geometry", "suffix-contradiction", "absent"],
)
def test_ambiguous_split_streams_are_refused(streams: list[SortedStream], shank: int | None, message: str) -> None:
    with pytest.raises(StreamResolutionError, match=message):
        resolve_split_stream(streams, "ProbeA", shank)


def test_manifest_rows_take_the_split_stream_shape(tmp_path: Path) -> None:
    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A, "ProbeB": SERIAL_B})
    resolved = resolve_manifest_streams(_manifest([("ProbeA", 0), ("ProbeB", 3)]), default_locator(tmp_path))
    assert resolved.errors == {}
    columns = ["probe_name", "ephys_collection", "logical_probe", "histology_shank", "ephys_shank"]
    assert resolved.df[columns].to_dict("records") == [
        {"probe_name": "ProbeA", "ephys_collection": "ProbeA-1", "logical_probe": "ProbeA", "histology_shank": 0,
         "ephys_shank": 0},
        {"probe_name": "ProbeB", "ephys_collection": "ProbeB-4", "logical_probe": "ProbeB", "histology_shank": 3,
         "ephys_shank": 0},
    ]  # fmt: skip
    assert [(r.ephys_collection, r.histology_shank, r.ephys_shank) for r in resolved.rows] == [
        ("ProbeA-1", 0, 0),
        ("ProbeB-4", 3, 0),
    ]


def test_resolving_a_resolved_manifest_changes_nothing(tmp_path: Path) -> None:
    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    once = resolve_manifest_streams(_manifest([("ProbeA", 1), ("ProbeA", 2)]), default_locator(tmp_path))
    csv = tmp_path / "resolved.csv"
    once.df.to_csv(csv, index=False)
    twice = resolve_manifest_streams(pd.read_csv(csv), default_locator(tmp_path))
    assert twice.resolutions == {} and twice.errors == {}
    assert [r.ephys_collection for r in twice.rows] == ["ProbeA-2", "ProbeA-3"]


def test_rows_without_a_mounted_sort_pass_through(tmp_path: Path) -> None:
    resolved = resolve_manifest_streams(_manifest([("ProbeA", 1)]), default_locator(tmp_path))
    assert resolved.resolutions == {} and resolved.errors == {}
    assert resolved.rows[0].ephys_collection == "ProbeA"


def test_standalone_runs_refuse_an_unresolvable_manifest(tmp_path: Path) -> None:
    _mount_quad_sort(tmp_path / "data", {"ProbeA": SERIAL_A})
    _manifest([("ProbeA", 1), ("ProbeA", "")]).to_csv(tmp_path / "manifest.csv", index=False)
    with pytest.raises(StreamResolutionError, match=r"row 2 \(20260826_ProbeA_Shank\): .*needs probe_shank"):
        read_resolved_manifest(tmp_path / "manifest.csv", tmp_path / "data")


def _discover_config(tmp_path: Path) -> SimpleNamespace:
    return SimpleNamespace(
        manifest_csv=tmp_path / "manifest.csv",
        results_root=tmp_path / "results",
        data_root=tmp_path / "data",
        skip_ephys=False,
        allow_partial=False,
    )


@pytest.fixture
def tracks_present(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(stages, "_track_annotation_present", lambda data_root, mr: (True, None))


def test_discover_hands_every_stage_the_resolved_stream(tmp_path: Path, tracks_present: None) -> None:
    config = _discover_config(tmp_path)
    _mount_quad_sort(config.data_root, {"ProbeA": SERIAL_A})
    _manifest([("ProbeA", shank) for shank in range(4)]).to_csv(config.manifest_csv, index=False)

    written = stages.stage_discover(config)

    filtered = pd.read_csv(config.results_root / "manifest.csv")
    assert list(filtered["ephys_collection"]) == [f"ProbeA-{n}" for n in range(1, 5)]
    assert list(filtered["histology_shank"]) == [0, 1, 2, 3]
    assert set(filtered["ephys_shank"]) == {0}
    assert sorted(json.loads(p.read_text())["ephys_collection"] for p in written) == [
        f"ProbeA-{n}" for n in range(1, 5)
    ]
    coverage = RunCoverage.read(config.results_root / COVERAGE_FILENAME)
    assert coverage.rows_viable == 4
    assert "ProbeA shank 2 resolved to per-shank stream ProbeA-3, serial 23328424092" in coverage.rows[2].notes


def test_discover_refuses_a_row_its_streams_cannot_resolve(tmp_path: Path, tracks_present: None) -> None:
    config = _discover_config(tmp_path)
    root = _mount_quad_sort(config.data_root, {"ProbeA": SERIAL_A})
    # The sort lost the third shank's stream.
    for analyzer in (root / "postprocessed").glob("*ProbeA-3_*"):
        (analyzer / "recording_info" / ".zattrs").unlink()
        (analyzer / "recording_info").rmdir()
        analyzer.rmdir()
    _manifest([("ProbeA", 1), ("ProbeA", 2)]).to_csv(config.manifest_csv, index=False)

    with pytest.raises(CoverageError):
        stages.stage_discover(config)

    coverage = RunCoverage.read(config.results_root / COVERAGE_FILENAME)
    assert [row.histology_track_id for row in coverage.dropped] == ["20260826_ProbeA_Shank3"]
    assert "no stream of ProbeA records shank 2" in str(coverage.dropped[0].reason)
