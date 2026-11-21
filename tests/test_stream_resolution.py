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
        (QUAD_A, None, "provide histology_shank"),
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


def test_standalone_structure_validation_uses_the_same_resolved_shanks(tmp_path: Path) -> None:
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A, "ProbeB": SERIAL_B})
    manifest = tmp_path / "manifest.csv"
    _ephys_only_manifest().to_csv(manifest, index=False)
    config = PipelineConfig(
        data_root=tmp_path,
        results_root=tmp_path / "results",
        neuroglancer_file=Path("ng.json"),
        manifest_csv=manifest,
    )
    validator = PipelineValidator(config, skip_resource_checks=True)
    validator.validate_manifest_structure()
    errors = [result for result in validator.results if not result.passed and result.severity == "error"]
    assert not errors, [result.message for result in errors]


@pytest.mark.parametrize("layout", ["shared", "split", "mixed"])
@pytest.mark.parametrize("optional_columns", [False, True])
def test_minimal_manifest_survives_discovery_and_pack(tmp_path: Path, layout: str, optional_columns: bool) -> None:
    """A physical-shank-only CSV must work alone as well as beside split streams."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig, ProcessResult
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    config = _discover_config(tmp_path)
    root = _mount_quad_sort(config.data_root, {"ProbeB": SERIAL_B})
    _analyzer(root / "postprocessed", "ProbeA", [_probe(SERIAL_A, list(range(4)))])
    probes = {"shared": ["ProbeA"], "split": ["ProbeB"], "mixed": ["ProbeA", "ProbeB"]}[layout]
    frame = _manifest([(probe, shank) for probe in probes for shank in range(4)]).rename(
        columns={"probe_name": "ephys_collection", "probe_shank": "histology_shank", "probe_id": "histology_track_id"}
    )
    frame = frame.drop(columns=["surface_finding", "annotation_format"])
    if optional_columns:
        for name in ("logical_probe", "ephys_shank", "annotation_format"):
            frame[name] = ""
    annotations = config.data_root / "annotations" / "ng.json"
    annotations.parent.mkdir()
    annotations.write_text(
        json.dumps(
            {
                "dimensions": {"z": [1e-6, "m"], "y": [1e-6, "m"], "x": [1e-6, "m"]},
                "layers": [
                    {
                        "type": "annotation",
                        "name": name,
                        "annotations": [
                            {"type": "point", "id": "a", "point": [1, 2, 3, 0]},
                            {"type": "point", "id": "b", "point": [4, 5, 6, 0]},
                        ],
                    }
                    for name in frame["histology_track_id"]
                ],
            }
        )
    )
    frame.to_csv(config.manifest_csv, index=False)
    validator = PipelineValidator(
        PipelineConfig(
            data_root=config.data_root,
            results_root=config.results_root,
            neuroglancer_file="annotations/ng.json",
            manifest_csv=config.manifest_csv,
        ),
        skip_resource_checks=True,
    )
    validator.validate_manifest_structure()
    errors = [r.message for r in validator.results if not r.passed and r.severity == "error"]
    assert errors == []

    written = stages.stage_discover(config)
    filtered = read_resolved_manifest(config.results_root / "manifest.csv", config.data_root)
    rows = [stages.ManifestRow.from_series(row) for _, row in filtered.iterrows()]
    assert all(row.annotation_format == "json" for row in rows)
    assert len(written) == (1 if "ProbeA" in probes else 0) + (4 if "ProbeB" in probes else 0)
    expected = [
        (probe if probe == "ProbeA" else f"ProbeB-{shank + 1}", shank, shank if probe == "ProbeA" else 0)
        for probe in probes
        for shank in range(4)
    ]
    assert [(r.ephys_collection, r.histology_shank, r.ephys_shank) for r in rows] == expected

    # Exercise the packer's actual row-to-pick mapping with the saved discovery CSV.
    mouse_root = config.results_root / "863973"
    for row in rows:
        folder = mouse_root / row.recording_id / str(row.ephys_collection)
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "channels.localCoordinates.npy").touch()
    results = [ProcessResult(probe_id=r.probe_id, recording_id=r.recording_id, wrote_files=True) for r in rows]
    entries = _build_probes(
        rows,
        results,
        SimpleNamespace(tracks_root=mouse_root / "tracks"),
        mouse_root,
        SimpleNamespace(skip_ephys=False, emit_qc=True),
    )[RECORDING]
    picks = [
        (name, pick.histology_shank, pick.ephys_shank) for name, entry in entries.items() for pick in entry.xyz_picks
    ]
    assert sorted(picks) == sorted(expected)


def test_explicit_local_shank_collisions_are_still_rejected(tmp_path: Path) -> None:
    """Defaulting omitted values must not overwrite supplied local shank mappings."""
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    frame = _manifest([("ProbeA", 0), ("ProbeA", 1)]).rename(columns={"probe_shank": "histology_shank"})
    frame["ephys_shank"] = 0
    manifest = tmp_path / "manifest.csv"
    frame.to_csv(manifest, index=False)
    validator = PipelineValidator(
        PipelineConfig(
            data_root=tmp_path, results_root=tmp_path / "out", neuroglancer_file="ng.json", manifest_csv=manifest
        ),
        skip_resource_checks=True,
    )
    validator.validate_manifest_structure()
    assert any(not r.passed and "ephys shanks" in r.message for r in validator.results)


def test_rows_without_a_mounted_sort_pass_through(tmp_path: Path) -> None:
    resolved = resolve_manifest_streams(_manifest([("ProbeA", 1)]), default_locator(tmp_path))
    assert resolved.resolutions == {} and resolved.errors == {}
    assert resolved.rows[0].ephys_collection == "ProbeA"


def test_standalone_runs_refuse_an_unresolvable_manifest(tmp_path: Path) -> None:
    _mount_quad_sort(tmp_path / "data", {"ProbeA": SERIAL_A})
    _manifest([("ProbeA", 1), ("ProbeA", "")]).to_csv(tmp_path / "manifest.csv", index=False)
    with pytest.raises(StreamResolutionError, match=r"row 2 \(20260826_ProbeA_Shank\): .*provide histology_shank"):
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


def _ephys_only_manifest() -> pd.DataFrame:
    """Use Anna's logical collection plus ephys_shank-only CSV layout."""
    frame = _manifest([(probe, shank) for probe in ("ProbeA", "ProbeB") for shank in range(4)])
    return frame.rename(
        columns={"probe_name": "ephys_collection", "probe_shank": "ephys_shank", "probe_id": "histology_track_id"}
    )


def test_discover_normalizes_ephys_only_shanks_and_records_interpretation(tmp_path: Path, tracks_present: None) -> None:
    config = _discover_config(tmp_path)
    _mount_quad_sort(config.data_root, {"ProbeA": SERIAL_A, "ProbeB": SERIAL_B})
    _ephys_only_manifest().to_csv(config.manifest_csv, index=False)

    written = stages.stage_discover(config)

    filtered = pd.read_csv(config.results_root / "manifest.csv")
    assert len(written) == 8
    assert list(filtered["ephys_collection"]) == [f"{probe}-{n}" for probe in ("ProbeA", "ProbeB") for n in range(1, 5)]
    assert list(filtered["histology_shank"]) == [0, 1, 2, 3] * 2
    assert set(filtered["ephys_shank"]) == {0}
    coverage = RunCoverage.read(config.results_root / COVERAGE_FILENAME)
    assert coverage.rows_viable == 8
    for row in coverage.rows:
        assert any(
            "interpreted input ephys_shank=" in note and "verifying recorded geometry" in note for note in row.notes
        )
    resolved = resolve_manifest_streams(filtered, default_locator(config.data_root))
    assert not resolved.resolutions and not resolved.errors


@pytest.mark.parametrize("collection", ["ProbeA", "ProbeA-2"])
def test_ephys_only_shank_keeps_its_meaning_for_exact_collections(tmp_path: Path, collection: str) -> None:
    root = _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    if collection == "ProbeA":
        _analyzer(root / "postprocessed", "ProbeA", [_probe(SERIAL_A, [0, 1, 2, 3])])
    frame = _ephys_only_manifest().iloc[[1]].copy()
    frame["ephys_collection"] = collection
    resolved = resolve_manifest_streams(frame, default_locator(tmp_path))
    assert not resolved.resolutions and not resolved.errors
    row = resolved.rows[0]
    assert (row.ephys_collection, row.ephys_shank, row.histology_shank) == (collection, 1, None)


def test_explicit_physical_shank_takes_precedence_over_ephys_shank(tmp_path: Path) -> None:
    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    frame = _ephys_only_manifest().iloc[[0]].copy()
    frame["histology_shank"] = 2
    resolved = resolve_manifest_streams(frame, default_locator(tmp_path))
    assert not resolved.errors
    assert resolved.rows[0].ephys_collection == "ProbeA-3"
    assert not resolved.resolutions[0].used_ephys_shank_fallback


@pytest.mark.parametrize(
    "streams",
    [
        _streams(("ProbeA-1", SERIAL_A, set())),
        QUAD_A[:2] + _streams(("ProbeA-3", SERIAL_B, {2})),
        _streams(("ProbeA-1", SERIAL_A, {1})),
        QUAD_A[1:],
    ],
    ids=["unreadable-geometry", "different-probes", "contradictory-suffix", "missing-shank"],
)
def test_ephys_only_fallback_requires_verified_geometry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, streams: list[SortedStream]
) -> None:
    from aind_ibl_ephys_alignment_preprocessing import stream_resolution

    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    monkeypatch.setattr(stream_resolution, "sorted_streams", lambda _: streams)
    resolved = resolve_manifest_streams(_ephys_only_manifest().iloc[[0]], default_locator(tmp_path))
    assert not resolved.resolutions
    assert "cannot interpret ephys_shank as histology_shank" in resolved.errors[0]
    assert resolved.rows[0].ephys_collection == "ProbeA"


@pytest.mark.parametrize("column", ["histology_shank", "probe_shank", "ephys_shank"])
def test_ephys_only_fallback_does_not_hide_invalid_shank_values(tmp_path: Path, column: str) -> None:
    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    frame = _ephys_only_manifest().iloc[[1]].copy()
    frame[column] = "invalid" if column != "ephys_shank" else 1.5
    resolved = resolve_manifest_streams(frame, default_locator(tmp_path))
    assert not resolved.resolutions and 0 in resolved.errors


def test_repeated_analyzers_for_one_collection_still_identify_one_stream(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from aind_ibl_ephys_alignment_preprocessing import stream_resolution

    _mount_quad_sort(tmp_path, {"ProbeA": SERIAL_A})
    monkeypatch.setattr(stream_resolution, "sorted_streams", lambda _: QUAD_A * 2)
    resolved = resolve_manifest_streams(_ephys_only_manifest().iloc[[2]], default_locator(tmp_path))
    assert not resolved.errors
    assert resolved.rows[0].ephys_collection == "ProbeA-3"
    assert resolved.resolutions[0].used_ephys_shank_fallback
