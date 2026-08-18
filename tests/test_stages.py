"""Tests for the fan-out pipeline stage entry points.

Only ``stage_discover`` and its viability/naming helpers are exercised here --
they depend on nothing heavier than pandas + ``role_dispatch`` (the NG and
sorting checks are monkeypatched). The histology/ephys/pack stages need
ants/iblatlas/SimpleITK and real assets, so they are covered by integration
runs rather than these smoke tests.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from aind_code_ocean_pipeline_utils.role_dispatch import find_stream_config

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.coverage import (
    COVERAGE_FILENAME,
    CoverageError,
    RunCoverage,
)
from aind_ibl_ephys_alignment_preprocessing.datapackage import (
    _find_mouse_output_trees,
    merge_pipeline_outputs,
)
from aind_ibl_ephys_alignment_preprocessing.ephys import find_sorted_session_dirs
from aind_ibl_ephys_alignment_preprocessing.stages import (
    EPHYS_STREAM_MARKER,
    _ephys_unit_name,
    _track_annotation_present,
    stage_discover,
)
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow

_MANIFEST_HEADER = "mouseid,sorted_recording,probe_file,probe_id,probe_name,surface_finding"


def _write_manifest(path: Path, rows: list[str]) -> None:
    """Write a minimal manifest CSV with the given data rows."""
    path.write_text("\n".join([_MANIFEST_HEADER, *rows]) + "\n")


def _config(tmp_path: Path, *, skip_ephys: bool = False, allow_partial: bool = False) -> SimpleNamespace:
    """A stand-in PipelineConfig with the attributes stage_discover reads."""
    return SimpleNamespace(
        manifest_csv=tmp_path / "manifest.csv",
        results_root=tmp_path / "results",
        data_root=tmp_path / "data",
        skip_ephys=skip_ephys,
        allow_partial=allow_partial,
    )


@pytest.fixture
def all_viable(monkeypatch: pytest.MonkeyPatch) -> None:
    """Force every probe viable so emission logic can be tested in isolation."""
    monkeypatch.setattr(stages, "_probe_viability", lambda config, mr, sorted_dirs=None: (True, None, []))


def test_public_stage_functions_exist() -> None:
    """Every pipeline stage entry point is importable from the module."""
    for name in (
        "stage_discover",
        "stage_histology",
        "stage_ephys_launch",
        "stage_ephys",
        "stage_ephys_collect",
        "stage_pack",
    ):
        assert callable(getattr(stages, name))


def test_ephys_unit_name_uses_collection_and_falls_back() -> None:
    """The unit name joins recording + collection, with an 'all' fallback."""
    assert _ephys_unit_name("rec_2025-10-08", "ProbeA") == "rec_2025-10-08__ProbeA"
    assert _ephys_unit_name("rec_2025-10-08", None) == "rec_2025-10-08__all"


def test_stage_discover_writes_one_config_per_unique_unit(tmp_path: Path, all_viable: None) -> None:
    """One ephys config per unique ``(recording, collection)``; duplicates collapse."""
    config = _config(tmp_path)
    # Two collections of one recording + one collection of another + a duplicate
    # of the first row -> three unique fan-out units.
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track2,ProbeB,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track3,ProbeA,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
        ],
    )

    written = stage_discover(config)

    assert len(written) == 3
    assert all(p.name == "config.json" for p in written)


def test_stage_discover_writes_filtered_manifest(tmp_path: Path, all_viable: None) -> None:
    """A filtered ``manifest.csv`` for histology/pack is written to results."""
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track3,ProbeB,",
        ],
    )

    stage_discover(config)

    out_manifest = config.results_root / "manifest.csv"
    assert out_manifest.is_file()
    assert out_manifest.read_text().strip().count("\n") == 2  # header + 2 rows


def test_stage_discover_config_payload_is_complete(tmp_path: Path, all_viable: None) -> None:
    """Each written config carries the marker and every field a worker needs."""
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,surface_791094"],
    )

    written = stage_discover(config)
    payload = json.loads(written[0].read_text())

    assert payload[EPHYS_STREAM_MARKER]
    assert payload["mouseid"] == "791094"
    assert payload["sorted_recording"] == "ecephys_791094_2025-10-08_sorted_x"
    assert payload["recording_id"] == "ecephys_791094_2025-10-08"
    assert payload["ephys_collection"] == "ProbeA"
    assert payload["surface_finding"] == "surface_791094"


def test_stage_discover_output_is_worker_discoverable(tmp_path: Path, all_viable: None) -> None:
    """A worker's find_stream_config locates exactly one config it wrote."""
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-09_sorted_y,ng,Track3,ProbeA,"],
    )

    stage_discover(config)

    _, cfg = find_stream_config(config.results_root, schema_marker=EPHYS_STREAM_MARKER)
    assert cfg["recording_id"] == "ecephys_791094_2025-10-09"
    assert cfg["ephys_collection"] == "ProbeA"


def test_stage_discover_handles_missing_surface_finding(tmp_path: Path, all_viable: None) -> None:
    """A blank surface_finding column serializes as null, not the string 'nan'."""
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,"],
    )

    written = stage_discover(config)
    payload = json.loads(written[0].read_text())

    assert payload["surface_finding"] is None


def test_stage_discover_filters_nonviable_probes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With partial coverage allowed, non-viable probes drop from manifest and configs."""
    config = _config(tmp_path, allow_partial=True)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track2,ProbeB,",
        ],
    )
    # ProbeB is not viable (e.g. failed sorting or missing track).
    monkeypatch.setattr(
        stages,
        "_probe_viability",
        lambda config, mr, sorted_dirs=None: (mr.ephys_collection != "ProbeB", "skipped for test", []),
    )

    written = stage_discover(config)

    # Only ProbeA's unit is emitted.
    assert len(written) == 1
    assert json.loads(written[0].read_text())["ephys_collection"] == "ProbeA"
    # ...and the filtered manifest keeps only ProbeA.
    filtered = (config.results_root / "manifest.csv").read_text()
    assert "ProbeA" in filtered
    assert "ProbeB" not in filtered


def test_stage_discover_skip_ephys_emits_no_configs(tmp_path: Path, all_viable: None) -> None:
    """With ephys disabled, discover writes the manifest but no fan-out configs."""
    config = _config(tmp_path, skip_ephys=True)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,"],
    )

    written = stage_discover(config)

    assert written == []
    assert (config.results_root / "manifest.csv").is_file()


# ---------------------------------------------------------------------------
# Coverage: what the manifest asked for vs what discover emitted
# ---------------------------------------------------------------------------

_SHANK_HEADER = "mouseid,sorted_recording,probe_file,probe_id,probe_name,ephys_shank,surface_finding"


def _write_shank_manifest(path: Path, rows: list[str]) -> None:
    """Write a manifest whose rows carry an explicit ``ephys_shank``."""
    path.write_text("\n".join([_SHANK_HEADER, *rows]) + "\n")


def _read_coverage(config: SimpleNamespace) -> RunCoverage:
    """Load the coverage record discover wrote."""
    return RunCoverage.read(config.results_root / COVERAGE_FILENAME)


def _mount_sorted_for_coverage(root: Path, *, input_recording: str, collection: str) -> None:
    """Mount a spike-sorted asset that names *input_recording* and sorted *collection*."""
    (root / "spikesorted").mkdir(parents=True, exist_ok=True)
    (root / "data_description.json").write_text(json.dumps({"input_data_name": input_recording}))
    (root / "postprocessed" / f"experiment1_Neuropix.{collection}-AP_recording1").mkdir(parents=True, exist_ok=True)


@pytest.fixture
def tracks_present(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pass the histology half of the gate so the sorting half can be tested for real."""
    monkeypatch.setattr(stages, "_track_annotation_present", lambda data_root, mr: (True, None))


def test_stage_discover_refuses_when_a_pinned_session_has_no_sorting(tmp_path: Path, tracks_present: None) -> None:
    """The 791094 regression: two pinned sessions, one resolvable, and the run stops.

    Before the coverage record this passed silently -- the unresolvable session's
    rows were dropped, the filtered manifest and the fan-out set shrank together,
    and every downstream check compared the two shrunken sets and agreed.
    """
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track2,ProbeB,",
        ],
    )
    # Only 2025-10-08's sorting is mounted; 2025-10-09's is absent entirely.
    _mount_sorted_for_coverage(
        config.data_root / "sorted_x", input_recording="ecephys_791094_2025-10-08", collection="ProbeA"
    )

    with pytest.raises(CoverageError) as excinfo:
        stage_discover(config)

    assert "ecephys_791094_2025-10-09" in str(excinfo.value)
    # The record that explains the refusal must survive it.
    coverage = _read_coverage(config)
    assert coverage.rows_requested == 2
    assert coverage.rows_viable == 1
    assert [row.recording_id for row in coverage.dropped] == ["ecephys_791094_2025-10-09"]
    assert "no spike-sorting output" in str(coverage.dropped[0].reason)


def test_stage_discover_coverage_records_every_row_and_its_unit(tmp_path: Path, all_viable: None) -> None:
    """Every manifest row appears in the record, naming the unit it fed."""
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track2,ProbeB,",
        ],
    )

    stage_discover(config)
    coverage = _read_coverage(config)

    assert coverage.mouse_id == "791094"
    assert coverage.rows_requested == coverage.rows_viable == 2
    assert coverage.dropped == []
    assert coverage.units_emitted == [
        "ecephys_791094_2025-10-08__ProbeA",
        "ecephys_791094_2025-10-09__ProbeB",
    ]
    assert [row.unit for row in coverage.rows] == coverage.units_emitted


def test_stage_discover_coverage_sees_a_dropped_shank(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A lost shank is a coverage gap even though its unit is still emitted.

    Fan-out units are ``(recording, collection)`` *deduplicated*, so a sibling
    shank keeps the unit alive and a unit-level diff reports full coverage. Only
    a row-keyed record shows the loss.
    """
    config = _config(tmp_path, allow_partial=True)
    _write_shank_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,0,",
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track2,ProbeA,1,",
        ],
    )
    monkeypatch.setattr(
        stages,
        "_probe_viability",
        lambda config, mr, sorted_dirs=None: (mr.ephys_shank != 1, "no track points for layer", []),
    )

    written = stage_discover(config)
    coverage = _read_coverage(config)

    # The unit survives on shank 0 alone -- which is exactly why unit counts lie.
    assert len(written) == 1
    assert coverage.units_emitted == ["ecephys_791094_2025-10-08__ProbeA"]
    # ...and the record still names the shank that was lost.
    assert coverage.rows_viable == 1
    assert [row.key for row in coverage.dropped] == ["ecephys_791094_2025-10-08/ProbeA/shank1"]


def test_stage_discover_allows_a_partial_run_when_asked(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``allow_partial`` downgrades the refusal to a warning, but still records it."""
    config = _config(tmp_path, allow_partial=True)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track2,ProbeB,",
        ],
    )
    monkeypatch.setattr(
        stages,
        "_probe_viability",
        lambda config, mr, sorted_dirs=None: (mr.ephys_collection == "ProbeA", "skipped for test", []),
    )

    written = stage_discover(config)

    assert len(written) == 1
    assert len(_read_coverage(config).dropped) == 1


def test_stage_discover_refuses_zero_viable_even_with_allow_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No viable probe is fatal regardless: a partial run never means "none"."""
    config = _config(tmp_path, allow_partial=True)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,"],
    )
    monkeypatch.setattr(
        stages,
        "_probe_viability",
        lambda config, mr, sorted_dirs=None: (False, "skipped for test", []),
    )

    with pytest.raises(CoverageError, match="no viable probes"):
        stage_discover(config)


def test_stage_discover_refuses_an_empty_manifest(tmp_path: Path, all_viable: None) -> None:
    """A header-only manifest is the same silent success, with no rows to blame.

    It has nothing to drop, so every gap-based check passes it. The first version
    of the gate skipped the zero-viable branch when ``rows_requested`` was 0 and
    let a manifest describing no work run green over nothing.
    """
    config = _config(tmp_path, allow_partial=True)
    _write_manifest(config.manifest_csv, [])

    with pytest.raises(CoverageError, match="lists no manifest rows"):
        stage_discover(config)

    # The record still lands, so the refusal is inspectable rather than log-only.
    assert _read_coverage(config).rows_requested == 0


def test_stage_discover_writes_coverage_when_ephys_is_skipped(tmp_path: Path, all_viable: None) -> None:
    """Histology-only runs get a record too -- they have no fan-out to notice a gap."""
    config = _config(tmp_path, skip_ephys=True)
    _write_manifest(
        config.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,"],
    )

    stage_discover(config)
    coverage = _read_coverage(config)

    assert coverage.skip_ephys is True
    assert coverage.units_emitted == []
    assert coverage.rows_viable == 1


def test_probe_viability_notes_a_name_matched_sorting(tmp_path: Path, tracks_present: None) -> None:
    """A sorting found only by directory name is recorded, not just logged.

    That fallback *is* the 791094 signature: the asset did not name the recording
    it was sorted from, so the structural match failed.
    """
    config = _config(tmp_path)
    # Mount a sorting whose data_description names a different recording, so the
    # structural lookup misses and only the name walk can find it.
    mount = config.data_root / "ecephys_791094_2025-10-08_sorted_x"
    _mount_sorted_for_coverage(mount, input_recording="ecephys_791094_1999-01-01", collection="ProbeA")

    row = ManifestRow.from_series(
        pd.Series(
            {
                "mouseid": "791094",
                "sorted_recording": "ecephys_791094_2025-10-08_sorted_x",
                "probe_file": "ng",
                "probe_id": "Track1",
                "probe_name": "ProbeA",
            },
            name=0,  # iterrows names each row with its index; from_series reads it
        )
    )
    ok, reason, notes = stages._probe_viability(config, row, find_sorted_session_dirs(config.data_root))

    assert ok and reason is None
    assert any("matched by directory name" in note for note in notes)


class _Row:
    """Minimal ManifestRow stand-in for the track-annotation check."""

    def __init__(self, probe_file: str, histology_track_id: str, annotation_format: str = "json") -> None:
        self.probe_file = probe_file
        self.histology_track_id = histology_track_id
        self.annotation_format = annotation_format


def test_track_annotation_present_missing_file(tmp_path: Path) -> None:
    """A missing annotation JSON is reported as not present."""
    ok, reason = _track_annotation_present(tmp_path, _Row("my_ng", "Track1"))
    assert not ok
    assert "annotation file not found" in str(reason)


def test_track_annotation_present_no_points(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An annotation whose layer has no points is reported as not present."""
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "my_ng.json").write_text("{}")
    monkeypatch.setattr("aind_s3_cache.json_utils.get_json", lambda p: {})
    monkeypatch.setattr(
        "aind_zarr_utils.neuroglancer.neuroglancer_annotations_to_indices",
        lambda data, layer_names=None: ({}, None),
    )

    ok, reason = _track_annotation_present(tmp_path, _Row("my_ng", "Track1"))
    assert not ok
    assert "no track points" in str(reason)


def test_track_annotation_present_ok(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A found annotation with layer points is reported present."""
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "my_ng.json").write_text("{}")
    monkeypatch.setattr("aind_s3_cache.json_utils.get_json", lambda p: {})
    monkeypatch.setattr(
        "aind_zarr_utils.neuroglancer.neuroglancer_annotations_to_indices",
        lambda data, layer_names=None: ({"Track1": [[1, 2, 3], [4, 5, 6]]}, None),
    )

    ok, reason = _track_annotation_present(tmp_path, _Row("my_ng", "Track1"))
    assert ok
    assert reason is None


# ---------------------------------------------------------------------------
# Fan-in merge (stage_pack)
# ---------------------------------------------------------------------------


def _write(path: Path, text: str = "x") -> None:
    """Create a file (and parents) with some content."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _histology_tree(root: Path, mouse_id: str) -> Path:
    """A minimal histology output subtree under ``root/<mouse_id>``."""
    mouse = root / mouse_id
    _write(mouse / "ccf_space_histology" / "histology_ccf.nrrd")
    _write(mouse / "track_data" / "bregma_xyz" / "probeA.npy")
    _write(mouse / "rec1" / "probeA" / "xyz_picks.json")  # histology's slice of the shared dir
    return mouse


def _ephys_tree(root: Path, mouse_id: str, probe: str) -> Path:
    """A minimal ephys output subtree under ``root/<mouse_id>`` for one probe."""
    mouse = root / mouse_id
    _write(mouse / "rec1" / probe / "spikes.times.npy")  # ephys's slice of the shared dir
    _write(mouse / "rec1" / probe / "channels.localCoordinates.npy")
    return mouse


def test_find_mouse_output_trees_across_depths(tmp_path: Path) -> None:
    """Mouse trees are found whether mounted as a direct child or one level deep."""
    data = tmp_path / "data"
    direct = _histology_tree(data / "histology_asset", "791094")  # depth 2
    indexed = _ephys_tree(data / "ephys_asset" / "0", "791094", "probeA")  # depth 3

    found = _find_mouse_output_trees(data, "791094")

    assert {p.resolve() for p in found} == {direct.resolve(), indexed.resolve()}


def test_find_mouse_output_trees_excludes_destination(tmp_path: Path) -> None:
    """The merge destination is never returned as one of its own sources."""
    data = tmp_path / "data"
    _histology_tree(data / "histology_asset", "791094")
    dest = _histology_tree(data / "results", "791094")  # lives under data_root

    found = _find_mouse_output_trees(data, "791094", exclude=dest)

    assert dest.resolve() not in {p.resolve() for p in found}


def test_merge_pipeline_outputs_unions_disjoint_trees(tmp_path: Path) -> None:
    """Histology + per-probe ephys trees union into one complete mouse tree."""
    data = tmp_path / "data"
    _histology_tree(data / "histology_asset", "791094")
    _ephys_tree(data / "ephys_asset_A" / "0", "791094", "probeA")
    _ephys_tree(data / "ephys_asset_B" / "0", "791094", "probeB")
    results = tmp_path / "results"

    merged = merge_pipeline_outputs(data, results, "791094")

    assert merged == results / "791094"
    # Histology outputs are present...
    assert (merged / "ccf_space_histology" / "histology_ccf.nrrd").is_file()
    assert (merged / "rec1" / "probeA" / "xyz_picks.json").is_file()
    # ...alongside both probes' ephys ALF in the shared per-probe dirs.
    assert (merged / "rec1" / "probeA" / "spikes.times.npy").is_file()
    assert (merged / "rec1" / "probeB" / "channels.localCoordinates.npy").is_file()


def test_merge_pipeline_outputs_warns_on_true_overlap(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """A genuine cross-node file overlap is surfaced as a warning (last wins)."""
    import logging

    data = tmp_path / "data"
    _write(data / "node_a" / "791094" / "rec1" / "probeA" / "spikes.times.npy", "a")
    _write(data / "node_b" / "791094" / "rec1" / "probeA" / "spikes.times.npy", "b")
    results = tmp_path / "results"

    with caplog.at_level(logging.WARNING):
        merged = merge_pipeline_outputs(data, results, "791094")

    assert any("overlap" in r.message for r in caplog.records)
    # Deterministic: sources are sorted, so node_b (last) wins.
    assert (merged / "rec1" / "probeA" / "spikes.times.npy").read_text() == "b"


def test_merge_pipeline_outputs_no_tree_raises(tmp_path: Path) -> None:
    """No mouse tree under the mount is a loud failure, not a silent empty pack."""
    with pytest.raises(FileNotFoundError, match="791094"):
        merge_pipeline_outputs(tmp_path / "data", tmp_path / "results", "791094")


def _stream_cfg(mouse_id: str, recording_id: str, collection: str | None) -> dict[str, object]:
    """A minimal ephys fan-out config as stage_discover would emit it."""
    return {
        "mouseid": mouse_id,
        "sorted_recording": f"{recording_id}_sorted",
        "recording_id": recording_id,
        "ephys_collection": collection,
        "surface_finding": None,
    }


def test_stage_ephys_namespaces_output_by_unit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Parallel ephys workers write disjoint top-level /results names, then pack unions them.

    Every worker writing a bare ``/results/<mouse_id>/`` would collide at the
    downstream Collect (``pack``) staging step. Namespacing by fan-out unit keeps
    the top-level names disjoint while staying mergeable by ``pack``.
    """
    config = SimpleNamespace(
        results_root=tmp_path / "results",
        data_root=tmp_path / "data",
        num_parallel_jobs=1,
    )

    def _fake_run(sorted_recording, recording_id, collection, surface, out, data_root, *, num_parallel_jobs):  # type: ignore[no-untyped-def]
        # Emulate a real per-probe ALF write into this unit's mouse tree.
        mouse_dir = out.tracks_root.parent
        _write(mouse_dir / "rec1" / str(collection) / "spikes.times.npy")

    monkeypatch.setattr(stages, "run_ephys_for_stream", _fake_run)

    stages.stage_ephys(config, stream_config=_stream_cfg("791094", "ecephys_791094_A", "ProbeA"))
    stages.stage_ephys(config, stream_config=_stream_cfg("791094", "ecephys_791094_A", "ProbeB"))

    # Disjoint top-level folders -> no Collect-stage "input file name collision".
    top_level = {p.name for p in config.results_root.iterdir()}
    assert top_level == {"ecephys_791094_A__ProbeA", "ecephys_791094_A__ProbeB"}
    assert (config.results_root / "ecephys_791094_A__ProbeA" / "791094").is_dir()

    # ...and pack's layout-agnostic merge still unions both units into one tree.
    merged = merge_pipeline_outputs(config.results_root, tmp_path / "packed", "791094")
    assert (merged / "rec1" / "ProbeA" / "spikes.times.npy").is_file()
    assert (merged / "rec1" / "ProbeB" / "spikes.times.npy").is_file()


# ------------------------------------------------------- stage_ephys_launch --


def _mount_sorted(data_root: Path, input_recording: str, *, mount_name: str = "sorted") -> None:
    """Mount a sorted asset under a fixed pipeline slot (name != asset name).

    The launcher must resolve it by content: a ``spikesorted`` child marks it as a
    sorted asset, and ``data_description.json``'s ``input_data_name`` identifies the
    raw recording it derives from (== ManifestRow.recording_id).
    """
    root = data_root / mount_name
    (root / "spikesorted").mkdir(parents=True)
    (root / "data_description.json").write_text(json.dumps({"input_data_name": input_recording}))


def test_stage_ephys_launch_scopes_to_the_mounted_sort(tmp_path: Path) -> None:
    """ephys-launch emits configs only for the sort whose sorted asset is mounted.

    The asset mounts under a fixed slot (``/data/sorted``), so scoping is by content
    (``input_data_name`` == recording_id), not by matching the asset name.
    """
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        [
            "786867,ecephys_A_sorted_x,ng,T1,ProbeA,",
            "786867,ecephys_A_sorted_x,ng,T2,ProbeB,",
            "786867,ecephys_B_sorted_y,ng,T3,ProbeA,",
        ],
    )
    # Recording A's sorted asset is mounted under the fixed slot /data/sorted.
    _mount_sorted(config.data_root, "ecephys_A")

    written = stages.stage_ephys_launch(config)

    # A's two collections emit; recording B (not mounted) is excluded.
    assert len(written) == 2
    for path in written:
        cfg = json.loads(path.read_text())
        assert cfg["recording_id"] == "ecephys_A"
        assert cfg["sorted_recording"] == "ecephys_A_sorted_x"


def test_stage_ephys_launch_raises_when_no_sorted_mounted(tmp_path: Path) -> None:
    """A missing sorted asset is a loud failure, not a silent empty fan-out."""
    config = _config(tmp_path)
    _write_manifest(config.manifest_csv, ["786867,ecephys_A_sorted_x,ng,T1,ProbeA,"])
    config.data_root.mkdir(parents=True)  # empty /data -- no spikesorted/ marker

    with pytest.raises(FileNotFoundError, match="spikesorted"):
        stages.stage_ephys_launch(config)


def test_stage_ephys_launch_raises_when_mounted_sort_absent_from_manifest(tmp_path: Path) -> None:
    """A mounted sort with no matching manifest rows fails loudly (0 configs)."""
    config = _config(tmp_path)
    _write_manifest(config.manifest_csv, ["786867,ecephys_A_sorted_x,ng,T1,ProbeA,"])
    # The mounted sort derives from a raw the (filtered) manifest never mentions.
    _mount_sorted(config.data_root, "ecephys_ZZZ")

    with pytest.raises(RuntimeError, match="0 fan-out configs"):
        stages.stage_ephys_launch(config)


# ------------------------------------------------------ stage_ephys_collect --


def test_stage_ephys_collect_merges_from_data_using_manifest_mouse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ephys-collect merges the /data trees into /results using the manifest's mouse id."""
    config = _config(tmp_path)
    _write_manifest(config.manifest_csv, ["786867,ecephys_A_sorted_x,ng,T1,ProbeA,"])

    calls: dict[str, object] = {}
    import aind_ibl_ephys_alignment_preprocessing.datapackage as dp

    monkeypatch.setattr(
        dp,
        "merge_pipeline_outputs",
        lambda src, dst, mouse: calls.update(src=src, dst=dst, mouse=mouse),
    )

    stages.stage_ephys_collect(config, merge_from=config.data_root)

    assert calls == {"src": config.data_root, "dst": config.results_root, "mouse": "786867"}


def test_coverage_units_match_the_emitted_config_names(tmp_path: Path, all_viable: None) -> None:
    """``units_emitted`` uses the same key the fan-out configs carry.

    That key is what the trigger's resolved-asset map is keyed by, and what the
    run record joins on. If the two ever diverge, the record still validates --
    it just answers nothing.
    """
    config = _config(tmp_path)
    _write_manifest(
        config.manifest_csv,
        [
            "791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,",
            "791094,ecephys_791094_2025-10-09_sorted_y,ng,Track2,ProbeB,",
        ],
    )

    written = stage_discover(config)
    coverage = _read_coverage(config)

    config_names = sorted(json.loads(path.read_text())["name"] for path in written)
    assert config_names == coverage.units_emitted


def test_pack_carries_the_run_record_forward(tmp_path: Path, all_viable: None) -> None:
    """End to end across the handoff: discover writes it, pack completes it.

    Exercises the real serialized file rather than a hand-built record, so a
    change on either side that breaks the join fails here. This step is what
    makes the record durable at all -- discover's copy lives in an intermediate
    the launcher deletes on success.
    """
    from aind_ibl_ephys_alignment_preprocessing.stages import _carry_coverage_forward

    discover = _config(tmp_path / "discover")
    discover.manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_manifest(
        discover.manifest_csv,
        ["791094,ecephys_791094_2025-10-08_sorted_x,ng,Track1,ProbeA,"],
    )
    stage_discover(discover)

    # pack runs against discover's *filtered* manifest, mounted under /data.
    pack = SimpleNamespace(
        manifest_csv=discover.results_root / "manifest.csv",
        results_root=tmp_path / "pack_results",
        data_root=tmp_path / "data",
        skip_ephys=False,
        allow_partial=False,
    )
    # Exactly the shape the trigger's _asset_resolution_arg emits.
    resolution = json.dumps(
        {
            "units": {
                "ecephys_791094_2025-10-08__ProbeA": {
                    "raw": {"id": "raw-id", "name": "ecephys_791094_2025-10-08"},
                    "sorted": {"id": "sorted-id", "name": "ecephys_791094_2025-10-08_sorted_x"},
                }
            },
            "smartspim": {"id": "spim-id", "name": "SmartSPIM_791094_stitched"},
        }
    )

    _carry_coverage_forward(pack, resolution)

    coverage = RunCoverage.read(pack.results_root / COVERAGE_FILENAME)
    assert set(coverage.unit_assets) == set(coverage.units_emitted)
    assert coverage.is_complete is True
    assert coverage.unit_assets["ecephys_791094_2025-10-08__ProbeA"].sorted.id == "sorted-id"
    assert coverage.smartspim.id == "spim-id"


def test_pack_survives_a_missing_run_record(tmp_path: Path) -> None:
    """The monolith has no discover stage; bookkeeping must not fail the run."""
    from aind_ibl_ephys_alignment_preprocessing.stages import _carry_coverage_forward

    config = _config(tmp_path)
    config.manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    config.manifest_csv.write_text("mouseid\n791094\n")

    _carry_coverage_forward(config, None)

    assert not (config.results_root / COVERAGE_FILENAME).exists()
