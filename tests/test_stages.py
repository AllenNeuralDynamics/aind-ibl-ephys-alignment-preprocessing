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

import pytest
from aind_code_ocean_pipeline_utils.role_dispatch import find_stream_config

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.stages import (
    EPHYS_STREAM_MARKER,
    _ephys_unit_name,
    _track_annotation_present,
    stage_discover,
)

_MANIFEST_HEADER = "mouseid,sorted_recording,probe_file,probe_id,probe_name,surface_finding"


def _write_manifest(path: Path, rows: list[str]) -> None:
    """Write a minimal manifest CSV with the given data rows."""
    path.write_text("\n".join([_MANIFEST_HEADER, *rows]) + "\n")


def _config(tmp_path: Path, *, skip_ephys: bool = False) -> SimpleNamespace:
    """A stand-in PipelineConfig with the attributes stage_discover reads."""
    return SimpleNamespace(
        manifest_csv=tmp_path / "manifest.csv",
        results_root=tmp_path / "results",
        data_root=tmp_path / "data",
        skip_ephys=skip_ephys,
    )


@pytest.fixture
def all_viable(monkeypatch: pytest.MonkeyPatch) -> None:
    """Force every probe viable so emission logic can be tested in isolation."""
    monkeypatch.setattr(stages, "_probe_viability", lambda config, mr: (True, None))


def test_public_stage_functions_exist() -> None:
    """The four pipeline stages are importable from the module."""
    for name in ("stage_discover", "stage_histology", "stage_ephys", "stage_pack"):
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
    """Non-viable probes are dropped from both the manifest and the ephys configs."""
    config = _config(tmp_path)
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
        lambda config, mr: (mr.ephys_collection != "ProbeB", "skipped for test"),
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
