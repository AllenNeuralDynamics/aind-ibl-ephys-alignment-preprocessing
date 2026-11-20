"""Run-level registration settings survive staged workflows and legacy CSV input."""

import json
from pathlib import Path

import pandas as pd
import pytest

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.discovery import manifest_per_brain_path
from aind_ibl_ephys_alignment_preprocessing.registration_settings import (
    SETTINGS_FILENAME,
    resolve_registration_settings,
    write_registration_settings,
)
from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

FIELDS = ("registration_asset", "template_to_ccf_asset")


def _config(tmp_path, field=None, legacy=None, configured=None):
    root = tmp_path / "data"
    root.mkdir(parents=True, exist_ok=True)
    header = "mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection"
    row = "853781,session_sorted,annotations,track,ProbeA"
    if legacy is not None:
        header += f",{field}"
        row += f",{legacy}"
    (root / "manifest.csv").write_text(f"{header}\n{row}\n")
    return PipelineConfig(
        data_root=root,
        results_root=tmp_path / "results",
        neuroglancer_file=Path("ng.json"),
        manifest_csv=Path("manifest.csv"),
        skip_ephys=True,
        **({field: configured} if configured is not None else {}),
    )


@pytest.mark.parametrize("field", FIELDS)
@pytest.mark.parametrize("source", ["run", "csv", "both"])
def test_run_and_legacy_values_resolve_to_the_same_directory(tmp_path, field, source):
    value = "registration/transforms"
    config = _config(
        tmp_path,
        field,
        legacy=value if source in ("csv", "both") else None,
        configured=value if source in ("run", "both") else None,
    )
    assert manifest_per_brain_path(config, field) == Path(value)


@pytest.mark.parametrize("field", FIELDS)
def test_run_and_csv_conflict_even_when_they_name_the_same_asset(tmp_path, field):
    config = _config(tmp_path, field, legacy="registration/a", configured="registration/b")
    with pytest.raises(ValueError, match=f"Conflicting {field}"):
        resolve_registration_settings(config)


def test_matching_relative_and_absolute_paths_are_accepted(tmp_path):
    config = _config(tmp_path, "registration_asset", legacy="registration/transforms/")
    config = config.model_copy(update={"registration_asset": config.data_root / "registration/transforms"})
    path = write_registration_settings(config, config.results_root)
    assert json.loads(path.read_text())["registration_asset"] == "registration/transforms"


@pytest.mark.parametrize("field", FIELDS)
@pytest.mark.parametrize("source", ["run", "csv", "prior_results"])
def test_discoverys_recorded_default_cannot_be_replaced_downstream(tmp_path, field, source):
    config = _config(tmp_path)
    write_registration_settings(config, config.manifest_csv.parent)
    if source == "run":
        config = config.model_copy(update={field: Path("different/transforms")})
        inherited = None
    elif source == "csv":
        frame = pd.read_csv(config.manifest_csv)
        frame[field] = "different/transforms"
        frame.to_csv(config.manifest_csv, index=False)
        inherited = None
    else:
        inherited = tmp_path / "prior"
        inherited.mkdir()
        (inherited / SETTINGS_FILENAME).write_text(json.dumps({"schema_version": 1, field: "different/transforms"}))
    with pytest.raises(ValueError, match=f"Conflicting {field}.*recorded pipeline default"):
        resolve_registration_settings(config, inherited_from=inherited)


def test_recorded_defaults_remain_usable_without_overrides(tmp_path):
    config = _config(tmp_path)
    write_registration_settings(config, config.manifest_csv.parent)
    assert resolve_registration_settings(config) == dict.fromkeys(FIELDS)


def test_discovery_keeps_run_settings_out_of_track_csv_and_stages_them_separately(tmp_path, monkeypatch):
    config = _config(tmp_path, "registration_asset", configured="registration/transforms")
    monkeypatch.setattr(stages, "_track_annotation_present", lambda *args: (True, None))
    stages.stage_discover(config)
    filtered = pd.read_csv(config.results_root / "manifest.csv")
    assert not set(FIELDS) & set(filtered.columns)
    saved = json.loads((config.results_root / SETTINGS_FILENAME).read_text())
    assert saved["registration_asset"] == "registration/transforms"

    # A downstream worker has different mount roots and receives no override CLI
    # flags. The settings accompanying discovery's manifest remain authoritative.
    downstream = PipelineConfig(
        data_root=tmp_path / "other_mount",
        results_root=tmp_path / "packed",
        neuroglancer_file=Path("ng.json"),
        manifest_csv=config.results_root / "manifest.csv",
    )
    assert manifest_per_brain_path(downstream, "registration_asset") == Path("registration/transforms")
    final_path = write_registration_settings(downstream, downstream.results_root / "853781")
    assert json.loads(final_path.read_text()) == saved
    changed = downstream.model_copy(update={"registration_asset": Path("different/transforms")})
    with pytest.raises(ValueError, match="Conflicting registration_asset"):
        resolve_registration_settings(changed)


@pytest.mark.parametrize("saved", [{"schema_version": 2}, {"schema_version": 1, "registration_asset": []}])
def test_invalid_discovery_settings_are_not_silently_ignored(tmp_path, saved):
    config = _config(tmp_path)
    (config.manifest_csv.parent / SETTINGS_FILENAME).write_text(json.dumps(saved))
    with pytest.raises(ValueError):
        resolve_registration_settings(config)


def test_pack_preserves_registration_settings_inherited_from_prior_results(tmp_path, monkeypatch):
    from aind_ibl_ephys_alignment_preprocessing import pipeline

    config = _config(tmp_path)
    source = tmp_path / "prior/853781"
    source.mkdir(parents=True)
    (source / "datapackage.json").write_text("{}")
    saved = {"schema_version": 1, "registration_asset": "reg/transforms", "template_to_ccf_asset": None}
    (source / SETTINGS_FILENAME).write_text(json.dumps(saved))

    def check_settings(actual_config):
        assert actual_config.registration_asset == Path("reg/transforms")
        assert manifest_per_brain_path(actual_config, "registration_asset") == Path("reg/transforms")
        return object()

    monkeypatch.setattr(pipeline, "find_asset_info", check_settings)
    monkeypatch.setattr(pipeline, "infer_process_results_from_outputs", lambda *a, **k: [])
    monkeypatch.setattr(pipeline, "build_datapackage", lambda *a: object())
    monkeypatch.setattr(pipeline, "producer_asset_overrides", lambda *a: {})
    monkeypatch.setattr(pipeline, "write_datapackage", lambda dp, root, **kw: root / "datapackage.json")
    path = stages.stage_pack(config, source_results=source)
    assert json.loads((path.parent / SETTINGS_FILENAME).read_text()) == saved


def test_standalone_cli_accepts_registration_settings(tmp_path, monkeypatch):
    import sys

    from aind_ibl_ephys_alignment_preprocessing import validation
    from aind_ibl_ephys_alignment_preprocessing.scripts import run

    seen = []

    class Validator:
        def __init__(self, config):
            seen.append(config)

        def validate_all(self):
            return []

        def print_summary(self, results):
            pass

        def has_errors(self, results):
            return False

    monkeypatch.setattr(validation, "PipelineValidator", Validator)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aind-ibl-preprocess",
            "--data-root",
            str(tmp_path),
            "--results-root",
            str(tmp_path / "out"),
            "--manifest",
            "manifest.csv",
            "--neuroglancer",
            "ng.json",
            "--validate-only",
            "--registration-asset",
            "reg/transforms",
            "--template-to-ccf-asset",
            "template/transforms",
        ],
    )
    with pytest.raises(SystemExit) as result:
        run.main()
    assert result.value.code == 0
    assert seen[0].registration_asset == Path("reg/transforms")
    assert seen[0].template_to_ccf_asset == Path("template/transforms")
