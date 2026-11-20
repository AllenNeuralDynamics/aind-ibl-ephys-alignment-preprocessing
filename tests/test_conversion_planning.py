"""Conversion decisions survive dispatch and control both extractors."""

import json
import shutil
from types import SimpleNamespace

import pandas as pd
import pytest
from aind_ephys_ibl_gui_conversion.conversion_plan import PLAN_FILENAME, inventory_digest

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.conversion_planning import (
    build_discovery_plans,
    carry_unit_diagnostics,
    read_plan_registry,
)
from aind_ibl_ephys_alignment_preprocessing.ephys_inventory import carry_inventory, discover_inventory

from .test_ephys_inventory import _config, _sources, _write


def test_canonical_unit_worker_and_pack_plans(tmp_path, monkeypatch):
    config = _config(tmp_path)
    sources = _sources(config)
    monkeypatch.setattr(stages, "_probe_viability", lambda *args: (True, None, []))
    (path,) = stages.stage_discover(config, inventory_sources=sources)
    inventory = json.loads((config.results_root / "ephys_inventory.json").read_text())
    registry = read_plan_registry(config.results_root, inventory)
    assert registry["execution_enabled"] is False
    assert set(registry["units"]) == {"rec__ProbeA"}
    plan = registry["units"]["rec__ProbeA"]
    assert plan["main_recording"]["experiment"] == "experiment7"
    assert plan["main_recording"]["reason"] == "longest_recorded_duration"

    flattened = tmp_path / "flattened"
    shutil.copytree(path.parent, flattened)
    cfg = stages.find_ephys_stream_config(flattened)
    assert cfg["conversion_plan"] == PLAN_FILENAME
    worker = tmp_path / "worker/rec__ProbeA"
    carry_unit_diagnostics(cfg, worker)
    assert json.loads((worker / PLAN_FILENAME).read_text()) == plan

    # Both saved artifacts reach the extraction wrapper together.
    calls = []
    monkeypatch.setattr(stages, "run_ephys_for_stream", lambda *a, **k: calls.append(k))
    stages.stage_ephys(config.model_copy(update={"results_root": tmp_path / "extract"}), stream_config=cfg)
    assert len(calls) == 1
    actual_plan, actual_inventory = calls[0]["conversion_decisions"]
    assert actual_plan == plan
    assert actual_inventory["scope"]["collection"] == "ProbeA"

    pack = config.model_copy(
        update={"manifest_csv": config.results_root / "manifest.csv", "results_root": tmp_path / "pack"}
    )
    carry_inventory(pack)
    assert json.loads((pack.results_root / "123" / PLAN_FILENAME).read_text()) == registry


def test_stale_plan_registry_and_unit_plan_are_refused(tmp_path):
    config = _config(tmp_path)
    inventory = discover_inventory(config, [], _sources(config))
    registry = json.loads((config.results_root / PLAN_FILENAME).read_text())
    registry["inventory_digest"] = "stale"
    _write(config.results_root / PLAN_FILENAME, registry)
    with pytest.raises(ValueError, match="do not match"):
        read_plan_registry(config.results_root, inventory)
    inv_path = tmp_path / "input/ephys_inventory.json"
    plan_path = tmp_path / "input" / PLAN_FILENAME
    _write(inv_path, inventory)
    _write(plan_path, {"inventory_digest": "stale"})
    with pytest.raises(ValueError, match="does not match"):
        carry_unit_diagnostics(
            {"_inventory_path": str(inv_path), "_conversion_plan_path": str(plan_path)}, tmp_path / "worker"
        )


def test_rows_in_one_unit_must_agree_on_overrides(tmp_path):
    config = _config(tmp_path)
    inventory = discover_inventory(config, [], _sources(config))
    frame = pd.read_csv(config.manifest_csv)
    frame["main_recording_id"] = "first"
    other = frame.copy()
    other["main_recording_id"] = "second"
    with pytest.raises(ValueError, match="Conflicting main_recording_id"):
        build_discovery_plans(inventory, pd.concat([frame, other]))


def test_override_survives_per_sort_launcher_with_relocated_outputs(tmp_path):
    config = _config(tmp_path)
    sources = _sources(config)
    raw = config.data_root / "raw/ecephys_compressed"
    original = next(raw.iterdir())
    longer = raw / original.name.replace("experiment7", "experiment8")
    shutil.copytree(original, longer)
    attrs = json.loads((longer / ".zattrs").read_text())
    attrs["annotations"]["experiment_name"] = "experiment8"
    _write(longer / ".zattrs", attrs)
    _write(longer / "traces_seg0/.zarray", {"shape": [300000, 2]})
    inventory = discover_inventory(config, [], sources)
    default = read_plan_registry(config.results_root, inventory)["units"]["rec__ProbeA"]
    assert default["main_recording"]["experiment"] == "experiment8"
    assert default["status"] == "blocked"  # no sorting for the longer block

    chosen = next(r for r in inventory["recordings"] if r["experiment"] == "experiment7")
    frame = pd.read_csv(config.manifest_csv)
    frame["main_recording_id"] = chosen["id"]
    frame.to_csv(config.manifest_csv, index=False)
    inventory = discover_inventory(config, [], sources)
    expected = read_plan_registry(config.results_root, inventory)["units"]["rec__ProbeA"]
    assert expected["main_recording"]["reason"] == "explicit_override"
    _write(config.data_root / "ephys_inventory.json", inventory)
    shutil.copy2(config.results_root / PLAN_FILENAME, config.data_root / PLAN_FILENAME)
    (config.data_root / "sorted/spikesorted").mkdir()
    _write(config.data_root / "sorted/data_description.json", {"source_data": ["rec"]})

    launcher = config.model_copy(update={"results_root": tmp_path / "launcher"})
    (path,) = stages.stage_ephys_launch(launcher)
    assert json.loads((path.parent / PLAN_FILENAME).read_text()) == expected
    # Older inventory-only intermediates can also resolve the CSV override.
    (config.data_root / PLAN_FILENAME).unlink()
    legacy = config.model_copy(update={"results_root": tmp_path / "legacy-launcher"})
    (path,) = stages.stage_ephys_launch(legacy)
    assert json.loads((path.parent / PLAN_FILENAME).read_text()) == expected


def test_rejected_manifest_rows_still_have_plans(tmp_path, monkeypatch):
    from aind_ibl_ephys_alignment_preprocessing.coverage import CoverageError

    config = _config(tmp_path)
    sources = _sources(config)
    monkeypatch.setattr(stages, "_probe_viability", lambda *a: (False, "no track", []))
    with pytest.raises(CoverageError):
        stages.stage_discover(config, inventory_sources=sources)
    inventory = json.loads((config.results_root / "ephys_inventory.json").read_text())
    registry = read_plan_registry(config.results_root, inventory)
    assert "rec__ProbeA" in registry["units"]
    assert registry["inventory_digest"] == inventory_digest(inventory)


@pytest.mark.parametrize("fail", [False, True])
def test_worker_shares_execution_and_records_outcome(tmp_path, monkeypatch, fail):
    import aind_ephys_ibl_gui_conversion.ephys as converter

    from aind_ibl_ephys_alignment_preprocessing import conversion_planning, ephys

    config = _config(tmp_path)
    (config.data_root / "sorted/spikesorted").mkdir(parents=True)
    (config.data_root / "raw/ecephys_compressed").mkdir(parents=True)
    statuses, calls = [], []
    execution = SimpleNamespace(
        check_request=lambda collection: None,
        stream=SimpleNamespace(output_folder=tmp_path / "output", blocks=[]),
        write_receipt=lambda destination, status, **kwargs: statuses.append(status),
    )
    monkeypatch.setattr(conversion_planning, "bind_extraction", lambda *args: execution)

    def extract(*args, **kwargs):
        calls.append(kwargs["execution"])
        if fail:
            raise RuntimeError("extraction failed")

    monkeypatch.setattr(converter, "extract_spikes", extract)
    monkeypatch.setattr(converter, "extract_continuous", extract)
    args = ("sort", "rec", "ProbeA", None, stages.prepare_result_dirs("123", config.results_root), config.data_root)
    if fail:
        with pytest.raises(RuntimeError, match="extraction failed"):
            ephys.run_ephys_for_stream(*args, conversion_decisions=({}, {}))
    else:
        ephys.run_ephys_for_stream(*args, conversion_decisions=({}, {}))
    assert calls == [execution, execution]
    assert statuses == ["started", "failed" if fail else "complete"]


def test_binder_uses_exact_asset_roots(tmp_path, monkeypatch):
    import aind_ephys_ibl_gui_conversion.planned_extraction as extraction

    from aind_ibl_ephys_alignment_preprocessing.conversion_planning import bind_extraction

    main, secondary, sorted_root = (tmp_path / name for name in ("raw", "surface", "sorted"))
    for root in (main, secondary, sorted_root):
        root.mkdir()
    inv = dict(assets=[dict(key="r", name="raw"), dict(key="s", name="surface")])
    plan = dict(
        plan_version="1.0.0",
        inventory_digest=inventory_digest(inv),
        status="ready",
        findings=[],
        main_recording=dict(id="main", asset="r"),
        analyzers=[dict(id="a", asset="sort")],
        sorting_family="f",
        sorting_families=[dict(id="f", analyzer_ids=["a"])],
        continuous_inputs=[
            dict(id="main", asset="r", selected=True, role="main"),
            dict(id="secondary", asset="s", selected=True, role="surface"),
        ],
    )
    captured = []
    monkeypatch.setattr(extraction, "prepare_extraction", lambda *args: captured.append(args))
    bind_extraction(plan, inv, tmp_path, sorted_root, main, tmp_path / "results")
    assert captured[0][2] == dict(r=main, s=secondary, sort=sorted_root)

    # Renamed fixed slots still bind by asset identity; main is not inferred
    # from the sorting's single-input resolver or from mount ordering.
    main.rename(tmp_path / "slot_two")
    secondary.rename(tmp_path / "slot_one")
    for slot, name in (("slot_two", "raw"), ("slot_one", "surface")):
        _write(tmp_path / slot / "data_description.json", dict(name=name))
        (tmp_path / slot / "ecephys_compressed").mkdir()
    bind_extraction(plan, inv, tmp_path, sorted_root, None, tmp_path / "results")
    assert captured[-1][2] == dict(r=tmp_path / "slot_two", s=tmp_path / "slot_one", sort=sorted_root)
