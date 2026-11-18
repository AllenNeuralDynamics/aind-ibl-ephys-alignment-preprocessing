"""Persist conversion decisions through discovery, fan-out, and pack."""

from __future__ import annotations

import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from aind_ephys_ibl_gui_conversion.conversion_plan import (
    PLAN_FILENAME,
    PLAN_VERSION,
    build_conversion_plan,
    inventory_digest,
)
from aind_ephys_ibl_gui_conversion.inventory import INVENTORY_FILENAME, inventory_slice, write_inventory

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow

logger = logging.getLogger(__name__)
OVERRIDE_COLUMNS = ("main_recording_id", "sorting_family_id")


def _one_value(frame: pd.DataFrame, column: str) -> str | None:
    """Allow blank row overrides while rejecting disagreements within a unit."""
    if column not in frame:
        return None
    values = {str(value).strip() for value in frame[column].dropna() if str(value).strip()}
    if len(values) > 1:
        raise ValueError(f"Conflicting {column} overrides for one conversion unit: {sorted(values)}")
    return next(iter(values), None)


def build_discovery_plans(inventory: dict[str, Any], frame: pd.DataFrame) -> dict[str, Any]:
    """Plan every requested unit before manifest viability can drop evidence."""
    grouped: dict[str, list[int]] = {}
    for index, (_, series) in enumerate(frame.iterrows()):
        row = ManifestRow.from_series(series)
        if row.ephys_collection is not None:
            unit = f"{row.recording_id}__{row.ephys_collection}"
            grouped.setdefault(unit, []).append(index)
    units = {}
    for unit, positions in grouped.items():
        rows = frame.iloc[positions]
        row = ManifestRow.from_series(rows.iloc[0])
        sort = _one_value(rows, "sorted_recording")
        overrides = {column: _one_value(rows, column) for column in OVERRIDE_COLUMNS}
        plan = build_conversion_plan(
            inventory,
            str(row.ephys_collection),
            sort,
            main_recording_id=overrides["main_recording_id"],
            sorting_family_id=overrides["sorting_family_id"],
            surface_finding=_one_value(rows, "surface_finding"),
        )
        plan["unit"] = unit
        units[unit] = plan
        main = plan["main_recording"] or {}
        logger.info(
            "[discover] conversion plan %s: main=%s duration_s=%s reason=%s analyzers=%d status=%s",
            unit,
            main.get("id"),
            main.get("duration_s"),
            main.get("reason"),
            len(plan["analyzers"]),
            plan["status"],
        )
        for finding in plan["findings"]:
            logger.warning("[discover] conversion plan %s: %s", unit, finding["code"])
    logger.info("[discover] conversion plans: %s", dict(Counter(p["status"] for p in units.values())))
    return dict(
        plan_version=PLAN_VERSION, inventory_digest=inventory_digest(inventory), execution_enabled=False, units=units
    )


def read_plan_registry(directory: Path, inventory: dict[str, Any]) -> dict[str, Any] | None:
    """Accept legacy inputs, but refuse decisions tied to a different inventory."""
    path = directory / PLAN_FILENAME
    if not path.is_file():
        return None
    plans: dict[str, Any] = json.loads(path.read_text())
    if plans.get("plan_version") != PLAN_VERSION:
        raise ValueError(f"Unsupported conversion plan version in {path}")
    if plans.get("inventory_digest") != inventory_digest(inventory):
        raise ValueError(f"Conversion plans do not match the discovery inventory in {directory}")
    return plans


def write_unit_plan(
    inventory: dict[str, Any],
    cfg: dict[str, Any],
    destination: Path,
    registry: dict[str, Any] | None = None,
) -> None:
    """Stage an existing decision unchanged, or plan legacy inventory inputs."""
    unit = cfg.get("name") or f"{cfg.get('recording_id', cfg.get('sorted_recording'))}__{cfg['ephys_collection']}"
    if registry is not None:
        if unit not in registry["units"]:
            raise ValueError(f"No conversion plan for emitted unit {unit}")
        plan = registry["units"][unit]
    else:
        plan = build_conversion_plan(
            inventory,
            cfg["ephys_collection"],
            cfg.get("sorted_recording"),
            surface_finding=cfg.get("surface_finding"),
        )
        plan["unit"] = unit
    scoped = inventory_slice(
        inventory, cfg.get("ephys_collection"), cfg.get("sorted_recording"), cfg.get("surface_finding")
    )
    _validate_unit_plan(plan, scoped, unit)
    write_inventory(plan, destination / PLAN_FILENAME)
    cfg["conversion_plan"] = PLAN_FILENAME


def _validate_unit_plan(plan: dict[str, Any], inventory: dict[str, Any], unit: str | None = None) -> None:
    """Reject stale, unsupported, or misrouted decisions before staging them."""
    if plan.get("inventory_digest") != inventory_digest(inventory):
        raise ValueError("Worker conversion plan does not match its inventory")
    if plan.get("plan_version") != PLAN_VERSION:
        raise ValueError("Unsupported worker conversion plan version")
    if unit is not None and plan.get("unit") != unit:
        raise ValueError(f"Conversion plan belongs to another unit, expected {unit}")


def carry_unit_diagnostics(cfg: dict[str, Any], destination: Path) -> tuple[dict[str, Any], dict[str, Any]] | None:
    """Validate and retain worker decisions, returning them for extraction."""
    inventory = None
    if cfg.get("_inventory_path"):
        inventory = json.loads(Path(cfg["_inventory_path"]).read_text())
        write_inventory(inventory, destination / INVENTORY_FILENAME)
    if cfg.get("_conversion_plan_path"):
        plan = json.loads(Path(cfg["_conversion_plan_path"]).read_text())
        if inventory is None:
            raise ValueError("Worker conversion plan does not match its inventory")
        _validate_unit_plan(plan, inventory, cfg.get("name"))
        write_inventory(plan, destination / PLAN_FILENAME)
        return plan, inventory
    return None


def bind_extraction(
    plan: dict[str, Any],
    inventory: dict[str, Any],
    data_root: Path,
    sorted_root: Path,
    raw_root: Path | None,
    results_folder: Path,
) -> Any:
    """Bind portable source keys to worker mounts before opening planned data."""
    from aind_ephys_ibl_gui_conversion.planned_extraction import prepare_extraction, validate_execution_plan

    from aind_ibl_ephys_alignment_preprocessing.ephys import (
        _find_marked_dirs,
        _read_data_description_name,
        find_session_dir,
    )

    validate_execution_plan(plan, inventory)
    roots: dict[str, Path] = {}
    sorted_keys = {r["asset"] for r in plan["analyzers"]}
    if len(sorted_keys) != 1:
        raise ValueError("Planned extraction requires one sorted source asset")
    roots[next(iter(sorted_keys))] = sorted_root
    references = [r for r in plan["continuous_inputs"] if r["selected"]]
    references += [r["paired_lfp"] for r in references if r.get("paired_lfp")]
    assets = {a["key"]: a for a in inventory["assets"]}
    raw_keys = {r["asset"] for r in references}
    mounted = {p.resolve() for p in _find_marked_dirs(data_root, ("ecephys_compressed", "ecephys")).values()}
    if raw_root is not None:
        mounted.add(raw_root.resolve())
    for key in raw_keys:
        asset = assets[key]
        matches: set[Path] = set()
        identities = {asset.get("name"), asset.get("id"), key, *asset.get("aliases", [])} - {None}
        for name in identities:
            found = find_session_dir(data_root, str(name))
            if found is not None:
                matches.add(found.resolve())
        matches.update(p for p in mounted if _read_data_description_name(p) in identities)
        # A sole raw mount is the worker's explicit raw slot. Multiple raw
        # sources require individual identity matches, never positional pairing.
        if not matches and len(raw_keys) == 1 and len(mounted) == 1:
            matches.update(mounted)
        if len(matches) != 1:
            raise ValueError(f"Cannot uniquely bind planned asset {key}: {sorted(matches)}")
        roots[key] = matches.pop()
    if len({roots[key].resolve() for key in raw_keys}) != len(raw_keys):
        raise ValueError("Different planned raw assets resolve to the same mount")
    return prepare_extraction(plan, inventory, roots, results_folder)
