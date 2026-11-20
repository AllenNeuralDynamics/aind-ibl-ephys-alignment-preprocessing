"""Resolve brain-wide registration choices without requiring track CSV columns."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

REGISTRATION_FIELDS = ("registration_asset", "template_to_ccf_asset")
SETTINGS_FILENAME = "registration_settings.json"


def _saved_settings(manifest: Path) -> dict[str, str | None]:
    """Read discovery's portable settings beside its staged manifest."""
    path = manifest.parent / SETTINGS_FILENAME
    if not path.is_file():
        return {}
    saved = json.loads(path.read_text())
    if not isinstance(saved, dict) or saved.get("schema_version") != 1:
        raise ValueError(f"Unsupported registration settings in {path}")
    result: dict[str, str | None] = {}
    for name in REGISTRATION_FIELDS:
        value = saved.get(name)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"Invalid {name} in {path}")
        result[name] = value
    return result


def resolve_registration_settings(
    config: PipelineConfig, *, inherited_from: Path | None = None
) -> dict[str, Path | None]:
    """Merge run settings, legacy CSV columns, and discovery settings; refuse conflicts."""
    legacy: dict[str, set[str]] = {name: set() for name in REGISTRATION_FIELDS}
    if config.manifest_csv.is_file():
        with config.manifest_csv.open(newline="", encoding="utf-8-sig") as source:
            for row in csv.DictReader(source):
                for name in REGISTRATION_FIELDS:
                    if value := str(row.get(name) or "").strip():
                        legacy[name].add(value)
    saved = _saved_settings(config.manifest_csv)
    inherited = _saved_settings(inherited_from / "manifest.csv") if inherited_from is not None else {}
    result: dict[str, Path | None] = {}
    for name in REGISTRATION_FIELDS:
        if len(legacy[name]) > 1:
            raise ValueError(
                f"{name} must name one directory per manifest; {config.manifest_csv} has {sorted(legacy[name])}"
            )
        configured = getattr(config, name, None)
        candidates = [
            Path(value)
            for value in (configured, *sorted(legacy[name]), saved.get(name), inherited.get(name))
            if value is not None
        ]
        identities = {(config.data_root / path).resolve() for path in candidates}
        if len(identities) > 1:
            raise ValueError(
                f"Conflicting {name} in run settings, manifest, or discovery settings: "
                f"{[str(path) for path in candidates]}; supply one registration directory"
            )
        result[name] = candidates[0] if candidates else None
    return result


def write_registration_settings(config: PipelineConfig, destination: Path) -> Path:
    """Record resolved choices for later stages and beside the final datapackage."""
    resolved = resolve_registration_settings(config)
    portable = {}
    for name, path in resolved.items():
        if path is not None and path.is_absolute() and path.is_relative_to(config.data_root):
            path = path.relative_to(config.data_root)
        portable[name] = path.as_posix() if path is not None else None
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / SETTINGS_FILENAME
    path.write_text(json.dumps({"schema_version": 1, **portable}, indent=2) + "\n")
    return path
