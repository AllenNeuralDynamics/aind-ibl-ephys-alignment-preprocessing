"""JSON Schema export helpers for the datapackage contract."""

from __future__ import annotations

import json
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any

from aind_ibl_ephys_alignment_preprocessing.datapackage import SCHEMA_VERSION, DataPackage

SCHEMA_NAME = "aind-ibl-ephys-alignment-datapackage"
SCHEMA_ID = f"https://schemas.allenneuraldynamics.org/{SCHEMA_NAME}/{SCHEMA_VERSION}/datapackage.schema.json"
JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
SCHEMA_RELATIVE_PARTS = ("schemas", SCHEMA_NAME, SCHEMA_VERSION, "datapackage.schema.json")


def datapackage_json_schema() -> dict[str, Any]:
    """Return the JSON Schema for the current datapackage contract."""
    schema = DataPackage.model_json_schema()
    schema["$schema"] = JSON_SCHEMA_DIALECT
    schema["$id"] = SCHEMA_ID
    return schema


def datapackage_schema_json() -> str:
    """Return a stable JSON representation of the datapackage schema."""
    return json.dumps(datapackage_json_schema(), indent=2, sort_keys=True) + "\n"


def source_tree_datapackage_schema_path() -> Path:
    """Return the source-tree path for the bundled datapackage schema artifact."""
    return Path(__file__).resolve().parent.joinpath(*SCHEMA_RELATIVE_PARTS)


def bundled_datapackage_schema_resource() -> Traversable:
    """Return the installed package resource for the bundled datapackage schema."""
    return files("aind_ibl_ephys_alignment_preprocessing").joinpath(*SCHEMA_RELATIVE_PARTS)
