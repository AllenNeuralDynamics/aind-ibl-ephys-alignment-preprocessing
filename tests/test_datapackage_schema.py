"""Tests for the exported datapackage JSON Schema contract."""

from __future__ import annotations

import json
from typing import Any

import pytest
from pydantic import ValidationError

from aind_ibl_ephys_alignment_preprocessing.datapackage import SCHEMA_VERSION, DataPackage, PathReference
from aind_ibl_ephys_alignment_preprocessing.datapackage_schema import (
    JSON_SCHEMA_DIALECT,
    SCHEMA_ID,
    bundled_datapackage_schema_resource,
    datapackage_json_schema,
    datapackage_schema_json,
)


def _defs(schema: dict[str, Any], name: str) -> dict[str, Any]:
    return schema["$defs"][name]


def _allows_null(field_schema: dict[str, Any]) -> bool:
    return any(option.get("type") == "null" for option in field_schema.get("anyOf", []))


def test_committed_datapackage_schema_matches_generated_schema():
    """The bundled JSON Schema artifact must stay in lockstep with the Pydantic model."""
    assert bundled_datapackage_schema_resource().read_text() == datapackage_schema_json()


def test_datapackage_schema_has_stable_identity_and_version():
    """The schema artifact names the contract independently from the Python package version."""
    schema = datapackage_json_schema()

    assert schema["$schema"] == JSON_SCHEMA_DIALECT
    assert schema["$id"] == SCHEMA_ID
    assert schema["properties"]["schema_version"]["const"] == SCHEMA_VERSION
    assert "schema_version" in schema["required"]


def test_datapackage_schema_encodes_nullable_qc_outputs():
    """QC-only CCF outputs can be null while GUI-required image-space picks remain required."""
    schema = datapackage_json_schema()

    xyz_picks = _defs(schema, "XyzPicks")
    assert "image_space" in xyz_picks["required"]
    assert _allows_null(xyz_picks["properties"]["ccf"])

    histology = _defs(schema, "HistologyPaths")
    assert _allows_null(histology["properties"]["ccf_space"])


def test_datapackage_schema_encodes_path_reference_shape():
    """A 3.x path reference explicitly carries both asset identity and a non-empty path."""
    path_reference = _defs(datapackage_json_schema(), "PathReference")

    assert set(path_reference["required"]) == {"asset", "path"}
    assert path_reference["properties"]["asset"]["anyOf"] == [{"type": "string"}, {"type": "null"}]
    assert path_reference["properties"]["path"]["minLength"] == 1

    assert PathReference.model_validate({"asset": None, "path": "local/file.txt"}) == PathReference(
        asset=None,
        path="local/file.txt",
    )
    with pytest.raises(ValidationError):
        PathReference.model_validate({"path": "missing-asset-key.txt"})
    with pytest.raises(ValidationError):
        PathReference.model_validate({"asset": None, "path": ""})


def test_bundled_schema_is_valid_json():
    """The committed schema artifact is plain JSON that non-Python consumers can load."""
    loaded = json.loads(bundled_datapackage_schema_resource().read_text())
    assert loaded == datapackage_json_schema()


def test_datapackage_rejects_unsupported_schema_version():
    """The Pydantic contract now matches the schema's exact 3.1.0 version."""
    with pytest.raises(ValidationError):
        DataPackage.model_validate(
            {
                "schema_version": "3.0.0",
                "mouse_id": "m1",
                "platform": "local",
                "external_assets": {},
                "transforms": {},
                "histology": {},
                "probes": {},
            }
        )
