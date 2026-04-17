"""Tests for the datapackage manifest builder."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from aind_ibl_ephys_alignment_preprocessing.manifest import (
    SCHEMA_VERSION,
    DataPackage,
    TransformPaths,
    _build_transforms,
)


def _fake_asset_info(reg_dir: Path) -> SimpleNamespace:
    # AssetInfo is a frozen dataclass with extra fields that _build_transforms
    # does not read; only registration_dir_path is consulted.
    return SimpleNamespace(registration_dir_path=reg_dir)


def _fake_config(template_to_ccf_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(template_to_ccf_dir=template_to_ccf_dir)


def test_schema_version_is_1_1_0():
    """The datapackage schema has been bumped to 1.1.0 (relative transforms)."""
    assert SCHEMA_VERSION == "1.1.0"


def test_transforms_are_relative_to_manifest_root(tmp_path):
    """_build_transforms should emit paths relative to the manifest root."""
    manifest_root = tmp_path / "results" / "mouse42"
    smartspim = tmp_path / "data" / "SmartSPIM_mouse42_123" / "image_atlas_alignment" / "Ex_561_Em_600"
    template_to_ccf = tmp_path / "data" / "spim_template_to_ccf"

    tx = _build_transforms(
        _fake_asset_info(smartspim),  # type: ignore[arg-type]
        _fake_config(template_to_ccf),  # type: ignore[arg-type]
        manifest_root,
    )

    # All four paths must resolve correctly when joined with manifest_root.
    for rel, expected in [
        (tx.image_to_template_affine, smartspim / "ls_to_template_SyN_0GenericAffine.mat"),
        (tx.image_to_template_warp, smartspim / "ls_to_template_SyN_1InverseWarp.nii.gz"),
        (tx.template_to_ccf_affine, template_to_ccf / "syn_0GenericAffine.mat"),
        (tx.template_to_ccf_warp, template_to_ccf / "syn_1InverseWarp.nii.gz"),
    ]:
        assert not Path(rel).is_absolute(), f"expected relative path, got {rel!r}"
        assert (manifest_root / rel).resolve() == expected.resolve()


def test_transforms_traverse_up_for_sibling_assets(tmp_path):
    """When the SmartSPIM asset is a sibling of the mouse root, path uses ``..``."""
    manifest_root = tmp_path / "results" / "mouse42"
    smartspim = tmp_path / "SmartSPIM_mouse42_123" / "image_atlas_alignment" / "Ex_561_Em_600"
    template_to_ccf = tmp_path / "spim_template_to_ccf"

    tx = _build_transforms(
        _fake_asset_info(smartspim),  # type: ignore[arg-type]
        _fake_config(template_to_ccf),  # type: ignore[arg-type]
        manifest_root,
    )

    assert tx.image_to_template_affine.startswith("../"), tx.image_to_template_affine
    assert tx.template_to_ccf_affine.startswith("../"), tx.template_to_ccf_affine


def test_transforms_posix_style_separators(tmp_path):
    """Paths in the manifest are POSIX-style regardless of host OS."""
    manifest_root = tmp_path / "results" / "mouse42"
    smartspim = tmp_path / "SmartSPIM_mouse42_123" / "image_atlas_alignment" / "Ex_561_Em_600"

    tx = _build_transforms(
        _fake_asset_info(smartspim),  # type: ignore[arg-type]
        _fake_config(tmp_path / "spim_template_to_ccf"),  # type: ignore[arg-type]
        manifest_root,
    )

    # No backslashes anywhere.
    for value in tx.model_dump().values():
        assert "\\" not in value, value


def test_transform_paths_fields_shape():
    """TransformPaths drops smartspim_asset/reference_asset; 4 fields remain."""
    assert set(TransformPaths.model_fields) == {
        "image_to_template_affine",
        "image_to_template_warp",
        "template_to_ccf_affine",
        "template_to_ccf_warp",
    }


def test_datapackage_round_trip(tmp_path):
    """Serialize + parse a minimal DataPackage to check schema is loadable."""
    from aind_ibl_ephys_alignment_preprocessing.manifest import (
        CcfSpaceHistology,
        HistologyPaths,
        ImageSpaceHistology,
        load_datapackage,
        write_datapackage,
    )

    dp = DataPackage(
        schema_version=SCHEMA_VERSION,
        mouse_id="mouse42",
        transforms=TransformPaths(
            image_to_template_affine="../SmartSPIM/image_atlas_alignment/Ex_561_Em_600/a.mat",
            image_to_template_warp="../SmartSPIM/image_atlas_alignment/Ex_561_Em_600/w.nii.gz",
            template_to_ccf_affine="../spim_template_to_ccf/a.mat",
            template_to_ccf_warp="../spim_template_to_ccf/w.nii.gz",
        ),
        histology=HistologyPaths(
            image_space=ImageSpaceHistology(
                registration="image_space_histology/histology_registration.nrrd",
                registration_pipeline="image_space_histology/histology_registration_pipeline.nrrd",
                ccf_template="image_space_histology/ccf_in_mouse.nrrd",
                labels="image_space_histology/labels_in_mouse.nrrd",
            ),
            ccf_space=CcfSpaceHistology(
                registration="ccf_space_histology/histology_registration.nrrd",
            ),
        ),
        probes={},
    )
    path = write_datapackage(dp, tmp_path)
    loaded = load_datapackage(path)
    assert loaded == dp
