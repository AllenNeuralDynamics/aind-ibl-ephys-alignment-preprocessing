"""The manifest's ``registration_asset`` column pins the registration directory.

Three campaign brains were CCF-registered on the wrong channel and re-registered
standalone. The column lets a run consume that without overwriting
``image_atlas_alignment/`` inside a published asset.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from aind_ibl_ephys_alignment_preprocessing.discovery import (
    REGISTRATION_TRANSFORMS,
    manifest_registration_override,
    resolve_registration_dir,
)
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, PipelineConfig

MANIFEST_HEADER = "mouseid,probe_file,sorted_recording,probe_id,probe_name"


def _config(tmp_path: Path, manifest_body: str | None) -> PipelineConfig:
    data_root = tmp_path / "data"
    data_root.mkdir(exist_ok=True)
    manifest = data_root / "manifest.csv"
    if manifest_body is not None:
        manifest.write_text(manifest_body)
    return PipelineConfig(
        data_root=data_root,
        results_root=tmp_path / "results",
        neuroglancer_file=Path("ng.json"),
        manifest_csv=Path("manifest.csv"),
    )


def _registration_dir(root: Path, name: str, *, complete: bool = True) -> Path:
    reg = root / name
    reg.mkdir(parents=True, exist_ok=True)
    for transform in REGISTRATION_TRANSFORMS if complete else REGISTRATION_TRANSFORMS[:1]:
        (reg / transform).write_text("x")
    return reg


# --- ManifestRow ----------------------------------------------------------


def test_manifest_row_defaults_registration_asset_to_none():
    """The eight correctly-registered brains have no such column."""
    import pandas as pd

    # name= mirrors iterrows(), which is how from_series is always called.
    row = ManifestRow.from_series(
        pd.Series(
            {"mouseid": "750108", "probe_file": "p", "sorted_recording": "s", "probe_id": "a", "probe_name": "b"},
            name=0,
        )
    )

    assert row.registration_asset is None


def test_manifest_row_parses_registration_asset():
    import pandas as pd

    row = ManifestRow.from_series(
        pd.Series(
            {
                "mouseid": "750108",
                "probe_file": "p",
                "sorted_recording": "s",
                "probe_id": "a",
                "probe_name": "b",
                "registration_asset": "SmartSPIM_750108_reg/ccf_Ex_639_Em_667",
            },
            name=0,
        )
    )

    assert row.registration_asset == Path("SmartSPIM_750108_reg/ccf_Ex_639_Em_667")


# --- reading the column off the manifest ----------------------------------


def test_no_column_means_no_override(tmp_path):
    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750108,p,s,a,b\n")

    assert manifest_registration_override(config) is None


def test_blank_column_means_no_override(tmp_path):
    """An empty cell must take the untouched path, not resolve to data_root."""
    config = _config(tmp_path, f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,\n")

    assert manifest_registration_override(config) is None


def test_column_is_read_per_brain_not_per_row(tmp_path):
    """Replicated across rows like mouseid; one value is not a conflict."""
    body = f"{MANIFEST_HEADER},registration_asset\n" + "".join(
        f"750108,p{i},s,a,b,reg/ccf_Ex_639_Em_667\n" for i in range(3)
    )
    config = _config(tmp_path, body)

    assert manifest_registration_override(config) == Path("reg/ccf_Ex_639_Em_667")


def test_disagreeing_rows_raise(tmp_path):
    """Two registrations for one brain is a manifest error, not a merge."""
    body = f"{MANIFEST_HEADER},registration_asset\n750108,p1,s,a,b,reg/ccf_A\n750108,p2,s,a,b,reg/ccf_B\n"
    config = _config(tmp_path, body)

    with pytest.raises(ValueError, match="must name one directory"):
        manifest_registration_override(config)


def test_missing_manifest_is_not_an_override(tmp_path):
    """The ephys stage runs with a placeholder manifest path."""
    config = _config(tmp_path, None)

    assert manifest_registration_override(config) is None


# --- resolving the directory ----------------------------------------------


def test_without_override_returns_in_asset_dir_and_no_channel(tmp_path):
    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750108,p,s,a,b\n")
    asset = config.data_root / "SmartSPIM_750108"

    reg_dir, channel = resolve_registration_dir(config, asset)

    assert reg_dir == asset / "image_atlas_alignment"
    # None means "keep whatever processing.json nominated".
    assert channel is None


def test_override_returns_dir_and_channel_from_its_name(tmp_path):
    config = _config(tmp_path, f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_639_Em_667\n")
    _registration_dir(config.data_root, "reg/ccf_Ex_639_Em_667")
    asset = config.data_root / "SmartSPIM_750108"

    reg_dir, channel = resolve_registration_dir(config, asset)

    assert reg_dir == config.data_root / "reg" / "ccf_Ex_639_Em_667"
    assert channel == "Ex_639_Em_667"


def test_override_accepts_the_pipelines_own_layout(tmp_path):
    """`image_atlas_alignment/<channel>` has no `ccf_` prefix to strip."""
    body = f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/image_atlas_alignment/Ex_561_Em_600\n"
    config = _config(tmp_path, body)
    _registration_dir(config.data_root, "reg/image_atlas_alignment/Ex_561_Em_600")

    _, channel = resolve_registration_dir(config, config.data_root / "SmartSPIM_750108")

    assert channel == "Ex_561_Em_600"


def test_missing_override_directory_raises_naming_the_asset(tmp_path):
    config = _config(tmp_path, f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_639_Em_667\n")

    with pytest.raises(FileNotFoundError, match="reg/ccf_Ex_639_Em_667"):
        resolve_registration_dir(config, config.data_root / "SmartSPIM_750108")


def test_incomplete_override_directory_raises_naming_the_transform(tmp_path):
    """A directory that exists but lacks a transform must not pass silently."""
    config = _config(tmp_path, f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_639_Em_667\n")
    _registration_dir(config.data_root, "reg/ccf_Ex_639_Em_667", complete=False)

    with pytest.raises(FileNotFoundError, match="ls_to_template_SyN_1InverseWarp.nii.gz"):
        resolve_registration_dir(config, config.data_root / "SmartSPIM_750108")


# --- datapackage provenance ------------------------------------------------


def _asset_info(asset_path: Path, reg_dir: Path):
    from aind_ibl_ephys_alignment_preprocessing.types import AssetInfo, PipelineRegistrationInfo, ZarrPaths

    return AssetInfo(
        asset_path=asset_path,
        asset_uri="s3://bucket/SmartSPIM_750108",
        zarr_volumes=ZarrPaths(registration="z.zarr", additional=[], metadata={}, processing={}),
        pipeline_registration_chains=PipelineRegistrationInfo(
            pt_tx_str=[], pt_tx_inverted=[], img_tx_str=[], img_tx_inverted=[]
        ),
        registration_dir_path=reg_dir,
    )


def test_in_asset_registration_stays_on_the_smartspim_asset(tmp_path):
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_external_assets, _build_transforms

    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750108,p,s,a,b\n")
    asset = config.data_root / "SmartSPIM_750108"
    info = _asset_info(asset, asset / "image_atlas_alignment" / "Ex_561_Em_600")

    assets = _build_external_assets(info, config)
    transforms = _build_transforms(info, config, tmp_path)

    assert "registration" not in assets
    assert transforms.image_to_template_affine.asset == "smartspim"
    assert transforms.image_to_template_affine.path == (
        "image_atlas_alignment/Ex_561_Em_600/ls_to_template_SyN_0GenericAffine.mat"
    )


def test_override_records_which_registration_was_used(tmp_path):
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_external_assets, _build_transforms

    config = _config(tmp_path, f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_639_Em_667\n")
    asset = config.data_root / "SmartSPIM_750108"
    info = _asset_info(asset, config.data_root / "reg" / "ccf_Ex_639_Em_667")

    assets = _build_external_assets(info, config)
    transforms = _build_transforms(info, config, tmp_path)

    assert assets["registration"].role == "registration_override"
    assert assets["registration"].name == "reg"
    # The refs must point at the override asset, not the stitched one -- this is
    # the whole point of the column over overwriting in place.
    assert transforms.image_to_template_affine.asset == "registration"
    assert transforms.image_to_template_affine.path == "ccf_Ex_639_Em_667/ls_to_template_SyN_0GenericAffine.mat"
    assert transforms.image_to_template_warp.asset == "registration"
    # template->CCF is unrelated to the override and must not move.
    assert transforms.template_to_ccf_affine.asset == "spim_template_to_ccf"


# --- the manifest CSV contract ---------------------------------------------


def test_registration_asset_is_optional():
    """A column added later must never invalidate manifests written before it."""
    from aind_ibl_ephys_alignment_preprocessing.types import OPTIONAL_MANIFEST_COLUMNS

    assert "registration_asset" in {c.name for c in OPTIONAL_MANIFEST_COLUMNS}


def test_manifest_without_optional_columns_validates(tmp_path):
    """The eight correctly-registered brains carry none of the optional columns."""
    from aind_ibl_ephys_alignment_preprocessing.types import OPTIONAL_MANIFEST_COLUMNS
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    config = _config(tmp_path, "mouseid,sorted_recording,probe_file,probe_id,probe_name\n750108,s,p,a,b\n")
    validator = PipelineValidator(config, skip_resource_checks=True)
    validator.validate_manifest_structure()

    failures = [r for r in validator.results if not r.passed and r.severity == "error"]
    assert failures == [], [r.message for r in failures]
    # and none of the optional names is ever reported as missing
    missing = " ".join(r.message for r in validator.results if "Missing required" in r.message)
    for column in OPTIONAL_MANIFEST_COLUMNS:
        assert column.name not in missing


def test_missing_required_column_still_fails(tmp_path):
    """Optionality must not have loosened the required set."""
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    config = _config(tmp_path, "mouseid,sorted_recording,probe_id,probe_name\n750108,s,a,b\n")
    validator = PipelineValidator(config, skip_resource_checks=True)
    validator.validate_manifest_structure()

    assert any(not r.passed and "probe_file" in r.message for r in validator.results)


def test_required_columns_accept_their_legacy_aliases(tmp_path):
    """probe_id/probe_name are the old spelling and must keep working."""
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    config = _config(tmp_path, "mouseid,sorted_recording,probe_file,probe_id,probe_name\n750108,s,p,a,b\n")
    validator = PipelineValidator(config, skip_resource_checks=True)
    validator.validate_manifest_structure()

    assert not [r for r in validator.results if not r.passed and r.severity == "error"]


def test_optional_columns_in_use_are_reported(tmp_path):
    from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

    body = f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_639_Em_667\n"
    config = _config(tmp_path, body)
    validator = PipelineValidator(config, skip_resource_checks=True)
    validator.validate_manifest_structure()

    reported = [r.message for r in validator.results if r.item == "optional_columns"]
    assert reported and "registration_asset" in reported[0]


def test_every_manifest_row_field_is_in_the_contract():
    """The contract must not drift from what the parser actually reads."""
    from aind_ibl_ephys_alignment_preprocessing.types import MANIFEST_COLUMNS

    declared = {n for c in MANIFEST_COLUMNS for n in c.accepted_names}
    # Fields ManifestRow derives rather than reads straight from a column.
    derived = {"probe_id", "probe_name", "probe_shank", "row_index"}
    parsed = {f.name for f in __import__("dataclasses").fields(ManifestRow)} - derived

    assert parsed <= declared, parsed - declared


# --- find_asset_info end to end --------------------------------------------
#
# The override's whole job happens inside find_asset_info, which was previously
# only smoke-tested. The zarr-utils calls are patched out so the wiring itself
# -- which directory, which channel, which URI reaches the pipeline chain -- is
# what gets exercised.

DEFAULT_CHANNEL = "Ex_561_Em_600"
OVERRIDE_CHANNEL = "Ex_639_Em_667"


def _stitched_asset(data_root: Path, *channels: str) -> Path:
    asset = data_root / "SmartSPIM_750108"
    omezarr = asset / "image_tile_fusing" / "OMEZarr"
    for channel in channels:
        (omezarr / f"{channel}.zarr").mkdir(parents=True, exist_ok=True)
    return asset


def _patch_zarr_utils(monkeypatch, asset: Path, default_channel: str) -> dict:
    """Patch discovery's zarr-utils calls; record what the pipeline chain saw."""
    from aind_ibl_ephys_alignment_preprocessing import discovery

    default_zarr = asset / "image_tile_fusing" / "OMEZarr" / f"{default_channel}.zarr"
    seen: dict = {}

    monkeypatch.setattr(discovery, "get_json", lambda *a, **k: {})
    monkeypatch.setattr(discovery, "get_image_sources", lambda *a, **k: {"layer": str(default_zarr)})
    monkeypatch.setattr(discovery, "as_pathlike", lambda uri: (None, None, str(uri)))
    monkeypatch.setattr(discovery, "_asset_from_zarr_pathlike", lambda p: asset.name)
    monkeypatch.setattr(
        discovery,
        "alignment_zarr_uri_and_metadata_from_zarr_or_asset_pathlike",
        lambda **k: (default_zarr.as_posix(), {"meta": 1}, {"processing": 1}),
    )

    def _chain(uri, processing, **kwargs):
        seen["uri"] = uri
        seen["processing"] = processing
        return [], [], [], []

    monkeypatch.setattr(discovery, "pipeline_transforms_local_paths", _chain)
    return seen


def test_find_asset_info_without_override_is_unchanged(tmp_path, monkeypatch):
    from aind_ibl_ephys_alignment_preprocessing.discovery import find_asset_info

    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750108,p,s,a,b\n")
    asset = _stitched_asset(config.data_root, DEFAULT_CHANNEL, OVERRIDE_CHANNEL)
    seen = _patch_zarr_utils(monkeypatch, asset, DEFAULT_CHANNEL)

    info = find_asset_info(config)

    assert info.registration_dir_path == asset / "image_atlas_alignment" / DEFAULT_CHANNEL
    assert info.zarr_volumes.registration.endswith(f"{DEFAULT_CHANNEL}.zarr")
    assert [Path(p).stem for p in info.zarr_volumes.additional] == [OVERRIDE_CHANNEL]
    assert seen["uri"].endswith(f"{DEFAULT_CHANNEL}.zarr")


def test_find_asset_info_override_moves_dir_channel_and_chain(tmp_path, monkeypatch):
    """All four things the override must change, in one assertion block."""
    from aind_ibl_ephys_alignment_preprocessing.discovery import find_asset_info

    body = f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_{OVERRIDE_CHANNEL}\n"
    config = _config(tmp_path, body)
    asset = _stitched_asset(config.data_root, DEFAULT_CHANNEL, OVERRIDE_CHANNEL)
    _registration_dir(config.data_root, f"reg/ccf_{OVERRIDE_CHANNEL}")
    seen = _patch_zarr_utils(monkeypatch, asset, DEFAULT_CHANNEL)

    info = find_asset_info(config)

    # 1. the transforms come from the override asset
    assert info.registration_dir_path == config.data_root / "reg" / f"ccf_{OVERRIDE_CHANNEL}"
    # 2. the histology volume follows the channel the transforms were computed from
    assert info.zarr_volumes.registration.endswith(f"{OVERRIDE_CHANNEL}.zarr")
    # 3. ...and that channel is no longer counted as an additional one
    assert [Path(p).stem for p in info.zarr_volumes.additional] == [DEFAULT_CHANNEL]
    # 4. the per-channel stitching chain follows the channel too...
    assert seen["uri"].endswith(f"{OVERRIDE_CHANNEL}.zarr")
    # ...but processing.json still comes from the stitched asset, not the override
    assert seen["processing"] == {"processing": 1}
    # the zarr itself never moves out of the stitched asset
    assert Path(info.zarr_volumes.registration).is_relative_to(asset)


def test_find_asset_info_rejects_a_channel_the_brain_does_not_have(tmp_path, monkeypatch):
    """Third failure mode: the column names a channel absent from the zarr set."""
    from aind_ibl_ephys_alignment_preprocessing.discovery import find_asset_info

    body = f"{MANIFEST_HEADER},registration_asset\n750108,p,s,a,b,reg/ccf_Ex_999_Em_999\n"
    config = _config(tmp_path, body)
    asset = _stitched_asset(config.data_root, DEFAULT_CHANNEL, OVERRIDE_CHANNEL)
    _registration_dir(config.data_root, "reg/ccf_Ex_999_Em_999")
    _patch_zarr_utils(monkeypatch, asset, DEFAULT_CHANNEL)

    with pytest.raises(FileNotFoundError) as excinfo:
        find_asset_info(config)

    message = str(excinfo.value)
    assert "Ex_999_Em_999" in message
    # and it says what the brain does have, so the fix is obvious
    assert DEFAULT_CHANNEL in message and OVERRIDE_CHANNEL in message
