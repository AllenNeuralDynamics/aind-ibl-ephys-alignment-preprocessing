"""The manifest's ``template_to_ccf_asset`` column replaces the SmartSPIM template's CCF registration.

A brain registered to another template (an iDISCO template, say) needs that template's template->CCF chain in place
of the one ``processing.json`` names.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from aind_ibl_ephys_alignment_preprocessing._constants import TEMPLATE_TO_CCF_TRANSFORMS
from aind_ibl_ephys_alignment_preprocessing.discovery import resolve_template_to_ccf_dir
from aind_ibl_ephys_alignment_preprocessing.types import (
    OPTIONAL_MANIFEST_COLUMNS,
    AssetInfo,
    ManifestRow,
    PipelineConfig,
    PipelineRegistrationInfo,
    ZarrPaths,
)

MANIFEST_HEADER = "mouseid,probe_file,sorted_recording,probe_id,probe_name"

POINT_CHAIN = [
    "/cache/smartspim/ls_to_template_SyN_0GenericAffine.mat",
    "/cache/smartspim/ls_to_template_SyN_1InverseWarp.nii.gz",
    "/cache/template/spim_template_to_ccf_syn_0GenericAffine_25.mat",
    "/cache/template/spim_template_to_ccf_syn_1InverseWarp_25.nii.gz",
]
IMAGE_CHAIN = [
    "/cache/template/spim_template_to_ccf_syn_1Warp_25.nii.gz",
    "/cache/template/spim_template_to_ccf_syn_0GenericAffine_25.mat",
    "/cache/smartspim/ls_to_template_SyN_1Warp.nii.gz",
    "/cache/smartspim/ls_to_template_SyN_0GenericAffine.mat",
]


def _config(tmp_path: Path, manifest_body: str) -> PipelineConfig:
    data_root = tmp_path / "data"
    data_root.mkdir(exist_ok=True)
    (data_root / "manifest.csv").write_text(manifest_body)
    return PipelineConfig(
        data_root=data_root,
        results_root=tmp_path / "results",
        neuroglancer_file=Path("ng.json"),
        manifest_csv=Path("manifest.csv"),
    )


def _template_dir(root: Path, *names: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for name in names:
        (root / name).write_text("x")
    return root


def _asset_info(tmp_path: Path, reg_dir: Path | None = None, template_dir: Path | None = None) -> AssetInfo:
    asset = tmp_path / "data" / "SmartSPIM_750107"
    return AssetInfo(
        asset_path=asset,
        asset_uri=None,
        zarr_volumes=ZarrPaths(registration="z.zarr", additional=[], metadata={}, processing={}),
        pipeline_registration_chains=PipelineRegistrationInfo(
            pt_tx_str=list(POINT_CHAIN),
            pt_tx_inverted=[True, False, True, False],
            img_tx_str=list(IMAGE_CHAIN),
            img_tx_inverted=[False] * 4,
        ),
        registration_dir_path=reg_dir or asset / "image_atlas_alignment" / "Ex_639_Em_667",
        template_to_ccf_dir=template_dir,
    )


def test_column_is_optional_and_parsed():
    import pandas as pd

    assert "template_to_ccf_asset" in {c.name for c in OPTIONAL_MANIFEST_COLUMNS}
    fields = {"mouseid": "750107", "probe_file": "p", "sorted_recording": "s", "probe_id": "a", "probe_name": "b"}
    assert ManifestRow.from_series(pd.Series(fields, name=0)).template_to_ccf_asset is None
    row = ManifestRow.from_series(pd.Series({**fields, "template_to_ccf_asset": "idisco/template_to_ccf"}, name=0))
    assert row.template_to_ccf_asset == Path("idisco/template_to_ccf")


def test_resolve_requires_both_point_transforms(tmp_path):
    config = _config(tmp_path, f"{MANIFEST_HEADER},template_to_ccf_asset\n750107,p,s,a,b,idisco/template_to_ccf\n")
    _template_dir(config.data_root / "idisco" / "template_to_ccf", TEMPLATE_TO_CCF_TRANSFORMS[0])
    with pytest.raises(FileNotFoundError, match="syn_1InverseWarp"):
        resolve_template_to_ccf_dir(config)
    _template_dir(config.data_root / "idisco" / "template_to_ccf", *TEMPLATE_TO_CCF_TRANSFORMS)
    assert resolve_template_to_ccf_dir(config) == config.data_root / "idisco" / "template_to_ccf"


def test_resolve_without_the_column_keeps_the_pipeline_chain(tmp_path):
    assert resolve_template_to_ccf_dir(_config(tmp_path, f"{MANIFEST_HEADER}\n750107,p,s,a,b\n")) is None


def test_point_chain_replaces_only_the_template_half(tmp_path):
    tmpl = tmp_path / "data" / "idisco" / "template_to_ccf"
    paths, inverted = _asset_info(tmp_path, template_dir=tmpl).point_chain()
    assert paths == POINT_CHAIN[:2] + [str(tmpl / "syn_0GenericAffine.mat"), str(tmpl / "syn_1InverseWarp.nii.gz")]
    assert inverted == [True, False, True, False]


def test_substitution_sees_through_the_s3_cache_prefix(tmp_path):
    """The template->CCF half is always downloaded, so the chain names cached copies."""
    cached = "/tmp/s3-ants-cache/aind-open-data__SmartSPIM-template_2024-05-16_11-26-14__fa7e696d0972__"
    info = _asset_info(tmp_path, template_dir=tmp_path / "t2c")
    chains = info.pipeline_registration_chains
    chains.pt_tx_str[2:] = [cached + Path(p).name for p in POINT_CHAIN[2:]]
    paths, _ = info.point_chain()
    assert paths[2:] == [
        str(tmp_path / "t2c" / "syn_0GenericAffine.mat"),
        str(tmp_path / "t2c" / "syn_1InverseWarp.nii.gz"),
    ]


def test_both_overrides_compose(tmp_path):
    reg = tmp_path / "data" / "reg" / "750107"
    tmpl = tmp_path / "data" / "idisco" / "template_to_ccf"
    paths, _ = _asset_info(tmp_path, reg, tmpl).point_chain()
    assert [Path(p).parent for p in paths] == [reg, reg, tmpl, tmpl]


def test_image_chain_substitutes_forward_warps(tmp_path):
    reg = _template_dir(
        tmp_path / "data" / "reg" / "750107", "ls_to_template_SyN_1Warp.nii.gz", "ls_to_template_SyN_0GenericAffine.mat"
    )
    tmpl = _template_dir(tmp_path / "data" / "idisco" / "t2c", "syn_1Warp.nii.gz", "syn_0GenericAffine.mat")
    paths, _ = _asset_info(tmp_path, reg, tmpl).image_chain()
    assert paths == [
        str(tmpl / "syn_1Warp.nii.gz"),
        str(tmpl / "syn_0GenericAffine.mat"),
        str(reg / "ls_to_template_SyN_1Warp.nii.gz"),
        str(reg / "ls_to_template_SyN_0GenericAffine.mat"),
    ]


def test_image_chain_reports_a_missing_forward_warp(tmp_path):
    tmpl = _template_dir(tmp_path / "data" / "idisco" / "t2c", *TEMPLATE_TO_CCF_TRANSFORMS)
    with pytest.raises(FileNotFoundError, match="syn_1Warp"):
        _asset_info(tmp_path, template_dir=tmpl).image_chain()


def test_datapackage_records_the_override_template(tmp_path):
    from aind_ibl_ephys_alignment_preprocessing.datapackage import (
        _build_external_assets,
        _build_transforms,
        producer_asset_overrides,
    )

    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750107,p,s,a,b\n")
    tmpl = config.data_root / "idisco_template_2026-10" / "template_to_ccf"
    info = _asset_info(tmp_path, template_dir=tmpl)

    transforms = _build_transforms(info, config, tmp_path)
    assets = _build_external_assets(info, config)

    assert transforms.template_to_ccf_affine.asset == "template_to_ccf"
    assert transforms.template_to_ccf_affine.path == "template_to_ccf/syn_0GenericAffine.mat"
    assert transforms.template_to_ccf_warp.path == "template_to_ccf/syn_1InverseWarp.nii.gz"
    assert assets["template_to_ccf"].name == "idisco_template_2026-10"
    assert assets["template_to_ccf"].role == "template_to_ccf_override"
    assert producer_asset_overrides(info, config)["template_to_ccf"] == config.data_root / "idisco_template_2026-10"


def test_datapackage_without_override_is_unchanged(tmp_path):
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_external_assets, _build_transforms

    config = _config(tmp_path, f"{MANIFEST_HEADER}\n750107,p,s,a,b\n")
    info = _asset_info(tmp_path)
    assert _build_transforms(info, config, tmp_path).template_to_ccf_affine.asset == "spim_template_to_ccf"
    assert "template_to_ccf" not in _build_external_assets(info, config)
