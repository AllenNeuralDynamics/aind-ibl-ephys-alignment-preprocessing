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


def test_schema_version_is_2_0_0():
    """Schema is 2.0.0 — probes are nested by ``recording_id`` then ``probe_name``."""
    assert SCHEMA_VERSION == "2.0.0"


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


def _fake_outputs(tmp_path: Path) -> SimpleNamespace:
    # _build_probes only reads tracks_root.parent (via row.gui_folder(outputs))
    # to build per-recording GUI folders. Make tracks_root a child of the
    # mouse root so gui_folder returns ``<mouse_root>/<recording_id>/<probe_name>``.
    return SimpleNamespace(tracks_root=tmp_path / "_unused_tracks")


def _fake_pipeline_config(skip_ephys: bool = False) -> SimpleNamespace:
    return SimpleNamespace(skip_ephys=skip_ephys)


def _fake_row(*, probe_id: str, probe_name: str, recording_id: str, probe_shank=None):
    """Build a stand-in for ManifestRow exposing only what _build_probes reads."""
    sorted_recording = f"{recording_id}_sorted_2025-01-01"
    return SimpleNamespace(
        probe_id=probe_id,
        probe_name=probe_name,
        recording_id=recording_id,
        sorted_recording=sorted_recording,
        probe_shank=probe_shank,
        gui_folder=lambda outputs, _rec=recording_id, _name=probe_name: (
            outputs.tracks_root.parent / _rec / _name
        ),
    )


def test_probes_nested_by_recording_then_name(tmp_path):
    """Probes are grouped under recording_id, then probe_name."""
    from aind_ibl_ephys_alignment_preprocessing.manifest import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [
        _fake_row(probe_id="pid-1", probe_name="probeA", recording_id="rec1"),
        _fake_row(probe_id="pid-2", probe_name="probeB", recording_id="rec1"),
        _fake_row(probe_id="pid-3", probe_name="probeC", recording_id="rec2"),
    ]
    results = [
        ProcessResult(probe_id=r.probe_id, recording_id=r.recording_id, wrote_files=True)
        for r in rows
    ]

    probes = _build_probes(
        rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config()
    )

    assert set(probes.keys()) == {"rec1", "rec2"}
    assert set(probes["rec1"].keys()) == {"probeA", "probeB"}
    assert set(probes["rec2"].keys()) == {"probeC"}


def test_same_probe_name_across_recordings_kept_distinct(tmp_path):
    """Same probe_name in two recordings stays distinct (the bug fixed in 2.0.0)."""
    from aind_ibl_ephys_alignment_preprocessing.manifest import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    # Same probe_name "45883-1" appears in rec1 and rec2 (re-inserted probe).
    rows = [
        _fake_row(probe_id="pid-1a", probe_name="45883-1", recording_id="rec1"),
        _fake_row(probe_id="pid-1b", probe_name="45883-1", recording_id="rec2"),
    ]
    results = [
        ProcessResult(probe_id=r.probe_id, recording_id=r.recording_id, wrote_files=True)
        for r in rows
    ]

    probes = _build_probes(
        rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config()
    )

    assert probes["rec1"]["45883-1"].probe_id == "pid-1a"
    assert probes["rec2"]["45883-1"].probe_id == "pid-1b"
    # Each entry has exactly one shank — they must NOT be merged into a
    # 2-shank probe under one recording (the pre-2.0.0 silent-merge bug).
    assert probes["rec1"]["45883-1"].num_shanks == 1
    assert probes["rec2"]["45883-1"].num_shanks == 1


def test_multi_shank_collapses_into_one_probe(tmp_path):
    """Rows differing only by probe_shank collapse into one ProbeEntry."""
    from aind_ibl_ephys_alignment_preprocessing.manifest import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [
        _fake_row(probe_id="pid", probe_name="probeM", recording_id="rec1", probe_shank=0),
        _fake_row(probe_id="pid", probe_name="probeM", recording_id="rec1", probe_shank=1),
    ]
    results = [
        ProcessResult(probe_id="pid", recording_id="rec1", wrote_files=True),
    ]

    probes = _build_probes(
        rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config()
    )

    assert list(probes["rec1"].keys()) == ["probeM"]
    entry = probes["rec1"]["probeM"]
    assert entry.num_shanks == 2
    assert {pk.shank for pk in entry.xyz_picks} == {1, 2}  # 1-based in datapackage
