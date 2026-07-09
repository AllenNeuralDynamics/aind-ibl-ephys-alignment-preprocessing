"""Tests for the datapackage builder."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from aind_ibl_ephys_alignment_preprocessing.datapackage import (
    SCHEMA_VERSION,
    DataPackage,
    ExternalAsset,
    PathReference,
    TransformPaths,
    _build_transforms,
)


def _fake_asset_info(reg_dir: Path) -> SimpleNamespace:
    # AssetInfo is a frozen dataclass with extra fields that _build_transforms
    # does not read; only asset_path/registration_dir_path are consulted.
    return SimpleNamespace(
        asset_path=reg_dir.parents[1],
        asset_uri=None,
        registration_dir_path=reg_dir,
    )


def _fake_config(template_to_ccf_dir: Path) -> SimpleNamespace:
    return SimpleNamespace(template_to_ccf_dir=template_to_ccf_dir, data_root=template_to_ccf_dir.parent)


def _local(path: str) -> PathReference:
    return PathReference(asset=None, path=path)


def _external(asset: str, path: str) -> PathReference:
    return PathReference(asset=asset, path=path)


def test_schema_version_is_3_0_0():
    """Schema is 3.0.0 — paths are reference objects, not locations."""
    assert SCHEMA_VERSION == "3.0.0"


def test_transforms_are_external_asset_references(tmp_path):
    """_build_transforms should emit paths inside external assets."""
    manifest_root = tmp_path / "results" / "mouse42"
    smartspim = tmp_path / "data" / "SmartSPIM_mouse42_123" / "image_atlas_alignment" / "Ex_561_Em_600"
    template_to_ccf = tmp_path / "data" / "spim_template_to_ccf"

    tx = _build_transforms(
        _fake_asset_info(smartspim),  # type: ignore[arg-type]
        _fake_config(template_to_ccf),  # type: ignore[arg-type]
        manifest_root,
    )

    # All four refs must resolve correctly through their asset roots.
    roots = {
        "smartspim": smartspim.parents[1],
        "spim_template_to_ccf": template_to_ccf,
    }
    for ref, expected in [
        (tx.image_to_template_affine, smartspim / "ls_to_template_SyN_0GenericAffine.mat"),
        (tx.image_to_template_warp, smartspim / "ls_to_template_SyN_1InverseWarp.nii.gz"),
        (tx.template_to_ccf_affine, template_to_ccf / "syn_0GenericAffine.mat"),
        (tx.template_to_ccf_warp, template_to_ccf / "syn_1InverseWarp.nii.gz"),
    ]:
        assert ref.asset is not None
        assert not Path(ref.path).is_absolute(), f"expected internal asset path, got {ref.path!r}"
        assert ".." not in Path(ref.path).parts
        assert (roots[ref.asset] / ref.path).resolve() == expected.resolve()


def test_transforms_do_not_traverse_to_sibling_assets(tmp_path):
    """Sibling assets are referenced by asset key, never by ``..`` traversal."""
    manifest_root = tmp_path / "results" / "mouse42"
    smartspim = tmp_path / "SmartSPIM_mouse42_123" / "image_atlas_alignment" / "Ex_561_Em_600"
    template_to_ccf = tmp_path / "spim_template_to_ccf"

    tx = _build_transforms(
        _fake_asset_info(smartspim),  # type: ignore[arg-type]
        _fake_config(template_to_ccf),  # type: ignore[arg-type]
        manifest_root,
    )

    assert tx.image_to_template_affine.asset == "smartspim"
    assert tx.template_to_ccf_affine.asset == "spim_template_to_ccf"
    for ref in tx.model_dump().values():
        assert ".." not in Path(ref["path"]).parts


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
        assert "\\" not in value["path"], value


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
    from aind_ibl_ephys_alignment_preprocessing.datapackage import (
        CcfSpaceHistology,
        HistologyPaths,
        ImageSpaceHistology,
        load_datapackage,
        write_datapackage,
    )

    dp = DataPackage(
        schema_version=SCHEMA_VERSION,
        mouse_id="mouse42",
        platform="local",
        external_assets={
            "smartspim": ExternalAsset(role="smartspim_registration", name="SmartSPIM"),
            "spim_template_to_ccf": ExternalAsset(role="template_to_ccf", name="spim_template_to_ccf"),
        },
        transforms=TransformPaths(
            image_to_template_affine=_external("smartspim", "image_atlas_alignment/Ex_561_Em_600/a.mat"),
            image_to_template_warp=_external("smartspim", "image_atlas_alignment/Ex_561_Em_600/w.nii.gz"),
            template_to_ccf_affine=_external("spim_template_to_ccf", "a.mat"),
            template_to_ccf_warp=_external("spim_template_to_ccf", "w.nii.gz"),
        ),
        histology=HistologyPaths(
            image_space=ImageSpaceHistology(
                registration=_local("image_space_histology/histology_registration.nrrd"),
                registration_pipeline=_local("image_space_histology/histology_registration_pipeline.nrrd"),
                ccf_template=_local("image_space_histology/ccf_in_mouse.nrrd"),
                labels=_local("image_space_histology/labels_in_mouse.nrrd"),
            ),
            ccf_space=CcfSpaceHistology(
                registration=_local("ccf_space_histology/histology_registration.nrrd"),
            ),
        ),
        probes={},
    )
    # validate=False: this exercises serialization round-trip, not on-disk
    # path existence (the fake paths above don't exist under tmp_path).
    path = write_datapackage(dp, tmp_path, validate=False)
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
    ephys_collection = probe_name
    return SimpleNamespace(
        probe_id=probe_id,
        probe_name=ephys_collection,
        histology_track_id=probe_id,
        logical_probe=ephys_collection,
        ephys_collection=ephys_collection,
        recording_id=recording_id,
        sorted_recording=sorted_recording,
        probe_shank=probe_shank,
        histology_shank=probe_shank,
        ephys_shank=probe_shank,
        gui_folder=lambda outputs, _rec=recording_id, _name=ephys_collection: outputs.tracks_root.parent / _rec / _name,
    )


def _fake_row_explicit(
    *,
    histology_track_id: str,
    logical_probe: str,
    ephys_collection: str,
    recording_id: str,
    histology_shank=None,
    ephys_shank=None,
):
    """Build a stand-in using the explicit manifest v2.1 vocabulary."""
    sorted_recording = f"{recording_id}_sorted_2025-01-01"
    return SimpleNamespace(
        probe_id=histology_track_id,
        probe_name=ephys_collection,
        histology_track_id=histology_track_id,
        logical_probe=logical_probe,
        ephys_collection=ephys_collection,
        recording_id=recording_id,
        sorted_recording=sorted_recording,
        probe_shank=None,
        histology_shank=histology_shank,
        ephys_shank=ephys_shank,
        gui_folder=lambda outputs, _rec=recording_id, _name=ephys_collection: outputs.tracks_root.parent / _rec / _name,
    )


def test_probes_nested_by_recording_then_name(tmp_path):
    """Probes are grouped under recording_id, then probe_name."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [
        _fake_row(probe_id="pid-1", probe_name="probeA", recording_id="rec1"),
        _fake_row(probe_id="pid-2", probe_name="probeB", recording_id="rec1"),
        _fake_row(probe_id="pid-3", probe_name="probeC", recording_id="rec2"),
    ]
    results = [ProcessResult(probe_id=r.probe_id, recording_id=r.recording_id, wrote_files=True) for r in rows]

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    assert set(probes.keys()) == {"rec1", "rec2"}
    assert set(probes["rec1"].keys()) == {"probeA", "probeB"}
    assert set(probes["rec2"].keys()) == {"probeC"}


def test_same_probe_name_across_recordings_kept_distinct(tmp_path):
    """Same probe_name in two recordings stays distinct (the bug fixed in 2.0.0)."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    # Same probe_name "45883-1" appears in rec1 and rec2 (re-inserted probe).
    rows = [
        _fake_row(probe_id="pid-1a", probe_name="45883-1", recording_id="rec1"),
        _fake_row(probe_id="pid-1b", probe_name="45883-1", recording_id="rec2"),
    ]
    results = [ProcessResult(probe_id=r.probe_id, recording_id=r.recording_id, wrote_files=True) for r in rows]

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    assert probes["rec1"]["45883-1"].probe_id == "pid-1a"
    assert probes["rec2"]["45883-1"].probe_id == "pid-1b"
    # Each entry has exactly one shank — they must NOT be merged into a
    # 2-shank probe under one recording (the pre-2.0.0 silent-merge bug).
    assert probes["rec1"]["45883-1"].num_shanks == 1
    assert probes["rec2"]["45883-1"].num_shanks == 1


def test_multi_shank_collapses_into_one_probe(tmp_path):
    """Rows differing only by probe_shank collapse into one ProbeEntry."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [
        _fake_row(probe_id="pid", probe_name="probeM", recording_id="rec1", probe_shank=0),
        _fake_row(probe_id="pid", probe_name="probeM", recording_id="rec1", probe_shank=1),
    ]
    results = [
        ProcessResult(probe_id="pid", recording_id="rec1", wrote_files=True),
    ]

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    assert list(probes["rec1"].keys()) == ["probeM"]
    entry = probes["rec1"]["probeM"]
    assert entry.num_shanks == 2
    assert {pk.shank for pk in entry.xyz_picks} == {1, 2}  # 1-based in datapackage
    assert entry.logical_probe == "probeM"
    assert entry.ephys_collection == "probeM"
    assert {pk.histology_shank for pk in entry.xyz_picks} == {0, 1}
    assert {pk.ephys_shank for pk in entry.xyz_picks} == {0, 1}


def test_split_stream_quadbase_maps_histology_to_ephys_shank(tmp_path):
    """A split quadbase stream can be one ephys shank of a logical probe."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [
        _fake_row_explicit(
            histology_track_id="track-probe0-shank3",
            logical_probe="probe0",
            ephys_collection="ProbeD",
            recording_id="rec1",
            histology_shank=3,
            ephys_shank=0,
        )
    ]
    results = [
        ProcessResult(
            probe_id="track-probe0-shank3",
            recording_id="rec1",
            wrote_files=True,
        )
    ]
    _touch(tmp_path / "rec1/ProbeD/channels.contactId.npy")

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    entry = probes["rec1"]["ProbeD"]
    assert entry.logical_probe == "probe0"
    assert entry.ephys_collection == "ProbeD"
    assert entry.num_shanks == 1
    assert entry.ephys == _local("rec1/ProbeD")
    assert entry.channel_table is not None
    assert entry.channel_table.contact_id == _local("rec1/ProbeD/channels.contactId.npy")
    assert entry.channel_table.shank_ind == _local("rec1/ProbeD/channels.shankInd.npy")
    assert entry.xyz_picks[0].ccf == _local("rec1/ProbeD/xyz_picks_shank4.json")
    assert entry.xyz_picks[0].histology_shank == 3
    assert entry.xyz_picks[0].ephys_shank == 0
    assert entry.xyz_picks[0].shank == 1


def test_channel_table_omits_missing_contact_id_for_existing_outputs(tmp_path):
    """Datapackage-only regeneration must tolerate older ephys outputs."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [_fake_row(probe_id="pid-1", probe_name="46101", recording_id="rec1")]
    results = [ProcessResult(probe_id="pid-1", recording_id="rec1", wrote_files=True)]

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    channel_table = probes["rec1"]["46101"].channel_table
    assert channel_table is not None
    assert channel_table.contact_id is None


def test_ephys_dir_is_probe_gui_folder_not_spikes_subdir(tmp_path):
    """ephys points at the probe gui_folder (where the ALF collection lives),
    not a nonexistent ``spikes`` subdirectory (regression for the GUI load bug)."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import _build_probes
    from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

    rows = [_fake_row(probe_id="pid-1", probe_name="46101", recording_id="rec1")]
    results = [ProcessResult(probe_id="pid-1", recording_id="rec1", wrote_files=True)]

    probes = _build_probes(rows, results, _fake_outputs(tmp_path), tmp_path, _fake_pipeline_config())

    ephys = probes["rec1"]["46101"].ephys
    assert ephys == _local("rec1/46101")
    assert ephys is not None
    assert not ephys.path.endswith("/spikes")
    # xyz_picks resolve to the same folder — ephys and picks are co-located.
    assert probes["rec1"]["46101"].xyz_picks[0].ccf == _local("rec1/46101/xyz_picks.json")


def test_infer_process_results_from_existing_outputs(tmp_path):
    """Datapackage-only regeneration can infer completed rows from xyz-picks."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import infer_process_results_from_outputs

    outputs = _fake_outputs(tmp_path)
    row_done = _fake_row(probe_id="pid-1", probe_name="probeA", recording_id="rec1")
    row_missing = _fake_row(probe_id="pid-2", probe_name="probeB", recording_id="rec1")

    done_dir = row_done.gui_folder(outputs)
    _touch(done_dir / "xyz_picks.json")
    _touch(done_dir / "xyz_picks_image_space.json")

    results = infer_process_results_from_outputs([row_done, row_missing], outputs)

    assert [(r.probe_id, r.wrote_files) for r in results] == [
        ("pid-1", True),
        ("pid-2", False),
    ]


def _touch(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x")


def _valid_datapackage_on_disk(root: Path) -> DataPackage:
    """Lay down a minimal valid output tree under *root* and return a matching
    DataPackage whose relative paths all resolve."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import (
        CcfSpaceHistology,
        ChannelTablePaths,
        HistologyPaths,
        ImageSpaceHistology,
        ProbeEntry,
        XyzPicks,
    )

    smartspim = root.parent / "SmartSPIM"
    template_to_ccf = root.parent / "spim_template_to_ccf"
    for p in [
        smartspim / "image_atlas_alignment/Ex_561_Em_600/ita.mat",
        smartspim / "image_atlas_alignment/Ex_561_Em_600/itw.nii.gz",
        template_to_ccf / "tca.mat",
        template_to_ccf / "tcw.nii.gz",
    ]:
        _touch(p)

    for rel in [
        "image_space_histology/histology_registration.nrrd",
        "image_space_histology/histology_registration_pipeline.nrrd",
        "image_space_histology/ccf_in_mouse.nrrd",
        "image_space_histology/labels_in_mouse.nrrd",
        "ccf_space_histology/histology_registration.nrrd",
        "rec1/46101/xyz_picks.json",
        "rec1/46101/xyz_picks_image_space.json",
        "rec1/46101/channels.localCoordinates.npy",
        "rec1/46101/channels.rawInd.npy",
        "rec1/46101/channels.contactId.npy",
        "rec1/46101/channels.shankInd.npy",
    ]:
        _touch(root / rel)

    return DataPackage(
        schema_version=SCHEMA_VERSION,
        mouse_id="m1",
        platform="local",
        external_assets={
            "smartspim": ExternalAsset(role="smartspim_registration", name="SmartSPIM"),
            "spim_template_to_ccf": ExternalAsset(role="template_to_ccf", name="spim_template_to_ccf"),
        },
        transforms=TransformPaths(
            image_to_template_affine=_external("smartspim", "image_atlas_alignment/Ex_561_Em_600/ita.mat"),
            image_to_template_warp=_external("smartspim", "image_atlas_alignment/Ex_561_Em_600/itw.nii.gz"),
            template_to_ccf_affine=_external("spim_template_to_ccf", "tca.mat"),
            template_to_ccf_warp=_external("spim_template_to_ccf", "tcw.nii.gz"),
        ),
        histology=HistologyPaths(
            image_space=ImageSpaceHistology(
                registration=_local("image_space_histology/histology_registration.nrrd"),
                registration_pipeline=_local("image_space_histology/histology_registration_pipeline.nrrd"),
                ccf_template=_local("image_space_histology/ccf_in_mouse.nrrd"),
                labels=_local("image_space_histology/labels_in_mouse.nrrd"),
            ),
            ccf_space=CcfSpaceHistology(
                registration=_local("ccf_space_histology/histology_registration.nrrd"),
            ),
        ),
        probes={
            "rec1": {
                "46101": ProbeEntry(
                    probe_id="pid-1",
                    logical_probe="46101",
                    ephys_collection="46101",
                    num_shanks=1,
                    ephys=_local("rec1/46101"),
                    channel_table=ChannelTablePaths(
                        local_coordinates=_local("rec1/46101/channels.localCoordinates.npy"),
                        raw_ind=_local("rec1/46101/channels.rawInd.npy"),
                        contact_id=_local("rec1/46101/channels.contactId.npy"),
                        shank_ind=_local("rec1/46101/channels.shankInd.npy"),
                    ),
                    xyz_picks=[
                        XyzPicks(
                            ccf=_local("rec1/46101/xyz_picks.json"),
                            image_space=_local("rec1/46101/xyz_picks_image_space.json"),
                            shank=None,
                        )
                    ],
                )
            }
        },
    )


def _asset_roots(root: Path) -> list[Path]:
    return [root.parent]


def test_write_datapackage_passes_when_all_paths_exist(tmp_path):
    """A datapackage whose paths all resolve writes and round-trips."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import (
        load_datapackage,
        write_datapackage,
    )

    dp = _valid_datapackage_on_disk(tmp_path)
    assert dp.missing_paths(tmp_path, asset_roots=_asset_roots(tmp_path)) == []
    path = write_datapackage(dp, tmp_path, asset_roots=_asset_roots(tmp_path))  # validate=True by default
    assert load_datapackage(path) == dp


def test_write_datapackage_raises_on_bad_ephys_dir(tmp_path):
    """The original bug: ephys points at a dir lacking the ALF channel file."""
    from aind_ibl_ephys_alignment_preprocessing.datapackage import (
        DataPackageError,
        write_datapackage,
    )

    dp = _valid_datapackage_on_disk(tmp_path)
    entry = dp.probes["rec1"]["46101"].model_copy(update={"ephys": _local("rec1/46101/spikes")})
    bad = dp.model_copy(update={"probes": {"rec1": {"46101": entry}}})

    assert "rec1/46101/spikes/ (ephys directory)" in bad.missing_paths(
        tmp_path,
        asset_roots=_asset_roots(tmp_path),
    )
    with pytest.raises(DataPackageError, match="not found"):
        write_datapackage(bad, tmp_path, asset_roots=_asset_roots(tmp_path))
    # Nothing was written.
    assert not (tmp_path / "datapackage.json").exists()


def test_missing_file_is_reported(tmp_path):
    """A missing referenced file (not just ephys) is flagged."""
    dp = _valid_datapackage_on_disk(tmp_path)
    (tmp_path / "image_space_histology/labels_in_mouse.nrrd").unlink()
    assert "image_space_histology/labels_in_mouse.nrrd" in dp.missing_paths(
        tmp_path,
        asset_roots=_asset_roots(tmp_path),
    )
