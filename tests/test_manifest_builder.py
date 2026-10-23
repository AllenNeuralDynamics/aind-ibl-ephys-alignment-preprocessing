"""Tests for building a manifest from a track-to-stream mapping."""

import json
from pathlib import Path

import pytest

from aind_ibl_ephys_alignment_preprocessing.manifest_builder import (
    ManifestBuildError,
    MouseMapping,
    SessionMapping,
    TrackAssignment,
    annotation_point_counts,
    build_from_files,
    build_manifest,
    local_streams,
    streams_in,
    unmapped_tracks,
)

SORT = "ecephys_1_2026-01-01_00-00-00_sorted_2026-01-02_00-00-00"
QUAD = {"46101-1", "46101-2", "46101-3", "46101-4"}


def _mapping(*tracks, annotation="tracks.json"):
    return MouseMapping("1", Path(annotation), (SessionMapping(SORT, tuple(tracks)),))


def _points(*names, n=5):
    return dict.fromkeys(names, n)


def _problems(mapping, streams, points):
    with pytest.raises(ManifestBuildError) as err:
        build_manifest(mapping, {SORT: streams}, points)
    return err.value.problems


def test_split_streams_share_a_probe_and_take_their_shank_from_the_stream():
    tracks = [TrackAssignment(f"t{i}", f"46101-{i}") for i in (1, 2, 3, 4)]
    df = build_manifest(_mapping(*tracks), {SORT: QUAD}, _points("t1", "t2", "t3", "t4"))
    assert df["logical_probe"].tolist() == ["46101"] * 4
    assert df["histology_shank"].tolist() == [0, 1, 2, 3]
    assert df["ephys_shank"].tolist() == [0] * 4


def test_a_single_stream_probe_leaves_shanks_blank():
    df = build_manifest(_mapping(TrackAssignment("bla", "46810")), {SORT: {"46810"}}, _points("bla"))
    row = df.iloc[0]
    assert (row["logical_probe"], row["ephys_collection"]) == ("46810", "46810")
    assert row[["histology_shank", "ephys_shank"]].isna().all()
    assert (row["mouseid"], row["probe_file"]) == ("1", "tracks")


def test_a_suffix_without_a_sibling_stream_is_the_probe_itself():
    df = build_manifest(_mapping(TrackAssignment("a", "ProbeA-2")), {SORT: {"ProbeA-2"}}, _points("a"))
    assert df.iloc[0]["logical_probe"] == "ProbeA-2"


def test_tracks_on_one_multi_shank_stream_take_the_given_shanks():
    tracks = [TrackAssignment("a", "45883", 0), TrackAssignment("b", "45883", 2)]
    df = build_manifest(_mapping(*tracks), {SORT: {"45883"}}, _points("a", "b"))
    assert df["histology_shank"].tolist() == [0, 2]
    assert df["ephys_shank"].tolist() == [0, 2]


def test_tracks_sharing_a_stream_without_distinct_shanks_are_refused():
    tracks = [TrackAssignment("a", "45883"), TrackAssignment("b", "45883")]
    assert _problems(_mapping(*tracks), {"45883"}, _points("a", "b")) == [
        f"{SORT}: tracks on stream 45883 need distinct shanks"
    ]


def test_a_split_stream_takes_one_track():
    tracks = [TrackAssignment("a", "46101-1"), TrackAssignment("b", "46101-1")]
    problems = _problems(_mapping(*tracks), {"46101-1", "46101-2"}, _points("a", "b"))
    assert f"{SORT}: split stream 46101-1 records one shank; give it one track and no shank" in problems


def test_every_problem_is_reported_at_once():
    tracks = [TrackAssignment("a", "46810"), TrackAssignment("ghost", "99999"), TrackAssignment("thin", "50213")]
    problems = _problems(_mapping(*tracks), {"46810", "50213", "50198"}, {"a": 5, "thin": 1})
    assert set(problems) == {
        "track ghost is not a layer of tracks.json",
        "track thin has 1 point(s)",
        f"{SORT}: stream 50198 has no track",
        f"{SORT}: stream 99999 is not in the sort",
    }


def test_a_track_assigned_in_two_sessions_is_refused():
    other = "ecephys_1_2026-01-03_00-00-00_sorted_2026-01-04_00-00-00"
    mapping = MouseMapping(
        "1",
        Path("tracks.json"),
        (
            SessionMapping(SORT, (TrackAssignment("a", "46810"),)),
            SessionMapping(other, (TrackAssignment("a", "46810"),)),
        ),
    )
    with pytest.raises(ManifestBuildError) as err:
        build_manifest(mapping, {SORT: {"46810"}, other: {"46810"}}, _points("a"))
    assert err.value.problems == ["track a is assigned 2 times"]


def test_streams_are_parsed_from_analyzer_names_across_experiments():
    names = [
        "experiment1_Record Node 111#Neuropix-PXI-100.46101-1_recording1.zarr",
        "experiment2_Record Node 111#Neuropix-PXI-100.46101-1_recording1.zarr",
        "experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1",
        "experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-LFP_recording1",
    ]
    assert streams_in(names) == {"46101-1", "ProbeA"}


def test_layers_are_counted_by_their_points():
    state = {
        "layers": [
            {"name": "Ex_639_Em_667.zarr", "type": "image"},
            {"name": "t1", "type": "annotation", "annotations": [{"type": "point"}] * 3},
            {"name": "empty", "type": "annotation", "annotations": []},
        ]
    }
    counts = annotation_point_counts(state)
    assert counts == {"t1": 3, "empty": 0}
    assert unmapped_tracks(_mapping(), counts) == ["t1"]


def test_a_mapping_file_builds_against_a_mounted_sort(tmp_path):
    post = tmp_path / "data" / "wrapper" / SORT / "postprocessed"
    for s in ("46810", "45883"):
        (post / f"experiment1_Record Node 111#Neuropix-PXI-100.{s}_recording1.zarr").mkdir(parents=True)
    layers = [{"name": n, "type": "annotation", "annotations": [{"type": "point"}] * 2} for n in ("bla", "pl", "x")]
    (tmp_path / "tracks.json").write_text(json.dumps({"layers": layers}))
    (tmp_path / "m.toml").write_text(
        f'mouseid = "1"\nannotation = "tracks.json"\n\n[[session]]\nsorted_recording = "{SORT}"\n\n'
        '[session.tracks]\nbla = "46810"\npl = { stream = "45883", shank = 2 }\n'
    )
    assert local_streams(tmp_path / "data", SORT) == {"46810", "45883"}
    df, unmapped = build_from_files(tmp_path / "m.toml", tmp_path / "data")
    assert df["histology_track_id"].tolist() == ["bla", "pl"]
    assert df.loc[1, "histology_shank"] == 1
    assert unmapped == ["x"]


def test_per_brain_overrides_are_written_only_when_set():
    tracks = (TrackAssignment("bla", "46810"),)
    plain = build_manifest(_mapping(*tracks), {SORT: {"46810"}}, _points("bla"))
    assert "registration_asset" not in plain and "template_to_ccf_asset" not in plain

    mapping = MouseMapping(
        "1",
        Path("tracks.json"),
        (SessionMapping(SORT, tracks),),
        registration_asset="reg/ccf_Ex_639_Em_667",
        template_to_ccf_asset="idisco/t2c",
    )
    df = build_manifest(mapping, {SORT: {"46810"}}, _points("bla"))
    assert df.columns[-2:].tolist() == ["registration_asset", "template_to_ccf_asset"]
    assert df.iloc[0]["template_to_ccf_asset"] == "idisco/t2c"
