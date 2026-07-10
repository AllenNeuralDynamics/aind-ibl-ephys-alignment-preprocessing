"""Regression tests for the Code Ocean asset resolver.

Fixtures are REAL Code Ocean assets pulled 2026-07-09 (mice 791094, 781370),
including real provenance.computation / external flags. These lock in the
ground-truth resolution the manifest encodes.
"""

from aind_ibl_ephys_alignment_preprocessing.co_asset_resolver import (
    CandidateAsset,
    _acquisition_from_uri,
    parse_pinned,
    resolve,
    sibling_captures,
    smartspim_acquisition_from_ng,
)

A = CandidateAsset

# --- 791094 real data -------------------------------------------------------

PINNED_791094 = [
    "ecephys_791094_2025-10-08_16-48-57_sorted_2026-04-24_13-43-00",  # -> 77902630
    "ecephys_791094_2025-10-09_14-26-31_sorted_2026-04-23_14-11-00",  # -> 9d137ff6 (truncated-name asset)
]
NG_ACQ_791094 = "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04"

RAW_791094 = [
    A("7d29045e", "ecephys_791094_2025-10-09_14-26-31", ("791094", "ecephys", "raw")),
    A("c06a9fae", "ecephys_791094_2025-10-08_16-48-57", ("791094", "ecephys", "raw")),
]
SPIM_791094 = [
    A("fa97a159", "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04", ("smartspim", "processed")),
    A("e9be1447", "SmartSPIM_791094_2025-11-12_16-36-24", ("791094", "SmartSPIM", "raw")),
]
TAGGED_791094 = [
    A(
        "9d137ff6",
        "ecephys_791094_2025-10-09_sorted_2026-04-23_14-11-00",
        ("ks4", "ecephys", "791094"),
        computation="559ae604",
        external=True,
    ),
    A(
        "757859f4",
        "ecephys_791094_2025-10-09_14-26-31_sorted_2026-04-24_05-44-57",
        ("derived", "ecephys", "791094"),
        computation="559ae604",
        external=False,
    ),
    A(
        "77902630",
        "ecephys_791094_2025-10-08_16-48-57_sorted_2026-04-24_13-43-00",
        ("ks4", "ecephys", "791094"),
        computation="7fa58705",
        external=True,
    ),
    A(
        "30ded4b9",
        "ecephys_791094_2025-10-08_16-48-57_sorted_2026-04-23_08-44-00",
        ("791094", "KS4", "SI", "sorted"),
        computation="7fa58705",
        external=False,
    ),
]


def _resolve_791094():
    return resolve("791094", PINNED_791094, RAW_791094, SPIM_791094, TAGGED_791094, NG_ACQ_791094)


def test_smartspim_from_ng_source():
    res = _resolve_791094()
    assert res.smartspim is not None
    assert res.smartspim.id == "fa97a159"


def test_raw_recordings_resolved():
    res = _resolve_791094()
    assert res.raw["791094_2025-10-08_16-48-57"].id == "c06a9fae"
    assert res.raw["791094_2025-10-09_14-26-31"].id == "7d29045e"


def test_pinned_sortings_match_manifest_ground_truth():
    # The whole point: pinned resolution reproduces the human's choices,
    # including the OLDER ks4 sorting for 2025-10-09 that "newest" gets wrong.
    res = _resolve_791094()
    assert res.sortings["791094_2025-10-08_16-48-57"].id == "77902630"
    assert res.sortings["791094_2025-10-09_14-26-31"].id == "9d137ff6"


def test_fuzzy_name_inconsistency_is_warned():
    res = _resolve_791094()
    assert any("fuzzy: name inconsistency" in w for w in res.warnings)


def test_sibling_captures_are_noted_not_attached():
    res = _resolve_791094()
    bundle_ids = {aid for aid, _ in res.data_assets()}
    # the internal-capture siblings must NOT be attached...
    assert "757859f4" not in bundle_ids
    assert "30ded4b9" not in bundle_ids
    # ...but must be surfaced as notes
    assert any("757859f4" in w for w in res.warnings)
    assert any("30ded4b9" in w for w in res.warnings)


def test_data_assets_bundle_shape():
    res = _resolve_791094()
    bundle = res.data_assets()
    ids = {aid for aid, _ in bundle}
    assert ids == {"fa97a159", "c06a9fae", "7d29045e", "77902630", "9d137ff6"}
    # mount defaults to the asset name
    assert all(mount for _, mount in bundle)


def test_prefer_external_capture():
    # 30ded4b9 is external=False, 77902630 external=True, same computation 7fa58705.
    # The pinned name for 2025-10-08 points at the external one; confirm it wins.
    res = _resolve_791094()
    assert res.sortings["791094_2025-10-08_16-48-57"].external is True


# --- 781370 real data (edge cases) -----------------------------------------

PINNED_781370 = ["ecephys_781370_2025-05-30_15-52-48_sorted_2025-10-01_07-00-38"]
NG_ACQ_781370 = "SmartSPIM_781370_2025-07-31_21-01-37_stitched_2025-09-06_14-34-07"
RAW_781370 = [A("d29c8b22", "ecephys_781370_2025-05-30_15-52-48", ("781370", "ecephys", "raw"))]
# duplicate-named SmartSPIM assets are real for 781370
SPIM_781370 = [
    A(
        "49c166a5",
        "SmartSPIM_781370_2025-07-31_21-01-37_stitched_2025-09-06_14-34-07",
        ("781370", "SmartSPIM", "derived"),
    ),
    A("51f3beb4", "SmartSPIM_781370_2025-07-31_21-01-37_stitched_2025-09-06_14-34-07", ("smartspim", "processed")),
    A("10db6bf1", "SmartSPIM_781370_2025-07-31_21-01-37", ("781370", "SmartSPIM", "raw")),
]
TAGGED_781370 = [
    A(
        "a487eb1a",
        "ecephys_781370_2025-05-30_15-52-48_sorted-curation-sprint_2026-05-13_02-43-03",
        ("derived", "781370", "ecephys"),
    ),
    A(
        "5e968806",
        "ecephys_781370_2025-05-30_15-52-48_preprocessed_2026-04-26_12-26-00",
        ("781370", "IBL", "preprocessed"),
    ),
    A("1e3581cd", "ecephys_781370_2025-05-30_15-52-48_sorted_2025-10-01_07-00-38", ("781370", "derived", "ecephys")),
]


def test_781370_single_sorting_and_excludes_derivatives():
    res = resolve("781370", PINNED_781370, RAW_781370, SPIM_781370, TAGGED_781370, NG_ACQ_781370)
    assert res.sortings["781370_2025-05-30_15-52-48"].id == "1e3581cd"
    assert res.raw["781370_2025-05-30_15-52-48"].id == "d29c8b22"
    # curation-sprint and preprocessed derivatives must not be selected
    ids = {aid for aid, _ in res.data_assets()}
    assert "a487eb1a" not in ids
    assert "5e968806" not in ids


def test_781370_duplicate_named_smartspim_warns_and_picks_one():
    res = resolve("781370", PINNED_781370, RAW_781370, SPIM_781370, TAGGED_781370, NG_ACQ_781370)
    assert res.smartspim is not None
    assert res.smartspim.id in {"49c166a5", "51f3beb4"}
    assert any("SmartSPIM" in w and "match" in w for w in res.warnings)


# --- unit tests on helpers --------------------------------------------------


def test_acquisition_from_uri_real_s3():
    uri = (
        "s3://aind-open-data/SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04"
        "/image_tile_fusing/OMEZarr/Ex_488_Em_525.zarr"
    )
    assert _acquisition_from_uri(uri) == "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04"


def test_acquisition_from_uri_none():
    assert _acquisition_from_uri(None) is None
    assert _acquisition_from_uri("s3://bucket/ecephys_x/y.zarr") is None


def test_smartspim_acquisition_from_ng_minimal():
    ng = {
        "layers": [
            {
                "type": "image",
                "name": "Ex_488",
                "source": (
                    "zarr://s3://aind-open-data/SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04"
                    "/image_tile_fusing/OMEZarr/Ex_488_Em_525.zarr"
                ),
            },
            {"type": "annotation", "name": "668_1_MD", "source": "local://annotations"},
        ]
    }
    assert smartspim_acquisition_from_ng(ng) == "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04"


def test_parse_pinned():
    assert parse_pinned("ecephys_791094_2025-10-08_16-48-57_sorted_2026-04-24_13-43-00") == (
        "791094_2025-10-08_16-48-57",
        "2026-04-24_13-43-00",
    )
    assert parse_pinned("not_a_sorting") is None


def test_sibling_captures_helper():
    chosen = TAGGED_791094[2]  # 77902630, computation 7fa58705
    sibs = sibling_captures(chosen, TAGGED_791094)
    assert [s.id for s in sibs] == ["30ded4b9"]
