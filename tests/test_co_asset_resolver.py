"""Regression tests for the Code Ocean asset resolver.

Fixtures are REAL Code Ocean assets pulled 2026-07-09 (mice 791094, 781370),
including real provenance.computation / external flags. These lock in the
ground-truth resolution the manifest encodes.
"""

from aind_ibl_ephys_alignment_preprocessing.co_asset_resolver import (
    CandidateAsset,
    _acquisition_from_uri,
    parse_pinned,
    raw_key_of,
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


# --- 750107 real data: corrected raw + duplicate names ----------------------
#
# Pulled 2026-08-11. The pinned sorting was produced from the *corrected* raw;
# two same-named uncorrected assets also exist. Name reconstruction attached an
# uncorrected one -- a different recording of the same session.

PINNED_750107 = ["ecephys_750107_2025-01-28_13-49-03_sorted_2026-07-30_11-52-00"]
NG_ACQ_750107 = "SmartSPIM_750107_2025-02-19_12-00-00"
KEY_750107 = "750107_2025-01-28_13-49-03"

# type=result, so absent from the type=Dataset raw search -- only reachable via tag:<mouseid>
CORRECTED_750107 = A("a31874d6", "ecephys_750107_2025-01-28_13-49-03_corrected", ("ecephys", "raw", "750107"))
# the two uncorrected assets share a name exactly
RAW_750107 = [
    A("d13f992d", "ecephys_750107_2025-01-28_13-49-03", ("750107", "ecephys", "raw")),
    A("3a13749e", "ecephys_750107_2025-01-28_13-49-03", ("750107", "ecephys", "raw")),
]
SORTING_750107 = A(
    "1b4f68df",
    "ecephys_750107_2025-01-28_13-49-03_sorted_2026-07-30_11-52-00",
    ("derived", "750107", "ecephys"),
    computation="741d3f21",
    external=True,
    source_assets=("a31874d6",),
)
SPIM_750107 = [A("aaaaaaaa", NG_ACQ_750107, ("750107", "SmartSPIM", "raw"))]
TAGGED_750107 = [SORTING_750107, CORRECTED_750107]


def _resolve_750107(raw=None, tagged=None):
    return resolve(
        "750107",
        PINNED_750107,
        RAW_750107 if raw is None else raw,
        SPIM_750107,
        TAGGED_750107 if tagged is None else tagged,
        NG_ACQ_750107,
    )


def test_750107_raw_follows_sorting_provenance_not_name():
    # The whole bug: the corrected raw is what the sorting was produced from.
    res = _resolve_750107()
    assert res.raw[KEY_750107].id == "a31874d6"
    assert res.raw[KEY_750107].name.endswith("_corrected")


def test_750107_corrected_raw_is_attached_though_absent_from_raw_search():
    # a31874d6 is type=result and appears only in the tag:<mouseid> pool.
    res = _resolve_750107()
    assert "a31874d6" in {aid for aid, _ in res.data_assets()}
    assert not any(a.id == "a31874d6" for a in RAW_750107)


def test_750107_uncorrected_duplicates_are_not_attached():
    res = _resolve_750107()
    ids = {aid for aid, _ in res.data_assets()}
    assert "d13f992d" not in ids
    assert "3a13749e" not in ids


def test_750107_provenance_path_is_silent():
    # Following provenance is the designed path, so it must not add noise.
    res = _resolve_750107()
    assert not any(KEY_750107 in w and "raw" in w.lower() for w in res.warnings)


def test_750107_name_fallback_warns_and_flags_ambiguity():
    # Strip the provenance: the resolver must fall back, say so, and report that
    # two assets share the key rather than silently letting the last one win.
    no_prov = SORTING_750107.__class__(
        id=SORTING_750107.id,
        name=SORTING_750107.name,
        tags=SORTING_750107.tags,
        computation=SORTING_750107.computation,
        external=SORTING_750107.external,
    )
    res = _resolve_750107(tagged=[no_prov, CORRECTED_750107])
    assert res.raw[KEY_750107].id in {"d13f992d", "3a13749e"}
    assert any("no provenance" in w for w in res.warnings)
    assert any("AMBIGUOUS" in w and "share this recording key" in w for w in res.warnings)


def test_raw_key_of_accepts_suffixed_upload_and_rejects_derived():
    assert raw_key_of("ecephys_750107_2025-01-28_13-49-03") == KEY_750107
    assert raw_key_of("ecephys_750107_2025-01-28_13-49-03_corrected") == KEY_750107
    assert raw_key_of("ecephys_750107_2025-01-28_13-49-03_sorted_2026-07-30_11-52-00") is None
    assert raw_key_of("ecephys_781370_2025-05-30_15-52-48_sorted-curation-sprint_2026-05-13_02-43-03") is None
    assert raw_key_of("ecephys_781370_2025-05-30_15-52-48_preprocessed_2026-04-26_12-26-00") is None
    assert raw_key_of("SmartSPIM_750107_2025-02-19_12-00-00") is None


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
