"""Tests for the postprocessed-sorting detection predicate and session lookup."""

import json
from pathlib import Path

import pytest

from aind_ibl_ephys_alignment_preprocessing.ephys import (
    find_raw_session_dir,
    find_session_dir,
    find_sorted_session_dir,
    has_sorting_output,
    read_sorted_input_recording,
    resolve_surface_finding,
)


def _mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _mount_sorted(root: Path, *, input_recording: str | None) -> Path:
    """Create a spike-sorted asset mount (spikesorted/ marker + data_description)."""
    _mkdir(root / "spikesorted")
    payload = {} if input_recording is None else {"input_data_name": input_recording}
    (root / "data_description.json").write_text(json.dumps(payload))
    return root


def _mount_raw(root: Path, *, name: str | None = None) -> Path:
    """Create a raw ecephys asset mount (ecephys_compressed/ marker + data_description)."""
    _mkdir(root / "ecephys_compressed")
    if name is not None:
        (root / "data_description.json").write_text(json.dumps({"name": name}))
    return root


def test_no_postprocessed_dir_returns_false(tmp_path):
    assert has_sorting_output(tmp_path, "ProbeA") is False


def test_matching_analyzer_returns_true(tmp_path):
    _mkdir(tmp_path / "postprocessed" / "experiment1_Record Node 104#Neuropix-PXI-100.ProbeA-AP_recording1")
    assert has_sorting_output(tmp_path, "ProbeA") is True


def test_matching_analyzer_zarr_returns_true(tmp_path):
    _mkdir(tmp_path / "postprocessed" / "experiment1_Neuropix.ProbeB-AP_recording1.zarr")
    assert has_sorting_output(tmp_path, "ProbeB") is True


def test_multi_shank_group_returns_true(tmp_path):
    _mkdir(tmp_path / "postprocessed" / "experiment1_Neuropix.ProbeM-AP_recording1_group0")
    _mkdir(tmp_path / "postprocessed" / "experiment1_Neuropix.ProbeM-AP_recording1_group1")
    assert has_sorting_output(tmp_path, "ProbeM") is True


def test_only_lfp_analyzer_returns_false(tmp_path):
    # An LFP-only postprocessed folder does not count as usable sorting output.
    _mkdir(tmp_path / "postprocessed" / "experiment1_Neuropix.ProbeC-LFP_recording1")
    assert has_sorting_output(tmp_path, "ProbeC") is False


def test_other_probe_only_returns_false(tmp_path):
    _mkdir(tmp_path / "postprocessed" / "experiment1_Neuropix.ProbeA-AP_recording1")
    assert has_sorting_output(tmp_path, "ProbeZ") is False


def test_has_sorting_output_none_folder_returns_false():
    # An unresolved session dir (find_session_dir returned None) is not usable output.
    assert has_sorting_output(None, "ProbeA") is False


# --------------------------------------------------------------- find_session_dir --

_SESSION = "ecephys_791093_2025-08-28_15-21-17"
_SORTED = "ecephys_791093_2025-08-28_15-21-17_sorted_2026-04-17_15-26-42"


def test_find_session_dir_flat_layout(tmp_path):
    _mkdir(tmp_path / _SORTED)
    assert find_session_dir(tmp_path, _SORTED) == tmp_path / _SORTED


def test_find_session_dir_nested_under_combined_asset(tmp_path):
    # A combined data asset wraps its members one level down.
    _mkdir(tmp_path / "ibl_combined" / _SORTED)
    assert find_session_dir(tmp_path, _SORTED) == tmp_path / "ibl_combined" / _SORTED


def test_find_session_dir_role_split_raw_and_sorted_unrelated_parents(tmp_path):
    # The whole point: raw and sorted under different parents are each found by name,
    # so extract_continuous no longer needs them to be siblings.
    _mkdir(tmp_path / "raw_combined" / _SESSION)
    _mkdir(tmp_path / "sorted_combined" / _SORTED)
    assert find_session_dir(tmp_path, _SESSION) == tmp_path / "raw_combined" / _SESSION
    assert find_session_dir(tmp_path, _SORTED) == tmp_path / "sorted_combined" / _SORTED


def test_find_session_dir_missing_returns_none(tmp_path):
    assert find_session_dir(tmp_path, _SORTED) is None


def test_find_session_dir_ambiguous_distinct_dirs_raises(tmp_path):
    _mkdir(tmp_path / "a" / _SESSION)
    _mkdir(tmp_path / "b" / _SESSION)
    with pytest.raises(ValueError, match="ambiguous"):
        find_session_dir(tmp_path, _SESSION)


def test_find_session_dir_symlinked_duplicate_is_not_ambiguous(tmp_path):
    # CO stages inputs via symlink chains: the same real dir reached two ways must
    # dedup to a single match, not read as an ambiguous collision.
    _mkdir(tmp_path / "real" / _SESSION)
    (tmp_path / "mirror").symlink_to(tmp_path / "real", target_is_directory=True)
    found = find_session_dir(tmp_path, _SESSION)
    assert found is not None
    assert found.resolve() == (tmp_path / "real" / _SESSION).resolve()


def test_find_session_dir_does_not_descend_into_a_match(tmp_path):
    # A same-named dir nested *inside* a match must not be walked (session trees are
    # huge); pruning also prevents a spurious ambiguity with the inner dir.
    _mkdir(tmp_path / _SESSION / _SESSION)
    assert find_session_dir(tmp_path, _SESSION) == tmp_path / _SESSION


def test_find_session_dir_is_depth_bounded(tmp_path):
    deep = tmp_path / "a" / "b" / "c" / "d" / _SESSION  # depth 5
    _mkdir(deep)
    assert find_session_dir(tmp_path, _SESSION, max_depth=2) is None
    assert find_session_dir(tmp_path, _SESSION, max_depth=8) == deep


# ---------------------------------------------------------- resolve_surface_finding --

_SURFACE = "ecephys_791093_2025-08-28_15-21-17_surface"


def test_resolve_surface_finding_prefers_direct_join(tmp_path):
    # Flat layout: the literal data_root/<name> wins (unchanged legacy behaviour).
    _mkdir(tmp_path / _SURFACE)
    assert resolve_surface_finding(tmp_path, _SURFACE) == tmp_path / _SURFACE


def test_resolve_surface_finding_walks_when_nested(tmp_path):
    # Combined-asset layout: the surface session moved under a wrapper; the direct
    # join misses, so it falls back to a name-walk and finds it.
    _mkdir(tmp_path / "ecephys_raw" / _SURFACE)
    assert resolve_surface_finding(tmp_path, _SURFACE) == tmp_path / "ecephys_raw" / _SURFACE


def test_resolve_surface_finding_missing_returns_direct_join(tmp_path):
    # Neither resolves: return the direct join so the caller/converter reports the
    # missing input, not this helper.
    assert resolve_surface_finding(tmp_path, _SURFACE) == tmp_path / _SURFACE


# ------------------------------------------ content-based resolution (fixed slots) --
# Code Ocean pipeline nodes mount inputs under fixed slot names, not the asset name,
# so the sorted/raw sessions are found by structure, not by matching the asset name.


def test_find_sorted_session_dir_by_marker_under_fixed_slot(tmp_path):
    # Mount name "sorted" != asset name; resolution is by the spikesorted/ child.
    _mount_sorted(tmp_path / "sorted", input_recording="ecephys_A")
    assert find_sorted_session_dir(tmp_path) == tmp_path / "sorted"


def test_find_sorted_session_dir_ignores_non_sorted_mounts(tmp_path):
    # discover + raw mounts have no spikesorted/ child and must be ignored.
    _mkdir(tmp_path / "discover")
    (tmp_path / "discover" / "manifest.csv").write_text("x")
    _mount_raw(tmp_path / "raw")
    _mount_sorted(tmp_path / "sorted", input_recording="ecephys_A")
    assert find_sorted_session_dir(tmp_path) == tmp_path / "sorted"


def test_find_sorted_session_dir_none_when_absent(tmp_path):
    _mkdir(tmp_path / "discover")
    assert find_sorted_session_dir(tmp_path) is None


def test_find_sorted_session_dir_ambiguous_raises(tmp_path):
    _mount_sorted(tmp_path / "sorted_a", input_recording="ecephys_A")
    _mount_sorted(tmp_path / "sorted_b", input_recording="ecephys_B")
    with pytest.raises(ValueError, match="ambiguous spike-sorted"):
        find_sorted_session_dir(tmp_path)


def test_read_sorted_input_recording(tmp_path):
    root = _mount_sorted(tmp_path / "sorted", input_recording="ecephys_786867_x")
    assert read_sorted_input_recording(root) == "ecephys_786867_x"


def test_read_sorted_input_recording_missing_field_returns_none(tmp_path):
    root = _mount_sorted(tmp_path / "sorted", input_recording=None)
    assert read_sorted_input_recording(root) is None


def test_find_raw_session_dir_single(tmp_path):
    _mount_raw(tmp_path / "raw")
    assert find_raw_session_dir(tmp_path) == tmp_path / "raw"


def test_find_raw_session_dir_none_when_absent(tmp_path):
    _mount_sorted(tmp_path / "sorted", input_recording="ecephys_A")
    assert find_raw_session_dir(tmp_path) is None


def test_find_raw_session_dir_disambiguates_surface_by_recording_id(tmp_path):
    # Main raw + a surface-finding raw both carry ecephys_compressed/; the one whose
    # data_description name matches recording_id is the main raw.
    _mount_raw(tmp_path / "raw", name="ecephys_786867_main")
    _mount_raw(tmp_path / "surface_0", name="ecephys_786867_surface")
    got = find_raw_session_dir(tmp_path, recording_id="ecephys_786867_main")
    assert got == tmp_path / "raw"


def test_find_raw_session_dir_ambiguous_without_match_raises(tmp_path):
    _mount_raw(tmp_path / "raw", name="ecephys_786867_main")
    _mount_raw(tmp_path / "surface_0", name="ecephys_786867_surface")
    with pytest.raises(ValueError, match="ambiguous raw ecephys"):
        find_raw_session_dir(tmp_path, recording_id="ecephys_does_not_match")


def _mount_nwb_output(sorted_root, recording_name="ecephys_771432_2025-03-07_18-22-06"):
    """Add zarr-backed NWB output to a sorted asset, as the 2026-07 sorts do.

    An NWB file is a *directory* when zarr-backed, and its processing module is
    literally named ``ecephys`` -- which is one of the raw-session markers.
    """
    for experiment in (1, 2, 3):
        nwb = sorted_root / "nwb" / f"{recording_name}_experiment{experiment}_recording1.nwb"
        (nwb / "processing" / "ecephys").mkdir(parents=True)


def test_nwb_internals_are_not_mistaken_for_raw_sessions(tmp_path):
    """771432: three NWB files inside the sorted asset made the raw lookup ambiguous."""
    _mount_raw(tmp_path / "raw", name="ecephys_771432_2025-05-07_18-22-06")
    sorted_root = _mount_sorted(tmp_path / "sorted", input_recording="ecephys_771432_2025-05-07_18-22-06")
    _mount_nwb_output(sorted_root)

    # Must resolve even without recording_id: only one real raw is mounted, so
    # there is nothing to disambiguate once the NWB internals stop matching.
    assert find_raw_session_dir(tmp_path) == tmp_path / "raw"


def test_nwb_internals_do_not_shadow_the_sorted_session(tmp_path):
    """The same pruning must not cost us the sorted asset it lives inside."""
    sorted_root = _mount_sorted(tmp_path / "sorted", input_recording="ecephys_A")
    _mount_nwb_output(sorted_root)

    assert find_sorted_session_dir(tmp_path) == sorted_root
