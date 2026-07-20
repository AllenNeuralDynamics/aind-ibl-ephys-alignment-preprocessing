"""Tests for the postprocessed-sorting detection predicate."""

from pathlib import Path

from aind_ibl_ephys_alignment_preprocessing.ephys import has_sorting_output


def _mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


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
