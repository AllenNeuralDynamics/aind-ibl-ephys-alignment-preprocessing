"""Timing records must be greppable, summable, and survive a failing step."""

from __future__ import annotations

import logging

import pytest

from aind_ibl_ephys_alignment_preprocessing._timing import format_timing, timed


def test_record_is_machine_readable():
    line = format_timing("histology.zarr_read", 12.3456, level=4, mvox=73)
    assert line == "[timing] step=histology.zarr_read seconds=12.346 level=4 mvox=73"


def test_timed_emits_one_record(caplog):
    with caplog.at_level(logging.INFO):
        with timed("histology.warp", interpolator="genericLabel"):
            pass

    lines = [r.message for r in caplog.records if r.message.startswith("[timing]")]
    assert len(lines) == 1
    assert "step=histology.warp" in lines[0]
    assert "interpolator=genericLabel" in lines[0]


def test_context_learned_inside_the_block_is_recorded(caplog):
    """A voxel count is only known after the read, but belongs on the read's record."""
    with caplog.at_level(logging.INFO):
        with timed("histology.zarr_read", level=3) as record:
            record["mvox"] = 646

    line = next(r.message for r in caplog.records if r.message.startswith("[timing]"))
    assert "level=3" in line and "mvox=646" in line


def test_a_failing_step_still_reports_its_duration(caplog):
    """The step that died after twenty minutes is the one worth timing."""
    with caplog.at_level(logging.INFO), pytest.raises(RuntimeError):
        with timed("histology.write", volume="labels_in_mouse.nrrd"):
            raise RuntimeError("boom")

    line = next(r.message for r in caplog.records if r.message.startswith("[timing]"))
    assert "step=histology.write" in line and "ok=0" in line
