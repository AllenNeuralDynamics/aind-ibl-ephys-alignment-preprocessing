"""Timing records must be greppable, summable, and survive a failing step."""

from __future__ import annotations

import asyncio
import logging
import re
import time

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


def _records(caplog):
    return [r.message for r in caplog.records if r.message.startswith("[timing]")]


def _field(line: str, key: str) -> float:
    return float(re.search(rf"\b{key}=([\d.]+)", line).group(1))


def test_concurrent_spans_overlap_and_must_not_be_summed(caplog):
    """Two steps that run together each span the whole period, so seconds double-count.

    This is the semantics the module documents; pinning it here so nobody
    "fixes" the numbers by summing them later.
    """

    async def step(name):
        with timed(name):
            await asyncio.to_thread(time.sleep, 0.15)

    async def both():
        async with asyncio.TaskGroup() as tg:
            tg.create_task(step("a"))
            tg.create_task(step("b"))

    with caplog.at_level(logging.INFO):
        started = time.perf_counter()
        asyncio.run(both())
        wall = time.perf_counter() - started

    spans = [_field(line, "seconds") for line in _records(caplog)]
    assert len(spans) == 2
    assert sum(spans) > wall * 1.5  # summing overcounts, by ~2x here
    assert max(spans) <= wall + 0.05  # any single span is bounded by the wall clock


def test_offsets_let_the_timeline_be_reconstructed(caplog):
    """t0/t1 are what make overlap visible when seconds alone cannot."""

    async def step(name):
        with timed(name):
            await asyncio.to_thread(time.sleep, 0.1)

    async def both():
        async with asyncio.TaskGroup() as tg:
            tg.create_task(step("a"))
            tg.create_task(step("b"))

    with caplog.at_level(logging.INFO):
        asyncio.run(both())

    a, b = _records(caplog)
    # Overlapping iff each starts before the other ends.
    assert _field(a, "t0") < _field(b, "t1")
    assert _field(b, "t0") < _field(a, "t1")
