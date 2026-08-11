"""Structured step timing, so a run's cost can be attributed rather than inferred.

The histology and ephys logs record what happened and in what order, never how
long it took, which leaves voxel counts as the only basis for deciding where a
run's time goes -- an inference that has been wrong in both directions.

Each timed step emits one line carrying a machine-readable ``step=<name>
seconds=<elapsed> t0=<start> t1=<end>``, where ``t0``/``t1`` are seconds since
this module was imported.

**The spans overlap, so do not sum them.** The stage runs channel reads, both
CCF warps and several writes concurrently, and a step's elapsed time is the
wall-clock span of a coroutine that was suspended across an ``await`` -- it
includes whatever else the loop and the thread pool were doing meanwhile. Two
concurrent 1-second writes each report ~1 s while together taking ~1 s. That is
why ``t0``/``t1`` are emitted: they let a timeline be reconstructed, overlap be
seen, and the critical path be found. Summing ``seconds`` gives a number with no
physical meaning.

What the records are good for: ranking steps within a run, comparing the same
step across runs (level 3 vs level 4 on the same volume), and spotting a step
whose span dwarfs everything else.

Deliberately dependency-free and cheap: a ``perf_counter`` pair around work that
takes seconds to minutes costs nothing measurable, so timing is always on and
there is no flag to forget.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any

logger = logging.getLogger(__name__)

#: Marker that makes timing lines greppable and parseable out of a mixed log.
_PREFIX = "[timing]"

#: Reference point for ``t0``/``t1``, so offsets are comparable within a run.
_ORIGIN = time.perf_counter()


def format_timing(step: str, seconds: float, **fields: Any) -> str:
    """Render one timing record.

    Parameters
    ----------
    step : str
        Dotted step name, e.g. ``"histology.zarr_read"``.
    seconds : float
        Elapsed wall-clock seconds.
    **fields
        Extra key/value context (level, shape, path, ...). Values are rendered
        with ``str`` and must not contain spaces if they are to parse cleanly.

    Returns
    -------
    str
        A line of the form ``[timing] step=<name> seconds=<elapsed> k=v ...``.
    """
    extra = "".join(f" {key}={value}" for key, value in fields.items())
    return f"{_PREFIX} step={step} seconds={seconds:.3f}{extra}"


@contextmanager
def timed(step: str, **fields: Any) -> Iterator[dict[str, Any]]:
    """Time the enclosed block and log one structured record on exit.

    The yielded dict is merged into the logged fields, so a caller can attach
    context it only learns inside the block (a shape, a voxel count, a level).

    Emits on failure too, tagged ``ok=0`` -- a step that died after twenty
    minutes is exactly the one worth knowing the duration of.

    Parameters
    ----------
    step : str
        Dotted step name.
    **fields
        Context known before the block runs.

    Yields
    ------
    dict[str, Any]
        Mutable field bag, merged into the emitted record.

    Examples
    --------
    >>> with timed("histology.zarr_read", level=4) as record:  # doctest: +SKIP
    ...     img = read(level)
    ...     record["mvox"] = img.size // 10**6
    """
    collected: dict[str, Any] = {}
    start = time.perf_counter()
    ok = True
    try:
        yield collected
    except BaseException:
        ok = False
        raise
    finally:
        end = time.perf_counter()
        merged = {**fields, **collected}
        if not ok:
            merged["ok"] = 0
        # Offsets, not just duration: concurrent steps overlap, so a timeline is
        # the only way to tell a slow step from a step that merely spanned one.
        merged["t0"] = f"{start - _ORIGIN:.3f}"
        merged["t1"] = f"{end - _ORIGIN:.3f}"
        logger.info(format_timing(step, end - start, **merged))
