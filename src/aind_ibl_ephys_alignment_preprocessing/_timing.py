"""Structured step timing, so a run's cost can be attributed rather than inferred.

The histology and ephys logs record what happened and in what order, never how
long it took, which leaves voxel counts as the only basis for deciding where a
run's time goes -- an inference that has been wrong in both directions.

Each timed step emits one line carrying a machine-readable ``step=<name>
seconds=<elapsed>``, so a run's steps can be summed and ranked with ``grep``
rather than eyeballed, and later folded into a datapackage's cost profile.

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
        elapsed = time.perf_counter() - start
        merged = {**fields, **collected}
        if not ok:
            merged["ok"] = 0
        logger.info(format_timing(step, elapsed, **merged))
