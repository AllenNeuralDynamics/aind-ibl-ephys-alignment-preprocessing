"""Concurrency primitives: IO lane throttling, semaphores, thread helpers."""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Callable
from contextlib import nullcontext
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any

from filelock import FileLock

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, OutputDirs

logger = logging.getLogger(__name__)


def _maybe_semaphore(limit: int | None) -> asyncio.BoundedSemaphore | nullcontext:  # type: ignore[type-arg]
    """Create a bounded semaphore if *limit* is set, else a no-op context manager."""
    return asyncio.BoundedSemaphore(limit) if limit is not None else nullcontext()


async def to_thread_logged(fn: Callable[..., Any], *a: Any, **kw: Any) -> Any:
    """Run *fn* in a thread with exception logging."""

    def _wrap() -> Any:
        try:
            return fn(*a, **kw)
        except Exception:
            logging.exception("Threaded call failed: %r%r%r", fn, a, kw)
            raise

    return await asyncio.to_thread(_wrap)


class IOLimits:
    """Per-path I/O lane throttling via semaphores.

    Parameters
    ----------
    scratch : int | None
        Concurrency limit for scratch path I/O.
    results : int | None
        Concurrency limit for results path I/O.
    data : int | None
        Concurrency limit for data path I/O.
    scratch_root : str
        Scratch root path string for lane matching.
    results_root : str
        Results root path string for lane matching.
    data_root : str
        Data root path string for lane matching.
    """

    def __init__(
        self,
        scratch: int | None = None,
        results: int | None = None,
        data: int | None = None,
        scratch_root: str = "/scratch",
        results_root: str = "/results",
        data_root: str = "/data",
    ) -> None:
        self._lanes: dict[str, asyncio.BoundedSemaphore | nullcontext] = {  # type: ignore[type-arg]
            scratch_root: _maybe_semaphore(scratch),
            results_root: _maybe_semaphore(results),
            data_root: _maybe_semaphore(data),
        }
        self._default_key = results_root

    def lane_for(self, path: str) -> asyncio.BoundedSemaphore | nullcontext:  # type: ignore[type-arg]
        """Return the semaphore/context for the lane matching *path*."""
        for root, sem in self._lanes.items():
            if path.startswith(root):
                return sem
        return self._lanes[self._default_key]


class Limits:
    """Concurrency limits for the async pipeline.

    Parameters
    ----------
    max_ephys : int | None
        Max concurrent ephys extractions.
    max_registration : int | None
        Max concurrent ANTs registration transforms.
    max_manifest_rows : int | None
        Max concurrent manifest row processing.
    max_scratch : int | None
        I/O concurrency for scratch.
    max_results : int | None
        I/O concurrency for results.
    max_data : int | None
        I/O concurrency for data.
    scratch_root : str
        Scratch root for IO lane matching.
    results_root : str
        Results root for IO lane matching.
    data_root : str
        Data root for IO lane matching.
    """

    ephys: asyncio.BoundedSemaphore | nullcontext  # type: ignore[type-arg]
    registration: asyncio.BoundedSemaphore | nullcontext  # type: ignore[type-arg]
    manifest_rows: asyncio.BoundedSemaphore | nullcontext  # type: ignore[type-arg]
    io: IOLimits
    max_ephys: int | None
    max_registration: int | None
    max_manifest_rows: int | None
    max_scratch: int | None
    max_results: int | None
    max_data: int | None

    def __init__(
        self,
        max_ephys: int | None = 2,
        max_registration: int | None = 1,
        max_manifest_rows: int | None = None,
        max_scratch: int | None = None,
        max_results: int | None = None,
        max_data: int | None = None,
        scratch_root: str = "/scratch",
        results_root: str = "/results",
        data_root: str = "/data",
    ) -> None:
        self.max_ephys = max_ephys
        self.max_registration = max_registration
        self.max_manifest_rows = max_manifest_rows
        self.max_scratch = max_scratch
        self.max_results = max_results
        self.max_data = max_data
        self.ephys = _maybe_semaphore(max_ephys)
        self.registration = _maybe_semaphore(max_registration)
        self.manifest_rows = _maybe_semaphore(max_manifest_rows)
        self.io = IOLimits(max_scratch, max_results, max_data, scratch_root, results_root, data_root)


async def io_to_thread_on(limits: Limits, target_path: str, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Run *fn* in a thread, throttled by the I/O lane for *target_path*."""
    async with limits.io.lane_for(str(target_path)):
        return await to_thread_logged(fn, *args, **kwargs)


def _run_ephys_sync(mr: ManifestRow, out: OutputDirs, data_root: Path) -> None:
    """Single recording ephys in a separate *process*.

    Idempotent via a disk marker and file lock.
    """
    from aind_ephys_ibl_gui_conversion.ephys import extract_continuous, extract_spikes

    sorted_rec = mr.sorted_recording
    recording_id = mr.recording_id
    results_folder = out.tracks_root.parent / recording_id
    results_folder.mkdir(parents=True, exist_ok=True)

    done = results_folder / ".ephys.done"
    lock = FileLock(str(results_folder.with_suffix(".lock")))
    with lock:
        if done.exists():
            logger.info("[Ephys %s] Skipping (already complete)", recording_id)
            return
        logger.info("[Ephys %s] Starting extraction (process pool)", recording_id)
        recording_folder = data_root / sorted_rec
        if mr.surface_finding is not None:
            extract_continuous(
                recording_folder,
                results_folder,
                probe_surface_finding=data_root / str(mr.surface_finding),
            )
        else:
            extract_continuous(recording_folder, results_folder)
        extract_spikes(recording_folder, results_folder)
        done.write_text("ok")
        logger.info("[Ephys %s] Completed", recording_id)


@dataclass
class EphysCoordinator:
    """Single-flight ephys deduplication across async tasks.

    Parameters
    ----------
    pool : ProcessPoolExecutor
        Process pool for ephys extraction.
    max_inflight : int
        Max concurrent ephys processes.
    """

    pool: Any  # ProcessPoolExecutor
    max_inflight: int = 4
    _tasks: dict[str, asyncio.Task[None]] = field(default_factory=dict, init=False)
    _lock: asyncio.Lock = field(default_factory=asyncio.Lock, init=False)
    _sem: asyncio.Semaphore = field(init=False)

    def __post_init__(self) -> None:
        self._sem = asyncio.Semaphore(self.max_inflight)

    async def ensure(self, key: str, run_sync: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        """Ensure exactly one process-pooled run per *key*. Others await the same Task."""
        async with self._lock:
            t = self._tasks.get(key)
            if t is None:

                async def _runner() -> None:
                    async with self._sem:
                        loop = asyncio.get_running_loop()
                        await loop.run_in_executor(self.pool, partial(run_sync, *args, **kwargs))

                t = asyncio.create_task(_runner(), name=f"ephys-{key}")

                def _cleanup(_: Any) -> None:
                    self._tasks.pop(key, None)

                t.add_done_callback(_cleanup)
                self._tasks[key] = t
        await t
