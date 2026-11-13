"""Resolve a manifest's physical probe onto the per-shank streams it was recorded as.

Quad-base probes record each shank as its own stream (``ProbeA-1`` .. ``ProbeA-4``),
so the ephys collection is not the probe a person names. A row giving the probe and
a shank is resolved to the one stream whose recorded contacts lie on that shank, and
rewritten into the split-stream shape :mod:`manifest_builder` produces. The shank is
read from the sort's own probe geometry (contact ids ``s<shank>e<electrode>``); the
stream-name suffix is only cross-checked against it, never trusted alone.
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Callable
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import pandas as pd

from aind_ibl_ephys_alignment_preprocessing.ephys import (
    _normalize_collection,
    analyzer_collection,
    find_session_dir,
    find_sorted_session_dirs,
)
from aind_ibl_ephys_alignment_preprocessing.manifest_builder import _SPLIT_STREAM_RE
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow

logger = logging.getLogger(__name__)

_CONTACT_SHANK_RE = re.compile(r"^s(\d+)e", re.IGNORECASE)


class StreamResolutionError(ValueError):
    """A row names a split probe but its streams do not identify one shank's stream."""


@dataclass(frozen=True, slots=True)
class SortedStream:
    """One spike-sorted collection and the probe geometry its analyzer recorded.

    Attributes
    ----------
    collection : str
        Ephys collection token (``ProbeA-2``).
    serial : str or None
        Probe serial number, when the geometry records one.
    shanks : frozenset of int
        Shanks the stream's contacts lie on; empty when the geometry is unreadable.
    """

    collection: str
    serial: str | None
    shanks: frozenset[int]


@dataclass(frozen=True, slots=True)
class StreamResolution:
    """Where a probe-and-shank row was resolved to."""

    probe: str
    shank: int
    collection: str
    serial: str | None

    def note(self) -> str:
        """Describe the resolution for logs and the coverage record."""
        serial = f", serial {self.serial}" if self.serial else ""
        return f"{self.probe} shank {self.shank} resolved to per-shank stream {self.collection}{serial}"


def _recording_info_attributes(analyzer: Path) -> dict[str, Any] | None:
    """Read a SortingAnalyzer's ``recording_info`` attributes (zarr v2 or v3)."""
    info = analyzer / "recording_info"
    try:
        if (info / ".zattrs").is_file():
            return dict(json.loads((info / ".zattrs").read_text()))
        if (info / "zarr.json").is_file():
            return dict(json.loads((info / "zarr.json").read_text()).get("attributes", {}))
        if (analyzer / ".zmetadata").is_file():
            metadata = json.loads((analyzer / ".zmetadata").read_text())["metadata"]
            return dict(metadata.get("recording_info/.zattrs", {}))
    except (OSError, ValueError, KeyError, TypeError) as error:
        logger.warning("[streams] unreadable recording_info in %s: %s", analyzer, error)
    return None


def _probe_geometry(attributes: dict[str, Any] | None) -> tuple[str | None, frozenset[int]]:
    """Return the serial and shanks of the single probe a per-stream analyzer records."""
    try:
        probes = attributes["probegroup"]["probes"] if attributes else []
    except (KeyError, TypeError):
        return None, frozenset()
    serials: set[str] = set()
    shanks: set[int] = set()
    for probe in probes:
        serial = (probe.get("annotations") or {}).get("serial_number")
        if serial:
            serials.add(str(serial))
        from_ids = {int(m.group(1)) for c in probe.get("contact_ids") or [] if (m := _CONTACT_SHANK_RE.match(str(c)))}
        if from_ids:
            shanks |= from_ids
        else:
            shanks |= {int(s) for s in probe.get("shank_ids") or [] if str(s).strip().isdigit()}
    return (serials.pop() if len(serials) == 1 else None), frozenset(shanks)


def sorted_streams(recording_folder: Path | None) -> list[SortedStream]:
    """List a sorted asset's spike-sorted collections with their recorded geometry.

    Parameters
    ----------
    recording_folder : Path or None
        The sorted-recording asset directory. ``None`` or a folder without
        ``postprocessed/`` yields no streams.

    Returns
    -------
    list of SortedStream
        One entry per non-LFP analyzer, in directory order.
    """
    postprocessed = recording_folder / "postprocessed" if recording_folder is not None else None
    if postprocessed is None or not postprocessed.is_dir():
        return []
    streams = []
    for analyzer in sorted(p for p in postprocessed.iterdir() if p.is_dir()):
        collection = analyzer_collection(analyzer.name)
        if collection is None:
            continue
        serial, shanks = _probe_geometry(_recording_info_attributes(analyzer))
        streams.append(SortedStream(collection, serial, shanks))
    return streams


def resolve_split_stream(
    streams: list[SortedStream], ephys_collection: str, shank: int | None
) -> StreamResolution | None:
    """Resolve a probe name to the per-shank stream recording ``shank``.

    Parameters
    ----------
    streams : list of SortedStream
        The recording's sorted collections (:func:`sorted_streams`).
    ephys_collection : str
        The row's collection, possibly a probe recorded as per-shank streams.
    shank : int or None
        The row's 0-based physical shank.

    Returns
    -------
    StreamResolution or None
        ``None`` when the collection matches a stream exactly or no per-shank
        streams of it exist; the row is then used as written.

    Raises
    ------
    StreamResolutionError
        If per-shank streams exist but do not identify exactly one stream: no
        shank given, streams from different probes, a stream spanning several
        shanks or none, a recorded shank that contradicts the stream-name
        suffix, or no stream on the requested shank.
    """
    wanted = _normalize_collection(ephys_collection)
    if wanted is None or any(s.collection == wanted for s in streams):
        return None
    split = [
        (int(m.group("position")), s)
        for s in streams
        if (m := _SPLIT_STREAM_RE.match(s.collection)) is not None and m.group("base") == wanted
    ]
    if not split:
        return None
    names = [s.collection for _, s in split]
    if shank is None:
        raise StreamResolutionError(f"{wanted} is recorded as per-shank streams {names}; the row needs probe_shank")
    serials = {s.serial for _, s in split if s.serial is not None}
    if len(serials) > 1:
        raise StreamResolutionError(f"streams {names} come from different probes (serials {sorted(serials)})")
    for position, stream in split:
        if len(stream.shanks) != 1:
            recorded = sorted(stream.shanks) or "no readable shank"
            raise StreamResolutionError(f"stream {stream.collection} records {recorded}; expected one shank")
        (recorded_shank,) = stream.shanks
        if recorded_shank != position - 1:
            raise StreamResolutionError(
                f"stream {stream.collection} records shank {recorded_shank}, contradicting its name"
            )
    matching = [s for _, s in split if s.shanks == {shank}]
    if not matching:
        raise StreamResolutionError(f"no stream of {wanted} records shank {shank}; streams are {names}")
    return StreamResolution(probe=wanted, shank=shank, collection=matching[0].collection, serial=matching[0].serial)


def resolved_row(mr: ManifestRow, resolution: StreamResolution) -> ManifestRow:
    """Rewrite a row into the split-stream shape: one-shank collection, physical probe kept."""
    return replace(
        mr,
        probe_name=resolution.collection,
        ephys_collection=resolution.collection,
        logical_probe=resolution.probe,
        histology_shank=resolution.shank,
        ephys_shank=0,
    )


@dataclass(frozen=True)
class ResolvedManifest:
    """A manifest with probe-and-shank rows resolved onto their per-shank streams.

    Attributes
    ----------
    df : pandas.DataFrame
        The manifest with resolved rows' columns rewritten; other rows untouched.
    rows : list of ManifestRow
        One row per ``df`` row, in order.
    resolutions : dict[int, StreamResolution]
        Position in ``df`` to its resolution.
    errors : dict[int, str]
        Position in ``df`` to why it could not be resolved.
    """

    df: pd.DataFrame
    rows: list[ManifestRow]
    resolutions: dict[int, StreamResolution]
    errors: dict[int, str]

    def raise_errors(self) -> None:
        """Raise one :class:`StreamResolutionError` naming every unresolved row."""
        if self.errors:
            lines = [f"row {position + 1} ({self.rows[position].probe_id}): {e}" for position, e in self.errors.items()]
            raise StreamResolutionError("manifest rows could not be resolved:\n" + "\n".join(lines))


def default_locator(data_root: Path) -> Callable[[ManifestRow], Path | None]:
    """Find a row's sorted asset by its pin, else by the recording it was sorted from."""
    by_recording: dict[str, Path] | None = None

    def locate(mr: ManifestRow) -> Path | None:
        nonlocal by_recording
        found = find_session_dir(data_root, str(mr.sorted_recording))
        if found is not None:
            return found
        if by_recording is None:
            by_recording = find_sorted_session_dirs(data_root)
        return by_recording.get(mr.recording_id)

    return locate


def read_resolved_manifest(manifest_csv: Path, data_root: Path, *, skip_ephys: bool = False) -> pd.DataFrame:
    """Read a manifest and resolve its per-shank streams, refusing any row that cannot be.

    For runs with no coverage record to carry a dropped row; :func:`stage_discover`
    drops and records them instead.

    Raises
    ------
    StreamResolutionError
        Naming every row whose per-shank streams do not identify one stream.
    """
    manifest_df = pd.read_csv(manifest_csv)
    if skip_ephys:
        return manifest_df
    resolved = resolve_manifest_streams(manifest_df, default_locator(data_root))
    resolved.raise_errors()
    return resolved.df


def resolve_manifest_streams(
    manifest_df: pd.DataFrame,
    locate: Callable[[ManifestRow], Path | None],
) -> ResolvedManifest:
    """Resolve every row naming a probe recorded as per-shank streams.

    Rows whose sorted asset is not mounted, or whose collection matches a stream
    exactly, pass through unchanged, so applying this to an already-resolved
    manifest is a no-op.

    Parameters
    ----------
    manifest_df : pandas.DataFrame
        The manifest as read.
    locate : callable
        Maps a row to its sorted asset directory, or ``None`` when not mounted
        (see :func:`default_locator`).

    Returns
    -------
    ResolvedManifest
        The rewritten manifest, its rows, and per-row resolutions and errors.
    """
    df = manifest_df.reset_index(drop=True).copy()
    rows = [ManifestRow.from_series(row) for _, row in df.iterrows()]
    streams_by_folder: dict[Path, list[SortedStream]] = {}
    resolutions: dict[int, StreamResolution] = {}
    errors: dict[int, str] = {}
    for position, mr in enumerate(rows):
        if mr.ephys_collection is None:
            continue
        folder = locate(mr)
        if folder is None:
            continue
        if folder not in streams_by_folder:
            streams_by_folder[folder] = sorted_streams(folder)
        try:
            resolution = resolve_split_stream(streams_by_folder[folder], str(mr.ephys_collection), mr.histology_shank)
        except StreamResolutionError as error:
            errors[position] = str(error)
            continue
        if resolution is None:
            continue
        resolutions[position] = resolution
        rows[position] = resolved_row(mr, resolution)
        logger.info("[streams] %s/%s: %s", mr.recording_id, mr.probe_id, resolution.note())
    if resolutions:
        # Explicit columns take precedence over the legacy aliases in ManifestRow.from_series, and every
        # row's derived values round-trip, so writing them for all rows changes only the resolved ones.
        df["ephys_collection"] = [r.ephys_collection for r in rows]
        df["logical_probe"] = [r.logical_probe for r in rows]
        for column in ("histology_shank", "ephys_shank"):
            df[column] = pd.array([getattr(r, column) for r in rows], dtype="Int64")
    return ResolvedManifest(df=df, rows=rows, resolutions=resolutions, errors=errors)
