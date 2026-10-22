"""Build a manifest CSV from the one thing a person has to decide: which track is which recording.

Track names follow no convention, so nothing is read from them. A mapping file pins each
session's sorted recording and assigns each annotated track to a stream of it; the remaining
columns are derived from the stream and the sort's own stream list, and the mapping is checked
for completeness against that list and the annotation file before anything is written.
"""

from __future__ import annotations

import json
import re
import tomllib
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from aind_ibl_ephys_alignment_preprocessing.ephys import analyzer_collection, find_session_dir

#: A split quad-base stream: one shank of a probe recorded as its own stream (``46101-3``).
_SPLIT_STREAM_RE = re.compile(r"^(?P<base>.+)-(?P<position>[1-9]\d*)$")

#: Columns written, in this order; ``registration_asset`` is added only when set.
OUTPUT_COLUMNS: tuple[str, ...] = (
    "mouseid",
    "sorted_recording",
    "surface_finding",
    "probe_file",
    "histology_track_id",
    "logical_probe",
    "ephys_collection",
    "histology_shank",
    "ephys_shank",
)

#: Fewest points that define a track.
MIN_TRACK_POINTS = 2


class ManifestBuildError(ValueError):
    """The mapping cannot produce a valid manifest; ``problems`` lists every reason."""

    def __init__(self, problems: list[str]) -> None:
        super().__init__("\n".join(problems))
        self.problems = problems


@dataclass(frozen=True)
class TrackAssignment:
    """One annotated track assigned to one stream.

    Parameters
    ----------
    track : str
        Neuroglancer layer name.
    stream : str
        Ephys collection, as parsed from the sort's stream names (``46101-1``, ``ProbeA``).
    shank : int | None
        0-based shank index within a single stream that records several shanks; ``None`` otherwise.
    """

    track: str
    stream: str
    shank: int | None = None


@dataclass(frozen=True)
class SessionMapping:
    """The tracks recorded in one sorted recording."""

    sorted_recording: str
    tracks: tuple[TrackAssignment, ...]
    surface_finding: str | None = None


@dataclass(frozen=True)
class MouseMapping:
    """Everything a person supplies for one mouse's manifest."""

    mouseid: str
    annotation: Path
    sessions: tuple[SessionMapping, ...]
    registration_asset: str | None = None

    @classmethod
    def from_toml(cls, path: Path) -> MouseMapping:
        """Read a mapping file; ``annotation`` resolves against the file's directory.

        A track's value is its stream, or ``{stream = ..., shank = N}`` for one shank of a
        single-stream multi-shank probe, where ``shank`` is a position: the first shank is 1.
        """
        data = tomllib.loads(path.read_text())
        sessions = tuple(_session_from_toml(s) for s in data.get("session", []))
        return cls(
            mouseid=str(data["mouseid"]),
            annotation=(path.parent / data["annotation"]).resolve(),
            sessions=sessions,
            registration_asset=data.get("registration_asset"),
        )


def _session_from_toml(raw: Mapping[str, Any]) -> SessionMapping:
    tracks = []
    for track, value in raw.get("tracks", {}).items():
        if isinstance(value, str):
            tracks.append(TrackAssignment(track, value))
        else:
            position = value.get("shank")
            tracks.append(TrackAssignment(track, str(value["stream"]), None if position is None else int(position) - 1))
    return SessionMapping(
        sorted_recording=str(raw["sorted_recording"]),
        tracks=tuple(tracks),
        surface_finding=raw.get("surface_finding"),
    )


def annotation_point_counts(ng_state: Mapping[str, Any]) -> dict[str, int]:
    """Count the point annotations in each annotation layer of a Neuroglancer state."""
    return {
        layer["name"]: sum(1 for a in layer.get("annotations", []) if a.get("type") == "point")
        for layer in ng_state.get("layers", [])
        if layer.get("type") == "annotation"
    }


def streams_in(analyzer_names: Iterable[str]) -> set[str]:
    """Return the spike-sorted collections named by ``postprocessed/`` entries, across experiments."""
    return {c for c in map(analyzer_collection, analyzer_names) if c is not None}


def local_streams(data_root: Path, sorted_recording: str) -> set[str] | None:
    """List a sort's collections from a mounted copy under *data_root*; ``None`` when not found."""
    folder = find_session_dir(data_root, sorted_recording)
    if folder is None or not (folder / "postprocessed").is_dir():
        return None
    return streams_in(p.name for p in (folder / "postprocessed").iterdir())


def s3_streams(sorted_recording: str, bucket: str = "aind-open-data") -> set[str]:
    """List a sort's collections from its ``postprocessed/`` prefix, one level deep, unsigned."""
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    client = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    names: list[str] = []
    pages = client.get_paginator("list_objects_v2").paginate(
        Bucket=bucket, Prefix=f"{sorted_recording}/postprocessed/", Delimiter="/"
    )
    for page in pages:
        names.extend(p["Prefix"].rstrip("/").rsplit("/", 1)[-1] for p in page.get("CommonPrefixes", []))
    return streams_in(names)


def _split_base(stream: str, session_streams: set[str]) -> str | None:
    """Return the probe a split stream belongs to, when the sort also holds a sibling shank."""
    match = _SPLIT_STREAM_RE.match(stream)
    if match is None:
        return None
    base = match.group("base")
    siblings = {s for s in session_streams if (m := _SPLIT_STREAM_RE.match(s)) and m.group("base") == base}
    return base if len(siblings) > 1 else None


def _session_rows(session: SessionMapping, session_streams: set[str], problems: list[str]) -> list[dict[str, Any]]:
    name = session.sorted_recording
    by_stream: dict[str, list[TrackAssignment]] = defaultdict(list)
    for t in session.tracks:
        by_stream[t.stream].append(t)
    for stream in sorted(session_streams - by_stream.keys()):
        problems.append(f"{name}: stream {stream} has no track")
    for stream in sorted(by_stream.keys() - session_streams):
        problems.append(f"{name}: stream {stream} is not in the sort")

    rows = []
    for stream, assigned in by_stream.items():
        base = _split_base(stream, session_streams)
        if base is not None:
            if len(assigned) > 1 or assigned[0].shank is not None:
                problems.append(f"{name}: split stream {stream} records one shank; give it one track and no shank")
                continue
            position = int(stream.rsplit("-", 1)[1])
            rows.append(_row(session, assigned[0], base, position - 1, 0))
        elif len(assigned) == 1 and assigned[0].shank is None:
            rows.append(_row(session, assigned[0], stream, None, None))
        else:
            shanks = [t.shank for t in assigned]
            if None in shanks or len(set(shanks)) != len(shanks):
                problems.append(f"{name}: tracks on stream {stream} need distinct shanks")
                continue
            rows.extend(_row(session, t, stream, t.shank, t.shank) for t in assigned)
    return rows


def _row(
    session: SessionMapping,
    track: TrackAssignment,
    logical_probe: str,
    histology_shank: int | None,
    ephys_shank: int | None,
) -> dict[str, Any]:
    return {
        "sorted_recording": session.sorted_recording,
        "surface_finding": session.surface_finding,
        "histology_track_id": track.track,
        "logical_probe": logical_probe,
        "ephys_collection": track.stream,
        "histology_shank": histology_shank,
        "ephys_shank": ephys_shank,
    }


def build_manifest(
    mapping: MouseMapping,
    streams: Mapping[str, set[str]],
    point_counts: Mapping[str, int],
) -> pd.DataFrame:
    """Derive and check the manifest for *mapping*.

    Parameters
    ----------
    mapping : MouseMapping
        The person's assignments.
    streams : Mapping[str, set[str]]
        Each sorted recording's spike-sorted collections.
    point_counts : Mapping[str, int]
        Points per annotation layer, from :func:`annotation_point_counts`.

    Returns
    -------
    pandas.DataFrame
        One row per track, in mapping order, with :data:`OUTPUT_COLUMNS`.

    Raises
    ------
    ManifestBuildError
        Listing every stream without a track, track without a stream or points, and reused track.
    """
    problems: list[str] = []
    uses = Counter(t.track for s in mapping.sessions for t in s.tracks)
    problems.extend(f"track {track} is assigned {n} times" for track, n in uses.items() if n > 1)
    for track in uses:
        n_points = point_counts.get(track)
        if n_points is None:
            problems.append(f"track {track} is not a layer of {mapping.annotation.name}")
        elif n_points < MIN_TRACK_POINTS:
            problems.append(f"track {track} has {n_points} point(s)")

    rows = []
    for session in mapping.sessions:
        if session.sorted_recording not in streams:
            problems.append(f"{session.sorted_recording}: no stream list")
            continue
        rows.extend(_session_rows(session, streams[session.sorted_recording], problems))
    if problems:
        raise ManifestBuildError(problems)

    order = {t.track: i for i, t in enumerate(t for s in mapping.sessions for t in s.tracks)}
    rows.sort(key=lambda r: order[r["histology_track_id"]])
    df = pd.DataFrame(rows).assign(mouseid=mapping.mouseid, probe_file=mapping.annotation.stem)
    columns: list[str] = list(OUTPUT_COLUMNS)
    if mapping.registration_asset:
        df["registration_asset"] = mapping.registration_asset
        columns.append("registration_asset")
    for shank in ("histology_shank", "ephys_shank"):
        df[shank] = df[shank].astype("Int64")
    return df[columns]


def unmapped_tracks(mapping: MouseMapping, point_counts: Mapping[str, int]) -> list[str]:
    """Annotation layers with points that no session assigns, for the person to confirm."""
    used = {t.track for s in mapping.sessions for t in s.tracks}
    return [name for name, n in point_counts.items() if n > 0 and name not in used]


def build_from_files(mapping_path: Path, data_root: Path | None = None) -> tuple[pd.DataFrame, list[str]]:
    """Build the manifest for a mapping file, listing streams under *data_root* or else from S3.

    Returns
    -------
    tuple[pandas.DataFrame, list[str]]
        The manifest, and the annotated layers it leaves out.
    """
    mapping = MouseMapping.from_toml(mapping_path)
    point_counts = annotation_point_counts(json.loads(mapping.annotation.read_text()))
    streams = {}
    for session in mapping.sessions:
        found = local_streams(data_root, session.sorted_recording) if data_root is not None else None
        streams[session.sorted_recording] = found if found is not None else s3_streams(session.sorted_recording)
    return build_manifest(mapping, streams, point_counts), unmapped_tracks(mapping, point_counts)
