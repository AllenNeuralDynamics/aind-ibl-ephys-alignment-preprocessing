"""What the manifest asked for versus what the run actually covered.

``discover`` is the only place that holds both sides of that comparison at once:
the manifest rows it read, and the fan-out units it chose to emit. Everything
downstream sees only the *surviving* work -- the filtered ``manifest.csv`` and
one ``stream_*/config.json`` per viable unit -- so a row dropped here leaves no
trace, and every later coverage check agrees with a set that has already shrunk.
A run that processes half a manifest then looks exactly like a run whose
manifest asked for half as much.

That is not hypothetical: 791094's ``2025-10-09`` session was silently dropped
from its 2026-08-11 preprocessing run. The trigger exited 0, the datapackage was
schema-valid over the surviving half, and nothing anywhere recorded that ten
probes had gone missing.

This module is the record that closes it. :func:`build_coverage` diffs the rows
read against the units emitted and keeps, for every dropped row, the reason
``_probe_viability`` gave -- reasons that previously existed only as log lines,
which is why diagnosing 791094 meant fetching a run log instead of reading a
file.

**Row granularity, not unit granularity.** A fan-out unit is
``(recording_id, ephys_collection)`` *deduplicated*, so a multi-shank probe
contributes several rows to one unit. If a single shank fails viability, a
surviving sibling row still emits the unit and a unit-level diff reports full
coverage while the filtered manifest is short one track. Rows are therefore
keyed by ``(recording_id, ephys_collection, ephys_shank)`` -- the same
producer-consumer identity the datapackage uses -- so a lost shank is visible.

**This is a run record, not part of the datapackage.** ``datapackage.json`` is
the alignment GUI's *consumption* contract: what to read and where. Coverage and
asset provenance answer questions about the run that produced it, which no
consumer reads -- folding them in would move a schema version, and therefore
every consumer's version gate, for reasons unrelated to anything they consume.
It ships as its own ``coverage.json`` beside the datapackage instead.

The file is written twice, progressively filled, because its two halves are
known in different places:

- ``discover`` writes the rows. It is the only stage holding both the manifest
  it read and the units it emitted.
- ``pack`` copies it forward and adds the resolved Code Ocean asset ids. Those
  are resolved by the trigger before any capsule launches and are invisible from
  inside one, so they arrive as a parameter. This second write is also what makes
  the record durable at all: discover's copy lives in an intermediate asset the
  trigger deletes on success.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow

logger = logging.getLogger(__name__)

#: Written to the root of ``discover``'s results beside the filtered manifest,
#: then again to ``pack``'s results beside the ``<mouse_id>/`` output tree.
COVERAGE_FILENAME = "coverage.json"

CoverageVersion = Literal["1.0.0"]
COVERAGE_VERSION: CoverageVersion = "1.0.0"


class CoverageError(RuntimeError):
    """Raised when a run would cover less than its manifest asked for.

    Deliberately a hard failure rather than a warning: the failure mode this
    guards against is not a crash but a *plausible success*, so the default has
    to be refusal or the silence recurs.
    """


class RowCoverage(BaseModel, frozen=True):
    """One manifest row's fate: judged viable and emitted, or dropped and why."""

    row_index: int | None = None
    recording_id: str
    ephys_collection: str | None = None
    ephys_shank: int | None = None
    histology_shank: int | None = None
    histology_track_id: str | None = None
    probe_file: str | None = None
    sorted_recording: str | None = None
    viable: bool
    #: Why the row was dropped; ``None`` when it was kept.
    reason: str | None = None
    #: Fan-out unit this row was folded into, when one was emitted for it.
    unit: str | None = None
    #: Non-fatal oddities worth recording -- e.g. a sorting matched by directory
    #: name because no mounted asset named the recording in its
    #: ``data_description``, which is exactly the 791094 signature.
    notes: list[str] = Field(default_factory=list)

    @property
    def key(self) -> str:
        """Canonical ``(recording, collection, shank)`` identity for this row."""
        shank = "-" if self.ephys_shank is None else str(self.ephys_shank)
        return f"{self.recording_id}/{self.ephys_collection or 'all'}/shank{shank}"

    def gap_message(self) -> str:
        """One line naming this row and why it was dropped."""
        track = f" track {self.histology_track_id}" if self.histology_track_id else ""
        return f"{self.key}{track}: {self.reason or 'dropped without a reason'}"


class AssetRef(BaseModel, frozen=True):
    """A Code Ocean data asset, by id and name.

    Both are optional because a resolution can be partial, and because the id is
    the durable half: asset *names* drift (published captures have shipped with
    the recording time missing), which is exactly why recording the id matters.
    """

    id: str | None = None
    name: str | None = None


class UnitAssets(BaseModel, frozen=True):
    """The assets one ephys stream was actually built from."""

    raw: AssetRef | None = None
    sorted: AssetRef | None = None


class RunCoverage(BaseModel, frozen=True):
    """What the manifest asked for, what ran, and what it ran on.

    Written by ``discover`` with the rows filled in, then rewritten by ``pack``
    with the asset ids added. One shape for both, so there is one file to find
    and one schema to read whichever stage's output you are holding.
    """

    coverage_version: CoverageVersion = COVERAGE_VERSION
    mouse_id: str
    #: The manifest ``discover`` read, as it was given (not the filtered copy).
    manifest: str
    skip_ephys: bool = False
    rows_requested: int
    rows_viable: int
    units_emitted: list[str] = Field(default_factory=list)
    rows: list[RowCoverage] = Field(default_factory=list)
    #: ``{fan-out unit name: the raw + sorted assets it was built from}``. Empty
    #: until ``pack`` merges the launcher's resolution in.
    unit_assets: dict[str, UnitAssets] = Field(default_factory=dict)
    #: The SmartSPIM asset the histology came from.
    smartspim: AssetRef | None = None

    @property
    def dropped(self) -> list[RowCoverage]:
        """Rows the viability gate refused."""
        return [row for row in self.rows if not row.viable]

    @property
    def is_complete(self) -> bool:
        """Whether the manifest asked for work and every requested row survived.

        A record with no rows at all is *not* complete. Equality alone reports an
        empty manifest as fully covered, which is true only vacuously and reads as
        reassurance about a run that did nothing.
        """
        return self.rows_requested > 0 and self.rows_requested == self.rows_viable

    def gaps(self) -> list[str]:
        """One human-readable line per dropped row, for a coverage refusal."""
        return [row.gap_message() for row in self.dropped]

    def with_assets(self, asset_resolution: str | Mapping[str, Any] | None) -> RunCoverage:
        """Return a copy carrying the launcher's resolved asset ids.

        *asset_resolution* is a JSON string or mapping shaped
        ``{"units": {<unit>: {"raw": {...}, "sorted": {...}}}, "smartspim": {...}}``.
        Blank or malformed input is warned about and ignored: provenance is worth
        recording, never worth failing a run that already did the work.
        """
        resolution = _parse_asset_resolution(asset_resolution)
        if not resolution:
            return self

        unit_assets = {
            str(unit): UnitAssets.model_validate(assets)
            for unit, assets in (resolution.get("units") or {}).items()
            if isinstance(assets, Mapping)
        }
        unknown = set(unit_assets) - set(self.units_emitted)
        if unknown:
            # The join key is the unit name discover emitted, not the ``stream_``
            # prefixed directory it wrote the config into. A mismatch yields a
            # record that validates and answers nothing, so say so loudly.
            logger.warning(
                "asset resolution names %d unit(s) discover never emitted: %s",
                len(unknown),
                ", ".join(sorted(unknown)),
            )
        smartspim = resolution.get("smartspim")
        return self.model_copy(
            update={
                "unit_assets": unit_assets,
                "smartspim": AssetRef.model_validate(smartspim) if isinstance(smartspim, Mapping) else None,
            }
        )

    def write(self, results_root: Path) -> Path:
        """Write the record to ``<results_root>/coverage.json`` and return its path."""
        path = Path(results_root) / COVERAGE_FILENAME
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2) + "\n")
        return path

    @classmethod
    def read(cls, path: Path) -> RunCoverage:
        """Load a coverage record from *path*."""
        return cls.model_validate(json.loads(Path(path).read_text()))

    @classmethod
    def find(cls, root: Path) -> RunCoverage | None:
        """Return the coverage record under *root*, or ``None`` if there is none.

        Looks beside *root* itself first (the producer's own results) and then
        one level down, which is where a captured ``discover`` asset lands when
        it is mounted under ``/data``. Returns ``None`` rather than raising: a
        monolith run, or a rerun against an older ``discover`` asset, simply has
        no record, and a missing record must not fail a run that already
        succeeded.

        Several candidates one level down is ambiguous and warned about, not
        silently resolved: only ``discover`` writes this file, but a rerun can
        mount a *previous* ``pack`` output whose root holds one too, and picking
        the alphabetically first would carry a stale row set forward as though it
        described this run.
        """
        root = Path(root)
        beside = root / COVERAGE_FILENAME
        nested = sorted(root.glob(f"*/{COVERAGE_FILENAME}"))
        if not beside.is_file() and len(nested) > 1:
            logger.warning(
                "%d coverage records under %s (%s); using %s -- the others are ignored",
                len(nested),
                root,
                ", ".join(str(path.parent.name) for path in nested),
                nested[0],
            )
        for candidate in (beside, *nested):
            if candidate.is_file():
                try:
                    return cls.read(candidate)
                except Exception as exc:  # a bad record must never fail a pack
                    logger.warning("ignoring unreadable coverage record %s: %s", candidate, exc)
        return None


def _parse_asset_resolution(asset_resolution: str | Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize the launcher's resolution parameter to a mapping (empty when absent)."""
    if asset_resolution is None:
        return {}
    if isinstance(asset_resolution, Mapping):
        return dict(asset_resolution)
    text = str(asset_resolution).strip()
    if not text:
        return {}
    try:
        parsed = json.loads(text)
    except ValueError as exc:
        logger.warning("ignoring unparsable asset resolution: %s", exc)
        return {}
    if not isinstance(parsed, dict):
        logger.warning("ignoring asset resolution: expected an object, got %s", type(parsed).__name__)
        return {}
    return parsed


def build_coverage(
    *,
    mouse_id: str,
    manifest_csv: Path,
    judged: list[tuple[ManifestRow, bool, str | None, list[str]]],
    unit_by_row: dict[int, str],
    skip_ephys: bool,
) -> RunCoverage:
    """Assemble the coverage record from ``discover``'s per-row verdicts.

    Parameters
    ----------
    mouse_id : str
        Mouse this manifest describes.
    manifest_csv : Path
        The manifest that was read, recorded as given.
    judged : list of (ManifestRow, bool, str or None, list of str)
        One ``(row, viable, reason, notes)`` tuple per manifest row, in file
        order. The reason is what ``_probe_viability`` returned.
    unit_by_row : dict[int, str]
        ``{index into judged: fan-out unit name}`` for rows that were folded
        into an emitted unit. Rows absent from this map contributed no unit --
        either because they were dropped, or because ephys is off.
    skip_ephys : bool
        Whether this run emits ephys units at all; recorded so an empty
        ``units_emitted`` is not mistaken for a coverage failure.

    Returns
    -------
    RunCoverage
        The record, ready to write. Asset ids are added later by ``pack``.
    """
    rows = [
        RowCoverage(
            row_index=mr.row_index,
            recording_id=mr.recording_id,
            ephys_collection=mr.ephys_collection,
            ephys_shank=mr.ephys_shank,
            histology_shank=mr.histology_shank,
            histology_track_id=mr.histology_track_id,
            probe_file=mr.probe_file,
            sorted_recording=mr.sorted_recording,
            viable=ok,
            reason=reason if not ok else None,
            unit=unit_by_row.get(index),
            notes=list(notes),
        )
        for index, (mr, ok, reason, notes) in enumerate(judged)
    ]
    return RunCoverage(
        mouse_id=mouse_id,
        manifest=str(manifest_csv),
        skip_ephys=skip_ephys,
        rows_requested=len(rows),
        rows_viable=sum(1 for row in rows if row.viable),
        units_emitted=sorted({unit for unit in unit_by_row.values()}),
        rows=rows,
    )
