"""Datapackage manifest for preprocessing pipeline outputs."""

from __future__ import annotations

import logging
import os
from collections import defaultdict
from pathlib import Path

from pydantic import BaseModel

from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    ManifestRow,
    OutputDirs,
    PipelineConfig,
    ProcessResult,
)

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "2.1.0"
# Why 2.1.0: 2.0.0 fixed cross-recording probe-name collisions by nesting
# probes under recording_id. 2.1.0 keeps that shape but makes the inner key
# explicitly be the ephys ALF collection, and adds logical_probe plus separate
# histology/ephys shank fields so split-stream quadbase probes and single-
# stream multi-shank probes are both representable without reinterpreting
# probe_name/probe_shank.


# ---------------------------------------------------------------------------
# Pydantic models mirroring the datapackage.json schema
# ---------------------------------------------------------------------------


class TransformPaths(BaseModel, frozen=True):
    """Paths to the ANTs transform chain files, relative to the manifest root.

    Paths are POSIX-style and may traverse up (``..``) to reach sibling assets
    (e.g. the SmartSPIM asset directory). Consumers resolve them against the
    directory containing ``datapackage.json``.
    """

    image_to_template_affine: str
    image_to_template_warp: str
    template_to_ccf_affine: str
    template_to_ccf_warp: str


class ImageSpaceHistology(BaseModel, frozen=True):
    """Paths to image-space histology volumes (relative to manifest)."""

    registration: str
    registration_pipeline: str
    ccf_template: str
    labels: str
    additional_channels: list[str] = []


class CcfSpaceHistology(BaseModel, frozen=True):
    """Paths to CCF-space histology volumes (relative to manifest)."""

    registration: str
    additional_channels: list[str] = []


class HistologyPaths(BaseModel, frozen=True):
    """Histology output paths grouped by coordinate space."""

    image_space: ImageSpaceHistology
    ccf_space: CcfSpaceHistology


class XyzPicks(BaseModel, frozen=True):
    """Paths to xyz-picks JSON files for one shank (or whole probe)."""

    ccf: str
    image_space: str
    histology_track_id: str | None = None
    histology_shank: int | None = None
    ephys_shank: int | None = None
    shank: int | None = None


class ChannelTablePaths(BaseModel, frozen=True):
    """Paths to the producer-owned ephys channel geometry table."""

    local_coordinates: str
    raw_ind: str
    shank_ind: str


class ProbeEntry(BaseModel, frozen=True):
    """Manifest entry for one ephys collection within a recording session.

    Uniquely identified by the ``(recording_id, ephys_collection)`` pair from
    the parent dict path. ``logical_probe`` may be shared by multiple ephys
    collections when the acquisition split one physical probe into streams.
    """

    probe_id: str
    logical_probe: str | None = None
    ephys_collection: str | None = None
    num_shanks: int
    ephys: str | None = None
    channel_table: ChannelTablePaths | None = None
    xyz_picks: list[XyzPicks]


class DataPackageError(Exception):
    """Raised when a datapackage references paths that are missing on disk."""


# File the GUI loads first from each ephys collection dir (load_data_local.py).
# Its presence is the load-bearing signal that ``ephys`` points at a real ALF
# collection rather than an empty/wrong directory.
_REQUIRED_EPHYS_FILE = "channels.localCoordinates.npy"


class DataPackage(BaseModel, frozen=True):
    """Top-level datapackage manifest for preprocessing outputs."""

    schema_version: str
    mouse_id: str
    transforms: TransformPaths
    histology: HistologyPaths
    # Nested ``recording_id -> ephys_collection -> ProbeEntry``. Per-shank
    # rows collapse into ``ProbeEntry.xyz_picks``.
    probes: dict[str, dict[str, ProbeEntry]]

    def referenced_files(self) -> list[str]:
        """Every relative *file* path this datapackage points at.

        Excludes ``ephys`` entries, which are directories (see
        :meth:`referenced_ephys_dirs`).
        """
        paths: list[str] = [
            self.transforms.image_to_template_affine,
            self.transforms.image_to_template_warp,
            self.transforms.template_to_ccf_affine,
            self.transforms.template_to_ccf_warp,
            self.histology.image_space.registration,
            self.histology.image_space.registration_pipeline,
            self.histology.image_space.ccf_template,
            self.histology.image_space.labels,
            *self.histology.image_space.additional_channels,
            self.histology.ccf_space.registration,
            *self.histology.ccf_space.additional_channels,
        ]
        for recording in self.probes.values():
            for probe in recording.values():
                if probe.channel_table is not None:
                    paths.append(probe.channel_table.local_coordinates)
                    paths.append(probe.channel_table.raw_ind)
                    paths.append(probe.channel_table.shank_ind)
                for xp in probe.xyz_picks:
                    paths.append(xp.ccf)
                    paths.append(xp.image_space)
        return paths

    def referenced_ephys_dirs(self) -> list[str]:
        """Relative ephys collection directories (one per probe with ephys)."""
        return [
            probe.ephys for recording in self.probes.values() for probe in recording.values() if probe.ephys is not None
        ]

    def missing_paths(self, root: Path) -> list[str]:
        """Return referenced paths that do not exist under *root*.

        *root* is the directory that will hold ``datapackage.json``; all stored
        paths are relative to it. Each ``ephys`` entry must be a directory
        containing ``channels.localCoordinates.npy`` (what the GUI loads first),
        so a wrong ephys directory is reported rather than passing silently.
        """
        root = Path(root)
        missing: list[str] = []
        for rel in self.referenced_files():
            if not (root / rel).exists():
                missing.append(rel)
        for rel in self.referenced_ephys_dirs():
            ephys = root / rel
            if not ephys.is_dir():
                missing.append(f"{rel}/ (ephys directory)")
            elif not (ephys / _REQUIRED_EPHYS_FILE).is_file():
                missing.append(f"{rel}/{_REQUIRED_EPHYS_FILE}")
        return sorted(missing)

    def validate_paths(self, root: Path) -> None:
        """Raise :class:`DataPackageError` if any referenced path is missing.

        Parameters
        ----------
        root : Path
            Directory the paths are relative to (where ``datapackage.json``
            will be written).
        """
        missing = self.missing_paths(root)
        if missing:
            listing = "\n  ".join(missing)
            raise DataPackageError(
                f"datapackage for mouse {self.mouse_id!r} references "
                f"{len(missing)} path(s) not found under {root}:\n  {listing}"
            )


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


def build_datapackage(
    mouse_id: str,
    manifest_rows: list[ManifestRow],
    results: list[ProcessResult],
    asset_info: AssetInfo,
    outputs: OutputDirs,
    config: PipelineConfig,
) -> DataPackage:
    """Construct a :class:`DataPackage` from pipeline state.

    Parameters
    ----------
    mouse_id : str
        Mouse identifier.
    manifest_rows : list[ManifestRow]
        All rows from the input manifest CSV.
    results : list[ProcessResult]
        Per-probe processing outcomes (parallel to *manifest_rows*).
    asset_info : AssetInfo
        Discovered SmartSPIM asset metadata.
    outputs : OutputDirs
        Output directory tree.
    config : PipelineConfig
        Pipeline configuration (for transform paths and flags).
    """
    # The manifest lives alongside the mouse output root.
    manifest_root = outputs.histology_img.parent

    transforms = _build_transforms(asset_info, config, manifest_root)
    histology = _build_histology(outputs, manifest_root)
    probes = _build_probes(manifest_rows, results, outputs, manifest_root, config)

    return DataPackage(
        schema_version=SCHEMA_VERSION,
        mouse_id=mouse_id,
        transforms=transforms,
        histology=histology,
        probes=probes,
    )


def _build_transforms(asset_info: AssetInfo, config: PipelineConfig, manifest_root: Path) -> TransformPaths:
    reg_dir = asset_info.registration_dir_path
    tmpl_dir = config.template_to_ccf_dir
    return TransformPaths(
        image_to_template_affine=_rel_up(reg_dir / "ls_to_template_SyN_0GenericAffine.mat", manifest_root),
        image_to_template_warp=_rel_up(reg_dir / "ls_to_template_SyN_1InverseWarp.nii.gz", manifest_root),
        template_to_ccf_affine=_rel_up(tmpl_dir / "syn_0GenericAffine.mat", manifest_root),
        template_to_ccf_warp=_rel_up(tmpl_dir / "syn_1InverseWarp.nii.gz", manifest_root),
    )


def _rel(path: Path, root: Path) -> str:
    """Return *path* relative to *root* as a POSIX string (no ``..``)."""
    return path.relative_to(root).as_posix()


def _rel_up(path: Path, root: Path) -> str:
    """Return *path* relative to *root*, allowing ``..`` traversal.

    Used for sibling-asset references (e.g. the SmartSPIM asset that lives
    next to the mouse output root rather than inside it).
    """
    return Path(os.path.relpath(path, root)).as_posix()


def _build_histology(outputs: OutputDirs, manifest_root: Path) -> HistologyPaths:
    img_dir = outputs.histology_img
    ccf_dir = outputs.histology_ccf

    # Image-space additional channels (Ex_\d+_Em_\d+.nrrd)
    img_additional = sorted(_rel(p, manifest_root) for p in img_dir.glob("Ex_*_Em_*.nrrd"))

    # CCF-space additional channels (histology_*.nrrd excluding registration)
    ccf_additional = sorted(
        _rel(p, manifest_root) for p in ccf_dir.glob("histology_*.nrrd") if p.name != "histology_registration.nrrd"
    )

    return HistologyPaths(
        image_space=ImageSpaceHistology(
            registration=_rel(img_dir / "histology_registration.nrrd", manifest_root),
            registration_pipeline=_rel(img_dir / "histology_registration_pipeline.nrrd", manifest_root),
            ccf_template=_rel(img_dir / "ccf_in_mouse.nrrd", manifest_root),
            labels=_rel(img_dir / "labels_in_mouse.nrrd", manifest_root),
            additional_channels=img_additional,
        ),
        ccf_space=CcfSpaceHistology(
            registration=_rel(ccf_dir / "histology_registration.nrrd", manifest_root),
            additional_channels=ccf_additional,
        ),
    )


def _build_probes(
    manifest_rows: list[ManifestRow],
    results: list[ProcessResult],
    outputs: OutputDirs,
    manifest_root: Path,
    config: PipelineConfig,
) -> dict[str, dict[str, ProbeEntry]]:
    # Build lookup of successful probe_ids
    successful = {r.probe_id for r in results if r.wrote_files}

    # Group rows by (recording_id, ephys_collection) — the ephys collection is
    # the ALF output folder created by aind-ephys-ibl-gui-conversion. Rows that
    # differ only in ephys_shank collapse into one entry's xyz_picks list. The
    # same collection label under two recording_ids stays distinct.
    groups: dict[tuple[str, str], list[ManifestRow]] = defaultdict(list)
    for row in manifest_rows:
        if _row_histology_track_id(row) in successful:
            groups[(row.recording_id, _row_ephys_collection(row))].append(row)

    probes: dict[str, dict[str, ProbeEntry]] = {}
    for (recording_id, ephys_collection), rows in groups.items():
        logical_probes = {_row_logical_probe(r) for r in rows}
        if len(logical_probes) != 1:
            raise ValueError(
                "Rows for one ephys_collection must share logical_probe: "
                f"{recording_id}/{ephys_collection} has {sorted(logical_probes)}"
            )
        ephys_shanks = {s for s in (_row_ephys_shank(r) for r in rows) if s is not None}
        num_shanks = len(ephys_shanks) if ephys_shanks else 1

        # Build xyz_picks list
        xyz_picks_list: list[XyzPicks] = []
        for row in rows:
            gui_folder = row.gui_folder(outputs)
            histology_shank = _row_histology_shank(row)
            ephys_shank = _row_ephys_shank(row)
            file_shank = histology_shank if histology_shank is not None else ephys_shank
            if file_shank is None:
                gui_ccf = "xyz_picks.json"
                gui_img = "xyz_picks_image_space.json"
                shank = None
            else:
                shank_id = int(file_shank) + 1
                gui_ccf = f"xyz_picks_shank{shank_id}.json"
                gui_img = f"xyz_picks_shank{shank_id}_image_space.json"
                shank = int(ephys_shank) + 1 if ephys_shank is not None else shank_id

            xyz_picks_list.append(
                XyzPicks(
                    ccf=_rel(gui_folder / gui_ccf, manifest_root),
                    image_space=_rel(gui_folder / gui_img, manifest_root),
                    histology_track_id=_row_histology_track_id(row),
                    histology_shank=histology_shank,
                    ephys_shank=ephys_shank,
                    shank=shank,
                )
            )

        # Ephys path (same for all shanks of an ephys collection). The
        # conversion writes the ALF collection (spikes/clusters/channels.*)
        # directly into gui_folder -- NOT a "spikes" subdirectory.
        first_row = rows[0]
        ephys_dir = first_row.gui_folder(outputs)
        ephys_path = _rel(ephys_dir, manifest_root) if not config.skip_ephys else None
        channel_table = (
            ChannelTablePaths(
                local_coordinates=_rel(
                    ephys_dir / "channels.localCoordinates.npy",
                    manifest_root,
                ),
                raw_ind=_rel(ephys_dir / "channels.rawInd.npy", manifest_root),
                shank_ind=_rel(ephys_dir / "channels.shankInd.npy", manifest_root),
            )
            if not config.skip_ephys
            else None
        )

        probes.setdefault(recording_id, {})[ephys_collection] = ProbeEntry(
            probe_id=_row_histology_track_id(first_row),
            logical_probe=_row_logical_probe(first_row),
            ephys_collection=ephys_collection,
            num_shanks=num_shanks,
            ephys=ephys_path,
            channel_table=channel_table,
            xyz_picks=xyz_picks_list,
        )

    return probes


def _row_histology_track_id(row: ManifestRow) -> str:
    return str(getattr(row, "histology_track_id", None) or row.probe_id)


def _row_logical_probe(row: ManifestRow) -> str:
    return str(getattr(row, "logical_probe", None) or getattr(row, "ephys_collection", None) or row.probe_name)


def _row_ephys_collection(row: ManifestRow) -> str:
    return str(getattr(row, "ephys_collection", None) or row.probe_name)


def _row_histology_shank(row: ManifestRow) -> int | None:
    value = getattr(row, "histology_shank", None)
    if value is None:
        value = getattr(row, "probe_shank", None)
    return int(value) if value is not None else None


def _row_ephys_shank(row: ManifestRow) -> int | None:
    value = getattr(row, "ephys_shank", None)
    if value is None:
        value = getattr(row, "probe_shank", None)
    return int(value) if value is not None else None


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def write_datapackage(dp: DataPackage, output_dir: Path, *, validate: bool = True) -> Path:
    """Write *dp* as ``datapackage.json`` in *output_dir*.

    Parameters
    ----------
    dp : DataPackage
        The datapackage to serialize.
    output_dir : Path
        Directory to write ``datapackage.json`` into. All paths in *dp* are
        relative to this directory.
    validate : bool
        When True (default), verify every path *dp* references exists under
        *output_dir* before writing, raising :class:`DataPackageError` on any
        miss. This turns a well-typed but dangling manifest (e.g. an ephys dir
        pointing at a nonexistent subdirectory) into a loud failure at write
        time instead of a runtime error in the GUI. Set False to serialize
        without touching the filesystem (round-trip tests, fixtures).

    Returns
    -------
    Path
        The path of the written file.
    """
    if validate:
        dp.validate_paths(output_dir)
    path = output_dir / "datapackage.json"
    path.write_text(dp.model_dump_json(indent=2, exclude_none=True))
    return path


def load_datapackage(path: Path) -> DataPackage:
    """Load and validate a ``datapackage.json`` file.

    Parameters
    ----------
    path : Path
        Path to the JSON file.

    Returns
    -------
    DataPackage
        Validated manifest.
    """
    return DataPackage.model_validate_json(path.read_text())
