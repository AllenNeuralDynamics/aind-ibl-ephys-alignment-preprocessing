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

SCHEMA_VERSION = "2.0.0"
# Why 2.0.0: probes were previously a flat ``dict[probe_name, ProbeEntry]``
# which silently merged rows from different recordings sharing a probe_name
# (e.g. the same physical probe re-inserted for a follow-up recording) into
# one entry, since ``probe_name`` alone is not a unique key. The unique key
# is the triplet ``(recording_id, probe_name, probe_shank)``. The fix nests
# probes by recording_id, then by probe_name, with shanks collapsed into the
# entry's ``xyz_picks`` list. The shape is incompatible with 1.x, so it's a
# major bump.


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
    shank: int | None = None


class ProbeEntry(BaseModel, frozen=True):
    """Manifest entry for a single probe within a recording session.

    Uniquely identified by the ``(recording_id, probe_name)`` pair from the
    parent dict path; ``recording_id`` and ``probe_name`` are therefore not
    repeated as fields here. Multi-shank probes collapse into a single entry
    with one ``XyzPicks`` per shank.
    """

    probe_id: str
    num_shanks: int
    ephys: str | None = None
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
    # Nested ``recording_id -> probe_name -> ProbeEntry``. Per-shank rows
    # collapse into ``ProbeEntry.xyz_picks``.
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

    # Group rows by (recording_id, probe_name) — the unique probe key.
    # Rows that differ only in probe_shank collapse into one entry's
    # xyz_picks list (multi-shank probe). The same probe_name appearing
    # under two recording_ids stays distinct.
    groups: dict[tuple[str, str], list[ManifestRow]] = defaultdict(list)
    for row in manifest_rows:
        if str(row.probe_id) in successful:
            groups[(row.recording_id, row.probe_name)].append(row)

    probes: dict[str, dict[str, ProbeEntry]] = {}
    for (recording_id, probe_name), rows in groups.items():
        has_shanks = any(r.probe_shank is not None for r in rows)
        num_shanks = len(rows) if has_shanks else 1

        # Build xyz_picks list
        xyz_picks_list: list[XyzPicks] = []
        for row in rows:
            gui_folder = row.gui_folder(outputs)
            if row.probe_shank is None:
                gui_ccf = "xyz_picks.json"
                gui_img = "xyz_picks_image_space.json"
                shank = None
            else:
                shank_id = int(row.probe_shank) + 1
                gui_ccf = f"xyz_picks_shank{shank_id}.json"
                gui_img = f"xyz_picks_shank{shank_id}_image_space.json"
                shank = shank_id

            xyz_picks_list.append(
                XyzPicks(
                    ccf=_rel(gui_folder / gui_ccf, manifest_root),
                    image_space=_rel(gui_folder / gui_img, manifest_root),
                    shank=shank,
                )
            )

        # Ephys path (same for all shanks of a probe). The conversion writes
        # the ALF collection (spikes/clusters/channels.*) directly into the
        # probe's gui_folder -- NOT a "spikes" subdirectory -- so ephys_dir is
        # the gui_folder itself, matching where xyz_picks.json also lives.
        first_row = rows[0]
        ephys_dir = first_row.gui_folder(outputs)
        ephys_path = _rel(ephys_dir, manifest_root) if not config.skip_ephys else None

        probes.setdefault(recording_id, {})[probe_name] = ProbeEntry(
            probe_id=first_row.probe_id,
            num_shanks=num_shanks,
            ephys=ephys_path,
            xyz_picks=xyz_picks_list,
        )

    return probes


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
