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

SCHEMA_VERSION = "1.1.0"


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
    """Manifest entry for a single probe."""

    probe_id: str
    recording_id: str
    num_shanks: int
    ephys: str | None = None
    xyz_picks: list[XyzPicks]


class DataPackage(BaseModel, frozen=True):
    """Top-level datapackage manifest for preprocessing outputs."""

    schema_version: str
    mouse_id: str
    transforms: TransformPaths
    histology: HistologyPaths
    probes: dict[str, ProbeEntry]


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
) -> dict[str, ProbeEntry]:
    # Build lookup of successful probe_ids
    successful = {r.probe_id for r in results if r.wrote_files}

    # Group rows by probe_name (multi-shank probes share probe_name)
    groups: dict[str, list[ManifestRow]] = defaultdict(list)
    for row in manifest_rows:
        if str(row.probe_id) in successful:
            groups[row.probe_name].append(row)

    probes: dict[str, ProbeEntry] = {}
    for probe_name, rows in groups.items():
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

        # Ephys path (same for all shanks of a probe)
        first_row = rows[0]
        ephys_dir = first_row.gui_folder(outputs) / "spikes"
        ephys_path = _rel(ephys_dir, manifest_root) if not config.skip_ephys else None

        probes[probe_name] = ProbeEntry(
            probe_id=first_row.probe_id,
            recording_id=first_row.recording_id,
            num_shanks=num_shanks,
            ephys=ephys_path,
            xyz_picks=xyz_picks_list,
        )

    return probes


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def write_datapackage(dp: DataPackage, output_dir: Path) -> Path:
    """Write *dp* as ``datapackage.json`` in *output_dir*.

    Returns
    -------
    Path
        The path of the written file.
    """
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
