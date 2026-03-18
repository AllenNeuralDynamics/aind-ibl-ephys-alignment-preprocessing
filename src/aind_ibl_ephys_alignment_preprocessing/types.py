"""Frozen dataclasses and Pydantic configuration for the preprocessing pipeline."""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import ants
import pandas as pd
from pydantic import BaseModel, model_validator

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PipelineConfig — environment-agnostic Pydantic model
# ---------------------------------------------------------------------------


class PipelineConfig(BaseModel, frozen=True):
    """Environment-agnostic pipeline configuration.

    Parameters
    ----------
    data_root : Path
        Root directory containing input data.
    results_root : Path
        Root directory for pipeline outputs.
    scratch_root : Path | None
        Scratch directory for temporary files.  When *None* a temporary
        directory is created at run time.
    neuroglancer_file : Path
        Path to the Neuroglancer JSON file.  Resolved against *data_root*
        when relative.
    manifest_csv : Path
        Path to the manifest CSV.  Resolved against *data_root* when
        relative.
    template_25 : Path
        SmartSPIM LCA template volume (25 µm).
    ccf_25 : Path
        Allen CCF average template volume (25 µm).
    ccf_labels_lateralized_25 : Path
        Lateralized CCF annotation labels (25 µm, compacted).
    ibl_atlas_histology_path : Path
        IBL atlas histology directory.
    ccf_labels_lateralized_25_unq_vals : Path
        Unique label values for the lateralized annotation.
    skip_ephys : bool
        If *True*, skip ephys extraction.
    desired_voxel_size_um : float
        Target voxel size in micrometers for multiscale level selection.
    """

    # I/O roots
    data_root: Path
    results_root: Path
    scratch_root: Path | None = None

    # Input files (resolved against data_root if relative)
    neuroglancer_file: Path
    manifest_csv: Path

    # Reference volumes (resolved against data_root if relative)
    template_25: Path = Path("smartspim_lca_template/smartspim_lca_template_25.nii.gz")
    ccf_25: Path = Path("allen_mouse_ccf/average_template/average_template_25.nii.gz")
    ccf_labels_lateralized_25: Path = Path(
        "allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_compact.nrrd"
    )
    ibl_atlas_histology_path: Path = Path("iblatlas_allenatlas/")
    ccf_labels_lateralized_25_unq_vals: Path = Path(
        "allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_unique_vals.npz"
    )

    # External reference data
    template_to_ccf_dir: Path = Path("spim_template_to_ccf/")

    # Processing options
    skip_ephys: bool = False
    desired_voxel_size_um: float = 25.0

    @model_validator(mode="after")
    def _resolve_relative_paths(self) -> PipelineConfig:
        """Make relative reference/input paths absolute against *data_root*."""
        path_fields = (
            "neuroglancer_file",
            "manifest_csv",
            "template_25",
            "ccf_25",
            "ccf_labels_lateralized_25",
            "ibl_atlas_histology_path",
            "ccf_labels_lateralized_25_unq_vals",
            "template_to_ccf_dir",
        )
        for name in path_fields:
            val = getattr(self, name)
            if not val.is_absolute():
                object.__setattr__(self, name, self.data_root / val)
        return self


# ---------------------------------------------------------------------------
# Lightweight frozen dataclasses (moved from capsule ibl_preprocess_types.py)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InputPaths:
    """Resolved input paths for the pipeline."""

    neuroglancer_file: Path
    manifest_csv: Path
    data_root: Path
    results_root: Path


@dataclass(frozen=True)
class ReferencePaths:
    """Paths to reference data volumes.

    Default values are *relative*; use :meth:`for_data_root` to resolve them
    against a concrete root directory.
    """

    template_25: Path = Path("smartspim_lca_template/smartspim_lca_template_25.nii.gz")
    ccf_25: Path = Path("allen_mouse_ccf/average_template/average_template_25.nii.gz")
    ccf_labels_lateralized_25: Path = Path(
        "allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_compact.nrrd"
    )
    ibl_atlas_histology_path: Path = Path("iblatlas_allenatlas/")
    ccf_labels_lateralized_25_unq_vals: Path = Path(
        "allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_unique_vals.npz"
    )

    @classmethod
    def for_data_root(cls, data_root: Path) -> ReferencePaths:
        """Create a :class:`ReferencePaths` with all paths resolved under *data_root*."""
        base = cls()
        return cls(
            template_25=data_root / base.template_25,
            ccf_25=data_root / base.ccf_25,
            ccf_labels_lateralized_25=data_root / base.ccf_labels_lateralized_25,
            ibl_atlas_histology_path=data_root / base.ibl_atlas_histology_path,
            ccf_labels_lateralized_25_unq_vals=data_root / base.ccf_labels_lateralized_25_unq_vals,
        )

    @classmethod
    def from_config(cls, config: PipelineConfig) -> ReferencePaths:
        """Create from an already-resolved :class:`PipelineConfig`."""
        return cls(
            template_25=config.template_25,
            ccf_25=config.ccf_25,
            ccf_labels_lateralized_25=config.ccf_labels_lateralized_25,
            ibl_atlas_histology_path=config.ibl_atlas_histology_path,
            ccf_labels_lateralized_25_unq_vals=config.ccf_labels_lateralized_25_unq_vals,
        )


@dataclass(frozen=True)
class ReferenceVolumes:
    """Loaded reference ANTs images."""

    ccf_25: ants.ANTsImage

    @classmethod
    def from_paths(cls, paths: ReferencePaths) -> ReferenceVolumes:
        """Load volumes synchronously."""
        ccf = ants.image_read(str(paths.ccf_25), pixeltype=None)
        return cls(ccf_25=ccf)

    @classmethod
    async def from_paths_async(cls, paths: ReferencePaths) -> ReferenceVolumes:
        """Load volumes in a background thread."""
        ccf = await asyncio.to_thread(ants.image_read, str(paths.ccf_25), pixeltype=None)
        return cls(ccf_25=ccf)


@dataclass(frozen=True)
class ZarrPaths:
    """OME-Zarr channel metadata."""

    registration: str
    additional: list[str]
    metadata: dict[str, Any]
    processing: dict[str, Any]


@dataclass(frozen=True)
class RegistrationInfo:
    """Registration directory layout."""

    registration_root: Path
    prep_image_folder: Path
    moved_image_folder: Path
    alignment_channel: str


@dataclass(frozen=True)
class PipelineRegistrationInfo:
    """ANTs transform chain paths and inversion flags."""

    pt_tx_str: list[str]
    pt_tx_inverted: list[bool]
    img_tx_str: list[str]
    img_tx_inverted: list[bool]


@dataclass(frozen=True)
class AssetInfo:
    """SmartSPIM asset discovery results."""

    asset_path: Path
    zarr_volumes: ZarrPaths
    pipeline_registration_chains: PipelineRegistrationInfo
    registration_dir_path: Path
    registration_in_ccf_precomputed: Path


@dataclass(frozen=True)
class OutputDirs:
    """Output directory tree for a single mouse."""

    histology_ccf: Path
    histology_img: Path
    tracks_root: Path
    spim: Path
    template: Path
    ccf: Path
    bregma_xyz: Path


@dataclass(frozen=True)
class ProcessResult:
    """Per-probe processing outcome."""

    probe_id: str
    recording_id: str
    wrote_files: bool
    skipped_reason: str | None = None


@dataclass(frozen=True)
class ManifestRow:
    """A single row from the manifest CSV.

    Parameters
    ----------
    probe_id : str
        Probe identifier.
    probe_name : str
        Subfolder name for GUI artifacts.
    probe_file : str
        Basename of the annotation file (no extension).
    sorted_recording : str
        Name of the spike-sorting folder.
    mouseid : str
        Mouse identifier.
    annotation_format : str
        Annotation file format (default ``"json"``).
    probe_shank : int | None
        0-based shank index.
    surface_finding : Path | None
        Optional surface-finding file path fragment.
    row_index : int | None
        Row index from the CSV for provenance.
    """

    probe_id: str
    probe_name: str
    probe_file: str
    sorted_recording: str
    mouseid: str
    annotation_format: str = "json"
    probe_shank: int | None = None
    surface_finding: Path | None = None
    row_index: int | None = None

    @property
    def recording_id(self) -> str:
        """Derive recording ID by stripping ``_sorted`` suffix."""
        return self.sorted_recording.split("_sorted")[0]

    def gui_folder(self, outputs: OutputDirs) -> Path:
        """Per-recording GUI output folder."""
        return outputs.tracks_root.parent / self.recording_id / self.probe_name

    @classmethod
    def from_series(cls, s: pd.Series[Any]) -> ManifestRow:
        """Construct from a pandas Series (one manifest row)."""

        def opt_int(x: Any) -> int | None:
            try:
                return int(x) if pd.notna(x) else None
            except Exception:
                return None

        def opt_path(x: Any) -> Path | None:
            return Path(str(x)) if pd.notna(x) and str(x) else None

        return cls(
            probe_id=str(s.get("probe_id")),
            probe_name=str(s.get("probe_name")),
            probe_file=str(s.get("probe_file")),
            sorted_recording=str(s.get("sorted_recording")),
            mouseid=str(s.get("mouseid")),
            annotation_format=str(s.get("annotation_format", "json")).lower(),
            probe_shank=opt_int(s.get("probe_shank")),
            surface_finding=opt_path(s.get("surface_finding")),
            row_index=int(s.name) if hasattr(s, "name") else None,  # type: ignore[call-overload]
        )
