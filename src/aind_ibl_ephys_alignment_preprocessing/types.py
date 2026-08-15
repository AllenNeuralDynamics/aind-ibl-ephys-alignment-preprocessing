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
    # Selects the multiscale level via `determine_desired_level`, which returns
    # the COARSEST level still at least as fine as this — i.e. it rounds toward
    # higher resolution, never below it. A 25.0 target (the atlas resolution)
    # therefore landed on SmartSPIM's 14.4 µm level, not its 28.8 µm one: 8x the
    # voxels, and that 8x multiplies through every volume the histology stage
    # touches (channel reads, both CCF warps, all six compressed NRRD writes,
    # and the peak memory holding them). 35.0 selects the ~30 µm tier on both
    # SmartSPIM base resolutions in use — 28.8 µm on a 1.8 µm base, 32 µm on a
    # 2.0 µm base — where 30.0 would fall back to 16 µm on the latter.
    # This is a display-resolution knob only: probe coordinates come from the
    # header-only anatomical stubs, which never read a multiscale level.
    desired_voxel_size_um: float = 35.0
    # Grid every image-space histology volume is resampled onto, isotropic, in µm.
    # Distinct from `desired_voxel_size_um` above, which only picks which stored
    # multiscale level to read: whether a volume *has* a level near the target is
    # an accident of when it was stitched, so reading alone leaves each mouse on
    # a different grid (14.4, 16/32 and 28.8/32 µm all occur across the current
    # eleven). That is a consistency-of-analysis problem before it is a speed one
    # — cross-animal comparisons end up differently sampled purely by stitching
    # vintage — and resampling is what removes it.
    #
    # Quoted in micrometres because that is the natural unit for these volumes,
    # but image spacing throughout the pipeline is MILLIMETRES: aind_zarr_utils
    # reads OME-Zarr with scale_unit="millimeter", so the pyramid's micrometre
    # scales arrive converted and the stored ANTs warps share that frame. The
    # conversion happens once, in resample_to_isotropic. Note the contrast with
    # `desired_voxel_size_um` above, which is compared against the raw OME-Zarr
    # metadata and so really is in micrometres — two knobs, two frames.
    #
    # A volume already coarser than this on every axis is left untouched; the
    # target is never used to manufacture resolution. Reaching a *common*
    # isotropic grid does mean interpolating up on a single axis for volumes
    # whose z is coarser than target (32 -> 30 is the worst case in hand, ~7%),
    # which is the price of every mouse sharing one grid.
    output_voxel_size_um: float = 30.0
    num_parallel_jobs: int = 4
    # QC/diagnostic outputs that the alignment GUI never reads: the Slicer FCSVs
    # (spim/template/ccf), the CCF/bregma xyz_picks, and the CCF-space histology
    # volumes. Off by default — producing them costs the ANTs point-warps and the
    # full-volume warps into CCF for no consumer. Turn on for QC/export runs.
    emit_qc: bool = False
    # Whether a run that covers less than its manifest asked for may proceed.
    # ``discover``'s viability gate drops a row with a warning and the run
    # finishes green over whatever is left, so the default has to be refusal --
    # see ``coverage.py`` for the session that was lost this way. Mirrors the
    # trigger capsule's ``--allow-partial`` so both layers enforce one policy.
    allow_partial: bool = False

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
    asset_uri: str | None
    zarr_volumes: ZarrPaths
    pipeline_registration_chains: PipelineRegistrationInfo
    registration_dir_path: Path


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
class ManifestColumn:
    """One column of the manifest CSV contract.

    Parameters
    ----------
    name : str
        Canonical column name.
    required : bool
        Whether a manifest without it (or any of its *aliases*) is invalid.
        Optional columns must never be demanded: every manifest predating a
        column has to keep working unchanged, which is the whole reason new
        columns are added optional.
    aliases : tuple[str, ...]
        Older names still accepted for the same field, most-preferred first.
        Satisfying any one of them satisfies the column.
    description : str
        What the column pins, for error messages and docs.
    """

    name: str
    required: bool
    aliases: tuple[str, ...] = ()
    description: str = ""

    @property
    def accepted_names(self) -> tuple[str, ...]:
        """Every column name that satisfies this field, canonical first."""
        return (self.name, *self.aliases)


#: The manifest CSV contract, in one place.
#:
#: Both the parser (:meth:`ManifestRow.from_series`) and the pre-flight
#: validator read this, so "which columns are required" cannot drift between
#: what a run accepts and what validation reports.
MANIFEST_COLUMNS: tuple[ManifestColumn, ...] = (
    ManifestColumn("mouseid", required=True, description="Mouse identifier; one per manifest."),
    ManifestColumn("sorted_recording", required=True, description="Spike-sorting folder name."),
    ManifestColumn("probe_file", required=True, description="Annotation file basename, no extension."),
    ManifestColumn(
        "histology_track_id",
        required=True,
        aliases=("probe_id",),
        description="Neuroglancer layer / histology track identifier.",
    ),
    ManifestColumn(
        "ephys_collection",
        required=True,
        aliases=("probe_name",),
        description="ALF/ephys output collection, from the Open Ephys stream.",
    ),
    ManifestColumn("annotation_format", required=False, description="Annotation format; defaults to json."),
    ManifestColumn(
        "logical_probe",
        required=False,
        description="Physical probe identity; split-stream collections share one.",
    ),
    ManifestColumn(
        "histology_shank",
        required=False,
        aliases=("probe_shank",),
        description="0-based physical/histology shank index.",
    ),
    ManifestColumn(
        "ephys_shank",
        required=False,
        aliases=("probe_shank",),
        description="0-based shank index within the ephys collection.",
    ),
    ManifestColumn(
        "surface_finding",
        required=False,
        description="Separate surface-finding recording to merge blocks from.",
    ),
    ManifestColumn(
        "registration_asset",
        required=False,
        description=(
            "Path to a registration directory holding the ls_to_template_SyN_* "
            "transforms to use instead of the stitched asset's own. Per-brain."
        ),
    ),
)

#: Columns whose absence invalidates a manifest, canonical name first.
REQUIRED_MANIFEST_COLUMNS: tuple[ManifestColumn, ...] = tuple(c for c in MANIFEST_COLUMNS if c.required)

#: Columns a manifest may omit entirely. Adding one here must never break a
#: manifest written before it existed.
OPTIONAL_MANIFEST_COLUMNS: tuple[ManifestColumn, ...] = tuple(c for c in MANIFEST_COLUMNS if not c.required)


@dataclass(frozen=True)
class ManifestRow:
    """A single row from the manifest CSV.

    Parameters
    ----------
    probe_id : str
        Legacy alias for ``histology_track_id``.
    probe_name : str
        Legacy alias for ``ephys_collection``.
    probe_file : str
        Basename of the annotation file (no extension).
    sorted_recording : str
        Name of the spike-sorting folder.
    mouseid : str
        Mouse identifier.
    annotation_format : str
        Annotation file format (default ``"json"``).
    probe_shank : int | None
        Legacy alias for both ``histology_shank`` and ``ephys_shank``.
    histology_track_id : str | None
        Neuroglancer layer / histology track identifier.
    logical_probe : str | None
        Physical/logical probe identity. Multiple ephys collections may share
        this value for split-stream probes.
    ephys_collection : str | None
        ALF/ephys output collection folder, derived from the Open Ephys stream.
    histology_shank : int | None
        0-based physical/histology shank index.
    ephys_shank : int | None
        0-based shank index within ``ephys_collection``.
    surface_finding : Path | None
        Optional surface-finding file path fragment.
    registration_asset : Path | None
        Optional path to a registration directory, relative to ``data_root``,
        holding the ``ls_to_template_SyN_*`` transforms to use instead of the
        stitched asset's own (e.g. ``SmartSPIM_750108_reg/ccf_Ex_639_Em_667``).
        Per-*brain*, not per-probe: like ``mouseid`` it is replicated across
        rows, and discovery reads the first non-empty value. Empty leaves the
        stitched asset's registration in use.
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
    histology_track_id: str | None = None
    logical_probe: str | None = None
    ephys_collection: str | None = None
    histology_shank: int | None = None
    ephys_shank: int | None = None
    surface_finding: Path | None = None
    registration_asset: Path | None = None
    row_index: int | None = None

    def __post_init__(self) -> None:
        """Populate new explicit fields from legacy manifest aliases."""
        histology_track_id = self.histology_track_id or self.probe_id
        ephys_collection = self.ephys_collection or self.probe_name
        logical_probe = self.logical_probe or ephys_collection
        histology_shank = self.histology_shank if self.histology_shank is not None else self.probe_shank
        ephys_shank = (
            self.ephys_shank
            if self.ephys_shank is not None
            else (self.probe_shank if self.probe_shank is not None else histology_shank)
        )

        object.__setattr__(self, "histology_track_id", histology_track_id)
        object.__setattr__(self, "ephys_collection", ephys_collection)
        object.__setattr__(self, "logical_probe", logical_probe)
        object.__setattr__(self, "histology_shank", histology_shank)
        object.__setattr__(self, "ephys_shank", ephys_shank)

    @property
    def recording_id(self) -> str:
        """Derive recording ID by stripping ``_sorted`` suffix."""
        return self.sorted_recording.split("_sorted")[0]

    def gui_folder(self, outputs: OutputDirs) -> Path:
        """Per-recording GUI output folder."""
        return outputs.tracks_root.parent / self.recording_id / str(self.ephys_collection)

    @classmethod
    def from_series(cls, s: pd.Series[Any]) -> ManifestRow:
        """Construct from a pandas Series (one manifest row)."""

        def opt_int(x: Any) -> int | None:
            try:
                return int(x) if pd.notna(x) else None
            except Exception:
                return None

        def opt_int_from(*names: str) -> int | None:
            for name in names:
                value = opt_int(s.get(name))
                if value is not None:
                    return value
            return None

        def opt_str_from(*names: str, default: str | None = None) -> str:
            for name in names:
                value = s.get(name)
                if pd.notna(value) and str(value) != "":
                    return str(value)
            if default is not None:
                return default
            return ""

        def opt_path(x: Any) -> Path | None:
            return Path(str(x)) if pd.notna(x) and str(x) else None

        histology_track_id = opt_str_from("histology_track_id", "probe_id")
        ephys_collection = opt_str_from("ephys_collection", "probe_name")
        logical_probe = opt_str_from("logical_probe", default=ephys_collection)
        histology_shank = opt_int_from("histology_shank", "probe_shank")
        ephys_shank = opt_int_from("ephys_shank", "probe_shank", "histology_shank")

        return cls(
            probe_id=histology_track_id,
            probe_name=ephys_collection,
            probe_file=str(s.get("probe_file")),
            sorted_recording=str(s.get("sorted_recording")),
            mouseid=str(s.get("mouseid")),
            annotation_format=str(s.get("annotation_format", "json")).lower(),
            probe_shank=opt_int(s.get("probe_shank")),
            histology_track_id=histology_track_id,
            logical_probe=logical_probe,
            ephys_collection=ephys_collection,
            histology_shank=histology_shank,
            ephys_shank=ephys_shank,
            surface_finding=opt_path(s.get("surface_finding")),
            registration_asset=opt_path(s.get("registration_asset")),
            row_index=int(s.name) if hasattr(s, "name") else None,  # type: ignore[call-overload]
        )
