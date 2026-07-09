"""Datapackage for preprocessing pipeline outputs."""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from pydantic import BaseModel

from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    ManifestRow,
    OutputDirs,
    PipelineConfig,
    ProcessResult,
)

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "3.0.0"
# Why 3.0.0: 2.x stored filesystem locations as strings. 3.0.0 changes every
# datapackage path into a reference object ``{asset, path}``, where
# ``asset=None`` means datapackage-local and external assets are resolved via
# the ``external_assets`` registry plus consumer/producers runtime policy.


# ---------------------------------------------------------------------------
# Pydantic models mirroring the datapackage.json schema
# ---------------------------------------------------------------------------


class PathReference(BaseModel, frozen=True):
    """A path inside either the datapackage asset or an external asset."""

    asset: str | None = None
    path: str


class ExternalAsset(BaseModel, frozen=True):
    """External asset descriptor used to resolve path references."""

    role: str
    name: str
    id: str | None = None
    uri: str | None = None
    checksum: str | None = None
    provenance: dict[str, Any] | None = None


class TransformPaths(BaseModel, frozen=True):
    """References to the ANTs transform chain files."""

    image_to_template_affine: PathReference
    image_to_template_warp: PathReference
    template_to_ccf_affine: PathReference
    template_to_ccf_warp: PathReference


class ImageSpaceHistology(BaseModel, frozen=True):
    """References to image-space histology volumes."""

    registration: PathReference
    registration_pipeline: PathReference
    ccf_template: PathReference
    labels: PathReference
    additional_channels: list[PathReference] = []


class CcfSpaceHistology(BaseModel, frozen=True):
    """References to CCF-space histology volumes."""

    registration: PathReference
    additional_channels: list[PathReference] = []


class HistologyPaths(BaseModel, frozen=True):
    """Histology output paths grouped by coordinate space."""

    image_space: ImageSpaceHistology
    ccf_space: CcfSpaceHistology


class XyzPicks(BaseModel, frozen=True):
    """References to xyz-picks JSON files for one shank (or whole probe)."""

    ccf: PathReference
    image_space: PathReference
    histology_track_id: str | None = None
    histology_shank: int | None = None
    ephys_shank: int | None = None
    shank: int | None = None


class ChannelTablePaths(BaseModel, frozen=True):
    """References to the producer-owned ephys channel geometry table."""

    local_coordinates: PathReference
    raw_ind: PathReference
    contact_id: PathReference | None = None
    shank_ind: PathReference


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
    ephys: PathReference | None = None
    channel_table: ChannelTablePaths | None = None
    xyz_picks: list[XyzPicks]


class DataPackageError(Exception):
    """Raised when a datapackage references paths that are missing on disk."""


# File the GUI loads first from each ephys collection dir (load_data_local.py).
# Its presence is the load-bearing signal that ``ephys`` points at a real ALF
# collection rather than an empty/wrong directory.
_REQUIRED_EPHYS_FILE = "channels.localCoordinates.npy"


class ReferenceResolver:
    """Resolve datapackage path references using local runtime policy."""

    def __init__(
        self,
        *,
        datapackage_dir: Path,
        external_assets: Mapping[str, ExternalAsset],
        asset_roots: Iterable[Path] | None = None,
        asset_overrides: Mapping[str, Path] | None = None,
    ) -> None:
        self.datapackage_dir = Path(datapackage_dir)
        self.external_assets = external_assets
        self.asset_roots = [Path(p) for p in asset_roots or []]
        self.asset_overrides = {str(k): Path(v) for k, v in (asset_overrides or {}).items()}

    def resolve(self, ref: PathReference) -> Path | None:
        """Return a local path for *ref*, or ``None`` if its asset is unknown."""
        if ref.asset is None:
            return self.datapackage_dir / ref.path

        asset = self.external_assets.get(ref.asset)
        if asset is None:
            return None

        for key in (ref.asset, asset.name):
            if key and key in self.asset_overrides:
                return self.asset_overrides[key] / ref.path

        for root in self.asset_roots:
            for key in (asset.name, asset.id):
                if key:
                    candidate = root / key / ref.path
                    if candidate.exists():
                        return candidate
        return None


class DataPackage(BaseModel, frozen=True):
    """Top-level datapackage for preprocessing outputs."""

    schema_version: str
    mouse_id: str
    platform: str
    external_assets: dict[str, ExternalAsset]
    transforms: TransformPaths
    histology: HistologyPaths
    # Nested ``recording_id -> ephys_collection -> ProbeEntry``. Per-shank
    # rows collapse into ``ProbeEntry.xyz_picks``.
    probes: dict[str, dict[str, ProbeEntry]]

    def referenced_files(self) -> list[PathReference]:
        """Every *file* reference this datapackage points at.

        Excludes ``ephys`` entries, which are directories (see
        :meth:`referenced_ephys_dirs`).
        """
        paths: list[PathReference] = [
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
                    if probe.channel_table.contact_id is not None:
                        paths.append(probe.channel_table.contact_id)
                    paths.append(probe.channel_table.shank_ind)
                for xp in probe.xyz_picks:
                    paths.append(xp.ccf)
                    paths.append(xp.image_space)
        return paths

    def referenced_ephys_dirs(self) -> list[PathReference]:
        """Ephys collection directory references (one per probe with ephys)."""
        return [
            probe.ephys for recording in self.probes.values() for probe in recording.values() if probe.ephys is not None
        ]

    def missing_paths(
        self,
        root: Path,
        *,
        asset_roots: Iterable[Path] | None = None,
        asset_overrides: Mapping[str, Path] | None = None,
    ) -> list[str]:
        """Return referenced paths that do not exist.

        *root* is the directory that will hold ``datapackage.json``; all stored
        local paths are relative to it. External paths are resolved via
        ``external_assets`` plus *asset_roots*/*asset_overrides*. Each ``ephys``
        entry must be a directory containing ``channels.localCoordinates.npy``
        (what the GUI loads first), so a wrong ephys directory is reported
        rather than passing silently.
        """
        root = Path(root)
        missing: list[str] = []
        resolver = ReferenceResolver(
            datapackage_dir=root,
            external_assets=self.external_assets,
            asset_roots=list(asset_roots or []),
            asset_overrides=dict(asset_overrides or {}),
        )
        for ref in self.referenced_files():
            resolved = resolver.resolve(ref)
            if resolved is None or not resolved.exists():
                missing.append(_display_ref(ref))
        for ref in self.referenced_ephys_dirs():
            ephys = resolver.resolve(ref)
            if ephys is None or not ephys.is_dir():
                missing.append(f"{_display_ref(ref)}/ (ephys directory)")
            elif not (ephys / _REQUIRED_EPHYS_FILE).is_file():
                missing.append(f"{_display_ref(ref)}/{_REQUIRED_EPHYS_FILE}")
        return sorted(missing)

    def validate_paths(
        self,
        root: Path,
        *,
        asset_roots: Iterable[Path] | None = None,
        asset_overrides: Mapping[str, Path] | None = None,
    ) -> None:
        """Raise :class:`DataPackageError` if any referenced path is missing.

        Parameters
        ----------
        root : Path
            Directory where ``datapackage.json`` will be written.
        """
        missing = self.missing_paths(root, asset_roots=asset_roots, asset_overrides=asset_overrides)
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

    external_assets = _build_external_assets(asset_info, config)
    transforms = _build_transforms(asset_info, config, manifest_root)
    histology = _build_histology(outputs, manifest_root)
    probes = _build_probes(manifest_rows, results, outputs, manifest_root, config)

    return DataPackage(
        schema_version=SCHEMA_VERSION,
        mouse_id=mouse_id,
        platform=_detect_platform(config),
        external_assets=external_assets,
        transforms=transforms,
        histology=histology,
        probes=probes,
    )


def _build_transforms(asset_info: AssetInfo, config: PipelineConfig, manifest_root: Path) -> TransformPaths:
    reg_dir = asset_info.registration_dir_path
    tmpl_dir = config.template_to_ccf_dir
    return TransformPaths(
        image_to_template_affine=_external_ref(
            "smartspim",
            reg_dir / "ls_to_template_SyN_0GenericAffine.mat",
            asset_info.asset_path,
        ),
        image_to_template_warp=_external_ref(
            "smartspim",
            reg_dir / "ls_to_template_SyN_1InverseWarp.nii.gz",
            asset_info.asset_path,
        ),
        template_to_ccf_affine=_external_ref(
            "spim_template_to_ccf",
            tmpl_dir / "syn_0GenericAffine.mat",
            tmpl_dir,
        ),
        template_to_ccf_warp=_external_ref(
            "spim_template_to_ccf",
            tmpl_dir / "syn_1InverseWarp.nii.gz",
            tmpl_dir,
        ),
    )


def _build_external_assets(asset_info: AssetInfo, config: PipelineConfig) -> dict[str, ExternalAsset]:
    platform = _detect_platform(config)
    return {
        "smartspim": ExternalAsset(
            role="smartspim_registration",
            name=asset_info.asset_path.name,
            uri=asset_info.asset_uri,
            provenance={"mounted_at": str(asset_info.asset_path), "on": platform},
        ),
        "spim_template_to_ccf": ExternalAsset(
            role="template_to_ccf",
            name=config.template_to_ccf_dir.name,
            provenance={"mounted_at": str(config.template_to_ccf_dir), "on": platform},
        ),
    }


def _detect_platform(config: PipelineConfig) -> str:
    """Best-effort platform descriptor for error messages and provenance."""
    if config.data_root == Path("/data") or config.data_root.as_posix().startswith("/data/"):
        return "code_ocean"
    return "local"


def _local_ref(path: Path, root: Path) -> PathReference:
    """Return a datapackage-local path reference."""
    return PathReference(asset=None, path=path.relative_to(root).as_posix())


def _external_ref(asset: str, path: Path, asset_root: Path) -> PathReference:
    """Return a path reference inside an external asset root."""
    return PathReference(asset=asset, path=path.relative_to(asset_root).as_posix())


def _display_ref(ref: PathReference) -> str:
    """Human-readable reference for validation errors."""
    return ref.path if ref.asset is None else f"{ref.asset}:{ref.path}"


def producer_asset_overrides(asset_info: AssetInfo, config: PipelineConfig) -> dict[str, Path]:
    """Local asset locations available during producer-side validation."""
    return {
        "smartspim": asset_info.asset_path,
        asset_info.asset_path.name: asset_info.asset_path,
        "spim_template_to_ccf": config.template_to_ccf_dir,
        config.template_to_ccf_dir.name: config.template_to_ccf_dir,
    }


def _build_histology(outputs: OutputDirs, manifest_root: Path) -> HistologyPaths:
    img_dir = outputs.histology_img
    ccf_dir = outputs.histology_ccf

    # Image-space additional channels (Ex_\d+_Em_\d+.nrrd)
    img_additional = sorted(
        (_local_ref(p, manifest_root) for p in img_dir.glob("Ex_*_Em_*.nrrd")),
        key=lambda ref: ref.path,
    )

    # CCF-space additional channels (histology_*.nrrd excluding registration)
    ccf_additional = sorted(
        (
            _local_ref(p, manifest_root)
            for p in ccf_dir.glob("histology_*.nrrd")
            if p.name != "histology_registration.nrrd"
        ),
        key=lambda ref: ref.path,
    )

    return HistologyPaths(
        image_space=ImageSpaceHistology(
            registration=_local_ref(img_dir / "histology_registration.nrrd", manifest_root),
            registration_pipeline=_local_ref(img_dir / "histology_registration_pipeline.nrrd", manifest_root),
            ccf_template=_local_ref(img_dir / "ccf_in_mouse.nrrd", manifest_root),
            labels=_local_ref(img_dir / "labels_in_mouse.nrrd", manifest_root),
            additional_channels=img_additional,
        ),
        ccf_space=CcfSpaceHistology(
            registration=_local_ref(ccf_dir / "histology_registration.nrrd", manifest_root),
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
            gui_ccf, gui_img, shank = _xyz_pick_names(histology_shank=histology_shank, ephys_shank=ephys_shank)

            xyz_picks_list.append(
                XyzPicks(
                    ccf=_local_ref(gui_folder / gui_ccf, manifest_root),
                    image_space=_local_ref(gui_folder / gui_img, manifest_root),
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
        ephys_path = _local_ref(ephys_dir, manifest_root) if not config.skip_ephys else None
        contact_id_path = ephys_dir / "channels.contactId.npy"
        channel_table = (
            ChannelTablePaths(
                local_coordinates=_local_ref(
                    ephys_dir / "channels.localCoordinates.npy",
                    manifest_root,
                ),
                raw_ind=_local_ref(ephys_dir / "channels.rawInd.npy", manifest_root),
                contact_id=_local_ref(contact_id_path, manifest_root) if contact_id_path.is_file() else None,
                shank_ind=_local_ref(ephys_dir / "channels.shankInd.npy", manifest_root),
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


def _xyz_pick_names(
    *,
    histology_shank: int | None,
    ephys_shank: int | None,
) -> tuple[str, str, int | None]:
    """Return CCF/image-space xyz-pick filenames and legacy 1-based shank."""
    file_shank = histology_shank if histology_shank is not None else ephys_shank
    if file_shank is None:
        return "xyz_picks.json", "xyz_picks_image_space.json", None
    shank_id = int(file_shank) + 1
    shank = int(ephys_shank) + 1 if ephys_shank is not None else shank_id
    return f"xyz_picks_shank{shank_id}.json", f"xyz_picks_shank{shank_id}_image_space.json", shank


def infer_process_results_from_outputs(manifest_rows: list[ManifestRow], outputs: OutputDirs) -> list[ProcessResult]:
    """Infer per-row success from existing xyz-pick outputs.

    This is intentionally lightweight: it does not rerun histology or ephys
    work. It is used to regenerate ``datapackage.json`` after schema/contract
    changes when the expensive preprocessing outputs already exist.
    """
    results: list[ProcessResult] = []
    for row in manifest_rows:
        gui_folder = row.gui_folder(outputs)
        ccf_name, img_name, _ = _xyz_pick_names(
            histology_shank=_row_histology_shank(row),
            ephys_shank=_row_ephys_shank(row),
        )
        wrote_files = (gui_folder / ccf_name).is_file() and (gui_folder / img_name).is_file()
        results.append(
            ProcessResult(
                probe_id=_row_histology_track_id(row),
                recording_id=row.recording_id,
                wrote_files=wrote_files,
                skipped_reason=None if wrote_files else "existing xyz-picks outputs not found",
            )
        )
    return results


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def write_datapackage(
    dp: DataPackage,
    output_dir: Path,
    *,
    validate: bool = True,
    asset_roots: Iterable[Path] | None = None,
    asset_overrides: Mapping[str, Path] | None = None,
) -> Path:
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
        *output_dir* or the configured external asset locations before
        writing, raising :class:`DataPackageError` on any miss. This turns a
        well-typed but dangling manifest (e.g. an ephys dir pointing at a
        nonexistent subdirectory) into a loud failure at write time instead of
        a runtime error in the GUI. Set False to serialize without touching the
        filesystem (round-trip tests, fixtures).

    Returns
    -------
    Path
        The path of the written file.
    """
    if validate:
        dp.validate_paths(output_dir, asset_roots=asset_roots, asset_overrides=asset_overrides)
    path = output_dir / "datapackage.json"
    path.write_text(dp.model_dump_json(indent=2))
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
