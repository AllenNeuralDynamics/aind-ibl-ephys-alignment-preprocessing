"""Asset discovery, zarr level selection, and output directory preparation."""

from __future__ import annotations

import logging
from pathlib import Path, PurePath
from typing import TYPE_CHECKING

import numpy as np
from aind_s3_cache.json_utils import get_json
from aind_s3_cache.uri_utils import as_pathlike
from aind_zarr_utils.neuroglancer import get_image_sources
from aind_zarr_utils.pipeline_transformed import (
    _asset_from_zarr_pathlike,
    alignment_zarr_uri_and_metadata_from_zarr_or_asset_pathlike,
    pipeline_transforms_local_paths,
)

from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    OutputDirs,
    PipelineRegistrationInfo,
    ZarrPaths,
)

if TYPE_CHECKING:
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

logger = logging.getLogger(__name__)

#: Transforms a registration directory must supply, whichever asset it lives in.
REGISTRATION_TRANSFORMS = (
    "ls_to_template_SyN_0GenericAffine.mat",
    "ls_to_template_SyN_1InverseWarp.nii.gz",
)


def manifest_registration_override(config: PipelineConfig) -> Path | None:
    """Return the manifest's ``registration_asset`` path fragment, if any.

    The column is per-*brain*: like ``mouseid`` it is replicated across rows, so
    the first non-empty value wins. Rows disagreeing is a manifest error, not
    something to resolve silently, so it raises.

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration with resolved paths.

    Returns
    -------
    Path or None
        The path fragment, relative to ``config.data_root``, or ``None`` when
        the column is absent/empty (the overwhelming majority of manifests).

    Raises
    ------
    ValueError
        If rows name different registration assets.
    """
    import pandas as pd

    if not config.manifest_csv.exists():
        return None
    df = pd.read_csv(config.manifest_csv)
    if "registration_asset" not in df.columns:
        return None
    values = {str(v).strip() for v in df["registration_asset"] if pd.notna(v) and str(v).strip()}
    if not values:
        return None
    if len(values) > 1:
        raise ValueError(
            f"registration_asset must name one directory per manifest; {config.manifest_csv} has {sorted(values)}"
        )
    return Path(values.pop())


def resolve_registration_dir(config: PipelineConfig, asset_path: Path) -> tuple[Path, str | None]:
    """Resolve the registration directory and the channel it was computed from.

    Without an override this is the stitched asset's own
    ``image_atlas_alignment/<channel>/``, and the channel is left to the caller
    (``processing.json`` names it). With one, the manifest points straight at a
    directory and the channel comes from its *name* -- ``processing.json`` names
    the channel that was registered in error, which is the whole reason an
    override exists, so it must not be consulted.

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration with resolved paths.
    asset_path : Path
        The stitched SmartSPIM asset root.

    Returns
    -------
    tuple[Path, str or None]
        ``(registration_dir, channel)``. *channel* is ``None`` when there is no
        override, meaning "keep whatever the asset's own metadata says".

    Raises
    ------
    FileNotFoundError
        If the override names a directory that does not exist, or one missing
        either transform.
    """
    override = manifest_registration_override(config)
    if override is None:
        return asset_path / "image_atlas_alignment", None

    reg_dir = config.data_root / override
    if not reg_dir.is_dir():
        raise FileNotFoundError(
            f"registration_asset {override} not found under {config.data_root} (resolved to {reg_dir})"
        )
    missing = [name for name in REGISTRATION_TRANSFORMS if not (reg_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"registration_asset {override} is missing {', '.join(missing)} in {reg_dir}")
    # ``ccf_<channel>`` is the standalone capsule's layout; the pipeline's own is
    # ``image_atlas_alignment/<channel>``, where the name is already the channel.
    channel = reg_dir.name.removeprefix("ccf_")
    logger.info("using registration override %s (channel %s)", reg_dir, channel)
    return reg_dir, channel


def find_asset_info(config: PipelineConfig) -> AssetInfo:
    """Discover SmartSPIM asset structure from the Neuroglancer JSON.

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration with resolved paths.

    Returns
    -------
    AssetInfo
        Discovered asset metadata (paths, zarr channels, transforms).

    Raises
    ------
    ValueError
        If no image sources are found in the Neuroglancer file.
    FileNotFoundError
        If the inferred asset path does not exist.
    """
    ng_data = get_json(str(config.neuroglancer_file))
    sources = get_image_sources(ng_data, remove_zarr_protocol=True)
    a_zarr_uri = next(iter(sources.values()), None)
    if a_zarr_uri is None:
        raise ValueError("No image sources found in neuroglancer data")
    _, _, a_zarr_pathlike = as_pathlike(a_zarr_uri)
    asset_pathlike = _asset_from_zarr_pathlike(a_zarr_pathlike)
    asset_path = config.data_root / asset_pathlike
    if not asset_path.exists():
        raise FileNotFoundError(f"Asset path not found: {asset_path}")
    asset_path_str = str(asset_path)
    zarr_path = asset_path / "image_tile_fusing" / "OMEZarr"
    image_channel_zarrs = [p for p in zarr_path.iterdir() if p.is_dir() and p.suffix == ".zarr"]
    alignment_zarr_uri, metadata, processing_data = alignment_zarr_uri_and_metadata_from_zarr_or_asset_pathlike(
        asset_uri=asset_path_str,
    )

    # An override moves the transforms out of the stitched asset, and with them
    # the channel they were computed from. The zarr itself never moves: only
    # the registration directory lives elsewhere.
    registration_dir_path, override_channel = resolve_registration_dir(config, asset_path)
    if override_channel is not None:
        override_zarr = zarr_path / f"{override_channel}.zarr"
        if not override_zarr.is_dir():
            raise FileNotFoundError(
                f"registration_asset names channel {override_channel!r}, but {override_zarr} does not exist; "
                f"available: {sorted(p.stem for p in image_channel_zarrs)}"
            )
        alignment_zarr_uri = override_zarr.as_posix()

    other_channels = list({p.as_posix() for p in image_channel_zarrs} - {alignment_zarr_uri})
    zarr_paths = ZarrPaths(
        registration=alignment_zarr_uri,
        additional=other_channels,
        metadata=metadata,
        processing=processing_data,
    )

    pt_tx_str, pt_tx_inverted, img_tx_str, img_tx_inverted = pipeline_transforms_local_paths(
        alignment_zarr_uri,
        processing_data,
        anonymous=True,
    )
    pipeline_reg_info = PipelineRegistrationInfo(
        pt_tx_str=pt_tx_str,
        pt_tx_inverted=pt_tx_inverted,
        img_tx_str=img_tx_str,
        img_tx_inverted=img_tx_inverted,
    )

    if override_channel is None:
        # No override: the registration lives beside the channel the asset's own
        # metadata nominated, so the stem is the channel.
        registration_dir_path = registration_dir_path / Path(alignment_zarr_uri).stem
    return AssetInfo(
        asset_path=asset_path,
        asset_uri=_infer_asset_uri(a_zarr_uri, asset_pathlike),
        zarr_volumes=zarr_paths,
        pipeline_registration_chains=pipeline_reg_info,
        registration_dir_path=registration_dir_path,
    )


def _infer_asset_uri(source_uri: str, asset_pathlike: str | PurePath) -> str | None:
    """Infer the asset-level URI from a source URI inside that asset.

    Neuroglancer sources often point inside the stitched SmartSPIM asset, for
    example ``s3://bucket/<asset>/image_tile_fusing/...``. The datapackage
    registry should name the asset that contains ``image_atlas_alignment/``,
    not the nested zarr. If the source URI does not contain the inferred asset
    pathlike verbatim, keep the durable URI unknown.
    """
    if "://" not in source_uri:
        return None
    marker = str(asset_pathlike).strip("/")
    if not marker:
        return None
    idx = source_uri.find(marker)
    if idx < 0:
        return None
    return source_uri[: idx + len(marker)]


def prepare_result_dirs(mouse_id: str, results_root: Path) -> OutputDirs:
    """Create the output directory tree for a single mouse.

    Parameters
    ----------
    mouse_id : str
        Mouse identifier.
    results_root : Path
        Root results directory.

    Returns
    -------
    OutputDirs
        Freshly-created directory tree.
    """
    histology_ccf = results_root / mouse_id / "ccf_space_histology"
    histology_img = results_root / mouse_id / "image_space_histology"
    tracks_root = results_root / mouse_id / "track_data"
    spim = tracks_root / "spim"
    template = tracks_root / "template"
    ccf = tracks_root / "ccf"
    bregma = tracks_root / "bregma_xyz"
    for d in (histology_ccf, histology_img, spim, template, ccf, bregma):
        d.mkdir(parents=True, exist_ok=True)
    return OutputDirs(histology_ccf, histology_img, tracks_root, spim, template, ccf, bregma)


def determine_desired_level(zarr_metadata: dict, desired_voxel_size_um: float) -> int:  # type: ignore[type-arg]
    """Select the highest-resolution multiscale level not exceeding *desired_voxel_size_um*.

    Parameters
    ----------
    zarr_metadata : dict
        OME-Zarr metadata containing ``coordinateTransformations``.
    desired_voxel_size_um : float
        Target voxel size in micrometers. Required rather than defaulted: the one
        value that governs a run lives on
        :attr:`~aind_ibl_ephys_alignment_preprocessing.types.PipelineConfig.desired_voxel_size_um`,
        and a second default here would be dead code that still reads like the
        answer.

    Returns
    -------
    int
        Zero-based multiscale level index.

    Notes
    -----
    The comparison is against ``min(z, y, x)`` of each level, so on an
    anisotropic pyramid the target acts as a floor on the *finest* axis rather
    than a target for the grid: a volume whose z is finer than its xy is selected
    on z and processed with xy coarser than asked for.
    """
    scales = np.array([np.array(x[0]["scale"][2:]).min() for x in zarr_metadata["coordinateTransformations"]])
    level: int = int(
        np.maximum(
            np.searchsorted(scales, desired_voxel_size_um, side="right") - 1,
            0,
        )
    )
    return level
