"""Asset discovery, zarr level selection, and output directory preparation."""

from __future__ import annotations

import logging
from pathlib import Path
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

    alignment_zarr_path = Path(alignment_zarr_uri)
    registration_channel_stem = alignment_zarr_path.stem
    registration_dir_path = asset_path / "image_atlas_alignment" / f"{registration_channel_stem}"
    registration_in_ccf_precomputed = registration_dir_path / "moved_ls_to_ccf.nii.gz"
    return AssetInfo(
        asset_path=asset_path,
        zarr_volumes=zarr_paths,
        pipeline_registration_chains=pipeline_reg_info,
        registration_dir_path=registration_dir_path,
        registration_in_ccf_precomputed=registration_in_ccf_precomputed,
    )


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


def determine_desired_level(zarr_metadata: dict, desired_voxel_size_um: float = 25.0) -> int:  # type: ignore[type-arg]
    """Select the highest-resolution multiscale level not exceeding *desired_voxel_size_um*.

    Parameters
    ----------
    zarr_metadata : dict
        OME-Zarr metadata containing ``coordinateTransformations``.
    desired_voxel_size_um : float
        Target voxel size in micrometers.

    Returns
    -------
    int
        Zero-based multiscale level index.
    """
    scales = np.array([np.array(x[0]["scale"][2:]).min() for x in zarr_metadata["coordinateTransformations"]])
    level: int = int(
        np.maximum(
            np.searchsorted(scales, desired_voxel_size_um, side="right") - 1,
            0,
        )
    )
    return level
