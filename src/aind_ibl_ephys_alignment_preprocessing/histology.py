"""Volume processing functions for histology data (synchronous).

Handles CCF reorientation, registration channel export, additional channel
transforms, and CCF-to-image-space inverse transforms.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import ants
import SimpleITK as sitk
from aind_registration_utils.annotations import expand_compacted_image
from ants.core import ANTsImage
from ants.utils import to_sitk

from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION
from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    OutputDirs,
    ReferencePaths,
    ReferenceVolumes,
)

try:
    import numpy as np
except ImportError:  # pragma: no cover
    pass

logger = logging.getLogger(__name__)


def convert_img_direction_and_write(
    img: sitk.Image,
    output_path: Path,
    direction: str = _BLESSED_DIRECTION,
) -> None:
    """Convert an image to the specified DICOM orientation and write as compressed NRRD.

    Parameters
    ----------
    img : sitk.Image
        Input SimpleITK image.
    output_path : Path
        Destination NRRD file path.
    direction : str
        Target DICOM orientation code (default ``"IRP"``).
    """
    img_oriented = sitk.DICOMOrient(img, direction)
    sitk.WriteImage(img_oriented, str(output_path), useCompression=True)


def copy_registration_channel_ccf_reorient(
    asset_info: AssetInfo,
    outputs: OutputDirs,
) -> None:
    """Copy precomputed CCF registration to results as reoriented NRRD.

    Parameters
    ----------
    asset_info : AssetInfo
        Asset metadata with precomputed registration path.
    outputs : OutputDirs
        Output directory tree.
    """
    if not asset_info.registration_in_ccf_precomputed.exists():
        raise FileNotFoundError(
            f"Precomputed registration in CCF not found: {asset_info.registration_in_ccf_precomputed}"
        )
    ccf_img = sitk.ReadImage(str(asset_info.registration_in_ccf_precomputed))
    img_in_ccf_dst = outputs.histology_ccf / "histology_registration.nrrd"
    convert_img_direction_and_write(ccf_img, img_in_ccf_dst)


def write_registration_channel_images(
    asset_info: AssetInfo,
    outputs: OutputDirs,
    *,
    level: int = 3,
    opened_zarr: tuple[Any, dict[str, Any]] | None = None,
) -> tuple[Path, Path]:
    """Write registration-channel outputs to CCF and image space.

    Parameters
    ----------
    asset_info : AssetInfo
        Asset metadata.
    outputs : OutputDirs
        Output directory tree.
    level : int
        Multiscale zarr level.
    opened_zarr : tuple | None
        Pre-opened zarr node and metadata.

    Returns
    -------
    tuple[Path, Path]
        ``(raw_img_path, pipeline_img_path)`` of written NRRD files.
    """
    from aind_zarr_utils.pipeline_transformed import base_and_pipeline_zarr_to_sitk
    from aind_zarr_utils.zarr import _open_zarr

    reg_zarr = asset_info.zarr_volumes.registration
    if opened_zarr is None:
        zarr_node, zarr_metadata = _open_zarr(reg_zarr)
    else:
        zarr_node, zarr_metadata = opened_zarr

    metadata = asset_info.zarr_volumes.metadata
    processing = asset_info.zarr_volumes.processing
    raw_img, pipeline_raw_img = base_and_pipeline_zarr_to_sitk(
        reg_zarr,
        metadata,
        processing,
        level=level,
        opened_zarr=(zarr_node, zarr_metadata),
    )
    raw_img_dst = outputs.histology_img / "histology_registration.nrrd"
    convert_img_direction_and_write(raw_img, raw_img_dst)
    del raw_img
    bugged_img_dst = outputs.histology_img / "histology_registration_pipeline.nrrd"
    convert_img_direction_and_write(pipeline_raw_img, bugged_img_dst)
    return raw_img_dst, bugged_img_dst


def process_additional_channels_pipeline(
    pipeline_histology_space_img: ANTsImage,
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    outputs: OutputDirs,
    scratch_root: Path,
    level: int = 3,
) -> None:
    """Process non-alignment OME-Zarr channels into CCF and image space.

    Parameters
    ----------
    pipeline_histology_space_img : ANTsImage
        Pipeline-space histology image for spatial metadata copying.
    asset_info : AssetInfo
        Asset metadata.
    refs : ReferenceVolumes
        Reference volumes (CCF template).
    outputs : OutputDirs
        Output directory tree.
    scratch_root : Path
        Temporary directory for intermediate files.
    level : int
        Multiscale zarr level.
    """
    from aind_zarr_utils.zarr import zarr_to_sitk

    for zarr_path in asset_info.zarr_volumes.additional:
        ch_str = Path(zarr_path).stem
        img_raw = zarr_to_sitk(zarr_path, asset_info.zarr_volumes.metadata, level=level)
        channel_dst = outputs.histology_img / f"{ch_str}.nrrd"
        convert_img_direction_and_write(img_raw, channel_dst)
        del img_raw

        ants_hist_img = ants.image_read(str(channel_dst), pixeltype=None)  # type: ignore[no-untyped-call]
        ants.copy_image_info(pipeline_histology_space_img, ants_hist_img)

        ch_in_ccf = ants.apply_transforms(
            refs.ccf_25,
            ants_hist_img,
            asset_info.pipeline_registration_chains.img_tx_str,
            whichtoinvert=asset_info.pipeline_registration_chains.img_tx_inverted,
        )
        ch_in_ccf_dst = outputs.histology_ccf / f"histology_{ch_str}.nrrd"
        ch_in_ccf_tmp_dst = scratch_root / f"histology-{ch_str}-ccf.nrrd"
        ants.image_write(ch_in_ccf, str(ch_in_ccf_tmp_dst))
        del ch_in_ccf
        try:
            compress_reorient_nrrd_file(
                ch_in_ccf_tmp_dst,
                ch_in_ccf_dst,
                force_orientation=_BLESSED_DIRECTION,
            )
        finally:
            ch_in_ccf_tmp_dst.unlink(missing_ok=True)


def apply_ccf_inverse_tx_then_fix_domain(
    ccf_space_img_moving: ANTsImage,
    pipeline_space_fixed_img: ANTsImage,
    correct_hist_domain_img: ANTsImage,
    asset_info: AssetInfo,
    **kwargs: Any,
) -> ANTsImage:
    """Apply inverse pipeline (CCF -> histology) transform then repair image domain.

    Parameters
    ----------
    ccf_space_img_moving : ANTsImage
        Image in CCF/template space to move into histology space.
    pipeline_space_fixed_img : ANTsImage
        Image in the pipeline's (buggy) histology space.
    correct_hist_domain_img : ANTsImage
        Reference histology image with correct spacing/origin/direction.
    asset_info : AssetInfo
        Pipeline registration chain paths.
    **kwargs
        Forwarded to ``ants.apply_transforms``.

    Returns
    -------
    ANTsImage
        Transformed image with corrected spatial domain.
    """
    pt_tx_str = asset_info.pipeline_registration_chains.pt_tx_str
    pt_tx_inverted = asset_info.pipeline_registration_chains.pt_tx_inverted
    ccf_space_img_in_hist_space: ANTsImage = ants.apply_transforms(
        fixed=pipeline_space_fixed_img,
        moving=ccf_space_img_moving,
        transformlist=pt_tx_str,
        whichtoinvert=pt_tx_inverted,
        **kwargs,
    )
    ccf_space_img_in_hist_space.set_spacing(correct_hist_domain_img.spacing)
    ccf_space_img_in_hist_space.set_origin(correct_hist_domain_img.origin)
    ccf_space_img_in_hist_space.set_direction(correct_hist_domain_img.direction)
    return ccf_space_img_in_hist_space


def compress_reorient_nrrd_file(
    input_path: Path,
    output_path: Path,
    force_orientation: str | None = None,
) -> None:
    """Re-compress and optionally reorient an NRRD file.

    Parameters
    ----------
    input_path : Path
        Source NRRD file.
    output_path : Path
        Destination NRRD file.
    force_orientation : str | None
        If set, reorient to this DICOM code before writing.
    """
    img = sitk.ReadImage(str(input_path))
    orientation_code = sitk.DICOMOrientImageFilter.GetOrientationFromDirectionCosines(img.GetDirection())
    if force_orientation is not None and orientation_code != force_orientation:
        logger.info("Reorienting %s from %s to %s", input_path, orientation_code, force_orientation)
        out_img = sitk.DICOMOrient(img, force_orientation)
    else:
        out_img = img
    temp_output_path = output_path.with_suffix(".temp.nrrd")
    sitk.WriteImage(out_img, str(temp_output_path), useCompression=True)
    temp_output_path.replace(output_path)


def transform_ccf_to_image_space(
    asset_info: AssetInfo,
    refs: ReferenceVolumes,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
) -> None:
    """Transform CCF template into native image space.

    Parameters
    ----------
    asset_info : AssetInfo
        Pipeline registration chain paths.
    refs : ReferenceVolumes
        Reference volumes.
    raw_hist_img : ANTsImage
        Histology image with correct spatial domain.
    pipeline_hist_domain_img : ANTsImage
        Histology image in pipeline (buggy) domain.
    outputs : OutputDirs
        Output directory tree.
    """
    ccf_in_hist_img = apply_ccf_inverse_tx_then_fix_domain(
        refs.ccf_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
    )
    ccf_in_hist_img_path = outputs.histology_img / "ccf_in_mouse.nrrd"
    ccf_in_hist_sitk = to_sitk(ccf_in_hist_img)
    del ccf_in_hist_img
    convert_img_direction_and_write(ccf_in_hist_sitk, ccf_in_hist_img_path)


def transform_ccf_labels_to_image_space(
    asset_info: AssetInfo,
    ref_paths: ReferencePaths,
    raw_hist_img: ANTsImage,
    pipeline_hist_domain_img: ANTsImage,
    outputs: OutputDirs,
) -> None:
    """Transform lateralized CCF labels into native image space.

    Parameters
    ----------
    asset_info : AssetInfo
        Pipeline registration chain paths.
    ref_paths : ReferencePaths
        Reference data paths.
    raw_hist_img : ANTsImage
        Histology image with correct spatial domain.
    pipeline_hist_domain_img : ANTsImage
        Histology image in pipeline (buggy) domain.
    outputs : OutputDirs
        Output directory tree.
    """
    ccf_labels_lateralized_25 = ants.image_read(
        str(ref_paths.ccf_labels_lateralized_25),
        pixeltype=None,  # type: ignore[no-untyped-call]
    )
    unq_vals = np.load(str(ref_paths.ccf_labels_lateralized_25_unq_vals))["unique_labels"]
    ccf_labels_in_hist_img = apply_ccf_inverse_tx_then_fix_domain(
        ccf_labels_lateralized_25,
        pipeline_space_fixed_img=pipeline_hist_domain_img,
        correct_hist_domain_img=raw_hist_img,
        asset_info=asset_info,
        interpolator="genericLabel",
    )
    del ccf_labels_lateralized_25
    ccf_labels_sitk = to_sitk(ccf_labels_in_hist_img)
    del ccf_labels_in_hist_img
    ccf_labels_expanded = expand_compacted_image(ccf_labels_sitk, unq_vals)
    ccf_labels_in_hist_img_path = outputs.histology_img / "labels_in_mouse.nrrd"
    convert_img_direction_and_write(ccf_labels_expanded, ccf_labels_in_hist_img_path)
