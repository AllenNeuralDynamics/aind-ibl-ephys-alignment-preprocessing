"""Per-probe manifest row processing (synchronous)."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import SimpleITK as sitk
from aind_anatomical_utils.coordinate_systems import convert_coordinate_system
from aind_ephys_ibl_gui_conversion.histology import create_slicer_fcsv
from aind_registration_utils.ants import apply_ants_transforms_to_point_arr
from aind_s3_cache.json_utils import get_json
from aind_zarr_utils.neuroglancer import neuroglancer_annotations_to_anatomical
from iblatlas.atlas import AllenAtlas

from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    ManifestRow,
    OutputDirs,
    ProcessResult,
)

logger = logging.getLogger(__name__)


def process_manifest_row(
    row: ManifestRow,
    asset_info: AssetInfo,
    hist_stub: sitk.Image,
    hist_stub_buggy: sitk.Image,
    ibl_atlas: AllenAtlas,
    outputs: OutputDirs,
    data_root: Path,
) -> ProcessResult:
    """End-to-end processing for a single manifest row.

    Steps
    -----
    1. Locate annotation, load NG points.
    2. Convert to image-space physical (LPS) coordinates.
    3. Write FCSV for SPIM/template/CCF.
    4. Convert to IBL xyz-picks and write JSONs (with shank handling).
    5. Copy xyz-picks into GUI sorting folder.

    Parameters
    ----------
    row : ManifestRow
        Manifest row describing the probe.
    asset_info : AssetInfo
        Asset metadata.
    hist_stub : sitk.Image
        Histology stub image with correct domain.
    hist_stub_buggy : sitk.Image
        Histology stub image in pipeline (buggy) domain.
    ibl_atlas : AllenAtlas
        Allen atlas instance for CCF-to-bregma conversion.
    outputs : OutputDirs
        Output directory tree.
    data_root : Path
        Root directory for input data.

    Returns
    -------
    ProcessResult
        Summary of success/skip for the row.
    """
    # 1) Locate annotation (JSON default)
    ext = "json" if row.annotation_format == "json" else None
    if ext is None:
        return ProcessResult(row.probe_id, row.recording_id, False, "Only JSON annotations supported")
    pattern = f"*/{row.probe_file}.{ext}"

    ann_path = next(data_root.glob(pattern), None)
    probe_id = str(row.probe_id)
    if ann_path is None:
        return ProcessResult(probe_id, str(row.sorted_recording), False, f"Annotation not found: {pattern}")

    # 2) Load NG points in histology space
    ng_data = get_json(str(ann_path))
    anno_zarr = asset_info.zarr_volumes.registration
    metadata = asset_info.zarr_volumes.metadata
    probe_pt_dict, _ = neuroglancer_annotations_to_anatomical(
        ng_data,
        anno_zarr,
        metadata,
        layer_names=[probe_id],
        stub_image=hist_stub,
    )
    probe_pts = probe_pt_dict.get(probe_id, None)
    if probe_pts is None:
        return ProcessResult(probe_id, str(row.sorted_recording), False, f"Probe points not found: {probe_id}")

    # 3) Write SPIM FCSV
    create_slicer_fcsv(str(outputs.spim / f"{probe_id}.fcsv"), probe_pts, direction="LPS")

    # 4) Image -> Template (points) via buggy pipeline transform
    probe_pt_dict_buggy, _ = neuroglancer_annotations_to_anatomical(
        ng_data,
        anno_zarr,
        metadata,
        layer_names=[probe_id],
        stub_image=hist_stub_buggy,
    )
    probe_pts_buggy = probe_pt_dict_buggy[probe_id]
    tx_list_pt_template = [
        str(asset_info.registration_dir_path / "ls_to_template_SyN_0GenericAffine.mat"),
        str(asset_info.registration_dir_path / "ls_to_template_SyN_1InverseWarp.nii.gz"),
    ]
    tx_list_pt_template_invert = [True, False]
    pts_template = apply_ants_transforms_to_point_arr(
        probe_pts_buggy,
        tx_list_pt_template,
        whichtoinvert=tx_list_pt_template_invert,
    )
    create_slicer_fcsv(str(outputs.template / f"{probe_id}.fcsv"), pts_template, direction="LPS")

    # 5) Template -> CCF (points)
    pts_ccf = apply_ants_transforms_to_point_arr(
        probe_pts_buggy,
        asset_info.pipeline_registration_chains.pt_tx_str,
        whichtoinvert=asset_info.pipeline_registration_chains.pt_tx_inverted,
    )
    create_slicer_fcsv(str(outputs.ccf / f"{row.probe_id}.fcsv"), pts_ccf, direction="LPS")

    # 6) IBL xyz-picks (um) from CCF (ML/AP/DV with signed flips)
    ccf_mlapdv_um = convert_coordinate_system(1000.0 * pts_ccf, src_coord="LPS", dst_coord="RPI")
    bregma_mlapdv_um = 1_000_000.0 * ibl_atlas.ccf2xyz(ccf_mlapdv_um, ccf_order="mlapdv")

    # Image-space xyz-picks (um)
    xyz_img = 1000.0 * convert_coordinate_system(probe_pts, src_coord="LPS", dst_coord="RAS")

    xyz_picks_image = {"xyz_picks": xyz_img.tolist()}
    xyz_picks_ccf = {"xyz_picks": bregma_mlapdv_um.tolist()}

    gui_folder = row.gui_folder(outputs)
    if row.probe_shank is None:
        img_name = f"{row.probe_id}_image_space.json"
        ccf_name = f"{row.probe_id}_ccf.json"
        gui_img = "xyz_picks_image_space.json"
        gui_ccf = "xyz_picks.json"
    else:
        shank_id = int(row.probe_shank) + 1
        img_name = f"{row.probe_id}_shank{shank_id}_image_space.json"
        ccf_name = f"{row.probe_id}_shank{shank_id}_ccf.json"
        gui_img = f"xyz_picks_shank{shank_id}_image_space.json"
        gui_ccf = f"xyz_picks_shank{shank_id}.json"

    # 7) Write bregma_xyz JSONs (global per-mouse)
    (outputs.bregma_xyz / img_name).write_text(json.dumps(xyz_picks_image))
    (outputs.bregma_xyz / ccf_name).write_text(json.dumps(xyz_picks_ccf))

    # 8) Per-recording GUI artifacts
    gui_folder.mkdir(parents=True, exist_ok=True)
    (gui_folder / gui_img).write_text(json.dumps(xyz_picks_image))
    (gui_folder / gui_ccf).write_text(json.dumps(xyz_picks_ccf))

    return ProcessResult(
        probe_id=str(row.probe_id),
        recording_id=row.recording_id,
        wrote_files=True,
        skipped_reason=None,
    )
