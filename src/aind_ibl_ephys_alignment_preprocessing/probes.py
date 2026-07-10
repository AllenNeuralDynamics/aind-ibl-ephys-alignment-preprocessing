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
    *,
    emit_qc: bool = False,
) -> ProcessResult:
    """End-to-end processing for a single manifest row.

    Always writes the image-space ``xyz_picks`` (the only track output the
    alignment GUI reads). The QC/diagnostic outputs — SPIM/template/CCF Slicer
    FCSVs and the CCF/bregma ``xyz_picks`` — are produced only when *emit_qc* is
    True; they require the ANTs point-warps (via ``hist_stub_buggy``) and the
    ``ibl_atlas``, and nothing in the alignment workflow consumes them.

    Parameters
    ----------
    row : ManifestRow
        Manifest row describing the probe.
    asset_info : AssetInfo
        Asset metadata.
    hist_stub : sitk.Image
        Histology stub image with correct domain.
    hist_stub_buggy : sitk.Image
        Histology stub image in pipeline (buggy) domain. Only used when
        *emit_qc* is True.
    ibl_atlas : AllenAtlas
        Allen atlas instance for CCF-to-bregma conversion. Only used when
        *emit_qc* is True.
    outputs : OutputDirs
        Output directory tree.
    data_root : Path
        Root directory for input data.
    emit_qc : bool
        Produce the GUI-unused QC outputs (FCSVs + CCF/bregma picks).

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
    probe_id = str(row.histology_track_id)
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

    # Image-space xyz-picks (um) — the ONLY track output the alignment GUI reads.
    # No ANTs warp: just a coordinate-system flip of the stub-domain NG points.
    xyz_img = 1000.0 * convert_coordinate_system(probe_pts, src_coord="LPS", dst_coord="RAS")
    xyz_picks_image = {"xyz_picks": xyz_img.tolist()}

    gui_folder = row.gui_folder(outputs)
    shank_suffix = "" if row.histology_shank is None else f"_shank{int(row.histology_shank) + 1}"
    img_name = f"{probe_id}{shank_suffix}_image_space.json"
    gui_img = f"xyz_picks{shank_suffix}_image_space.json"

    (outputs.bregma_xyz / img_name).write_text(json.dumps(xyz_picks_image))
    gui_folder.mkdir(parents=True, exist_ok=True)
    (gui_folder / gui_img).write_text(json.dumps(xyz_picks_image))

    if emit_qc:
        _write_qc_probe_outputs(
            asset_info=asset_info,
            ng_data=ng_data,
            probe_id=probe_id,
            probe_pts=probe_pts,
            hist_stub_buggy=hist_stub_buggy,
            ibl_atlas=ibl_atlas,
            outputs=outputs,
            gui_folder=gui_folder,
            shank_suffix=shank_suffix,
        )

    return ProcessResult(
        probe_id=probe_id,
        recording_id=row.recording_id,
        wrote_files=True,
        skipped_reason=None,
    )


def _write_qc_probe_outputs(
    *,
    asset_info: AssetInfo,
    ng_data: dict,  # type: ignore[type-arg]
    probe_id: str,
    probe_pts: object,
    hist_stub_buggy: sitk.Image,
    ibl_atlas: AllenAtlas,
    outputs: OutputDirs,
    gui_folder: Path,
    shank_suffix: str,
) -> None:
    """Write GUI-unused QC outputs: SPIM/template/CCF FCSVs + CCF/bregma xyz-picks.

    These require the two ANTs point-warps (via ``hist_stub_buggy``) and the
    ``ibl_atlas``; nothing in the alignment workflow reads them.
    """
    # SPIM FCSV (no warp)
    create_slicer_fcsv(str(outputs.spim / f"{probe_id}.fcsv"), probe_pts, direction="LPS")

    # Image -> Template (points) via buggy pipeline transform
    probe_pt_dict_buggy, _ = neuroglancer_annotations_to_anatomical(
        ng_data,
        asset_info.zarr_volumes.registration,
        asset_info.zarr_volumes.metadata,
        layer_names=[probe_id],
        stub_image=hist_stub_buggy,
    )
    probe_pts_buggy = probe_pt_dict_buggy[probe_id]
    tx_list_pt_template = [
        str(asset_info.registration_dir_path / "ls_to_template_SyN_0GenericAffine.mat"),
        str(asset_info.registration_dir_path / "ls_to_template_SyN_1InverseWarp.nii.gz"),
    ]
    pts_template = apply_ants_transforms_to_point_arr(
        probe_pts_buggy,
        tx_list_pt_template,
        whichtoinvert=[True, False],
    )
    create_slicer_fcsv(str(outputs.template / f"{probe_id}.fcsv"), pts_template, direction="LPS")

    # Template -> CCF (points)
    pts_ccf = apply_ants_transforms_to_point_arr(
        probe_pts_buggy,
        asset_info.pipeline_registration_chains.pt_tx_str,
        whichtoinvert=asset_info.pipeline_registration_chains.pt_tx_inverted,
    )
    create_slicer_fcsv(str(outputs.ccf / f"{probe_id}.fcsv"), pts_ccf, direction="LPS")

    # IBL CCF/bregma xyz-picks (um), ML/AP/DV with signed flips
    ccf_mlapdv_um = convert_coordinate_system(1000.0 * pts_ccf, src_coord="LPS", dst_coord="RPI")
    bregma_mlapdv_um = 1_000_000.0 * ibl_atlas.ccf2xyz(ccf_mlapdv_um, ccf_order="mlapdv")
    xyz_picks_ccf = {"xyz_picks": bregma_mlapdv_um.tolist()}

    ccf_name = f"{probe_id}{shank_suffix}_ccf.json"
    gui_ccf = f"xyz_picks{shank_suffix}.json"
    (outputs.bregma_xyz / ccf_name).write_text(json.dumps(xyz_picks_ccf))
    (gui_folder / gui_ccf).write_text(json.dumps(xyz_picks_ccf))
