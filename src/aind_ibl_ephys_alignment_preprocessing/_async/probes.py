"""Async per-probe manifest row processing."""

from __future__ import annotations

import asyncio
import json
import logging
from pathlib import Path
from typing import Any

import SimpleITK as sitk
from aind_anatomical_utils.coordinate_systems import convert_coordinate_system
from aind_ephys_ibl_gui_conversion.histology import create_slicer_fcsv
from aind_registration_utils.ants import apply_ants_transforms_to_point_arr
from aind_zarr_utils.neuroglancer import neuroglancer_annotations_to_anatomical
from iblatlas.atlas import AllenAtlas

from aind_ibl_ephys_alignment_preprocessing._async.concurrency import (
    Limits,
    io_to_thread_on,
    to_thread_logged,
)
from aind_ibl_ephys_alignment_preprocessing.types import (
    AssetInfo,
    ManifestRow,
    OutputDirs,
    ProcessResult,
)

logger = logging.getLogger(__name__)


async def read_json_in_thread(path: Path, limits: Limits) -> Any:
    """Read a JSON file in a background thread."""

    def _read() -> Any:
        with open(path) as f:
            return json.load(f)

    return await io_to_thread_on(limits, str(path), _read)


async def process_manifest_row_async(
    row: ManifestRow,
    asset_info: AssetInfo,
    hist_stub: sitk.Image,
    hist_stub_buggy: sitk.Image,
    ibl_atlas: AllenAtlas,
    outputs: OutputDirs,
    limits: Limits,
    data_root: Path,
) -> ProcessResult:
    """Async end-to-end processing for a single manifest row.

    Parameters
    ----------
    row : ManifestRow
        Manifest row describing the probe.
    asset_info : AssetInfo
        Asset metadata.
    hist_stub : sitk.Image
        Histology stub with correct domain.
    hist_stub_buggy : sitk.Image
        Histology stub in pipeline (buggy) domain.
    ibl_atlas : AllenAtlas
        Allen atlas for CCF-to-bregma conversion.
    outputs : OutputDirs
        Output directory tree.
    limits : Limits
        Concurrency limits.
    data_root : Path
        Root directory for input data.

    Returns
    -------
    ProcessResult
        Summary of success/skip.
    """
    logger.info("[Probe %s] Starting processing for recording: %s", row.probe_id, row.recording_id)

    ext = "json" if row.annotation_format == "json" else None
    if ext is None:
        return ProcessResult(row.probe_id, row.recording_id, False, "Only JSON annotations supported")
    pattern = f"*/{row.probe_file}.{ext}"
    ann_path = next(data_root.glob(pattern), None)
    probe_id = str(row.probe_id)
    if ann_path is None:
        return ProcessResult(probe_id, str(row.sorted_recording), False, f"Annotation not found: {pattern}")

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

    p_img = outputs.bregma_xyz / img_name
    p_ccf = outputs.bregma_xyz / ccf_name
    if p_img.exists() and p_ccf.exists():
        return ProcessResult(probe_id, str(row.sorted_recording), True, "Already processed")

    # Load NG points
    ng_data = await read_json_in_thread(ann_path, limits)
    anno_zarr = asset_info.zarr_volumes.registration
    metadata = asset_info.zarr_volumes.metadata
    async with asyncio.TaskGroup() as tg:
        probe_pt_dict_task = tg.create_task(
            to_thread_logged(
                neuroglancer_annotations_to_anatomical,
                ng_data,
                anno_zarr,
                metadata,
                layer_names=[probe_id],
                stub_image=hist_stub,
            ),
            name=f"ng-to-anat-{probe_id}",
        )
        probe_pt_dict_buggy_task = tg.create_task(
            to_thread_logged(
                neuroglancer_annotations_to_anatomical,
                ng_data,
                anno_zarr,
                metadata,
                layer_names=[probe_id],
                stub_image=hist_stub_buggy,
            ),
            name=f"ng-to-anat-buggy-{probe_id}",
        )
    probe_pt_dict, _ = probe_pt_dict_task.result()
    probe_pt_dict_buggy, _ = probe_pt_dict_buggy_task.result()
    probe_pts = probe_pt_dict.get(probe_id, None)
    if probe_pts is None:
        return ProcessResult(probe_id, str(row.sorted_recording), False, f"Probe points not found: {probe_id}")
    probe_pts_buggy = probe_pt_dict_buggy[probe_id]

    # Write SPIM FCSV
    await io_to_thread_on(
        limits,
        str(outputs.spim),
        create_slicer_fcsv,
        str(outputs.spim / f"{probe_id}.fcsv"),
        probe_pts,
        direction="LPS",
    )

    # Image -> Template
    tx_list_pt_template = [
        str(asset_info.registration_dir_path / "ls_to_template_SyN_0GenericAffine.mat"),
        str(asset_info.registration_dir_path / "ls_to_template_SyN_1InverseWarp.nii.gz"),
    ]
    tx_list_pt_template_invert = [True, False]
    pts_template = await to_thread_logged(
        apply_ants_transforms_to_point_arr,
        probe_pts_buggy,
        tx_list_pt_template,
        whichtoinvert=tx_list_pt_template_invert,
    )
    await io_to_thread_on(
        limits,
        str(outputs.template),
        create_slicer_fcsv,
        str(outputs.template / f"{probe_id}.fcsv"),
        pts_template,
        direction="LPS",
    )

    # Template -> CCF
    pts_ccf = await to_thread_logged(
        apply_ants_transforms_to_point_arr,
        probe_pts_buggy,
        asset_info.pipeline_registration_chains.pt_tx_str,
        whichtoinvert=asset_info.pipeline_registration_chains.pt_tx_inverted,
    )
    await io_to_thread_on(
        limits,
        str(outputs.ccf),
        create_slicer_fcsv,
        str(outputs.ccf / f"{row.probe_id}.fcsv"),
        pts_ccf,
        direction="LPS",
    )

    # IBL xyz-picks
    ccf_mlapdv_um = convert_coordinate_system(1000.0 * pts_ccf, src_coord="LPS", dst_coord="RPI")
    bregma_mlapdv_um = 1_000_000.0 * ibl_atlas.ccf2xyz(ccf_mlapdv_um, ccf_order="mlapdv")
    xyz_img = 1000.0 * convert_coordinate_system(probe_pts, src_coord="LPS", dst_coord="RAS")

    xyz_picks_image = {"xyz_picks": xyz_img.tolist()}
    xyz_picks_ccf = {"xyz_picks": bregma_mlapdv_um.tolist()}

    # Write bregma_xyz JSONs
    img_json_str = json.dumps(xyz_picks_image)
    ccf_json_str = json.dumps(xyz_picks_ccf)
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            io_to_thread_on(limits, str(p_img), p_img.write_text, img_json_str), name=f"write-bregma-img-{probe_id}"
        )
        tg.create_task(
            io_to_thread_on(limits, str(p_ccf), p_ccf.write_text, ccf_json_str), name=f"write-bregma-ccf-{probe_id}"
        )

    # Per-recording GUI artifacts
    await io_to_thread_on(limits, str(gui_folder), gui_folder.mkdir, parents=True, exist_ok=True)
    async with asyncio.TaskGroup() as tg:
        tg.create_task(
            io_to_thread_on(limits, str(gui_folder), (gui_folder / gui_img).write_text, img_json_str),
            name=f"write-gui-img-{probe_id}",
        )
        tg.create_task(
            io_to_thread_on(limits, str(gui_folder), (gui_folder / gui_ccf).write_text, ccf_json_str),
            name=f"write-gui-ccf-{probe_id}",
        )

    logger.info("[Probe %s] Completed: wrote xyz_picks to 4 locations", row.probe_id)
    return ProcessResult(
        probe_id=str(row.probe_id),
        recording_id=row.recording_id,
        wrote_files=True,
        skipped_reason=None,
    )


async def process_manifest_row_limit_async(
    row: ManifestRow,
    asset_info: AssetInfo,
    hist_stub: sitk.Image,
    hist_stub_buggy: sitk.Image,
    ibl_atlas: AllenAtlas,
    outputs: OutputDirs,
    limits: Limits,
    data_root: Path,
) -> ProcessResult:
    """Rate-limited wrapper around :func:`process_manifest_row_async`."""
    async with limits.manifest_rows:
        return await process_manifest_row_async(
            row, asset_info, hist_stub, hist_stub_buggy, ibl_atlas, outputs, limits, data_root
        )


async def process_manifest_row_safe_async(
    row: ManifestRow,
    asset_info: AssetInfo,
    hist_stub: sitk.Image,
    hist_stub_buggy: sitk.Image,
    ibl_atlas: AllenAtlas,
    outputs: OutputDirs,
    limits: Limits,
    data_root: Path,
) -> ProcessResult:
    """Error-catching wrapper that converts exceptions to :class:`ProcessResult`."""
    try:
        return await process_manifest_row_limit_async(
            row, asset_info, hist_stub, hist_stub_buggy, ibl_atlas, outputs, limits, data_root
        )
    except Exception as e:
        logger.exception("Row failed: probe=%s recording=%s", row.probe_id, row.recording_id)
        return ProcessResult(
            probe_id=str(row.probe_id),
            recording_id=row.recording_id,
            wrote_files=False,
            skipped_reason=f"{type(e).__name__}: {e}",
        )
