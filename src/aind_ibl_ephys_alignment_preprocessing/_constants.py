"""Internal constants for the preprocessing pipeline."""

from __future__ import annotations

_BLESSED_DIRECTION: str = "IRP"
"""DICOM orientation code expected by the IBL ephys alignment GUI."""

REGISTRATION_TRANSFORMS: tuple[str, ...] = (
    "ls_to_template_SyN_0GenericAffine.mat",
    "ls_to_template_SyN_1InverseWarp.nii.gz",
)
"""Transforms a registration directory must supply, whichever asset it lives in.

These are the individual->template half of the *point* chain, which is also what
resamples CCF volumes into image space. The forward warp is deliberately absent:
only the QC channels->CCF warp needs it, and standalone registrations are not
required to ship it.
"""

TEMPLATE_TO_CCF_TRANSFORMS: tuple[str, ...] = ("syn_0GenericAffine.mat", "syn_1InverseWarp.nii.gz")
"""Transforms a template->CCF directory must supply: the template->CCF half of the point chain."""

TEMPLATE_TO_CCF_FORWARD_WARP: str = "syn_1Warp.nii.gz"
"""The template->CCF forward warp, needed only by the QC channels->CCF warp."""

IMAGE_CHAIN_FORWARD_WARP: str = "ls_to_template_SyN_1Warp.nii.gz"
"""The individual->template forward warp, needed only by the QC channels->CCF warp."""

TEMPLATE_TO_CCF_BY_PIPELINE_NAME: dict[str, str] = {
    "spim_template_to_ccf_syn_0GenericAffine_25.mat": "syn_0GenericAffine.mat",
    "spim_template_to_ccf_syn_1InverseWarp_25.nii.gz": "syn_1InverseWarp.nii.gz",
    "spim_template_to_ccf_syn_1Warp_25.nii.gz": "syn_1Warp.nii.gz",
}
"""The pipeline's template->CCF file names and what a template->CCF directory calls the same transform."""
