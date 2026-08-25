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
