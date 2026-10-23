"""Templates and transform sidecars small enough to write in a test."""

from __future__ import annotations

from pathlib import Path

import SimpleITK as sitk
from aind_ants_transform_sidecar import Domain, SynTriplet, TransformSidecarV1, dump_package

from aind_ibl_ephys_alignment_preprocessing.registration_frame import REGISTRATION_SIDECAR_NAMES, _domain_from_file

#: An iDISCO-like template: same grid as the stand-in pipeline template, millimetres away.
OTHER_TEMPLATE_ORIGIN = (6.0, -4.0, 3.0)


def write_template(path: Path, origin: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> Domain:
    """Write a tiny 25 um template volume and return its domain."""
    path.parent.mkdir(parents=True, exist_ok=True)
    img = sitk.Image((40, 50, 30), sitk.sitkUInt8)
    img.SetSpacing((0.025,) * 3)
    img.SetOrigin(origin)
    sitk.WriteImage(img, str(path))
    return _domain_from_file(path)


def write_sidecar(
    directory: Path,
    fixed: Domain,
    moving: Domain | None = None,
    name: str = REGISTRATION_SIDECAR_NAMES[0],
    prefix: str = "ls_to_template_SyN",
) -> Path:
    """Write a SyN sidecar; *moving* defaults to *fixed* where the frame is not under test."""
    directory.mkdir(parents=True, exist_ok=True)
    transform = SynTriplet(
        affine=f"{prefix}_0GenericAffine.mat",
        warp=f"{prefix}_1Warp.nii.gz",
        inverse_warp=f"{prefix}_1InverseWarp.nii.gz",
    )
    sidecar = TransformSidecarV1(transform=transform, fixed_domain=fixed, moving_domain=moving or fixed)
    path = directory / name
    path.write_text(dump_package(sidecar))
    return path
