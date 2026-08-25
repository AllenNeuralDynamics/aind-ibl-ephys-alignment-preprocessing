"""Decide which physical frame an image-to-template transform expects its input in.

The SmartSPIM pipeline registers an anatomical volume whose origin comes from an
undocumented anchoring overlay, so points must be re-gridded onto that anchored
geometry before the ANTs chain -- the behaviour ``aind_zarr_utils`` exists to
mimic. Off-pipeline registrations ship a sidecar declaring the domain they were
actually computed in; there the re-grid is wrong, and the declared domain is what
must be honored. Getting this backwards is silent: the points still land
somewhere, and the region names they pick up still read plausibly.

Presence of a sidecar selects the branch, but the *domain it declares* is what is
used, so a future off-pipeline capsule that adopts the anchored geometry stays
correct without a change here.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import SimpleITK as sitk
from aind_ants_transform_sidecar import Domain, load_package
from aind_registration_utils.domains import ImageDomainAxisAligned, ImageHeader

logger = logging.getLogger(__name__)

#: Sidecar filenames looked for beside the image-to-template affine, in order.
#: Deliberately narrower than the alignment GUI's list, which also tries
#: ``<affine stem>.json``: a hard failure on a mismatched domain must not be
#: reachable by an unrelated JSON that happens to sit next to the transform.
REGISTRATION_SIDECAR_NAMES = (
    "ls_to_template_transform_information.json",
    "transform_information.json",
)

#: Domain agreement is checked in voxels, not in floating-point epsilon. The two
#: candidate frames are ~12 mm apart, while a genuine match is off by well under a
#: voxel: the sidecar's grid is the volume the registration ran on, the image here
#: has been resampled to another resolution, and ``bbox`` is voxel-*center*, so
#: even identical physical extents differ by half the spacing difference per side.
#: Two voxels of the coarser grid separates those scales by a factor of a hundred.
DOMAIN_TOLERANCE_VOXELS = 2.0

_AXES = ("L", "P", "S")


@dataclass(frozen=True)
class RegistrationFrame:
    """Which frame the registration's input points must be expressed in.

    Attributes
    ----------
    regrid_to_pipeline : bool
        Re-grid points onto the pipeline's anchored geometry before the ANTs
        chain. False when a sidecar documents the transform's own domain.
    reason : str
        Human-readable justification, recorded in logs.
    sidecar_path : Path or None
        The sidecar that decided it, when there was one.
    """

    regrid_to_pipeline: bool
    reason: str
    sidecar_path: Path | None = None


def find_registration_sidecar(registration_dir: Path) -> Path | None:
    """Return the transform sidecar in *registration_dir*, or None.

    Parameters
    ----------
    registration_dir : Path
        Directory holding the image-to-template transforms.

    Returns
    -------
    Path or None
        First matching sidecar, or None when the directory has none.
    """
    for name in REGISTRATION_SIDECAR_NAMES:
        candidate = registration_dir / name
        if candidate.is_file():
            return candidate
    return None


def resolve_registration_frame(registration_dir: Path, moving_image: sitk.Image) -> RegistrationFrame:
    """Decide whether points need the pipeline re-grid before the ANTs chain.

    Parameters
    ----------
    registration_dir : Path
        Directory holding the image-to-template transforms.
    moving_image : sitk.Image
        The anatomical image the points are expressed in. Only its header is
        read, so a stub is enough.

    Returns
    -------
    RegistrationFrame
        The decision and why it was made.

    Raises
    ------
    ValueError
        If a sidecar is present but declares no ``moving_domain``, or declares
        one that disagrees with *moving_image*. Neither has a safe recovery:
        falling back to the re-grid would apply the pipeline's compensation to a
        transform that is not the pipeline's, which is the original defect.
    """
    sidecar_path = find_registration_sidecar(registration_dir)
    if sidecar_path is None:
        return RegistrationFrame(
            regrid_to_pipeline=True,
            reason=f"no transform sidecar in {registration_dir}; assuming a pipeline registration",
        )

    sidecar = load_package(sidecar_path.read_text())
    if sidecar.moving_domain is None:
        raise ValueError(
            f"{sidecar_path} declares no moving_domain. A sidecar means the transform is off-pipeline, "
            "so the pipeline re-grid would be wrong, but without a domain there is nothing to honor "
            "instead. Regenerate the sidecar."
        )

    image_domain = _domain_from_image(moving_image)
    deltas = _bbox_deltas(sidecar.moving_domain, image_domain)
    tolerance = DOMAIN_TOLERANCE_VOXELS * max(max(sidecar.moving_domain.spacing_LPS), max(image_domain.spacing_LPS))
    if max(abs(d) for d in deltas.values()) > tolerance:
        raise ValueError(
            f"{sidecar_path} moving_domain disagrees with the anatomical image it is paired with.\n"
            f"  sidecar bbox: {_format_bbox(sidecar.moving_domain)}\n"
            f"  image bbox:   {_format_bbox(image_domain)}\n"
            f"  deltas (mm):  {_format_deltas(deltas)}\n"
            f"  tolerance:    {tolerance:.6g} mm ({DOMAIN_TOLERANCE_VOXELS} voxels)\n"
            "Either the transform was computed on a different volume, or the sidecar is stale."
        )

    reason = (
        f"{sidecar_path.name} declares the transform's own domain "
        f"(max bbox delta {max(abs(d) for d in deltas.values()):.6g} mm, tolerance {tolerance:.6g} mm); "
        "honoring it instead of re-gridding"
    )
    logger.info("Registration frame: %s", reason)
    return RegistrationFrame(regrid_to_pipeline=False, reason=reason, sidecar_path=sidecar_path)


def _domain_from_image(image: sitk.Image) -> Domain:
    """Describe *image*'s physical extent in the sidecar's own vocabulary."""
    return ImageDomainAxisAligned.from_header(ImageHeader.from_sitk(image)).to_sidecar()


def _bbox_deltas(left: Domain, right: Domain) -> dict[str, float]:
    """Signed per-bound differences between two domains' voxel-center bboxes."""
    deltas: dict[str, float] = {}
    for axis in _AXES:
        lo_l, hi_l = getattr(left.bbox, axis)
        lo_r, hi_r = getattr(right.bbox, axis)
        deltas[f"{axis}min"] = lo_l - lo_r
        deltas[f"{axis}max"] = hi_l - hi_r
    return deltas


def _format_bbox(domain: Domain) -> str:
    """Render a domain's bbox as ``L [lo, hi] P [...] S [...]``."""
    return " ".join(
        f"{axis} [{getattr(domain.bbox, axis)[0]:.4g}, {getattr(domain.bbox, axis)[1]:.4g}]" for axis in _AXES
    )


def _format_deltas(deltas: dict[str, float]) -> str:
    """Render per-bound deltas in a stable order."""
    return " ".join(f"{name}={value:+.4g}" for name, value in deltas.items())
