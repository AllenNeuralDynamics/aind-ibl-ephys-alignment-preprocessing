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

#: Domain agreement is checked in voxels, not in floating-point epsilon, because
#: the sidecar's grid is the volume the registration ran on and the image here is
#: a different resolution of the same acquisition.
#:
#: Resampling moves the *extent*, never the *placement*: a zarr stub's default
#: origin is ``(0, 0, 0)`` at every level (``aind_zarr_utils.zarr``, origin_type
#: "none"), so one bbox bound per axis is exactly 0 whatever level was read, and
#: only the far bound carries the level. The pipeline-anchored frame moves that
#: zero corner by ~12 mm, which is the thing being detected. So what varies with
#: resampling is not what discriminates.
#:
#: The far bound's disagreement is computable rather than empirical. Both grids
#: pin voxel 0's *centre* at the origin, so for ``N * s_fine == n * s_coarse`` the
#: bounds are ``(N-1) * s_fine`` and ``(n-1) * s_coarse`` -- a difference of
#: exactly ``s_coarse - s_fine``, plus up to one coarse voxel when the downsample
#: does not divide evenly. That bound is just under two coarse voxels, so two is
#: the threshold with no headroom (1.07x on the campaign's own grids); four gives
#: 2.1x while still sitting ~190x below the millimetres it must catch.
#:
#: Converting to voxel *edges* does not remove this -- centre-pinned grids are
#: neither centre- nor corner-aligned across levels. An exact check means building
#: the comparison stub at the sidecar's own level, which is only possible when the
#: registration read a raw pyramid level.
#:
#: Comparing all six bounds rather than just the zero corner keeps this
#: independent of that origin convention, and still catches a sidecar paired with
#: an altogether different volume.
DOMAIN_TOLERANCE_VOXELS = 4.0

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


def resolve_registration_frame(
    registration_dir: Path, moving_image: sitk.Image, size_ijk: tuple[int, int, int]
) -> RegistrationFrame:
    """Decide whether points need the pipeline re-grid before the ANTs chain.

    Parameters
    ----------
    registration_dir : Path
        Directory holding the image-to-template transforms.
    moving_image : sitk.Image
        The anatomical image the points are expressed in. Only origin, spacing
        and direction are read from it.
    size_ijk : tuple[int, int, int]
        Voxel counts in SimpleITK (x, y, z) order. Required rather than taken
        from *moving_image*, because the anatomical images here are header-only
        stubs: ``AnatomicalHeader.as_sitk_stub`` returns a 1x1x1 image, whose own
        ``GetSize`` would collapse the domain to a single point and fail every
        comparison. ``base_and_pipeline_anatomical_stub`` returns this alongside
        the stubs.

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

    image_domain = _domain_from_image(moving_image, size_ijk)
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


def _domain_from_image(image: sitk.Image, size_ijk: tuple[int, int, int]) -> Domain:
    """Describe *image*'s physical extent in the sidecar's own vocabulary."""
    return ImageDomainAxisAligned.from_header(ImageHeader.from_sitk(image, size_ijk)).to_sidecar()


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
