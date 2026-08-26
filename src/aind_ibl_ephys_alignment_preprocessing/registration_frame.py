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

#: Floor for the plausibility check, in mm. The frame is chosen by which candidate
#: is *nearer*, so this decides nothing -- it only rejects a sidecar that belongs to
#: some other volume. Set well above any legitimate resampling disagreement (tens of
#: microns) and well below the millimetres that separate the two frames.
DOMAIN_SANITY_MM = 1.0

#: Voxel-scaled companion to :data:`DOMAIN_SANITY_MM`, so a registration run at a
#: very coarse resolution is not judged against a floor finer than its own grid.
#: The larger of the two applies.
DOMAIN_TOLERANCE_VOXELS = 4.0

_AXES = ("L", "P", "S")


@dataclass(frozen=True)
class RegistrationFrame:
    """Which frame the registration's input points must be expressed in.

    Attributes
    ----------
    regrid_to_pipeline : bool
        Re-grid points onto the pipeline's anchored geometry before the ANTs
        chain. False when a sidecar shows the transform was trained in the
        volume's own frame.
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
    registration_dir: Path,
    native_image: sitk.Image,
    pipeline_image: sitk.Image,
    size_ijk: tuple[int, int, int],
) -> RegistrationFrame:
    """Decide whether points need the pipeline re-grid before the ANTs chain.

    The choice is made *relatively* -- which of the two candidate frames the
    sidecar's declared domain is nearer to -- so no threshold decides it and a
    mis-sized tolerance cannot pick the wrong frame. The absolute check that
    follows only asks whether the winner is plausible at all, which is what
    catches a sidecar paired with a different volume.

    Parameters
    ----------
    registration_dir : Path
        Directory holding the image-to-template transforms.
    native_image, pipeline_image : sitk.Image
        The anatomical geometry as the volume's own header gives it, and as the
        pipeline's anchoring overlay gives it. Only origin, spacing and
        direction are read.
    size_ijk : tuple[int, int, int]
        Voxel counts in SimpleITK (x, y, z) order, shared by both. Required
        rather than read from the images: these are header-only stubs, and
        ``AnatomicalHeader.as_sitk_stub`` returns a 1x1x1 image whose own
        ``GetSize`` would collapse each domain to a point.

    Returns
    -------
    RegistrationFrame
        The decision and why it was made.

    Raises
    ------
    ValueError
        If a sidecar declares no ``moving_domain``, or declares one that matches
        neither candidate. Neither has a safe recovery: falling back to the
        re-grid would apply the pipeline's compensation to a transform that is
        not the pipeline's, which is the original defect.
    """
    sidecar_path = find_registration_sidecar(registration_dir)
    if sidecar_path is None:
        return RegistrationFrame(
            regrid_to_pipeline=True,
            reason=f"no transform sidecar in {registration_dir}; assuming a pipeline registration",
        )

    sidecar = load_package(sidecar_path.read_text())
    moving = sidecar.moving_domain
    if moving is None:
        raise ValueError(
            f"{sidecar_path} declares no moving_domain. A sidecar means the transform is off-pipeline, "
            "so the pipeline re-grid would be wrong, but without a domain there is nothing to honor "
            "instead. Regenerate the sidecar."
        )

    candidates = {
        "native": _domain_from_image(native_image, size_ijk),
        "pipeline": _domain_from_image(pipeline_image, size_ijk),
    }
    distances = {name: _bbox_distance(moving, d) for name, d in candidates.items()}
    winner = min(distances, key=lambda name: distances[name])

    # Both frames describe the same acquisition, so an exact tie means the
    # anchoring overlay was a no-op and the branches are indistinguishable.
    tolerance = _sanity_tolerance(moving, candidates[winner])
    if distances[winner] > tolerance:
        raise ValueError(
            f"{sidecar_path} moving_domain matches neither frame of the anatomical image it is paired with.\n"
            f"  sidecar bbox:  {_format_bbox(moving)}\n"
            f"  native bbox:   {_format_bbox(candidates['native'])}  (off by {distances['native']:.4g} mm)\n"
            f"  pipeline bbox: {_format_bbox(candidates['pipeline'])}  (off by {distances['pipeline']:.4g} mm)\n"
            f"  plausible within {tolerance:.4g} mm\n"
            "Either the transform was computed on a different volume, or the sidecar is stale."
        )

    reason = (
        f"{sidecar_path.name} moving_domain is nearest the {winner} frame "
        f"(native {distances['native']:.4g} mm, pipeline {distances['pipeline']:.4g} mm; "
        f"plausible within {tolerance:.4g} mm)"
    )
    logger.info("Registration frame: %s", reason)
    return RegistrationFrame(
        regrid_to_pipeline=winner == "pipeline",
        reason=reason,
        sidecar_path=sidecar_path,
    )


def _domain_from_image(image: sitk.Image, size_ijk: tuple[int, int, int]) -> Domain:
    """Describe *image*'s physical extent in the sidecar's own vocabulary."""
    return ImageDomainAxisAligned.from_header(ImageHeader.from_sitk(image, size_ijk)).to_sidecar()


def _bbox_distance(left: Domain, right: Domain) -> float:
    """Largest disagreement between two domains' voxel-center bboxes, in mm.

    Shape is deliberately not compared: the sidecar's grid is the resolution the
    registration ran at, which legitimately differs from any volume here.
    """
    return max(abs(v) for v in _bbox_deltas(left, right).values())


def _sanity_tolerance(left: Domain, right: Domain) -> float:
    """How far apart two domains of the same volume may legitimately sit, in mm.

    Both grids pin voxel 0's *centre* at the origin, so for ``N * s_fine ==
    n * s_coarse`` the far bounds differ by exactly ``s_coarse - s_fine``, plus up
    to one coarse voxel when the downsample does not divide evenly -- just under
    two coarse voxels. This is no longer what picks the frame, so it is set well
    above that bound; it exists to reject a sidecar belonging to another volume.
    """
    coarsest = max(max(left.spacing_LPS), max(right.spacing_LPS))
    return max(DOMAIN_SANITY_MM, DOMAIN_TOLERANCE_VOXELS * coarsest)


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
