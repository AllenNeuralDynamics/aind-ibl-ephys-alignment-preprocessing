"""The sidecar decides the frame; disagreement is fatal, not recoverable.

The two candidate frames sit ~12 mm apart while a genuine match is off by a
fraction of a voxel, so these tests pin both ends of that gap: a resampled grid
still agrees, and a pipeline-anchored origin does not.
"""

from __future__ import annotations

import pytest

sitk = pytest.importorskip("SimpleITK")

from aind_ants_transform_sidecar import (  # noqa: E402
    BBox,
    Domain,
    SynTriplet,
    TransformSidecarV1,
    dump_package,
)

from aind_ibl_ephys_alignment_preprocessing.registration_frame import (  # noqa: E402
    REGISTRATION_SIDECAR_NAMES,
    resolve_registration_frame,
)

#: The SPIM-native box 776259's override was registered in, in LPS mm.
NATIVE_BBOX = {"L": (-13.579, 0.0), "P": (0.0, 16.114), "S": (-7.664, 0.0)}

#: The same brain placed as the pipeline anchors it. Illustrative rather than
#: measured -- what matters is the ~12 mm displacement, and that no bound is 0.
PIPELINE_BBOX = {"L": (-1.511, 12.068), "P": (-1.5, 14.614), "S": (-9.164, 1.5)}


def _stub(bbox=NATIVE_BBOX, spacing=0.0144):
    """A 1x1x1 header carrier plus its real size -- exactly what production passes.

    ``AnatomicalHeader.as_sitk_stub`` allocates no pixels, so the size travels
    beside the image rather than inside it.
    """
    size = tuple(max(2, int(round((hi - lo) / spacing)) + 1) for lo, hi in bbox.values())
    img = sitk.Image((1, 1, 1), sitk.sitkUInt8)
    img.SetSpacing((spacing,) * 3)
    img.SetOrigin(tuple(lo for lo, _ in bbox.values()))
    return img, size


def _domain_of(stub):
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import _domain_from_image

    img, size = stub
    return _domain_from_image(img, size)


def _write_sidecar(directory, domain, name=REGISTRATION_SIDECAR_NAMES[0]):
    transform = SynTriplet(
        affine="ls_to_template_SyN_0GenericAffine.mat",
        warp="ls_to_template_SyN_1Warp.nii.gz",
        inverse_warp="ls_to_template_SyN_1InverseWarp.nii.gz",
    )
    kwargs = {} if domain is None else {"fixed_domain": domain, "moving_domain": domain}
    sidecar = TransformSidecarV1(transform=transform, **kwargs)
    path = directory / name
    path.write_text(dump_package(sidecar))
    return path


def test_no_sidecar_keeps_the_pipeline_regrid(tmp_path):
    frame = resolve_registration_frame(tmp_path, *_stub())
    assert frame.regrid_to_pipeline
    assert frame.sidecar_path is None
    assert "no transform sidecar" in frame.reason


def test_matching_sidecar_is_honored(tmp_path):
    stub = _stub()
    path = _write_sidecar(tmp_path, _domain_of(stub))
    frame = resolve_registration_frame(tmp_path, *stub)
    assert not frame.regrid_to_pipeline
    assert frame.sidecar_path == path


def test_resampled_grid_still_agrees(tmp_path):
    """A different voxel size over the same box is a match, not a mismatch."""
    _write_sidecar(tmp_path, _domain_of(_stub(spacing=0.032)))
    frame = resolve_registration_frame(tmp_path, *_stub(spacing=0.0144))
    assert not frame.regrid_to_pipeline


def test_pipeline_anchored_sidecar_is_fatal(tmp_path):
    """The failure this whole module exists to catch, at its real magnitude."""
    _write_sidecar(tmp_path, _domain_of(_stub(bbox=PIPELINE_BBOX)))
    with pytest.raises(ValueError, match="disagrees with the anatomical image"):
        resolve_registration_frame(tmp_path, *_stub(bbox=NATIVE_BBOX))


def test_mismatch_message_names_both_boxes_and_the_deltas(tmp_path):
    _write_sidecar(tmp_path, _domain_of(_stub(bbox=PIPELINE_BBOX)))
    with pytest.raises(ValueError) as excinfo:
        resolve_registration_frame(tmp_path, *_stub(bbox=NATIVE_BBOX))
    message = str(excinfo.value)
    assert "sidecar bbox" in message and "image bbox" in message
    assert "Lmin=" in message and "tolerance" in message


def test_sidecar_without_a_domain_is_fatal(tmp_path):
    """Presence says off-pipeline, so the re-grid is wrong -- and there is
    nothing to honor instead. Falling back would be the original defect."""
    _write_sidecar(tmp_path, None)
    with pytest.raises(ValueError, match="declares no moving_domain"):
        resolve_registration_frame(tmp_path, *_stub())


def test_unrelated_json_beside_the_transform_is_ignored(tmp_path):
    """A hard failure must not be reachable by a file that is not a sidecar."""
    (tmp_path / "ls_to_template_SyN_0GenericAffine.json").write_text('{"note": "not a sidecar"}')
    assert resolve_registration_frame(tmp_path, *_stub()).regrid_to_pipeline


def test_alternate_sidecar_name_is_found(tmp_path):
    stub = _stub()
    _write_sidecar(tmp_path, _domain_of(stub), name=REGISTRATION_SIDECAR_NAMES[1])
    assert not resolve_registration_frame(tmp_path, *stub).regrid_to_pipeline


def test_bbox_uses_a_domain_the_sidecar_can_round_trip():
    """Guards the vocabulary: our Domain must compare against a parsed one."""
    domain = _domain_of(_stub())
    assert isinstance(domain, Domain)
    assert isinstance(domain.bbox, BBox)


def test_an_eightfold_downsample_still_agrees():
    """The campaign's real case: sidecar at pyramid level 3, stub at level 0.

    Both grids pin voxel 0's centre at the origin, so their far bounds differ by
    ``s_coarse - s_fine`` -- 14 um for 776259's 16 um sidecar. Pinned here because
    the tolerance is sized against that number, not measured against an asset.
    """
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import DOMAIN_TOLERANCE_VOXELS

    s_coarse = 0.016
    s_fine = s_coarse / 8
    worst_case = (s_coarse - s_fine) + s_coarse  # + up to one coarse voxel of truncation
    assert DOMAIN_TOLERANCE_VOXELS * s_coarse > worst_case * 2, "tolerance has lost its headroom"


def test_a_millimetre_scale_offset_is_never_within_tolerance(tmp_path):
    """The margin that matters: the signal is ~200x the threshold."""
    shifted = {axis: (lo + 1.0, hi + 1.0) for axis, (lo, hi) in NATIVE_BBOX.items()}
    _write_sidecar(tmp_path, _domain_of(_stub(bbox=shifted)))
    with pytest.raises(ValueError, match="disagrees with the anatomical image"):
        resolve_registration_frame(tmp_path, *_stub(bbox=NATIVE_BBOX))


def test_the_size_comes_from_the_caller_not_the_stub(tmp_path):
    """Regression: the orchestrators pass a 1x1x1 header carrier.

    ``base_and_pipeline_anatomical_stub`` returns images with no pixel buffer, so
    reading ``GetSize()`` off one collapses the domain to a single voxel at the
    origin and every real asset fails the comparison. Caught pre-flighting 776259
    against its actual sidecar, not by a run.
    """
    img, size = _stub()
    assert img.GetSize() == (1, 1, 1), "the fixture must reproduce the header-only stub"

    honest = _domain_from_image_public(img, size)
    collapsed = _domain_from_image_public(img, (1, 1, 1))
    assert honest.bbox.L != collapsed.bbox.L

    _write_sidecar(tmp_path, honest)
    assert not resolve_registration_frame(tmp_path, img, size).regrid_to_pipeline


def _domain_from_image_public(img, size):
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import _domain_from_image

    return _domain_from_image(img, size)
