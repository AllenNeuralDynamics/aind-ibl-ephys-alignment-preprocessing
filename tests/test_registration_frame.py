"""The sidecar decides the frame; disagreement is fatal, not recoverable.

The two candidate frames sit ~12 mm apart while a genuine match is off by a
fraction of a voxel, so these tests pin both ends of that gap: a resampled grid
still agrees, and a pipeline-anchored origin does not.
"""

from __future__ import annotations

import pytest

sitk = pytest.importorskip("SimpleITK")

from aind_ants_transform_sidecar import (  # noqa: E402
    SynTriplet,
    TransformSidecarV1,
    dump_package,
)

from aind_ibl_ephys_alignment_preprocessing.registration_frame import (  # noqa: E402
    REGISTRATION_SIDECAR_NAMES,
    resolve_registration_frame,
)

#: The SPIM-native box 776259's override was registered in, in LPS mm. One bound
#: per axis is exactly 0: the zarr stub's default origin, at every pyramid level.
NATIVE_BBOX = {"L": (-13.579, 0.0), "P": (0.0, 16.114), "S": (-7.664, 0.0)}

#: The pipeline's anchoring overlay changes the *origin* and nothing else -- same
#: voxel grid, same spacing -- so the anchored box is the native one translated.
#: The offset is illustrative; its magnitude (millimetres) is what matters.
PIPELINE_OFFSET = (12.068, -1.5, -1.5)
PIPELINE_BBOX = {
    axis: (lo + shift, hi + shift) for (axis, (lo, hi)), shift in zip(NATIVE_BBOX.items(), PIPELINE_OFFSET, strict=True)
}


def _header(bbox, spacing):
    """A 1x1x1 header carrier -- exactly what ``as_sitk_stub`` hands production."""
    img = sitk.Image((1, 1, 1), sitk.sitkUInt8)
    img.SetSpacing((spacing,) * 3)
    img.SetOrigin(tuple(lo for lo, _ in bbox.values()))
    return img


def _pair(native=NATIVE_BBOX, pipeline=PIPELINE_BBOX, spacing=0.0144):
    """The (native, pipeline, size) triple the orchestrators resolve a frame from."""
    size = tuple(max(2, int(round((hi - lo) / spacing)) + 1) for lo, hi in native.values())
    return _header(native, spacing), _header(pipeline, spacing), size


def _domain_of(bbox=NATIVE_BBOX, spacing=0.0144):
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import _domain_from_image

    size = tuple(max(2, int(round((hi - lo) / spacing)) + 1) for lo, hi in bbox.values())
    return _domain_from_image(_header(bbox, spacing), size)


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
    frame = resolve_registration_frame(tmp_path, *_pair())
    assert frame.regrid_to_pipeline
    assert frame.sidecar_path is None
    assert "no transform sidecar" in frame.reason


def test_a_native_frame_sidecar_is_honored(tmp_path):
    path = _write_sidecar(tmp_path, _domain_of(NATIVE_BBOX))
    frame = resolve_registration_frame(tmp_path, *_pair())
    assert not frame.regrid_to_pipeline
    assert frame.sidecar_path == path


def test_a_pipeline_frame_sidecar_keeps_the_regrid(tmp_path):
    """A capsule that adopted the anchored geometry stays correct without a change here."""
    _write_sidecar(tmp_path, _domain_of(PIPELINE_BBOX))
    assert resolve_registration_frame(tmp_path, *_pair()).regrid_to_pipeline


def test_a_resampled_grid_does_not_change_the_answer(tmp_path):
    """The decision is relative, so a different resolution cannot flip it."""
    _write_sidecar(tmp_path, _domain_of(NATIVE_BBOX, spacing=0.032))
    assert not resolve_registration_frame(tmp_path, *_pair(spacing=0.0018)).regrid_to_pipeline


def test_a_sidecar_for_another_volume_is_fatal(tmp_path):
    """Nearest-of-two always names a winner; the sanity bound is what rejects one."""
    far = {axis: (lo + 40.0, hi + 40.0) for axis, (lo, hi) in NATIVE_BBOX.items()}
    _write_sidecar(tmp_path, _domain_of(far))
    with pytest.raises(ValueError, match="matches neither frame"):
        resolve_registration_frame(tmp_path, *_pair())


def test_the_mismatch_message_names_both_candidates(tmp_path):
    far = {axis: (lo + 40.0, hi + 40.0) for axis, (lo, hi) in NATIVE_BBOX.items()}
    _write_sidecar(tmp_path, _domain_of(far))
    with pytest.raises(ValueError) as excinfo:
        resolve_registration_frame(tmp_path, *_pair())
    message = str(excinfo.value)
    assert "native bbox" in message and "pipeline bbox" in message and "plausible within" in message


def test_sidecar_without_a_domain_is_fatal(tmp_path):
    _write_sidecar(tmp_path, None)
    with pytest.raises(ValueError, match="declares no moving_domain"):
        resolve_registration_frame(tmp_path, *_pair())


def test_unrelated_json_beside_the_transform_is_ignored(tmp_path):
    (tmp_path / "ls_to_template_SyN_0GenericAffine.json").write_text('{"note": "not a sidecar"}')
    assert resolve_registration_frame(tmp_path, *_pair()).regrid_to_pipeline


def test_alternate_sidecar_name_is_found(tmp_path):
    _write_sidecar(tmp_path, _domain_of(NATIVE_BBOX), name=REGISTRATION_SIDECAR_NAMES[1])
    assert not resolve_registration_frame(tmp_path, *_pair()).regrid_to_pipeline


def test_the_size_comes_from_the_caller_not_the_stub(tmp_path):
    """Regression: the orchestrators pass 1x1x1 header carriers with no pixels.

    Reading ``GetSize()`` off one collapses both domains to a point at their own
    origins, which flips the comparison. Caught pre-flighting 776259 offline.
    """
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import _domain_from_image

    native, pipeline, size = _pair()
    assert native.GetSize() == (1, 1, 1), "the fixture must reproduce the header-only stub"
    assert _domain_from_image(native, size).bbox.L != _domain_from_image(native, (1, 1, 1)).bbox.L

    _write_sidecar(tmp_path, _domain_of(NATIVE_BBOX))
    assert not resolve_registration_frame(tmp_path, native, pipeline, size).regrid_to_pipeline


def test_the_sanity_bound_never_decides_the_frame():
    """It sits far above any real disagreement and far below the frames' separation."""
    from aind_ibl_ephys_alignment_preprocessing.registration_frame import _bbox_distance, _sanity_tolerance

    native, pipeline = _domain_of(NATIVE_BBOX, 0.0018), _domain_of(PIPELINE_BBOX, 0.0018)
    sidecar = _domain_of(NATIVE_BBOX, 0.0144)  # the 8x pyramid step, as 776259 has
    bound = _sanity_tolerance(sidecar, native)
    assert _bbox_distance(sidecar, native) < bound < _bbox_distance(sidecar, pipeline)
