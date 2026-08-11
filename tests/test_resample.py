"""The output grid must be the same for every mouse, without inventing signal.

Spacing here is **millimetres**, as ``aind_zarr_utils`` produces it. Fixtures
written in micrometres pass against code that makes the same mistake, so they are
deliberately in the pipeline's real units.
"""

from __future__ import annotations

import numpy as np
import pytest

sitk = pytest.importorskip("SimpleITK")

from aind_ibl_ephys_alignment_preprocessing._async.histology import (  # noqa: E402
    resample_to_isotropic,
)

UM = 1.0 / 1000.0  # spacing is in mm; targets are quoted in um


def _volume(spacing, size=(64, 64, 64), value=1000):
    img = sitk.GetImageFromArray(np.full(size[::-1], value, dtype=np.uint16))
    img.SetSpacing(spacing)
    img.SetOrigin((3.0, -5.0, 7.0))
    return img


@pytest.mark.parametrize(
    "spacing",
    [
        (14.4 * UM, 14.4 * UM, 16.0 * UM),  # pre-September stitching
        (28.8 * UM, 28.8 * UM, 32.0 * UM),  # September onward
        (32.0 * UM, 32.0 * UM, 16.0 * UM),  # anisotropic the other way
    ],
)
def test_every_vintage_lands_on_one_grid(spacing):
    """The point of the resample: stitching vintage must stop setting the grid."""
    out = resample_to_isotropic(_volume(spacing), 30.0, "test")

    assert out.GetSpacing() == pytest.approx((30 * UM, 30 * UM, 30 * UM))


def test_physical_extent_and_origin_are_preserved():
    """Resampling changes sampling, not where the data sits in the world."""
    img = _volume((14.4 * UM, 14.4 * UM, 16.0 * UM))
    out = resample_to_isotropic(img, 30.0, "test")

    assert out.GetOrigin() == img.GetOrigin()
    assert out.GetDirection() == img.GetDirection()
    for axis in range(3):
        before = img.GetSize()[axis] * img.GetSpacing()[axis]
        after = out.GetSize()[axis] * out.GetSpacing()[axis]
        assert after == pytest.approx(before, rel=0.02)


def test_already_coarser_volumes_are_left_alone():
    """The target must never manufacture resolution."""
    img = _volume((50.0 * UM, 50.0 * UM, 50.0 * UM))

    out = resample_to_isotropic(img, 30.0, "test")

    assert out is img


def test_dtype_survives_the_float_round_trip():
    """Smoothing needs floats; the channels must not become float on disk."""
    img = _volume((14.4 * UM, 14.4 * UM, 16.0 * UM))

    out = resample_to_isotropic(img, 30.0, "test")

    assert out.GetPixelID() == img.GetPixelID() == sitk.sitkUInt16


def test_downsampling_is_anti_aliased():
    """sitk.Resample point-samples; without a low-pass, fine detail aliases."""
    # A single-voxel-period checkerboard: pure content above the new Nyquist.
    size = 96
    arr = np.indices((size, size, size)).sum(axis=0) % 2
    img = sitk.GetImageFromArray((arr * 4000).astype(np.uint16))
    img.SetSpacing((10.0 * UM, 10.0 * UM, 10.0 * UM))

    out = resample_to_isotropic(img, 30.0, "test")
    values = sitk.GetArrayFromImage(out)

    # Smoothed first, the checkerboard averages to its mean rather than
    # point-sampling to whichever phase each output voxel happens to land on.
    interior = values[2:-2, 2:-2, 2:-2]
    assert interior.std() < 200, f"aliasing survived: std={interior.std():.0f}"
    assert 1500 < interior.mean() < 2500


def test_a_target_in_the_wrong_unit_fails_loudly():
    """A um target used as mm is a 1000x error that collapses the volume in silence."""
    from aind_ibl_ephys_alignment_preprocessing._async import histology

    img = _volume((28.8 * UM, 28.8 * UM, 32.0 * UM), size=(640, 462, 247))

    # 30_000 um == 30 mm: what the code did before the conversion was added.
    with pytest.raises(ValueError, match="different units"):
        histology.resample_to_isotropic(img, 30_000.0, "registration")


def test_real_volume_keeps_its_extent():
    """The regression itself: a full volume must not become a single voxel."""
    img = _volume((28.8 * UM, 28.8 * UM, 32.0 * UM), size=(640, 462, 247))

    out = resample_to_isotropic(img, 30.0, "registration")

    assert out.GetSize() == (614, 444, 263)
    for axis in range(3):
        before = img.GetSize()[axis] * img.GetSpacing()[axis]
        after = out.GetSize()[axis] * out.GetSpacing()[axis]
        assert after == pytest.approx(before, rel=0.02)


# --- base and pipeline images must stay one voxel grid ----------------------


def _pair():
    """A base/pipeline pair: same voxels, two headers, slightly different spacing.

    The pipeline spacing is the corrected level-0 spacing scaled by 2**level and
    need not equal the base's, which comes straight from the Zarr scales.
    """
    arr = np.zeros((60, 80, 100), dtype=np.uint16)
    base = sitk.GetImageFromArray(arr)
    base.SetSpacing((28.8 * UM, 28.8 * UM, 32.0 * UM))
    base.SetOrigin((11.82, -1.5, 1.5))
    pipeline = sitk.GetImageFromArray(arr)
    pipeline.SetSpacing((29.1 * UM, 29.1 * UM, 32.3 * UM))
    pipeline.SetOrigin((11.90, -1.6, 1.4))
    return base, pipeline


def test_pipeline_grid_stays_the_same_size_as_the_base():
    """The CCF warp lands on the pipeline grid and is stamped with the base header."""
    from aind_ibl_ephys_alignment_preprocessing._async.histology import (
        _mirror_pipeline_grid,
    )

    base, pipeline = _pair()
    assert base.GetSize() == pipeline.GetSize()

    resampled_base = resample_to_isotropic(base, 30.0, "base")
    resampled_pipeline = _mirror_pipeline_grid(pipeline, base, resampled_base)

    assert resampled_pipeline.GetSize() == resampled_base.GetSize()


def test_pipeline_header_arithmetic_is_self_consistent():
    """Index i after sits where index i*factor sat before, in the pipeline frame.

    Weak by construction: the assertion is algebraically what
    ``_mirror_pipeline_grid`` sets, so it catches a typo and nothing more. The
    real check is
    :func:`test_mirrored_pixels_equal_resampling_the_pipeline_image_directly`.
    """
    from aind_ibl_ephys_alignment_preprocessing._async.histology import (
        _mirror_pipeline_grid,
    )

    base, pipeline = _pair()
    resampled_base = resample_to_isotropic(base, 30.0, "base")
    out = _mirror_pipeline_grid(pipeline, base, resampled_base)

    factor = [a / b for b, a in zip(base.GetSpacing(), resampled_base.GetSpacing(), strict=True)]
    for index in [(0, 0, 0), (10, 20, 5), (37.5, 12.25, 8.0)]:
        equivalent = [i * f for i, f in zip(index, factor, strict=True)]
        assert out.TransformContinuousIndexToPhysicalPoint(index) == pytest.approx(
            pipeline.TransformContinuousIndexToPhysicalPoint(equivalent)
        )


def test_an_untouched_base_leaves_the_pipeline_image_alone():
    """Nothing to mirror when the base was already coarse enough."""
    from aind_ibl_ephys_alignment_preprocessing._async.histology import (
        _mirror_pipeline_grid,
    )

    base, pipeline = _pair()
    base.SetSpacing((50.0 * UM,) * 3)

    unchanged = resample_to_isotropic(base, 30.0, "base")
    assert unchanged is base
    assert _mirror_pipeline_grid(pipeline, base, unchanged) is pipeline


def test_mirrored_pixels_equal_resampling_the_pipeline_image_directly():
    """Copying the base's resampled voxels must equal resampling the pipeline image.

    The independent check: rather than asserting the header arithmetic against
    itself, resample the *pipeline* image in its own physical frame onto the
    mirrored grid and require the pixels to match. If the index correspondence
    were off by any amount, distinctive content would disagree.

    Two things have to be matched for the comparison to mean anything, and
    getting either wrong produces a large, convincing, wrong answer:

    - the **anti-alias filter**. Its sigma is physical, so the same filter on the
      same voxel array needs rescaling by the spacing ratio between the frames.
      Omitting it compares smoothed against unsmoothed.
    - the **output dtype**. The mirrored image is cast back to uint16; leaving
      the comparison in float shows a uniform sub-1 difference that is just
      rounding.

    That the filter must be rescaled is itself the argument for doing the
    resample once: the filter's job is to remove content above the *output*
    Nyquist, which the resample factor sets, and that factor is identical for
    both frames. Filtering per-frame in physical units would apply two slightly
    different filters to one array.
    """
    from aind_ibl_ephys_alignment_preprocessing._async.histology import (
        _mirror_pipeline_grid,
    )

    rng = np.random.default_rng(0)
    arr = rng.integers(0, 4000, (60, 80, 100)).astype(np.uint16)
    direction = (0.0, 0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0)

    base = sitk.GetImageFromArray(arr)
    base.SetSpacing((28.8 * UM, 28.8 * UM, 32.0 * UM))
    base.SetOrigin((11.82, -1.5, 1.5))
    base.SetDirection(direction)
    pipeline = sitk.GetImageFromArray(arr)
    pipeline.SetSpacing((29.1 * UM, 29.1 * UM, 32.3 * UM))
    pipeline.SetOrigin((11.90, -1.6, 1.4))
    pipeline.SetDirection(direction)

    mirrored = _mirror_pipeline_grid(pipeline, base, resample_to_isotropic(base, 30.0, "base"))

    target_mm = 30.0 * UM
    work = sitk.Cast(pipeline, sitk.sitkFloat32)
    for axis, (base_sp, pipe_sp) in enumerate(zip(base.GetSpacing(), pipeline.GetSpacing(), strict=True)):
        if target_mm > base_sp:  # the axes the base's filter touched
            work = sitk.RecursiveGaussian(work, sigma=((target_mm / 2.0) / base_sp) * pipe_sp, direction=axis)
    direct = sitk.Resample(
        work,
        mirrored.GetSize(),
        sitk.Transform(),
        sitk.sitkLinear,
        mirrored.GetOrigin(),
        mirrored.GetSpacing(),
        mirrored.GetDirection(),
        0.0,
        sitk.sitkFloat32,
    )

    got = sitk.GetArrayFromImage(mirrored).astype(np.int64)
    want = sitk.GetArrayFromImage(sitk.Cast(direct, sitk.sitkUInt16)).astype(np.int64)
    assert np.array_equal(got, want)
