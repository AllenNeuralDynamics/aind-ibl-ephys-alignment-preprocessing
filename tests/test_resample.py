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
