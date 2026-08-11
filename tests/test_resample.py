"""The output grid must be the same for every mouse, without inventing signal."""

from __future__ import annotations

import numpy as np
import pytest

sitk = pytest.importorskip("SimpleITK")

from aind_ibl_ephys_alignment_preprocessing._async.histology import (  # noqa: E402
    resample_to_isotropic,
)


def _volume(spacing, size=(64, 64, 64), value=1000):
    img = sitk.GetImageFromArray(np.full(size[::-1], value, dtype=np.uint16))
    img.SetSpacing(spacing)
    img.SetOrigin((3.0, -5.0, 7.0))
    return img


@pytest.mark.parametrize(
    "spacing",
    [
        (14.4, 14.4, 16.0),  # pre-September stitching
        (28.8, 28.8, 32.0),  # September onward
        (32.0, 32.0, 16.0),  # anisotropic the other way
    ],
)
def test_every_vintage_lands_on_one_grid(spacing):
    """The point of the resample: stitching vintage must stop setting the grid."""
    out = resample_to_isotropic(_volume(spacing), 30.0, "test")

    assert out.GetSpacing() == pytest.approx((30.0, 30.0, 30.0))


def test_physical_extent_and_origin_are_preserved():
    """Resampling changes sampling, not where the data sits in the world."""
    img = _volume((14.4, 14.4, 16.0))
    out = resample_to_isotropic(img, 30.0, "test")

    assert out.GetOrigin() == img.GetOrigin()
    assert out.GetDirection() == img.GetDirection()
    for axis in range(3):
        before = img.GetSize()[axis] * img.GetSpacing()[axis]
        after = out.GetSize()[axis] * out.GetSpacing()[axis]
        assert after == pytest.approx(before, rel=0.02)


def test_already_coarser_volumes_are_left_alone():
    """The target must never manufacture resolution."""
    img = _volume((50.0, 50.0, 50.0))

    out = resample_to_isotropic(img, 30.0, "test")

    assert out is img


def test_dtype_survives_the_float_round_trip():
    """Smoothing needs floats; the channels must not become float on disk."""
    img = _volume((14.4, 14.4, 16.0))

    out = resample_to_isotropic(img, 30.0, "test")

    assert out.GetPixelID() == img.GetPixelID() == sitk.sitkUInt16


def test_downsampling_is_anti_aliased():
    """sitk.Resample point-samples; without a low-pass, fine detail aliases."""
    # A single-voxel-period checkerboard: pure content above the new Nyquist.
    size = 96
    arr = np.indices((size, size, size)).sum(axis=0) % 2
    img = sitk.GetImageFromArray((arr * 4000).astype(np.uint16))
    img.SetSpacing((10.0, 10.0, 10.0))

    out = resample_to_isotropic(img, 30.0, "test")
    values = sitk.GetArrayFromImage(out)

    # Smoothed first, the checkerboard averages to its mean rather than
    # point-sampling to whichever phase each output voxel happens to land on.
    interior = values[2:-2, 2:-2, 2:-2]
    assert interior.std() < 200, f"aliasing survived: std={interior.std():.0f}"
    assert 1500 < interior.mean() < 2500
