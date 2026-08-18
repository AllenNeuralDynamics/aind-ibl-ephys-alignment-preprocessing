"""Header-built ANTs domains must warp identically to write-and-read-back ones.

The domains used to come from writing a volume and reading it with
``ants.image_read``. Building them from a header instead is only safe if the warp
cannot tell the difference, so that is what is asserted -- byte-identical output,
not "close enough".
"""

from __future__ import annotations

import numpy as np
import pytest

sitk = pytest.importorskip("SimpleITK")
ants = pytest.importorskip("ants")

from aind_anatomical_utils.anatomical_volume import AnatomicalHeader  # noqa: E402

from aind_ibl_ephys_alignment_preprocessing.histology import (  # noqa: E402
    ants_domain_stub,
    ants_warp_domain,
)

UM = 1.0 / 1000.0


def _volume(size=(30, 32, 34), spacing=(30 * UM,) * 3, fill=None):
    rng = np.random.default_rng(0)
    arr = rng.random(size[::-1]).astype("float32") if fill is None else np.full(size[::-1], fill, "float32")
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    img.SetOrigin((11.82, -1.5, 1.5))
    img.SetDirection((0.0, 0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0))
    return img


def _np_dtype(img):
    """The numpy dtype of a SimpleITK image, without copying it."""
    return sitk.GetArrayViewFromImage(img).dtype


def _write_read_back(img, tmp_path, name):
    """Today's path: reorient, write, read back with ANTs."""
    from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION

    path = tmp_path / name
    sitk.WriteImage(sitk.DICOMOrient(img, _BLESSED_DIRECTION), str(path), useCompression=True)
    return ants.image_read(str(path), pixeltype=None)


@pytest.mark.parametrize("interpolator", ["linear", "genericLabel"])
def test_warp_output_is_identical_to_the_read_back_domain(tmp_path, interpolator):
    """The whole justification for dropping the write and the read-back."""
    fixed_volume = _volume()
    moving = ants.from_numpy(
        np.random.default_rng(1).random((20, 22, 24)).astype("float32"),
        spacing=(25 * UM,) * 3,
        origin=(11.9, -1.6, 1.4),
    )

    old = _write_read_back(fixed_volume, tmp_path, "pipeline.nrrd")
    new = ants_warp_domain(AnatomicalHeader.from_sitk(fixed_volume), "pipeline", _np_dtype(fixed_volume))

    assert new.shape == old.shape
    assert new.pixeltype == old.pixeltype
    a = ants.apply_transforms(fixed=old, moving=moving, transformlist=[], interpolator=interpolator)
    b = ants.apply_transforms(fixed=new, moving=moving, transformlist=[], interpolator=interpolator)
    assert np.array_equal(a.numpy(), b.numpy())


def test_domain_geometry_matches_the_read_back_image(tmp_path):
    """Geometry, not just output: spacing/origin/direction must agree exactly."""
    volume = _volume()
    old = _write_read_back(volume, tmp_path, "pipeline.nrrd")
    new = ants_warp_domain(AnatomicalHeader.from_sitk(volume), "pipeline", _np_dtype(volume))

    assert new.spacing == pytest.approx(old.spacing)
    assert new.origin == pytest.approx(old.origin)
    assert np.allclose(new.direction, old.direction)


def test_stub_carries_the_same_geometry_at_one_voxel(tmp_path):
    """`correct_hist_domain_img` is read only for spacing/origin/direction."""
    volume = _volume()
    old = _write_read_back(volume, tmp_path, "raw.nrrd")
    stub = ants_domain_stub(AnatomicalHeader.from_sitk(volume), "raw")

    assert stub.shape == (1, 1, 1)
    assert stub.spacing == pytest.approx(old.spacing)
    assert stub.origin == pytest.approx(old.origin)
    assert np.allclose(stub.direction, old.direction)


def test_domain_is_uint8_zeros_not_a_copy_of_the_data():
    """The saving is that content is never materialised."""
    volume = _volume(fill=1234.0)
    domain = ants_warp_domain(AnatomicalHeader.from_sitk(volume), "pipeline", _np_dtype(volume))

    assert domain.numpy().max() == 0


def test_a_too_narrow_domain_wraps_the_moving_data():
    """Why np_eltype is required, and why the failure is worse than clipping.

    `apply_transforms` gives the output the fixed image's pixel type and offers
    no separate output-type argument, so a domain narrower than the data is
    destructive. It does not saturate -- it **wraps**: label 1340 in a uint8
    domain comes back as 60, which is not an obviously-broken value but a
    different, entirely plausible structure id. Nothing raises.
    """
    header = AnatomicalHeader.from_sitk(_volume())
    moving = ants.from_numpy(
        np.full((20, 22, 24), 1340.0, dtype="float32"),
        spacing=(25 * UM,) * 3,
        origin=(11.9, -1.6, 1.4),
    )

    narrow = ants.apply_transforms(fixed=ants_warp_domain(header, "x", np.uint8), moving=moving, transformlist=[])
    wide = ants.apply_transforms(fixed=ants_warp_domain(header, "x", np.uint16), moving=moving, transformlist=[])

    assert narrow.numpy().max() == 1340 % 256 == 60  # wrapped into a valid-looking id
    assert wide.numpy().max() == 1340
