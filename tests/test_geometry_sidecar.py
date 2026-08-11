"""The geometry sidecar must describe the volume on disk, exactly.

Its whole purpose is to replace opening a volume, so a consumer using it must
land on the same physical points as one that opened the file.
"""

from __future__ import annotations

import numpy as np
import pytest

sitk = pytest.importorskip("SimpleITK")

from aind_ibl_ephys_alignment_preprocessing._async.histology import (  # noqa: E402
    image_geometry,
)

# Deliberately not axis-aligned-trivial: an identity direction would pass even
# if rows and columns were transposed.
IRP_DIRECTION = (0.0, 0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0)


def _volume():
    img = sitk.GetImageFromArray(np.zeros((7, 9, 11), dtype=np.uint16))
    img.SetSpacing((30.0, 30.0, 30.0))
    img.SetOrigin((11.82, -1.5, 1.5))
    img.SetDirection(IRP_DIRECTION)
    return img


def _rebuild(payload) -> sitk.Image:
    """Reconstruct a geometry-only image the way a consumer would."""
    h = payload["header"]
    img = sitk.Image(*h["size_ijk"], sitk.sitkUInt8)
    img.SetSpacing(tuple(h["spacing"]))
    img.SetOrigin(tuple(h["origin"]))
    img.SetDirection(tuple(h["direction"]))
    return img


def test_size_is_index_order_not_array_shape():
    """The reversal is the silent failure this sidecar could cause."""
    img = _volume()
    payload = image_geometry(img)

    assert tuple(payload["header"]["size_ijk"]) == img.GetSize() == (11, 9, 7)
    assert tuple(payload["header"]["size_ijk"]) != sitk.GetArrayFromImage(img).shape


def test_direction_columns_are_the_index_axes_in_lps():
    """Row-major flattening with columns as index-axis vectors -- not the transpose."""
    payload = image_geometry(_volume())
    direction = np.array(payload["header"]["direction"]).reshape(3, 3)

    # i, j, k unit vectors expressed in LPS.
    assert np.allclose(direction[:, 0], [0, 0, -1])
    assert np.allclose(direction[:, 1], [0, -1, 0])
    assert np.allclose(direction[:, 2], [1, 0, 0])


@pytest.mark.parametrize(
    "index",
    [(0, 0, 0), (10, 8, 6), (5.5, 2.25, 3.75), (10.0, 0.0, 0.0)],
)
def test_reconstructed_geometry_maps_points_identically(index):
    """The point of the whole exercise: same physical answer, no volume opened."""
    img = _volume()
    rebuilt = _rebuild(image_geometry(img))

    assert rebuilt.TransformContinuousIndexToPhysicalPoint(index) == pytest.approx(
        img.TransformContinuousIndexToPhysicalPoint(index)
    )


def test_sidecar_describes_the_file_as_written_not_as_held(tmp_path):
    """Reorientation changes origin and direction; the sidecar must follow it."""
    from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION

    img = _volume()
    oriented = sitk.DICOMOrient(img, _BLESSED_DIRECTION)
    path = tmp_path / "histology_registration_pipeline.nrrd"
    sitk.WriteImage(oriented, str(path), useCompression=True)

    payload = image_geometry(oriented)
    on_disk = sitk.ReadImage(str(path))

    assert tuple(payload["header"]["origin"]) == pytest.approx(on_disk.GetOrigin())
    assert tuple(payload["header"]["direction"]) == pytest.approx(on_disk.GetDirection())
    assert tuple(payload["header"]["size_ijk"]) == on_disk.GetSize()
    # And capturing before the conversion would have been wrong:
    assert tuple(payload["header"]["origin"]) != pytest.approx(img.GetOrigin())


def test_payload_declares_its_conventions():
    """Numbers alone cannot say LPS or micrometres; ITK is unit-agnostic."""
    payload = image_geometry(_volume())

    assert payload["space"] == "left-posterior-superior"
    assert payload["units"] == "micrometer"
    assert payload["schema"].startswith("anatomical-header/")
