"""Tests for the CCF to IBL bregma-relative conversion."""

import numpy as np

from aind_ibl_ephys_alignment_preprocessing.bregma import BREGMA_CCF_MLAPDV_UM, ccf_to_bregma_um


def test_bregma_maps_to_origin() -> None:
    """Bregma itself lands at the origin."""
    np.testing.assert_allclose(ccf_to_bregma_um(BREGMA_CCF_MLAPDV_UM), 0.0)


def test_axes_flip_ap_and_dv_only() -> None:
    """Moving right, posterior and ventral in CCF is right, posterior and ventral in IBL's signs."""
    offsets = np.eye(3) * 100.0
    np.testing.assert_allclose(ccf_to_bregma_um(BREGMA_CCF_MLAPDV_UM + offsets), np.diag([100.0, -100.0, -100.0]))


def test_matches_iblatlas_ccf2xyz_values() -> None:
    """Pinned against ``AllenAtlas(25).ccf2xyz(..., ccf_order='mlapdv') * 1e6`` from iblatlas 1.2.1."""
    ccf = np.array([[0.0, 0.0, 0.0], [11400.0, 13200.0, 8000.0], [5000.0, 7000.0, 3500.0]])
    expected = np.array([[-5739.0, 5400.0, 332.0], [5661.0, -7800.0, -7668.0], [-739.0, -1600.0, -3168.0]])
    np.testing.assert_allclose(ccf_to_bregma_um(ccf), expected)
