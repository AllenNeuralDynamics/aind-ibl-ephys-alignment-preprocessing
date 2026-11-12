"""CCF coordinates to the IBL bregma-relative frame, in physical units with no atlas volume."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

#: Bregma in CCF (ML, AP, DV) micrometres from the volume's front-left-top corner: IBL's
#: ``ALLEN_CCF_LANDMARKS_MLAPDV_UM['bregma']``.
BREGMA_CCF_MLAPDV_UM = np.array([5739.0, 5400.0, 332.0])

#: IBL's ML, AP, DV point right, anterior and dorsal; CCF's ML, AP, DV point right, posterior and ventral.
_CCF_TO_IBL_SIGN = np.array([1.0, -1.0, -1.0])


def ccf_to_bregma_um(ccf_mlapdv_um: NDArray[Any]) -> NDArray[np.float64]:
    """Map CCF (ML, AP, DV) micrometres to IBL (ML, AP, DV) micrometres relative to bregma.

    Matches ``iblatlas`` ``AllenAtlas.ccf2xyz`` times 1e6 at any resolution: that call converts to
    voxel indices and back, so the resolution cancels and only the bregma offset and the AP and DV
    flips remain.

    Parameters
    ----------
    ccf_mlapdv_um : NDArray
        ``(..., 3)`` CCF coordinates in micrometres, ML right, AP posterior, DV ventral.

    Returns
    -------
    NDArray[np.float64]
        ``(..., 3)`` bregma-relative coordinates in micrometres, ML right, AP anterior, DV dorsal.
    """
    bregma_um: NDArray[np.float64] = (
        np.asarray(ccf_mlapdv_um, dtype=np.float64) - BREGMA_CCF_MLAPDV_UM
    ) * _CCF_TO_IBL_SIGN
    return bregma_um
