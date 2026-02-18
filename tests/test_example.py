"""Example test module."""

import aind_ibl_ephys_alignment_preprocessing


def test_version():
    """Test that version is defined."""
    assert aind_ibl_ephys_alignment_preprocessing.__version__ is not None
    assert isinstance(aind_ibl_ephys_alignment_preprocessing.__version__, str)
