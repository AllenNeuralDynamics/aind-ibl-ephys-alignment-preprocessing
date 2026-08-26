"""The sync and async histology paths must stay the same pipeline.

They are two spellings of one algorithm, and they drift silently: the async side
gained resampling, dtype narrowing and the geometry sidecar while the sync side
kept none of them, which left sync emitting a datapackage that referenced a
sidecar it never wrote.
"""

from __future__ import annotations

import inspect

from aind_ibl_ephys_alignment_preprocessing import histology
from aind_ibl_ephys_alignment_preprocessing._async import histology as ahistology

#: Steps whose absence on either side is a behaviour difference, not a style one.
SHARED_STEPS = [
    "resample_to_isotropic",
    "_mirror_pipeline_grid",
    "_narrow_dtype",
    "image_geometry",
]


def test_the_shared_helpers_are_one_implementation():
    """Not copies: the async module imports the sync module's functions."""
    for name in SHARED_STEPS:
        assert getattr(ahistology, name) is getattr(histology, name), name


def _calls(fn) -> str:
    return inspect.getsource(fn)


def test_both_registration_writers_resample_and_mirror():
    """The base/pipeline pairing has to hold on whichever path ran."""
    for writer in (
        histology.write_registration_channel_images,
        ahistology.write_registration_channel_images_async,
    ):
        source = _calls(writer)
        assert "resample_to_isotropic" in source, writer.__name__
        assert "_mirror_pipeline_grid" in source, writer.__name__
        assert "output_voxel_size_um" in source, writer.__name__


def test_both_registration_writers_emit_the_geometry_sidecar():
    """`_build_histology` always references it, so both paths must write it."""
    for writer in (
        histology.write_registration_channel_images,
        ahistology.write_registration_channel_images_async,
    ):
        source = _calls(writer)
        assert ".json" in source, writer.__name__


def test_both_channel_paths_resample():
    for fn in (
        histology.process_additional_channels_pipeline,
        ahistology.process_additional_channel_pipeline_async,
    ):
        assert "resample_to_isotropic" in _calls(fn), fn.__name__


def test_both_ccf_paths_narrow_their_dtypes():
    for fn in (
        histology.transform_ccf_to_image_space,
        ahistology.transform_ccf_to_image_space_async,
        histology.transform_ccf_labels_to_image_space,
        ahistology.transform_ccf_labels_to_image_space_async,
    ):
        assert "_narrow_dtype" in _calls(fn), fn.__name__


def test_both_ccf_paths_branch_on_the_registration_frame():
    """The re-grid is a compensation for the pipeline's anchoring, not a step.

    Applying it to an off-pipeline transform is the defect that put 776259's
    channels outside the atlas, so neither spelling may hard-code it.
    """
    for fn in (
        histology.apply_ccf_inverse_tx_then_fix_domain,
        ahistology.apply_ccf_inverse_tx_then_fix_domain_async,
    ):
        source = _calls(fn)
        assert "frame.regrid_to_pipeline" in source, fn.__name__
        # The header repair is meaningless once the resample ran on the real grid.
        assert source.count("set_origin") == 1, fn.__name__
        # ants.apply_transforms takes the output size from `fixed`, so a header
        # stub yields a single-voxel volume with no error. Both paths refuse one.
        assert "header stub" in source, fn.__name__


def test_both_ccf_paths_honor_a_registration_override():
    """``pipeline_registration_chains`` names the registration being replaced."""
    for fn in (
        histology.apply_ccf_inverse_tx_then_fix_domain,
        ahistology.apply_ccf_inverse_tx_then_fix_domain_async,
    ):
        source = _calls(fn)
        assert "point_chain()" in source, fn.__name__
        assert "pipeline_registration_chains" not in source, fn.__name__


def test_both_probe_paths_branch_and_honor_the_override():
    from aind_ibl_ephys_alignment_preprocessing import probes
    from aind_ibl_ephys_alignment_preprocessing._async import probes as aprobes

    for fn in (probes._write_qc_probe_outputs, aprobes._write_qc_probe_outputs_async):
        source = _calls(fn)
        assert "frame.regrid_to_pipeline" in source, fn.__name__
        assert "point_chain()" in source, fn.__name__
        assert "pipeline_registration_chains" not in source, fn.__name__
