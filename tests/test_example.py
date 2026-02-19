"""Smoke tests: verify every public module is importable."""

import aind_ibl_ephys_alignment_preprocessing


def test_version():
    """Test that version is defined."""
    assert aind_ibl_ephys_alignment_preprocessing.__version__ is not None
    assert isinstance(aind_ibl_ephys_alignment_preprocessing.__version__, str)


def test_types_importable():
    """Smoke test: types module imports."""
    from aind_ibl_ephys_alignment_preprocessing.types import (
        AssetInfo,
        ManifestRow,
        OutputDirs,
        PipelineConfig,
        ProcessResult,
        ReferencePaths,
        ReferenceVolumes,
    )

    assert PipelineConfig is not None
    assert ProcessResult is not None
    assert ManifestRow is not None
    assert OutputDirs is not None
    assert AssetInfo is not None
    assert ReferencePaths is not None
    assert ReferenceVolumes is not None


def test_constants_importable():
    """Smoke test: constants module imports."""
    from aind_ibl_ephys_alignment_preprocessing._constants import _BLESSED_DIRECTION

    assert _BLESSED_DIRECTION == "IRP"


def test_validation_importable():
    """Smoke test: validation module imports."""
    from aind_ibl_ephys_alignment_preprocessing.validation import (
        PipelineValidator,
        ValidationResult,
    )

    assert PipelineValidator is not None
    assert ValidationResult is not None


def test_discovery_importable():
    """Smoke test: discovery module imports."""
    from aind_ibl_ephys_alignment_preprocessing.discovery import (
        determine_desired_level,
        find_asset_info,
        prepare_result_dirs,
    )

    assert callable(find_asset_info)
    assert callable(prepare_result_dirs)
    assert callable(determine_desired_level)


def test_histology_importable():
    """Smoke test: histology module imports."""
    from aind_ibl_ephys_alignment_preprocessing.histology import (
        apply_ccf_inverse_tx_then_fix_domain,
        compress_reorient_nrrd_file,
        convert_img_direction_and_write,
        copy_registration_channel_ccf_reorient,
        process_additional_channels_pipeline,
        transform_ccf_labels_to_image_space,
        transform_ccf_to_image_space,
        write_registration_channel_images,
    )

    assert callable(convert_img_direction_and_write)
    assert callable(copy_registration_channel_ccf_reorient)
    assert callable(write_registration_channel_images)
    assert callable(process_additional_channels_pipeline)
    assert callable(apply_ccf_inverse_tx_then_fix_domain)
    assert callable(compress_reorient_nrrd_file)
    assert callable(transform_ccf_to_image_space)
    assert callable(transform_ccf_labels_to_image_space)


def test_probes_importable():
    """Smoke test: probes module imports."""
    from aind_ibl_ephys_alignment_preprocessing.probes import process_manifest_row

    assert callable(process_manifest_row)


def test_ephys_importable():
    """Smoke test: ephys module imports."""
    from aind_ibl_ephys_alignment_preprocessing.ephys import run_ephys_for_recording

    assert callable(run_ephys_for_recording)


def test_pipeline_importable():
    """Smoke test: pipeline module imports."""
    from aind_ibl_ephys_alignment_preprocessing.pipeline import run_pipeline

    assert callable(run_pipeline)


def test_async_concurrency_importable():
    """Smoke test: async concurrency module imports."""
    from aind_ibl_ephys_alignment_preprocessing._async.concurrency import (
        EphysCoordinator,
        IOLimits,
        Limits,
        io_to_thread_on,
        to_thread_logged,
    )

    assert IOLimits is not None
    assert Limits is not None
    assert EphysCoordinator is not None
    assert callable(to_thread_logged)
    assert callable(io_to_thread_on)


def test_async_pipeline_importable():
    """Smoke test: async pipeline module imports."""
    from aind_ibl_ephys_alignment_preprocessing._async.pipeline import run_pipeline_async

    assert callable(run_pipeline_async)


def test_scripts_importable():
    """Smoke test: scripts module imports."""
    from aind_ibl_ephys_alignment_preprocessing.scripts.run import main

    assert callable(main)


def test_pipeline_config_basic():
    """Verify PipelineConfig resolves relative paths."""
    from pathlib import Path

    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

    config = PipelineConfig(
        data_root=Path("/test/data"),
        results_root=Path("/test/results"),
        neuroglancer_file=Path("ng_file.json"),
        manifest_csv=Path("manifest.csv"),
    )
    assert config.neuroglancer_file == Path("/test/data/ng_file.json")
    assert config.manifest_csv == Path("/test/data/manifest.csv")
    assert config.template_25 == Path("/test/data/smartspim_lca_template/smartspim_lca_template_25.nii.gz")
    assert config.scratch_root is None
    assert config.skip_ephys is False


def test_reference_paths_for_data_root():
    """Verify ReferencePaths.for_data_root resolves correctly."""
    from pathlib import Path

    from aind_ibl_ephys_alignment_preprocessing.types import ReferencePaths

    rp = ReferencePaths.for_data_root(Path("/mydata"))
    assert rp.template_25 == Path("/mydata/smartspim_lca_template/smartspim_lca_template_25.nii.gz")
    assert rp.ibl_atlas_histology_path == Path("/mydata/iblatlas_allenatlas")
