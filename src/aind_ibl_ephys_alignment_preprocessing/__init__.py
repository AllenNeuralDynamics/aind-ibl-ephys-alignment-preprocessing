"""Library to prepare histology and ephys for the IBL ephys alignment GUI."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("aind-ibl-ephys-alignment-preprocessing")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

# Public API re-exports (lazy — no heavy imports at package level)
__all__ = [
    "__version__",
    "PipelineConfig",
    "ProcessResult",
    "run_pipeline",
    "run_pipeline_async",
]


def __getattr__(name: str) -> object:
    """Lazy attribute access to avoid importing heavy dependencies at import time."""
    if name == "PipelineConfig":
        from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

        return PipelineConfig
    if name == "ProcessResult":
        from aind_ibl_ephys_alignment_preprocessing.types import ProcessResult

        return ProcessResult
    if name == "run_pipeline":
        from aind_ibl_ephys_alignment_preprocessing.pipeline import run_pipeline

        return run_pipeline
    if name == "run_pipeline_async":
        from aind_ibl_ephys_alignment_preprocessing._async.pipeline import run_pipeline_async

        return run_pipeline_async
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
