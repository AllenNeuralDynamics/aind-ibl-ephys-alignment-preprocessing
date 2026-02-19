"""Standalone CLI entry point for ``aind-ibl-preprocess``."""

from __future__ import annotations

import argparse
import asyncio
import sys
from pathlib import Path


def main() -> None:
    """Parse args, build :class:`PipelineConfig`, and run the pipeline."""
    parser = argparse.ArgumentParser(
        prog="aind-ibl-preprocess",
        description="SmartSPIM histology + ephys preprocessing for the IBL alignment GUI.",
    )
    parser.add_argument("--data-root", type=Path, required=True, help="Root of input data.")
    parser.add_argument("--results-root", type=Path, required=True, help="Root of output results.")
    parser.add_argument("--scratch-root", type=Path, default=None, help="Scratch directory (default: tempdir).")
    parser.add_argument("--neuroglancer", required=True, help="Neuroglancer JSON file path.")
    parser.add_argument("--manifest", required=True, help="Manifest CSV file path.")
    parser.add_argument("--skip-ephys", action="store_true", help="Skip ephys extraction.")
    parser.add_argument("--validate-only", action="store_true", help="Run validation only, then exit.")
    parser.add_argument("--async", dest="run_async", action="store_true", help="Run pipeline asynchronously.")
    args = parser.parse_args()

    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

    config = PipelineConfig(
        data_root=args.data_root,
        results_root=args.results_root,
        scratch_root=args.scratch_root,
        neuroglancer_file=Path(args.neuroglancer),
        manifest_csv=Path(args.manifest),
        skip_ephys=args.skip_ephys,
    )

    if args.validate_only:
        from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator

        validator = PipelineValidator(config)
        results = validator.validate_all()
        validator.print_summary(results)
        sys.exit(0 if not validator.has_errors(results) else 1)

    if args.run_async:
        from aind_ibl_ephys_alignment_preprocessing._async.pipeline import run_pipeline_async

        asyncio.run(run_pipeline_async(config))
    else:
        from aind_ibl_ephys_alignment_preprocessing.pipeline import run_pipeline

        run_pipeline(config)


if __name__ == "__main__":
    main()
