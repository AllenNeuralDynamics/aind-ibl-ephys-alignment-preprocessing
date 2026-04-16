# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install dependencies
uv sync

# Run all checks (formatting, linting, type checking, tests)
./scripts/run_linters_and_checks.sh -c

# Run tests
uv run pytest

# Run a single test file
uv run pytest tests/test_example.py

# Run a single test by name
uv run pytest -k "test_name"

# Formatting
uv run ruff format

# Linting
uv run ruff check
uv run ruff check --fix

# Type checking
uv run mypy

# Spell checking
uv run codespell --check-filenames

# Build docs
uv run sphinx-build -b html docs/source/ docs/build/html

```

Always use `uv run` to execute commands, `uv add` to add dependencies, and `uv sync` to set up the environment. Never use bare `pip` or `python`.

## Architecture

This is a Python package using a `src/` layout. Source code lives in `src/aind_ibl_ephys_alignment_preprocessing/`, tests in `tests/`.

- Build system: hatchling
- Formatting/linting: ruff (line length 120, numpy docstring convention)
- Testing: pytest with coverage reporting
- Type checking: mypy (strict mode)
- Versioning: commitizen (semantic versioning via conventional commits)
- Documentation: Sphinx with furo theme (`docs/source/`)
- CLI entry points: defined in `[project.scripts]` in pyproject.toml, implemented in `src/aind_ibl_ephys_alignment_preprocessing/scripts/`

### Overview

Preprocessing pipeline that prepares SmartSPIM histology and Neuropixels ephys
data for the IBL ephys alignment GUI. Registers histology volumes to the Allen
CCF, converts Neuroglancer probe-track annotations into atlas coordinates, and
optionally extracts spike-sorted ephys into IBL ALF format.

### Module layout

```
src/aind_ibl_ephys_alignment_preprocessing/
  __init__.py          # Lazy-loading public API (__version__, PipelineConfig, ProcessResult, run_pipeline, run_pipeline_async)
  types.py             # All Pydantic models and frozen dataclasses (PipelineConfig, ManifestRow, OutputDirs, etc.)
  _constants.py        # _BLESSED_DIRECTION = "IRP" (DICOM orientation for IBL GUI)
  discovery.py         # Asset discovery from Neuroglancer JSON, zarr level selection, output dir creation
  validation.py        # Pre-flight validation (PipelineValidator)
  histology.py         # Sync volume processing: CCF reorientation, channel transforms, inverse transforms
  probes.py            # Per-probe annotation loading and coordinate conversion
  ephys.py             # Sync ephys extraction wrapper (delegates to aind-ephys-ibl-gui-conversion)
  pipeline.py          # Sync orchestrator: run_pipeline()
  manifest.py          # Output datapackage.json generation (DataPackage model, schema v1.0.0)
  scripts/run.py       # CLI entry point: aind-ibl-preprocess
  _async/
    concurrency.py     # IO lane throttling (per-path semaphores), EphysCoordinator, thread helpers
    histology.py       # Async volume processing wrappers
    probes.py          # Async probe processing
    ephys.py           # Async ephys coordination (single-flight dedup via ProcessPoolExecutor)
    pipeline.py        # Async orchestrator: run_pipeline_async() with TaskGroups
```

### Key patterns

- **Lazy loading**: `__init__.py` uses `__getattr__` to defer heavy imports
  (ants, pandas, etc.) until accessed.
- **Path resolution**: `PipelineConfig` is a frozen Pydantic v2 model with a
  `@model_validator(mode="after")` that resolves all relative paths against
  `data_root`.
- **Frozen data structures**: All dataclasses and Pydantic models are frozen.
  Use `dataclasses.replace()` or `.model_copy(update=...)` to create modified
  copies.
- **Async concurrency**: IO lanes throttle per-path (scratch/results/data).
  Registration is serialized (`max_registration=1`), ephys limited to 2
  concurrent. `EphysCoordinator` deduplicates runs per recording key.
- **Transform chain**: Neuroglancer pixels -> SPIM (LPS) -> Template -> CCF
  (LPS/um) -> Bregma (IBL). Each stage writes intermediate coordinate files.

### Public API

```python
from aind_ibl_ephys_alignment_preprocessing import (
    PipelineConfig,
    ProcessResult,
    run_pipeline,          # sync
    run_pipeline_async,    # async with TaskGroups
)
```

### CLI

```bash
aind-ibl-preprocess \
    --data-root /path/to/data \
    --results-root /path/to/results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv \
    [--scratch-root /tmp] [--skip-ephys] [--validate-only] [--async]
```

## Code style

- **Line length**: 120
- **Target**: Python 3.11+
- **Docstrings**: NumPy convention
- **Type checking**: mypy strict. Histology modules (`histology.py`,
  `_async/histology.py`) have `disallow_untyped_calls = false` due to
  SimpleITK/ANTs lacking stubs.
- **Ruff rules**: `Q, RUF100, C90, I, F, E, W, D, UP, ANN, PYI`.
  Tests skip `ANN` and `D`.
- **McCabe complexity**: max 14

## Testing

Tests are currently smoke tests (import verification, path resolution logic).
No mocking, no integration tests with real data. Coverage threshold is 0%.

## Dependencies with missing type stubs

These are `ignore_missing_imports = true` in mypy: `ants`, `SimpleITK`,
`aind_anatomical_utils`, `aind_ephys_ibl_gui_conversion`,
`aind_registration_utils`, `aind_s3_cache`, `aind_zarr_utils`, `iblatlas`,
`filelock`.

## Coordinate conventions

- **IBL**: x=ML (right+), y=AP (anterior+), z=DV (dorsal+), origin=bregma, meters
- **LPS**: Left-Posterior-Superior (ANTs/SimpleITK convention)
- **`_BLESSED_DIRECTION`**: `"IRP"` -- the DICOM orientation code the IBL GUI expects
