# Run preprocessing locally

This guide is for running the Python package directly. For the usual Code
Ocean workflow, use the [trigger instructions](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-trigger-capsule).
Both routes use the same [manifest](../README.md#manifest-csv).

## Install the package

From PyPI:

```bash
uv venv
uv pip install aind-ibl-ephys-alignment-preprocessing
```

From source:

```bash
git clone https://github.com/AllenNeuralDynamics/aind-ibl-ephys-alignment-preprocessing.git
cd aind-ibl-ephys-alignment-preprocessing
uv sync
```

## Process a mouse

Prepare the [input folders](input_data.md) and [manifest](../README.md#manifest-csv)
first. From the checkout prepared with `uv sync`, run:

```bash
uv run aind-ibl-preprocess \
    --data-root /path/to/data \
    --results-root /path/to/results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv
```

### CLI options

| Flag | Description |
|------|-------------|
| `--data-root PATH` | **(required)** Root directory containing all input data |
| `--results-root PATH` | **(required)** Root directory for output |
| `--neuroglancer PATH` | **(required)** Neuroglancer JSON file (absolute, or relative to `data_root`) |
| `--manifest PATH` | **(required)** Manifest CSV file (absolute, or relative to `data_root`) |
| `--registration-asset PATH` | Optional separate brain registration directory, relative to `data_root` or absolute. |
| `--template-to-ccf-asset PATH` | Optional registration from an alternative template to the Allen CCF. |
| `--scratch-root PATH` | Scratch directory for temporary files (default: system temp) |
| `--skip-ephys` | Skip electrophysiology extraction |
| `--validate-only` | Run pre-flight validation checks and exit |
| `--datapackage-only` | Regenerate `datapackage.json` from existing outputs without rerunning histology or ephys preprocessing |
| `--source-results PATH` | Existing results asset/root or mouse output directory to copy into `results_root` before `--datapackage-only` regeneration |
| `--async` | Run the pipeline asynchronously with concurrency |

### Validation

Before running the full pipeline, you can check that all required inputs are
present and correctly structured:

```bash
uv run aind-ibl-preprocess \
    --data-root /path/to/data \
    --results-root /path/to/results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv \
    --validate-only
```

This checks file existence, manifest structure, reference data availability,
disk space, and available RAM. Passing these checks does not establish that
sorting outputs exist for every shank or that every image registration is correct.

### Regenerating only `datapackage.json`

If preprocessing has already produced the histology/ephys outputs, you can
rewrite just the datapackage metadata:

```bash
uv run aind-ibl-preprocess \
    --data-root /path/to/data \
    --results-root /path/to/results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv \
    --datapackage-only
```

This path rediscovers asset metadata, infers successful manifest rows from
existing `xyz_picks*.json` outputs, validates referenced files, and writes a
fresh `datapackage.json`. It does not rerun ANTs transforms, zarr reads, probe
track conversion, or ephys extraction.

On Code Ocean, previous `/results` assets are immutable when attached back to a
new capsule. Attach the old preprocessed results asset under `/data`, then copy
it into the new run's `/results` before rewriting the datapackage:

```bash
uv run aind-ibl-preprocess \
    --data-root /data \
    --results-root /results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv \
    --datapackage-only \
    --source-results /data/<old_preprocessed_results_asset>
```

`--source-results` may point either at the prior results asset root, containing
`<mouseid>/datapackage.json`, or directly at the `<mouseid>` output directory.

## Python API

The package exposes a small public API:

```python
from aind_ibl_ephys_alignment_preprocessing import (
    PipelineConfig,
    run_pipeline,
    run_pipeline_async,
)

config = PipelineConfig(
    data_root="/path/to/data",
    results_root="/path/to/results",
    neuroglancer_file="neuroglancer.json",
    manifest_csv="manifest.csv",
    skip_ephys=False,            # set True to skip ephys extraction
    desired_voxel_size_um=25.0,  # multiscale zarr level selection
    num_parallel_jobs=4,         # parallel workers for ephys
)

# Synchronous
results = run_pipeline(config)

# Asynchronous (structured concurrency with TaskGroups)
import asyncio
results = asyncio.run(run_pipeline_async(config))
```

Call either the synchronous or asynchronous function for a run. Each returned
`ProcessResult` identifies a probe and recording, its output files, and an
optional skip reason.
