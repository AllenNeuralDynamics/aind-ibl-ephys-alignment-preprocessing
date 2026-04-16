# AIND IBL Ephys Alignment Preprocessing

![CI](https://github.com/AllenNeuralDynamics/aind-ibl-ephys-alignment-preprocessing/actions/workflows/ci-call.yml/badge.svg)
[![PyPI - Version](https://img.shields.io/pypi/v/aind-ibl-ephys-alignment-preprocessing)](https://pypi.org/project/aind-ibl-ephys-alignment-preprocessing/)
[![semantic-release: angular](https://img.shields.io/badge/semantic--release-angular-e10079?logo=semantic-release)](https://github.com/semantic-release/semantic-release)
[![License](https://img.shields.io/badge/license-MIT-brightgreen)](LICENSE)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Copier](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-grayscale-inverted-border.json)](https://github.com/copier-org/copier)

Preprocessing pipeline that prepares SmartSPIM histology and Neuropixels
electrophysiology data for the
[IBL ephys alignment GUI](https://github.com/AllenNeuralDynamics/ibl-ephys-alignment-gui).
It registers histology volumes to the Allen Common Coordinate Framework (CCF),
converts Neuroglancer probe-track annotations into atlas coordinates, and
optionally extracts spike-sorted ephys data into
[IBL ALF format](https://int-brain-lab.github.io/ONE/alf_intro.html).

## Installation

From PyPI:

```bash
pip install aind-ibl-ephys-alignment-preprocessing
```

From source:

```bash
git clone https://github.com/AllenNeuralDynamics/aind-ibl-ephys-alignment-preprocessing.git
cd aind-ibl-ephys-alignment-preprocessing
pip install .
```

## Quick start

The package provides a CLI entry point, `aind-ibl-preprocess`:

```bash
aind-ibl-preprocess \
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
| `--scratch-root PATH` | Scratch directory for temporary files (default: system temp) |
| `--skip-ephys` | Skip electrophysiology extraction |
| `--validate-only` | Run pre-flight validation checks and exit |
| `--async` | Run the pipeline asynchronously with concurrency |

### Validation

Before running the full pipeline, you can check that all required inputs are
present and correctly structured:

```bash
aind-ibl-preprocess \
    --data-root /path/to/data \
    --results-root /path/to/results \
    --neuroglancer neuroglancer.json \
    --manifest manifest.csv \
    --validate-only
```

This checks file existence, manifest structure, reference data availability,
disk space (warns below 50 GB, errors below 10 GB), and available RAM
(minimum 8 GB).

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

Each element of `results` is a `ProcessResult` containing the probe ID,
recording ID, list of files written, and an optional skip reason.

## Input data requirements

### Data root layout

All input paths can be absolute or relative to `data_root`. The expected
directory structure is:

```
data_root/
|
|-- neuroglancer.json
|-- manifest.csv
|
|-- smartspim_lca_template/
|   +-- smartspim_lca_template_25.nii.gz
|
|-- allen_mouse_ccf/
|   +-- average_template/
|       +-- average_template_25.nii.gz
|
|-- allen_mouse_ccf_annotations_lateralized_compact/
|   |-- ccf_2017_annotation_25_lateralized_compact.nrrd
|   +-- ccf_2017_annotation_25_lateralized_unique_vals.npz
|
|-- iblatlas_allenatlas/
|   +-- ...
|
|-- spim_template_to_ccf/
|   |-- syn_0GenericAffine.mat
|   +-- syn_1InverseWarp.nii.gz
|
|-- <smartspim_asset>/
|   |-- image_tile_fusing/
|   |   +-- OMEZarr/
|   |       |-- <registration_channel>.zarr/
|   |       +-- <additional_channels>.zarr/
|   +-- image_atlas_alignment/
|       +-- <registration_channel_stem>/
|           |-- moved_ls_to_ccf.nii.gz
|           |-- ls_to_template_SyN_0GenericAffine.mat
|           +-- ls_to_template_SyN_1InverseWarp.nii.gz
|
|-- <sorted_recording>/
|   +-- ...  (spike sorting output)
|
+-- **/<probe_file>.json
```

Each of these is described in detail below.

### Manifest CSV

A CSV file describing which probes to process. Each row represents one probe
(or one shank of a multi-shank probe).

**Required columns:**

| Column | Description |
|--------|-------------|
| `mouseid` | Mouse identifier. All rows must reference the same mouse. |
| `sorted_recording` | Name of the spike-sorted recording folder under `data_root`. The recording ID is derived by stripping a `_sorted` suffix if present. |
| `probe_file` | Basename (without extension) of the Neuroglancer annotation file. Resolved via glob `*/<probe_file>.<annotation_format>` under `data_root`. |
| `probe_id` | Unique probe identifier. |
| `probe_name` | Subfolder name used for GUI output artifacts. |

**Optional columns:**

| Column | Default | Description |
|--------|---------|-------------|
| `annotation_format` | `json` | File extension for the annotation file (lowercase). |
| `probe_shank` | *null* | 0-based shank index for multi-shank probes. Leave empty for single-shank. |
| `surface_finding` | *null* | Path (relative to `data_root`) to a surface-finding file. |

**Constraints:**
- All rows must have the same `mouseid`.
- The tuple `(recording_id, probe_name, probe_shank)` must be unique across
  rows.
- For multi-shank probes, multiple rows can share the same `probe_name` but
  must differ in `probe_shank`.

**Example:**

```csv
mouseid,sorted_recording,probe_file,probe_id,probe_name,probe_shank
mouse001,2024-06-01_rec_sorted,track_annotations_probeA,A0001,probeA,
mouse001,2024-06-01_rec_sorted,track_annotations_probeB_shank0,B0001_s0,probeB,0
mouse001,2024-06-01_rec_sorted,track_annotations_probeB_shank1,B0001_s1,probeB,1
```

### Neuroglancer JSON

A Neuroglancer state JSON file that contains image source URIs pointing to the
SmartSPIM OME-Zarr volumes. The pipeline extracts the first image source URI
to locate the SmartSPIM asset directory and discover the registration channel
and any additional channels.

The asset directory must contain:
- `image_tile_fusing/OMEZarr/<channel>.zarr/` -- fused OME-Zarr volumes
- `image_atlas_alignment/<registration_channel_stem>/` -- ANTs registration
  outputs including `moved_ls_to_ccf.nii.gz` (the precomputed light-sheet to
  CCF registration)

### Reference volumes

These files are required for atlas registration. Default paths are relative to
`data_root` and can be overridden in `PipelineConfig`.

| File | Default path | Format | Description |
|------|--------------|--------|-------------|
| SmartSPIM LCA template | `smartspim_lca_template/smartspim_lca_template_25.nii.gz` | NIfTI (.nii.gz) | 25 um SmartSPIM template volume |
| CCF average template | `allen_mouse_ccf/average_template/average_template_25.nii.gz` | NIfTI (.nii.gz) | Allen CCF average template at 25 um |
| CCF lateralized labels | `allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_compact.nrrd` | NRRD | Lateralized CCF annotation labels at 25 um |
| CCF label lookup | `allen_mouse_ccf_annotations_lateralized_compact/ccf_2017_annotation_25_lateralized_unique_vals.npz` | NumPy (.npz) | Unique label values for the annotation volume |
| IBL atlas histology | `iblatlas_allenatlas/` | Directory | IBL atlas histology data used by `iblatlas.AllenAtlas` |
| Template-to-CCF transforms | `spim_template_to_ccf/` | Directory | ANTs SyN registration outputs (`syn_0GenericAffine.mat`, `syn_1InverseWarp.nii.gz`) |

### Probe annotation files

Per-probe Neuroglancer point annotation files. For each row in the manifest,
the pipeline searches `data_root` for a file matching the glob pattern
`*/<probe_file>.<annotation_format>` (e.g., `*/track_annotations_probeA.json`).

These files contain 3D point annotations marking the probe track in the
SmartSPIM image space. They are converted through the transform chain:

```
Neuroglancer pixels --> SPIM (LPS) --> Template --> CCF (LPS/um) --> Bregma (IBL)
```

### Electrophysiology data

Unless `--skip-ephys` is passed, the pipeline expects spike-sorted data at
`data_root/<sorted_recording>/` for each unique recording in the manifest. The
ephys extraction (delegated to
[aind-ephys-ibl-gui-conversion](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion))
converts this into IBL ALF format (e.g., `spikes.times.npy`,
`spikes.clusters.npy`, `channels.localCoordinates.npy`, etc.).

Ephys extraction runs once per unique recording ID (deduplicated across
manifest rows).

## Output structure

All outputs are written under `results_root/<mouseid>/`:

```
results_root/
|-- manifest.csv                              # Copy of input manifest
+-- <mouseid>/
    |-- ccf_space_histology/
    |   |-- histology_registration.nrrd       # Registration channel in CCF space
    |   +-- histology_<channel>.nrrd          # Additional channels in CCF space
    |
    |-- image_space_histology/
    |   |-- histology_registration.nrrd       # Registration channel in image space
    |   |-- histology_registration_pipeline.nrrd
    |   |-- ccf_in_mouse.nrrd                 # CCF template warped to image space
    |   +-- labels_in_mouse.nrrd              # CCF labels warped to image space
    |
    |-- track_data/
    |   |-- spim/<probe_id>.*                 # Track coordinates in SPIM space
    |   |-- template/<probe_id>.*             # Track coordinates in template space
    |   |-- ccf/<probe_id>.*                  # Track coordinates in CCF space
    |   |-- bregma_xyz/<probe_id>.*           # Track coordinates in IBL bregma space
    |   +-- datapackage.json                  # Machine-readable output manifest
    |
    +-- <recording_id>/
        +-- <probe_name>/
            |-- xyz_picks.json                # Probe track picks (CCF coordinates)
            |-- xyz_picks_image_space.json    # Probe track picks (image space)
            |-- xyz_picks_shank<N>.json       # Per-shank picks (multi-shank only)
            +-- spikes/                       # Ephys ALF output (if not skipped)
```

The `datapackage.json` file contains a structured manifest of all outputs,
including transform chain paths, histology volume paths, and per-probe
metadata. It can be loaded back with:

```python
from aind_ibl_ephys_alignment_preprocessing.manifest import load_datapackage

dp = load_datapackage("/path/to/datapackage.json")
```

## Pipeline overview

The pipeline performs the following steps:

1. **Asset discovery** -- Parse the Neuroglancer JSON to locate the SmartSPIM
   OME-Zarr volumes and ANTs registration outputs.
2. **Validation** -- Check that all required inputs exist and are correctly
   structured.
3. **Histology processing** -- Reorient the registration channel to CCF space,
   export additional channels, and compute inverse transforms (CCF template and
   labels warped back to image space).
4. **Probe processing** -- For each manifest row, load the Neuroglancer point
   annotations and convert them through the transform chain (SPIM -> template
   -> CCF -> bregma), writing coordinate files at each stage and producing
   `xyz_picks.json` files for the alignment GUI.
5. **Ephys extraction** (optional) -- Convert spike-sorted data into IBL ALF
   format.
6. **Manifest generation** -- Write `datapackage.json` summarizing all outputs.

## Adapting to your own data

If you are not using the AIND SmartSPIM pipeline but want to use this package
with your own histology and ephys data, the key requirements are:

1. **OME-Zarr volumes** -- Your fused histology images must be in OME-Zarr
   format, discoverable via a Neuroglancer JSON state file.
2. **ANTs registration outputs** -- You need ANTs SyN registration transforms
   mapping your histology images to a common template and from that template to
   the Allen CCF. The pipeline expects standard ANTs output files
   (`*_0GenericAffine.mat`, `*_1InverseWarp.nii.gz`).
3. **Reference volumes** -- The Allen CCF template, lateralized annotation
   labels, and the SmartSPIM LCA template at 25 um resolution. These can be
   overridden in `PipelineConfig` if your paths differ from the defaults.
4. **Probe annotations** -- Neuroglancer point annotation JSON files marking
   probe tracks in histology image space.
5. **Spike-sorted data** (optional) -- If you want ephys extraction, provide
   spike-sorted output compatible with
   [aind-ephys-ibl-gui-conversion](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion).

The coordinate convention used throughout is the IBL convention: x = ML
(right positive), y = AP (anterior positive), z = DV (dorsal positive), with
origin at bregma, in meters.

## Development

To develop the code, run:
```bash
uv sync
```

Please test your changes using the full linting and testing suite:

```bash
./scripts/run_linters_and_checks.sh -c
```

Or run individual commands:
```bash
uv run --frozen ruff format          # Code formatting
uv run --frozen ruff check           # Linting
uv run --frozen mypy                 # Type checking
uv run --frozen interrogate -v src   # Documentation coverage
uv run --frozen codespell --check-filenames  # Spell checking
uv run --frozen pytest --cov aind_ibl_ephys_alignment_preprocessing  # Tests with coverage
```

### Documentation
```bash
sphinx-build -b html docs/source/ docs/build/html
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
