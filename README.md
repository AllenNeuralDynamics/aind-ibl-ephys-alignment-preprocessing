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
| `--datapackage-only` | Regenerate `datapackage.json` from existing outputs without rerunning histology or ephys preprocessing |
| `--source-results PATH` | Existing results asset/root or mouse output directory to copy into `results_root` before `--datapackage-only` regeneration |
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

### Regenerating only `datapackage.json`

If preprocessing has already produced the histology/ephys outputs, you can
rewrite just the datapackage metadata:

```bash
aind-ibl-preprocess \
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
aind-ibl-preprocess \
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

A CSV file describing how histology tracks map onto ephys collections. Each
row represents one histology track (or one histology shank) mapped to one
ephys collection and ephys shank.

The contract is declared once in `types.MANIFEST_COLUMNS` and drives both the
parser and pre-flight validation, so what a run accepts cannot drift from what
validation demands. **Optional columns are never required**: a manifest written
before a column existed stays valid.

**Required columns:**

| Column | Description |
|--------|-------------|
| `mouseid` | Mouse identifier. All rows must reference the same mouse. |
| `sorted_recording` | Name of the spike-sorted recording folder under `data_root`. The recording ID is derived by stripping a `_sorted` suffix if present. |
| `probe_file` | Basename (without extension) of the Neuroglancer annotation file. Resolved via glob `*/<probe_file>.<annotation_format>` under `data_root`. |
| `histology_track_id` | Neuroglancer layer / histology track identifier. Legacy alias: `probe_id`. |
| `ephys_collection` | Ephys ALF output folder produced by `aind-ephys-ibl-gui-conversion` (for example `ProbeA`). Legacy alias: `probe_name`. |

**Optional columns:**

| Column | Default | Description |
|--------|---------|-------------|
| `annotation_format` | `json` | File extension for the annotation file (lowercase). |
| `logical_probe` | `ephys_collection` | Physical/logical probe identity. Split quadbase streams can share one `logical_probe` while using different `ephys_collection` values. |
| `histology_shank` | `probe_shank` | 0-based shank index in the physical/histology probe. |
| `ephys_shank` | `probe_shank`, else `histology_shank` | 0-based shank index within `ephys_collection`. For split quadbase streams this is usually `0` even when `histology_shank` is 0..3. |
| `probe_id` | *legacy* | Alias for `histology_track_id`. |
| `probe_name` | *legacy* | Alias for `ephys_collection`. |
| `probe_shank` | *legacy* | Alias for both `histology_shank` and `ephys_shank`. |
| `surface_finding` | *null* | Path (relative to `data_root`) to a surface-finding file. |
| `registration_asset` | *null* | Path relative to `data_root` to a replacement registration directory, e.g. `SmartSPIM_750108_reg/ccf_Ex_639_Em_667`. One value per brain; all nonblank rows must agree. See [Separate registration asset](#separate-registration-asset) for its contents and coordinate requirements. |

**Constraints:**
- All rows must have the same `mouseid`.
- The tuple `(mouseid, histology_track_id, histology_shank)` must be unique.
- The tuple `(recording_id, ephys_collection, histology_shank)` must be unique
  for GUI filename generation.
- The tuple `(recording_id, ephys_collection, ephys_shank)` must be unique for
  ephys shank mapping.
- For single-stream multi-shank probes, multiple rows share one
  `ephys_collection` and differ in `ephys_shank`.
- For split quadbase probes, multiple rows can share one `logical_probe` but
  use different `ephys_collection` values; each collection can have
  `ephys_shank=0`.

**Example:**

```csv
mouseid,sorted_recording,probe_file,histology_track_id,logical_probe,ephys_collection,histology_shank,ephys_shank
mouse001,2024-06-01_rec_sorted,track_annotations_probeA,A0001,probeA,ProbeA,,
mouse001,2024-06-01_rec_sorted,track_annotations_probeB_shank0,B0001,probeB,ProbeB,0,0
mouse001,2024-06-01_rec_sorted,track_annotations_probeB_shank1,B0001,probeB,ProbeB,1,1
mouse001,2024-06-01_rec_sorted,track_annotations_quad_shank3,Q0001,quad0,ProbeD,3,0
```

### Neuroglancer JSON

A Neuroglancer state JSON file that contains image source URIs pointing to the
SmartSPIM OME-Zarr volumes. The pipeline extracts the first image source URI
to locate the SmartSPIM asset directory and discover the registration channel
and any additional channels.

The asset directory must contain:
- `image_tile_fusing/OMEZarr/<channel>.zarr/` -- fused OME-Zarr volumes
- `image_atlas_alignment/<registration_channel_stem>/` -- ANTs registration
  outputs: the light-sheet to template affine and inverse warp (or supply a
  separate registration asset as described below)

### Separate registration asset

The optional manifest column `registration_asset` selects replacement
light-sheet-to-SmartSPIM-template transforms, for example after re-registering a
brain on a corrected channel. It points directly to the directory containing
the transforms, relative to `PipelineConfig.data_root`, rather than just to the
asset root. Omit the column or leave it blank to use the stitched asset's
registration. Repeat the same path across rows for the brain; blank rows are
ignored, but conflicting nonblank paths raise an error.

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,registration_asset
750108,recording_sorted,probeA,A0001,ProbeA,SmartSPIM_750108_reg/ccf_Ex_639_Em_667
```

Mount or place the asset so this layout exists under `data_root`:

```text
SmartSPIM_750108_reg/
└── ccf_Ex_639_Em_667/
    ├── ls_to_template_SyN_0GenericAffine.mat
    ├── ls_to_template_SyN_1InverseWarp.nii.gz
    └── ls_to_template_transform_information.json  # geometry sidecar; see below
```

The affine and inverse warp are required, with exactly these filenames. Supply
valid ANTs SyN outputs for this brain's selected channel registered to the
SmartSPIM LCA template used by preprocessing. Despite the `ccf_` directory
prefix, these files represent the light-sheet-to-template registration; the
separate template-to-CCF transforms are still required. A forward warp
(`ls_to_template_SyN_1Warp.nii.gz`) and a moved image are not required for the
normal preprocessing path.

The final directory name selects the registration channel: preprocessing strips
one leading `ccf_`, so `ccf_Ex_639_Em_667` and `Ex_639_Em_667` both select
`Ex_639_Em_667`. That channel must exist at
`<smartspim_asset>/image_tile_fusing/OMEZarr/Ex_639_Em_667.zarr/`. The stitched
asset still supplies the image volumes and metadata; the registration asset
only supplies the replacement transforms and their geometry. Its own
`processing.json` or copy of the OME-Zarr volumes is not required. Missing
directories, required transforms, or the selected channel raise an error.

**Coordinate frame:** For a standalone registration computed in the image's
native frame, include the transform producer's geometry sidecar beside the
transforms. Preprocessing looks for
`ls_to_template_transform_information.json` first, then
`transform_information.json`, and parses it with `aind-ants-transform-sidecar`.
It must declare a `moving_domain` describing the actual registration input's
physical geometry in LPS millimeters. Preprocessing compares its voxel-center
bounds with the selected volume's native and pipeline-anchored frames, and
uses the nearer plausible frame. A different registration resolution is
allowed; a missing `moving_domain` or bounds matching neither frame raises an
error. Keep the sidecar paired with the transforms that produced it.

Without a recognized sidecar, preprocessing assumes the standard SmartSPIM
pipeline's anchored frame, even for a separate asset. Therefore, a native-frame
standalone registration needs the sidecar to avoid incorrect coordinates;
placing it in a separate asset alone does not select the native frame.

**QC limitation:** Optional additional-channel CCF volumes (`emit_qc=True`) still
use the stitched asset's pipeline image-transform chain. Those QC volumes do
not use the replacement registration. Probe coordinates and CCF-to-image
resampling use the selected registration.

#### Generate the sidecar from ANTs outputs

Use the original **fixed image** (SmartSPIM template) and **moving image**
(histology channel) passed to ANTs, including any resampling or origin changes
made before registration. The affine and warp files alone do not record enough
information to reconstruct both input domains. Do not substitute the warped
output image for the moving input or assign an origin just to make validation
pass.

In an environment containing the preprocessing dependencies (`uv sync` in the
preprocessing library checkout), save the following as
`generate_registration_sidecar.py`:

```python
import sys
from pathlib import Path

import SimpleITK as sitk
from aind_ants_transform_sidecar import dump_package, load_package
from aind_registration_utils.domains import ImageHeader, make_syn_sidecar


def read_header(path):
    reader = sitk.ImageFileReader()
    reader.SetFileName(str(path))
    reader.ReadImageInformation()  # Read geometry without loading voxel data.
    if reader.GetDimension() != 3:
        raise ValueError("Registration inputs must be 3D")
    return ImageHeader(
        origin=reader.GetOrigin(),
        spacing=reader.GetSpacing(),
        direction=reader.GetDirection(),
        size=reader.GetSize(),
    )


registration_dir, fixed_image, moving_image = map(Path, sys.argv[1:])
sidecar = make_syn_sidecar(
    fixed_image_header=read_header(fixed_image),
    moving_image_header=read_header(moving_image),
    affine="ls_to_template_SyN_0GenericAffine.mat",
    warp="ls_to_template_SyN_1Warp.nii.gz",
    inverse_warp="ls_to_template_SyN_1InverseWarp.nii.gz",
)
text = dump_package(sidecar)
load_package(text)  # Validate the serialized package.
destination = registration_dir / "ls_to_template_transform_information.json"
destination.write_text(text + "\n")
print(destination)
```

Run it with the transform directory, fixed input, and moving input, in that order:

```bash
uv run python generate_registration_sidecar.py \
    /path/to/SmartSPIM_750108_reg/ccf_Ex_639_Em_667 \
    /path/to/fixed_template_used_by_ants.nii.gz \
    /path/to/moving_histology_used_by_ants.nii.gz
```

Image headers must express origin and spacing in **millimeters**, with
axis-aligned directions (axis permutations and flips are supported; oblique
directions are rejected). SimpleITK supplies LPS geometry; the helper computes
voxel-center bounds and spatial signatures and supplies both `fixed_domain`
and `moving_domain`, as required by the sidecar schema.

The transform filenames in the example are relative to the sidecar directory;
rename or copy differently named raw ANTs outputs to these expected names
before publishing the asset. The schema requires a forward-warp filename even
though normal preprocessing only needs the affine and inverse warp. The
helper does not check that the named transform files exist or inspect their
contents. Keep the actual forward warp when available; naming it in the sidecar
does not make it a required input for normal preprocessing.

Generate this JSON in a writable staging directory before uploading the asset,
then include it beside the affine and inverse warp. Code Ocean mounts under
`/data` are read-only. Reloading the JSON checks its schema; preprocessing
separately checks whether its moving domain matches the selected histology
volume.

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
    |-- ccf_space_histology/                  # QC only (emit_qc); GUI never reads it
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
        +-- <ephys_collection>/
            |-- xyz_picks.json                # Probe track picks (CCF coordinates)
            |-- xyz_picks_image_space.json    # Probe track picks (image space)
            |-- xyz_picks_shank<N>.json       # Per-shank picks (multi-shank only)
            |-- channels.localCoordinates.npy  # Ephys channel table
            |-- channels.rawInd.npy
            |-- channels.contactId.npy          # Optional; new ephys conversion outputs
            |-- channels.shankInd.npy
            +-- band_corr/                    # Full matrices + row_channels.json
```

The `datapackage.json` file contains a structured manifest of all outputs,
including transform chain paths, histology volume paths, and per-probe
metadata. It can be loaded back with:

```python
from aind_ibl_ephys_alignment_preprocessing.datapackage import load_datapackage

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
