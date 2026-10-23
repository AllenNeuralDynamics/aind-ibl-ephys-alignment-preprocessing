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
| `template_to_ccf_asset` | *null* | Path relative to `data_root` to the template-to-CCF registration of the template `registration_asset` was registered to, when that is not the SmartSPIM LCA template. One value per brain, like `registration_asset`. See [Separate template-to-CCF registration](#separate-template-to-ccf-registration). |

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

#### Building a manifest from a mapping

Most of a manifest repeats or follows from two choices: which sorted recording each session uses,
and which stream each annotated track was recorded on. `aind-ibl-build-manifest` takes only those,
in a TOML file, and derives the rest. Track names are never parsed, because annotators follow no
naming convention.

```toml
mouseid = "837772"
annotation = "837772-probe-annotations-CK-2026-09-29.json"  # relative to this file
# registration_asset = "..."                                  # optional, per brain
# template_to_ccf_asset = "..."                               # optional, per brain

[[session]]
sorted_recording = "ecephys_837772_2026-06-27_14-08-21_sorted_2026-08-27_18-59-49"
# surface_finding = "..."                                     # optional

[session.tracks]
"732_MD_1" = "46101-1"   # one shank of a split quad-base probe
"733_BLA" = "46810"      # a single-stream probe
"740_PL" = { stream = "45883", shank = 1 }  # one shank of a single-stream multi-shank probe; the first shank is 1
```

```bash
aind-ibl-build-manifest 837772.toml -o manifest_837772.csv [--data-root /data]
```

A stream is the collection parsed from the sort's `postprocessed/` analyzer names. They are listed
from a copy mounted under `--data-root`, or else unsigned from that sort's prefix in
`aind-open-data`. How the columns are derived:

- `probe_file` is the annotation file's name without its extension.
- **Split quad-base probe:** a stream `<probe>-<N>`, when the sort also holds another
  `<probe>-<M>`, gets `logical_probe = <probe>`, `histology_shank = N - 1` and `ephys_shank = 0`.
- **Any other stream** is its own `logical_probe`. Its shanks stay blank, or both take the given
  shank's index.

The build writes nothing if any of these fail, and it reports every failure:
- every stream in each sort has a track, and every assigned stream is in its sort;
- every track is a layer of the annotation file with at least two points;
- no track is assigned twice;
- a split stream has exactly one track, and the tracks sharing any other stream have distinct
  shanks.

Annotated layers left out of the mapping are printed as notes, so a missed track gets noticed.

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

Set the manifest's optional `registration_asset` column to the transform
directory relative to `data_root`:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,registration_asset
750108,recording_sorted,probeA,A0001,ProbeA,SmartSPIM_750108_reg/ccf_Ex_639_Em_667
```

Use one path per brain; conflicting nonblank values fail. An absent or blank
column uses the stitched asset's registration.

The directory must contain `ls_to_template_SyN_0GenericAffine.mat`,
`ls_to_template_SyN_1InverseWarp.nii.gz` and a geometry sidecar (below): ANTs
SyN outputs registering this brain's selected channel to the SmartSPIM LCA
template, or to another template whose CCF registration
`template_to_ccf_asset` (below) supplies.

The final directory name selects the channel after stripping a leading
`ccf_`. For this example, the stitched asset must contain
`image_tile_fusing/OMEZarr/Ex_639_Em_667.zarr/`. Images and metadata come
from the stitched asset; the override needs only transforms and their sidecar.

Every separate registration needs a geometry sidecar beside the transforms;
a directory without one fails. Preprocessing checks
`ls_to_template_transform_information.json`, then `transform_information.json`.
Its `moving_domain` must describe the actual registration input in LPS
millimeters and match the selected volume's native or pipeline-anchored frame.
Its `fixed_domain` must match the template the template-to-CCF chain starts
from: the SmartSPIM LCA template unless `template_to_ccf_asset` is set. Missing
or incompatible domains fail, since a wrong frame or template still puts the
probes somewhere plausible.

**QC:** Additional-channel CCF volumes (`emit_qc=True`) also need
`ls_to_template_SyN_1Warp.nii.gz` in the override directory; without it the
QC warp fails rather than falling back to the stitched asset's registration.

#### Generate the sidecar from ANTs outputs

Use the actual fixed template and moving histology images passed to ANTs,
including any preprocessing of their geometry. Transform files alone cannot
recover these domains. Headers must be 3D, axis-aligned, and in millimeters.

Run `uv sync` in the preprocessing library checkout, then save this as
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

Pass the transform directory, fixed image, and moving image in that order:

```bash
uv run python generate_registration_sidecar.py \
    /path/to/SmartSPIM_750108_reg/ccf_Ex_639_Em_667 \
    /path/to/fixed_template_used_by_ants.nii.gz \
    /path/to/moving_histology_used_by_ants.nii.gz
```

The script reads only image headers. Transform filenames are relative to the
sidecar; rename raw ANTs outputs to the expected names. The schema requires a
forward-warp filename, but normal preprocessing needs only the affine and
inverse warp. The helper does not verify transform files.

Generate the JSON beside the transforms before uploading the asset;
Code Ocean's `/data` mounts are read-only.

### Separate template-to-CCF registration

A brain registered to a template other than the SmartSPIM LCA template (an
iDISCO template, say) needs that template's CCF registration. Set the optional
`template_to_ccf_asset` column to its directory relative to `data_root`, one
path per brain like `registration_asset`:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,registration_asset,template_to_ccf_asset
750107,recording_sorted,probeA,A0001,ProbeA,idisco_reg/750107,AIND-MouseiDISCO-v1.0-to-CCFv3-r1_2026-10-22_17-34-50
```

The directory must contain `template_to_ccf_0GenericAffine.mat` and
`template_to_ccf_1InverseWarp.nii.gz`, ANTs SyN outputs with CCF fixed and the
template moving; `emit_qc=True` also needs `template_to_ccf_1Warp.nii.gz`. They replace the
pipeline's `spim_template_to_ccf_syn_*_25` transforms in the probe and QC
chains, and the datapackage refers to them under the `template_to_ccf` asset,
named after the directory's first component under `data_root`.

The column requires a `registration_asset` with a sidecar, whose `fixed_domain`
is checked against this template. If the directory ships its own sidecar
(`template_to_ccf_transform_information.json` or `transform_information.json`),
the registration's `fixed_domain` must match its `moving_domain`. Without one,
the registration's `fixed_domain` must *not* match the SmartSPIM LCA template;
to replace the LCA template's own CCF registration, ship the sidecar.

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

The `datapackage.json` file describes the mouse, recordings, ephys collections,
transform chain, histology, channel tables, and track mappings. It references
each ephys directory; individual spike/QC files and `band_corr/` contents are
loaded by convention rather than enumerated in the JSON. See
[the datapackage contract](docs/external_asset_references.md) for asset
descriptions versus runtime location resolution, Code Ocean/desktop examples,
shank joins, and validation limits. The model can be loaded back with:

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
6. **Manifest generation** -- Write `datapackage.json` describing explicit
   references and ephys collection directories; preserve their conventional
   contents, including QC files not individually listed in the JSON.

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

See [CONTRIBUTING.md](CONTRIBUTING.md) for developer setup, checks, and
repository conventions. Coding agents should also read [AGENTS.md](AGENTS.md).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
