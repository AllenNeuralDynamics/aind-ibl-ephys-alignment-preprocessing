# AIND IBL Ephys Alignment Preprocessing

![CI](https://github.com/AllenNeuralDynamics/aind-ibl-ephys-alignment-preprocessing/actions/workflows/ci-call.yml/badge.svg)
[![PyPI](https://img.shields.io/pypi/v/aind-ibl-ephys-alignment-preprocessing)](https://pypi.org/project/aind-ibl-ephys-alignment-preprocessing/)
[![License](https://img.shields.io/badge/license-MIT-brightgreen)](LICENSE)

This package prepares histology images, annotated probe tracks, and sorted
electrophysiology recordings for the
[IBL ephys alignment GUI](https://github.com/AllenNeuralDynamics/ibl-ephys-alignment-gui).
It combines the images, tracks, and recording data into one package for a
mouse. It uses existing image registrations, which describe how the mouse's
histology images align with a reference brain image.

## Choose how to run

On Code Ocean, use the shared [preprocessing trigger](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-trigger-capsule).
Its instructions cover first-time credentials setup, launching, and finding
the output. The trigger finds the input assets and launches the
[worker capsule](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-pipeline-capsule)
for each processing stage. You normally only need to prepare the annotations
asset described below.

To use this as a standalone Python package, follow [Run preprocessing locally](docs/running_locally.md)
for installation, command-line use, and the Python API. The input requirements
and manifest format are shared by both routes.

## Prepare the inputs

You need the mouse's processed SmartSPIM images and registration results,
its raw and sorted ephys recordings, a saved Neuroglancer state, and a
manifest CSV. The manifest connects the annotated tracks to the recordings
they belong to. You can prepare it in a spreadsheet using the examples below.

### Neuroglancer layers

The Neuroglancer JSON must identify the histology images and contain the
point annotations marking the probe tracks. Image layers supply the source
locations of the processed SmartSPIM OME-Zarr images. Named annotation layers
supply the points for each track. Keep the image sources and annotation
coordinates from the same imaging data.

Each `histology_track_id` in the CSV must match an annotation layer's name
exactly. Each `probe_file` must identify the JSON containing that layer.
Several rows can refer to one JSON, and a manifest can refer to several JSONs.
The CSV contains the names and mappings; the points themselves stay in the
Neuroglancer files.

## Manifest CSV

The manifest connects each track drawn in Neuroglancer to its ephys
recording. Make it in a spreadsheet and export it as CSV, with one row per
track. A four-shank probe with four annotated tracks needs four rows.

The examples below use `example_sorted` for the sorted asset and `tracks.json`
for the annotation file. Replace these with your own asset and filename, and
replace the track IDs with their exact names in Neuroglancer.

### Start with these five columns

| Column | Enter this value |
|---|---|
| `mouseid` | The mouse ID, the same on every row. |
| `sorted_recording` | The sorted recording folder name; on Code Ocean, use the sorted asset's exact name. |
| `probe_file` | The annotation JSON filename, such as `tracks.json`. The older form `tracks` also works. |
| `histology_track_id` | The exact track/layer name inside that JSON, such as `A_shank1`. |
| `ephys_collection` | The recorded probe name, such as `ProbeA`. For a split quadbase, use its common probe name as explained below. |

For a single-shank recording, those five columns are enough:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection
858398,example_sorted,tracks.json,A_track,ProbeA
```

### Four-shank NP2.0 recorded together

If all four shanks are recorded together as `ProbeA`, repeat `ProbeA` on all
four rows. Add `histology_shank` to identify the physical shank for each track.

CSV shank numbers start at **0**: the four shanks are numbered `0`, `1`, `2`,
and `3`. In this example, the track named `A_shank1` belongs to physical shank
`0`. Enter the number that matches your track; preprocessing does not read
shank numbers from track names.

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,histology_shank
858398,example_sorted,tracks.json,A_shank1,ProbeA,0
858398,example_sorted,tracks.json,A_shank2,ProbeA,1
858398,example_sorted,tracks.json,A_shank3,ProbeA,2
858398,example_sorted,tracks.json,A_shank4,ProbeA,3
```

Preprocessing uses the physical shank number to select the corresponding
data within `ProbeA`. You can omit `ephys_shank` and `logical_probe`.

### Quadbase shanks recorded separately

A stream is a group of recording channels stored together. A split quadbase
records its four shanks as separate streams, for example
`ProbeB-1`, `ProbeB-2`, `ProbeB-3`, and `ProbeB-4`. Enter the common probe name
`ProbeB` and the physical `histology_shank` for each track:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,histology_shank
858398,example_sorted,tracks.json,B_shank1,ProbeB,0
858398,example_sorted,tracks.json,B_shank2,ProbeB,1
858398,example_sorted,tracks.json,B_shank3,ProbeB,2
858398,example_sorted,tracks.json,B_shank4,ProbeB,3
```

Preprocessing checks which physical shank each stream recorded and fills in
the matching stream name: `ProbeB` with physical shank `2`, for example, becomes
`ProbeB-3`. That stream contains one shank, whose index within the stream is
`0`. Preprocessing fills in this distinction; leave `ephys_shank` out of this
form of the CSV. If the saved recording information cannot identify a unique
stream, the run reports which mapping it could not resolve.

You can combine single-shank, four-shank, and split quadbase rows in one
manifest. Use the six-column header from either four-shank example and leave
`histology_shank` blank for a single-shank row.

### Save the input files

Save the CSV with `manifest` in its filename, for example
`858398_manifest.csv`. Keep every annotation JSON named by `probe_file`
available with the other inputs. The [input folder guide](docs/input_data.md)
shows the layout for a direct run. Capsule-specific instructions are linked
under [Choose how to run](#choose-how-to-run).

Use the correct `sorted_recording` on each row when a mouse has several
recording sessions. Dates in annotation filenames or track names are labels;
the `sorted_recording` column selects the recording to process.

### Optional columns

`ephys_shank` selects a shank within an explicitly named recording stream.
It is useful when you need to specify that mapping yourself. `logical_probe`
groups separate streams belonging to one physical probe. `annotation_format`
defaults to `json`, including when its cells are blank. These are optional,
supported fields; the examples above do not need them. See the
[column reference](docs/manifest_reference.md) for their defaults and for
separate surface-finding recordings.

### Legacy columns

These columns are still supported for existing manifests, but deprecated
for authoring new ones:

| Legacy column | Use for new manifests |
|---|---|
| `probe_id` | `histology_track_id` |
| `probe_name` | `ephys_collection` |
| `probe_shank` | `histology_shank`; use `ephys_shank` too when an explicit stream mapping needs a different local shank number. |
| `registration_asset` | A registration setting for the whole run, outside the track CSV. |
| `template_to_ccf_asset` | A template-to-atlas registration setting for the whole run, outside the track CSV. |

The two registration names remain current run settings; only their placement
in the CSV is legacy. The [column reference](docs/manifest_reference.md)
explains legacy defaults and how explicit and inferred shank mappings interact.

### Building a manifest from a mapping

If you prefer to generate the CSV, the optional
[manifest builder](docs/building_a_manifest.md) takes a configuration file
mapping tracks to recordings and streams.

## Registration overrides

Most runs use the registrations already supplied with the imaging data. If
the brain was registered separately, supply its transform directory through
the run's `registration_asset` setting. If it targets a reference image other
than the standard SmartSPIM template, also set `template_to_ccf_asset` to
that template's registration to the Allen brain atlas.

The [registration guide](docs/registration_overrides.md) describes the files
needed for these cases. Registration is a choice for the whole brain, so new
manifests can leave it to the run settings. Existing CSV registration columns
remain supported.

## Run and check the output

Run preprocessing using your chosen execution method. The output is a mouse
directory containing the images, converted recording data, and track points
needed by the alignment GUI.

Check the run's warnings as well as its completion status. A run may finish
successfully when some recorded shanks lack sorting outputs. For example,
histology and raw-channel measurements can be available for all four shanks
while spikes are available for only two. An output asset cannot supply spikes
that were absent from its sorted input.

The output keeps a record of those gaps and the recordings selected for
conversion. The [inventory and conversion-plan guide](docs/ephys_inventory.md)
explains where to find that record if you need to investigate a warning.

### Output structure

The mouse's output directory contains `datapackage.json`, which tells the
alignment GUI where to find the histology images, track points, and recording
data. Keep the complete directory together when copying it to another machine.

```text
<mouseid>/
├── datapackage.json
├── image_space_histology/
├── track_data/
└── <recording_id>/
    └── <ephys_collection>/
        ├── xyz_picks*.json
        ├── channels.*.npy
        ├── spikes.*.npy
        └── band_corr/
```

This is a shortened view. The recording folders also contain measurements
used for the GUI's quality-control plots. Copying only the files named inside
`datapackage.json` can lose those measurements. See the
[output contract](docs/external_asset_references.md) for file and portability
details.

### Regenerating only `datapackage.json`

To repair the package description without repeating the conversion, follow
the [regeneration instructions](docs/running_locally.md#regenerating-only-datapackagejson).
This uses existing outputs; it does not recover missing sorting data.

## Development

For development commands, architecture, and validation, see
[CONTRIBUTING.md](CONTRIBUTING.md). This project uses the [MIT license](LICENSE).
