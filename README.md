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

You need the mouse's processed SmartSPIM imaging asset with its registration
results, the sorted ephys recording assets, and the Neuroglancer JSON files
containing your probe-track annotations. Keep the image layers in the saved
JSON: the trigger uses them to find the imaging asset. The annotations must
refer to the imaging data used for that mouse. Preprocessing also uses the
raw ephys data; the trigger finds it from the sorted recording's input
information.

Your annotations asset contains the annotation JSONs and a manifest CSV. The
CSV connects the tracks you drew to the recordings they belong to. You can
prepare it in a spreadsheet using the examples below.

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
| `sorted_recording` | The exact name of the sorted Code Ocean asset. The trigger also accepts its asset ID. |
| `probe_file` | The annotation JSON filename without `.json`. For `tracks.json`, enter `tracks`. |
| `histology_track_id` | The exact track/layer name inside that JSON, such as `A_shank1`. |
| `ephys_collection` | The recorded probe name, such as `ProbeA`. For a split quadbase, use its common probe name as explained below. |

For a single-shank recording, those five columns are enough:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection
858398,example_sorted,tracks,A_track,ProbeA
```

### Four-shank NP2.0 recorded together

If all four shanks are recorded together as `ProbeA`, repeat `ProbeA` on all
four rows. Add `histology_shank` to identify the physical shank for each track,
and `ephys_shank` to select that shank within the shared recording. Both
columns contain the same number for this recording layout.

CSV shank numbers start at **0**: the four shanks are numbered `0`, `1`, `2`,
and `3`. In this example, the track named `A_shank1` belongs to physical shank
`0`. Enter the number that matches your track; preprocessing does not read
shank numbers from track names.

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,histology_shank,ephys_shank
858398,example_sorted,tracks,A_shank1,ProbeA,0,0
858398,example_sorted,tracks,A_shank2,ProbeA,1,1
858398,example_sorted,tracks,A_shank3,ProbeA,2,2
858398,example_sorted,tracks,A_shank4,ProbeA,3,3
```

Here `ProbeA` contains four shanks, so `ephys_shank` distinguishes the data
for the four tracks. You can omit `logical_probe`.

### Quadbase shanks recorded separately

A stream is a group of recording channels stored together. A split quadbase
records its four shanks as separate streams, for example
`ProbeB-1`, `ProbeB-2`, `ProbeB-3`, and `ProbeB-4`. Enter the common probe name
`ProbeB` and the physical `histology_shank` for each track:

```csv
mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection,histology_shank
858398,example_sorted,tracks,B_shank1,ProbeB,0
858398,example_sorted,tracks,B_shank2,ProbeB,1
858398,example_sorted,tracks,B_shank3,ProbeB,2
858398,example_sorted,tracks,B_shank4,ProbeB,3
```

Preprocessing checks which physical shank each stream recorded and fills in
the matching stream name: `ProbeB` with physical shank `2`, for example, becomes
`ProbeB-3`. That stream contains one shank, whose index within the stream is
`0`. Preprocessing fills in this distinction; leave `ephys_shank` out of this
form of the CSV. If the saved recording information cannot identify a unique
stream, the run reports which mapping it could not resolve.

You can combine single-shank, four-shank, and split quadbase rows in one
manifest. Use the seven-column header from the four-shank example and leave
unused shank cells blank: both cells for a single-shank row, or just
`ephys_shank` for a quadbase row written in the form above.

### Save the annotations asset

Save the CSV with `manifest` in its filename, for example
`858398_manifest.csv`. Put it and every annotation JSON named by `probe_file`
in one internal Code Ocean asset, with exactly one manifest CSV. Give that
asset's ID to the trigger. Its [run instructions](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-trigger-capsule#run-a-mouse)
cover the launch.

Use the correct `sorted_recording` on each row when a mouse has several
recording sessions. Dates in annotation filenames or track names are labels;
the `sorted_recording` column selects the recording to process.

### Other columns and existing manifests

Existing manifests can continue to name exact streams and specify
`ephys_shank`. The [column reference](docs/manifest_reference.md) explains that
form, optional columns, and older names such as `probe_id`, `probe_name`, and
`probe_shank`. You do not need to add those columns to the examples above.

### Building a manifest from a mapping

If you prefer to generate the CSV, the optional
[manifest builder](docs/building_a_manifest.md) takes a configuration file
mapping tracks to recordings and streams.

## Registration overrides

Most runs use the registrations already supplied with the imaging data. If
the brain was registered separately, choose that registration through the
trigger's **Registration override** field. If it targets a reference image
other than the standard SmartSPIM template, also set **Template-to-CCF
override** to that template's registration to the Allen brain atlas. The
trigger attaches those assets for the workers.

The [registration guide](docs/registration_overrides.md) describes the files
needed for these cases. Registration is a choice for the whole brain, so new
manifests can leave it to the run settings. Existing CSV registration columns
remain supported.

## Run and check the output

Give the annotations asset's ID to the trigger and follow its
[run instructions](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-trigger-capsule#run-a-mouse).
The trigger processes the histology and recordings and creates one output data
asset. Use that output asset with the alignment GUI. For a direct worker run,
follow the [worker instructions](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-pipeline-capsule#standalone-use-the-monolith).

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
