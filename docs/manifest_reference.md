# Manifest column reference

For a new CSV, start with the [manifest guide](../README.md#manifest-csv).
It gives complete examples for single-shank recordings, four-shank recordings
in one stream, and quadbase recordings split across streams. This reference
covers explicit stream mappings, additional columns, and existing manifests.

## Find the recording's probe names

Use the names saved with the recording rather than the labels you chose for
the histology tracks. In a sorted asset, the directories under `postprocessed/`
include the stream name. For example, a directory named
`experiment1_Record Node 103#Neuropix-PXI-111.ProbeB-3_recording1.zarr` contains
stream `ProbeB-3`. If the same recording has `ProbeB-1` through `ProbeB-4`,
the README's simplified quadbase form uses their common name, `ProbeB`, and
a physical shank number.

## Physical shanks and shanks within a stream

`histology_shank` identifies a physical shank on the probe. `ephys_shank`
identifies a shank within the selected `ephys_collection`. Both are zero-based.
They have the same value when all of a probe's shanks share one collection.
They differ for a split quadbase: physical shank `2` can be recorded in
`ProbeB-3`, which contains just one shank, indexed locally as `0`.

The README's quadbase example names `ProbeB` and supplies `histology_shank`.
Discovery verifies the sorting's recorded probe geometry, resolves each row
to its separate stream, and sets that stream's local `ephys_shank` to `0`.
The saved discovery manifest contains the resolved values.

## Naming the exact quadbase streams

Existing manifests can name the separate streams directly. Set `ephys_shank`
to `0` on every row, since each stream contains one shank:

| `logical_probe` | `ephys_collection` | `histology_shank` | `ephys_shank` |
|---|---|---|---|
| `ProbeB` | `ProbeB-1` | 0 | 0 |
| `ProbeB` | `ProbeB-2` | 1 | 0 |
| `ProbeB` | `ProbeB-3` | 2 | 0 |
| `ProbeB` | `ProbeB-4` | 3 | 0 |

The five required columns from the README are still needed; the table shows
only the mapping columns. `logical_probe` groups these collections under the
physical probe name. It is optional and defaults to `ephys_collection`.

Exact collection names keep their existing meaning. In particular, setting
`ephys_collection=ProbeB-3` and `ephys_shank=2` asks for the third shank
within that collection; it does not mean the physical shank recorded there.

## Additional columns

| Column | Use |
|---|---|
| `logical_probe` | Group separate collections belonging to one physical probe. Defaults to `ephys_collection`. |
| `histology_shank` | Identify the physical shank for a track. Needed when a probe name must be resolved to one of several split streams. |
| `ephys_shank` | Select a shank within a collection. Include it for multiple tracks sharing a four-shank collection. |
| `annotation_format` | Annotation file extension, default `json`. With the default, `probe_file=tracks` selects `tracks.json`. |
| `surface_finding` | Identify a separate surface-finding recording to include. On Code Ocean, use its asset name; for local runs, use its location relative to `data_root`. |
| `main_recording_id` | Override the automatically selected main recording for a conversion unit. See the [inventory and conversion-plan contract](ephys_inventory.md). |
| `sorting_family_id` | Override the sorting-family selection for a conversion unit. See the same contract. |

## Older column names

| Older name | Current name |
|---|---|
| `probe_id` | `histology_track_id` |
| `probe_name` | `ephys_collection` |
| `probe_shank` | Compatibility alias for both `histology_shank` and `ephys_shank` |

Older names remain accepted. `probe_shank` is optional; new manifests should
use the explicit columns shown in the README examples. The parser fills an
omitted `ephys_shank` from `probe_shank`, then `histology_shank`. For tracks
sharing a four-shank collection, include `ephys_shank` explicitly as shown in
the guide so the CSV also satisfies the structural uniqueness checks.

Some existing quadbase CSVs name the common probe, such as `ProbeB`, and put
the physical shank in `ephys_shank` alone. Discovery accepts that layout only
when there is no exact collection named `ProbeB` and the recorded geometry
identifies a unique split stream. It then records the physical
`histology_shank`, resolved collection, and local `ephys_shank=0` in the saved
manifest. The log and `coverage.json` notes explain the interpretation.
Missing, contradictory, or unreadable geometry prevents that resolution.

## Registration settings

New manifests should leave registration choices to the run's
`registration_asset` and `template_to_ccf_asset` settings. The corresponding
CSV columns remain supported. Nonblank values must agree across rows and
with the run settings, including the transform subdirectory.

Discovery saves the resolved choices in `registration_settings.json` beside
its manifest. Keep that file with the manifest when staging a histology or
pack-only run. Pack retains it beside `datapackage.json`. See the
[registration guide](registration_overrides.md) for the required files and
geometry information.

## Structural checks

All rows must name the same mouse. The following combinations must be unique:

| Combination | Purpose |
|---|---|
| `(mouseid, histology_track_id, histology_shank)` | Histology track filenames |
| `(recording_id, ephys_collection, histology_shank)` | GUI track filenames |
| `(recording_id, ephys_collection, ephys_shank)` | Ephys shank mapping |

`recording_id` is derived by preprocessing; it is not a required CSV column.
The accepted column names are declared in `types.MANIFEST_COLUMNS` and shared
by the parser and preflight validation.
