# Datapackage And External Asset Contract

Status: current implementation. The producer writes schema 4.1.0; the GUI
supports input schema majors 3 and 4. Saved GUI alignment outputs use a
separate schema.

## Purpose And Layers

The same mouse package must work in Code Ocean, on a desktop, or on another
computer with a different storage layout. The datapackage describes data and
references assets; deployment configuration locates local copies. The package
need not contain all external dependencies.

| Layer | Responsibility |
|---|---|
| Datapackage structure | Mouse, recordings, collections, channel tables, histology picks, volumes, and transforms |
| Asset descriptions | Logical asset keys and available identity/provenance metadata |
| Location resolution | Consumer runtime configuration maps references onto the current machine |
| Ephys collection contents | Conversion and GUI share ALF filenames, QC directories, array ordering, and shank/row-map conventions |

The JSON schema covers structure and references, not an inventory of all
metrics or their array semantics. See the
[conversion collection contract](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion/blob/main/docs/shank_channel_metadata_spec.md)
and [GUI input guide](https://github.com/AllenNeuralDynamics/ibl-ephys-alignment-gui/blob/main/docs/input_datapackage_contract.md).

## Datapackage Structure

| Field | Meaning |
|---|---|
| `schema_version` | Input contract version, independent of package release versions |
| `mouse_id` | Subject identity |
| `platform` | Producer metadata; does not choose the GUI resolver policy |
| `external_assets` | Logical asset key to descriptor |
| `transforms` | Four image-to-template and template-to-CCF transform references |
| `histology` | Image-space volume/geometry references and optional CCF QC references |
| `probes` | `recording_id -> ephys_collection -> collection entry` |

A collection entry records `probe_id`, `logical_probe`, `ephys_collection`,
`num_shanks`, `ephys`, `channel_table`, and `xyz_picks`. `logical_probe` can
repeat across split streams. Stable collection identity is
`(recording_id, ephys_collection)`.

`ephys` references the **directory** containing ALF data and QC products.
`channel_table` explicitly references `local_coordinates`, `raw_ind`,
`shank_ind`, and optional `contact_id`. Individual spike, cluster, RMS, PSD,
line-amplitude, and matrix files are loaded by convention inside the collection.
They are not enumerated by the datapackage.

`xyz_picks` maps physical tracks to collection shanks. `histology_shank` is the
physical track index; `ephys_shank` is the index within the collection. For a
split stream these can be 3 and 0 respectively. Manifest rows for one
recording/collection become one entry with multiple picks. Legacy `shank` is
one-based; explicit shank fields are zero-based. The GUI normalizes sorted
unique pick values, preferring `ephys_shank` over `shank`; a one-shank
collection uses its first pick. Counts alone do not prove that channel shank
identity and histology picks agree.

Histology `emit_qc` controls optional CCF volumes, CCF/bregma picks, and Slicer
outputs. It does not disable electrophysiology RMS/PSD/correlation products.
The GUI uses image-space picks and recomputes CCF coordinates. Schema 4.0
replaced the duplicate pipeline volume with `registration_pipeline_geometry`;
4.1 removed the unused CCF registration reference.

## References And Asset Descriptions

Canonical local references have two fields:

```json
{"asset": null, "path": "image_space_histology/histology_registration.nrrd"}
```

They resolve relative to the directory containing `datapackage.json`, not the
working directory, `/data`, or a search root. External references name a
logical key in `external_assets` and a path inside that asset:

```json
{
  "asset": "smartspim",
  "path": "image_atlas_alignment/Ex_639_Em_680/ls_to_template_SyN_0GenericAffine.mat"
}
```

An illustrative registry fragment is:

```json
{
  "external_assets": {
    "smartspim": {
      "role": "smartspim_registration",
      "name": "SmartSPIM_mouse001_stitched",
      "id": null,
      "uri": null,
      "checksum": null,
      "provenance": {"mounted_at": "/data/SmartSPIM_mouse001_stitched"}
    }
  }
}
```

The producer writes these logical keys:

| Key | Role | Present |
|---|---|---|
| `smartspim` | `smartspim_registration` | Always: the stitched asset holding the volumes and, by default, the registration |
| `spim_template_to_ccf` | `template_to_ccf` | Always: the SmartSPIM LCA template's CCF registration |
| `registration` | `registration_override` | With the manifest's `registration_asset` |
| `template_to_ccf` | `template_to_ccf_override` | With the manifest's `template_to_ccf_asset`; `transforms.template_to_ccf_*` then name it instead of `spim_template_to_ccf` |

These fragments are not complete datapackages. Role/name are required; ID,
URI, checksum, and provenance are optional. A URI should describe the asset
actually containing the referenced files, not an earlier acquisition asset.
Use relative internal paths, without mount paths or `../../data` offsets
between independent assets. Producer platform and provenance mount locations
are descriptive and do not control resolution.

The GUI resolver locates already-present files. It does not fetch URIs or
verify checksums; recording those fields does not establish integrity.

## GUI Resolution And Configuration

`RootSearchResolver` resolves references as follows:

1. Local reference: append its path to the datapackage directory.
2. External reference: try an override by logical key, then recorded name.
   Overrides point directly to asset directories. A missing requested path
   raises an error without falling through to root search.
3. Visit configured roots in order. Within each root try recorded asset name,
   then asset ID, appending the internal path. First existing candidate wins.
4. Raise an asset-resolution error with identity and searched paths on a miss.

Roots contain asset directories; search is not recursive. The GUI does not
implicitly use `/data`, provenance mounts, or the initial folder chooser root.

| Setting | Precedence |
|---|---|
| Roots | Explicit `load_mouse_root(asset_roots=...)`, else nonempty `IBL_ASSET_ROOTS`, else config-file `asset_roots`, else empty |
| Overrides | Config-file mapping, overlaid by `IBL_ASSET_OVERRIDES`, then explicit arguments; later values win for identical keys |

`IBL_ASSET_CONFIG` selects the JSON file. `IBL_ASSET_ROOTS` uses `:` on
Linux/macOS and `;` on Windows. `IBL_ASSET_OVERRIDES` is the mapping itself,
without an outer `asset_overrides` key. Use absolute config paths; `~` inside
JSON is not expanded. Logical-key lookup precedes name-key lookup even if the
name-keyed value came from a higher-precedence configuration source.

## The Same Package On Two Machines

Suppose the registry names `SmartSPIM_mouse001_stitched` and
`spim_template_to_ccf`. On Code Ocean, attach the preprocessed output and the
external assets. Set `IBL_ASSET_ROOTS=/data` and open the directory containing
the JSON, for example `/data/preprocessed_asset/mouse001`. External paths
resolve under `/data/<recorded asset name>/<internal path>`.

On a desktop, copy the complete mouse directory to
`/home/alice/alignments/mouse001`, put external assets under
`/mnt/ibl-assets/<recorded asset name>`, and use:

```json
{"asset_roots": ["/mnt/ibl-assets"], "asset_overrides": {}}
```

Set `IBL_ASSET_CONFIG` to that JSON file's absolute path before launching.
For a renamed copy, add a name-keyed override such as
`"SmartSPIM_mouse001_stitched": "/mnt/histology/mouse001"`.
Name-based root search and name-keyed overrides suit multiple subjects;
a logical-key override such as `smartspim` applies across packages opened in
that process. Windows uses the same JSON structure with paths such as
`D:/ibl-assets`. Keep location configuration separate from the datapackage.

Moving the whole mouse tree preserves local references. External assets must
also be available and configured. Moving only the JSON or omitting ephys QC
subdirectories does not preserve a usable package.

## Validation And Compatibility

Producer `DataPackage.validate_paths()` resolves explicit file references and
requires ephys references to be directories containing
`channels.localCoordinates.npy`. It does not check every convention-loaded
metric or its shape. Producer validation accepts roots/overrides as arguments;
it does not inherit the GUI environment/config-file loader.

GUI schema validation checks structure. External resolution checks locations;
subsequent loading checks channel arrays and individual products. Missing ALF
objects can disable plots with warnings. Schema validity alone does not
establish data completeness or correct shank joins.

GUI bundled schemas are 3.0.0, 3.1.0, 3.2.0, 4.0.0, and 4.1.0. Another
version uses the newest bundled schema in the same supported major not newer
than the input, with the version constant relaxed. Breaking consumed fields
requires a major bump. Filename/array changes need coordinated review even
when the JSON shape stays unchanged.

Regenerate old string-path datapackages with `--datapackage-only`; see
[the regeneration instructions](https://github.com/AllenNeuralDynamics/aind-ibl-ephys-alignment-preprocessing/blob/main/README.md#regenerating-only-datapackagejson). For read-only prior
results on Code Ocean, use `--source-results` to copy them into writable
results before regeneration. Regenerating metadata does not repair missing or
inconsistent ephys artifacts.

Implementation authorities are `datapackage.py` and bundled schemas in this
package, and GUI `io/datapackage_loader.py`, `io/datapackage_schema.py`, and
`services/ephys_data.py`. Keep producer path-validation and GUI resolver tests
aligned when changing references.
