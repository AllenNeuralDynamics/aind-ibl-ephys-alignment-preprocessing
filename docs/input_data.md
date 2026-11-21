# Input folders for direct runs

The trigger finds and attaches the data assets needed by its worker runs. Use
this guide when running the library locally or preparing a direct worker run.
Direct worker setup is described in the [worker README](https://github.com/AllenNeuralDynamics/aind-ibl-preprocess-pipeline-capsule).

## Data root layout

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

Also provide the raw recording folders for ephys extraction. The imaging
asset needs fused OME-Zarr volumes under `image_tile_fusing/OMEZarr/` and its
image-to-template registration under `image_atlas_alignment/`, or a
[separate registration](registration_overrides.md).

## Reference volumes

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

## Probe annotation files

Per-probe Neuroglancer point annotation files. For each row in the manifest,
the pipeline searches `data_root` for a file matching the glob pattern
`*/<probe_file>.<annotation_format>` (e.g., `*/track_annotations_probeA.json`).

These files contain 3D point annotations marking the probe track in the
SmartSPIM image space. They are converted through the transform chain:

```
Neuroglancer pixels --> SPIM (LPS) --> Template --> CCF (LPS/um) --> Bregma (IBL)
```

## Electrophysiology data

Unless `--skip-ephys` is passed, the pipeline expects spike-sorted data at
`data_root/<sorted_recording>/` for each unique recording in the manifest. The
ephys extraction (delegated to
[aind-ephys-ibl-gui-conversion](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion))
converts this into IBL ALF format (e.g., `spikes.times.npy`,
`spikes.clusters.npy`, `channels.localCoordinates.npy`, etc.).

The staged Code Ocean workflow processes each recording/collection pair
separately. Repeated track rows do not require duplicate extraction of that
collection. Provide the raw recordings as well as sorting outputs: continuous
recording data is needed for the ephys quality-control calculations.
