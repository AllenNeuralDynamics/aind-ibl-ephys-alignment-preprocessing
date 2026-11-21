# Registration overrides

Use a registration override when a brain's transforms were generated separately
from its stitched imaging asset, or when the brain was registered to another
template. These are brain-wide run settings, independent of the track mapping
in the [manifest](../README.md#manifest-csv). Leave them unset to use the usual
stitched registration and SmartSPIM-template-to-CCF registration.

## Separate registration asset

Set the brain-wide **Registration override** App Panel field, or pass the
transform directory relative to `data_root` on the CLI:

```bash
aind-ibl-preprocess --data-root /data --results-root /results \
    --manifest manifest.csv --neuroglancer annotations.json \
    --registration-asset SmartSPIM_750108_reg/ccf_Ex_639_Em_667
```

The programmatic setting is `PipelineConfig(registration_asset=...)`.
Existing CSVs can still supply a `registration_asset` column. If both sources
are provided, they must identify the same directory. Leaving both unset uses
the stitched asset's registration.

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

### Generate the sidecar from ANTs outputs

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

## Separate template-to-CCF registration

A brain registered to a template other than the SmartSPIM LCA template (an
iDISCO template, say) needs that template's CCF registration. Set the
**Template-to-CCF override** App Panel field or `--template-to-ccf-asset` to its
directory relative to `data_root`:

```bash
aind-ibl-preprocess --data-root /data --results-root /results \
    --manifest manifest.csv --neuroglancer annotations.json \
    --registration-asset idisco_reg/750107 \
    --template-to-ccf-asset AIND-MouseiDISCO-v1.0-to-CCFv3-r1_2026-10-22_17-34-50
```

`PipelineConfig.template_to_ccf_asset` and the legacy CSV column of the same
name are also accepted, with the same conflict checks.

The directory must contain `template_to_ccf_0GenericAffine.mat` and
`template_to_ccf_1InverseWarp.nii.gz`, ANTs SyN outputs with CCF fixed and the
template moving; `emit_qc=True` also needs `template_to_ccf_1Warp.nii.gz`. They replace the
pipeline's `spim_template_to_ccf_syn_*_25` transforms in the probe and QC
chains, and the datapackage refers to them under the `template_to_ccf` asset,
named after the directory's first component under `data_root`.

The override requires a `registration_asset` with a sidecar, whose `fixed_domain`
is checked against this template. If the directory ships its own sidecar
(`template_to_ccf_transform_information.json` or `transform_information.json`),
the registration's `fixed_domain` must match its `moving_domain`. Without one,
the registration's `fixed_domain` must *not* match the SmartSPIM LCA template;
to replace the LCA template's own CCF registration, ship the sidecar.
