# Image-source fixtures

These are the image layers from real Neuroglancer states under
`/mnt/vast/scratch/ephys/persist/data/light_sheet/NeuroglancerAnnotations`.
Only layer `name`, `type`, and `source` are retained; annotation points and viewer
state are omitted.

- `nested_image_sources.json`: `717381/717381_Neuroglancer_annotation_GL.json`.
  Includes nested S3 prefixes, protocol wrappers, trailing slashes, and an object
  source with `url`.
- `transformed_image_sources.json`:
  `837772/837772-probe-annotations-CK-2026-09-29.json`. Includes plain string
  sources and a source object carrying a transform.

Tests also remove the modality prefix from these sources to cover the new asset
naming convention. Copies live in both repositories so their checks run
independently.
