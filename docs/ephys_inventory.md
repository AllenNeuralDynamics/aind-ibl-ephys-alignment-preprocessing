# Discovery inventory handoff

Discovery writes `ephys_inventory.json` before filtering manifest rows. The
versioned observation contract and metadata readers live in
[`aind-ephys-ibl-gui-conversion`](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion/blob/main/docs/ephys_inventory.md).
This artifact accompanies `coverage.json`; it does not change the GUI datapackage
schema. Conversion decisions are recorded separately from observations.

Discovery also writes `ephys_conversion_plan.json`. Its `units` map contains a
plan for every requested recording/collection before viability filtering. Each
plan binds to its inventory slice by SHA-256 hash and records the longest raw
block as main, its analyzer family, physical shank coverage, candidate QC inputs,
and selection findings. The registry also binds to the canonical inventory.
See the [conversion-plan contract](https://github.com/AllenNeuralDynamics/aind-ephys-ibl-gui-conversion/blob/main/docs/conversion_plan.md).

The ephys worker validates and executes these decisions through a shared
preflight context for spike and continuous extraction. Plans retain
`execution_enabled=false`; a separate collection-level `ephys_extraction.json`
records activation, outcome, and actual block/segment membership. Ready plans
execute; an incomplete plan executes only for known main shanks without sorting.
Other unresolved or blocked plans stop extraction before ALF arrays are written.

## Optional selection overrides

Existing manifests use duration-based defaults without additional columns.
Two optional CSV columns can override a unit's plan:

| Column | Value |
|---|---|
| `main_recording_id` | A raw recording `id` from its discovery inventory |
| `sorting_family_id` | A family `id` from the plan's `sorting_families` candidates |

Blank rows defer to another row's override for that unit, or to the default.
All nonblank values within one resolved recording/collection must agree;
conflicts stop discovery after saving the inventory. An explicit family still
has to belong to the selected main block and pass source/geometry checks.
Overrides determine the saved selection consumed by both extractors.

`stage_discover(..., inventory_sources=...)` accepts explicit asset descriptors.
Each descriptor has a stable `key`, a `raw` or `sorted` role, optional asset
`id`/`name`/`aliases`, and declared `source_keys` for a sorting. Deployment fields
are `local_mount` relative to `data_root`, or an explicitly resolved S3 `uri` and
`anonymous` access policy. Deployment fields are not serialized into the inventory.

Without descriptors, discovery inventories manifest-related local sorted/raw
mounts and any supplied surface-finding mount. Unknown asset IDs stay unknown.
The reader never constructs a bucket or searches Code Ocean from a recording name.
With `skip_ephys`, no inventory reads or fan-out inventory files are generated.

Unavailable sources remain in `assets` with `metadata_status=incomplete`.
Read failures have diagnostic codes and asset-relative metadata paths. Geometry
checks with insufficient identity remain unknown. Inventory findings remain
observations; the planner evaluates whether they prevent execution. Inventory
creation itself does not choose among several analyzers on one stream.

## Outputs and stage prerequisites

| Stage | Reads | Writes |
|---|---|---|
| Discover | Related mounted asset metadata under `data_root`, or bounded metadata from explicitly resolved S3 roots | `/results/ephys_inventory.json`; `/results/stream_<unit>/ephys_inventory.json`; relative `inventory` reference in each unit's `config.json` |
| Ephys launcher | Canonical inventory beside discovery's filtered `manifest.csv`, when present | A self-contained slice beside each emitted unit config |
| Ephys worker | Saved plan and inventory beside its staged unit config; selected mounted raw zarrs and analyzers | `/results/<unit>/ephys_inventory.json`, saved unit plan, and collection-level `ephys_extraction.json` plus ALF/QC arrays |
| Pack | Canonical inventory beside discovery's filtered `manifest.csv`, when present | `/results/<mouseid>/ephys_inventory.json`, beside `datapackage.json` |

Each inventory output now has an `ephys_conversion_plan.json` companion. The
discovery and final mouse copies contain the complete registry; each fan-out
directory and worker namespace contains only that unit's plan. Unit configs
reference it through the relative `conversion_plan` field. Workers check that
the plan version, unit, and inventory hash agree before retaining the copy.
Pack and the per-sort launcher read the canonical registry beside the filtered
manifest. A present registry must match its inventory; it is never silently
replaced with new defaults. An older inventory-only intermediate can generate
plans using the filtered manifest and its optional overrides.

Worker slices retain the candidate-reference closure as well as their declared
source assets. Ambiguous sources therefore remain inspectable after Flatten;
they do not become implicit selections. Worker output paths use the existing
per-unit namespace to avoid Collect collisions.

The complete discovery artifact must reach pack before intermediate cleanup.
Legacy discovery outputs without an inventory remain usable. A present but
unreadable canonical artifact raises instead of silently discarding diagnostics.

The S3 adapter uses signed access by default and unsigned access only when the
descriptor explicitly requests it. It lists immediate child prefixes and reads
JSON metadata and bounded byte ranges for timestamp scalars, never voltage
chunks or analyzer extensions. Timestamp ranges must return exactly the requested
bytes; a server that ignores the range is rejected. Sources without an
accessible mount or S3 location are reported unavailable; this release does not
add Code Ocean credentials to the worker to read internal raw assets remotely.

Compressed timestamp chunks use a separate bounded read capped at 8 MiB.
The converter enforces the same limit before decoding into a fixed-size
buffer. These reads inspect timestamps only; discovery still does not read
voltage data.

## Deployment

Release the conversion package with `inventory.py` first, then preprocessing,
then rebuild the worker capsule against those releases. Update the trigger only
after the worker accepts the `inventory-sources` named parameter. These source
changes do not themselves release packages, rebuild images, or launch runs.

## Extraction preconditions and compatibility

Workers bind source keys to mounted asset roots, then verify exact references,
segment lengths, sample rates, channel geometry and order, and analyzer unit
counts against discovery. A present plan never falls back to Neo discovery or
another experiment's sorting. Known missing main sorting shanks remain missing;
QC can still cover their raw channels. A single sorted segment can match one
unique raw segment after source, geometry, rate, and acquisition identity
checks. The plan records that segment and its local clock; the worker uses it
for spikes and main QC. Repeated lengths and disagreeing group mappings remain
blocked, as do multisegment analyzers. Other QC segments keep independent clocks.

The execution receipt lives inside the ordinary ephys collection and therefore
survives pack with its arrays. It records `started`, `complete`, or `failed`, the
original decision, and skipped short QC blocks. Saved worker diagnostics remain
available when preflight rejects the inputs. Worker configs without saved plans
use the legacy extractor with a warning; older direct library and monolithic
pipeline calls retain their existing behavior.
