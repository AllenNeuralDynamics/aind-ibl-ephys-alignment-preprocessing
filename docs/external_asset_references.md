# Design: durable external-asset references in the datapackage

Status: **implemented in the schema 3.0.0 feature branches** for reference
objects, external asset registry, runtime resolver roots/overrides, and
datapackage-only regeneration. Remaining future decision points are marked
⟨DECIDE⟩.

Scope: how the preprocessed datapackage (`datapackage.json`, produced by
`aind-ibl-ephys-alignment-preprocessing`) should reference files that live in
*other* data assets — the SmartSPIM histology asset and the shared
`spim_template_to_ccf` asset — so that the alignment GUI
(`ibl-ephys-alignment-gui`) can resolve them on **any** platform (Code Ocean,
desktop) and years from now.

Sibling contract doc: channel geometry —
`aind-ephys-ibl-gui-conversion/docs/shank_channel_metadata_spec.md`.

## 1. The problem (and why the obvious fixes are wrong)

The transform-chain files (LS→template, template→CCF) do not live in the
preprocessed asset; they live in separate assets. The current datapackage stores
them as `..`-relative paths from the mouse root:

```
"image_to_template_affine":
  "../../data/SmartSPIM_791094_..._stitched_.../image_atlas_alignment/Ex_639_Em_680/ls_to_template_SyN_0GenericAffine.mat"
```

At **read** time on Code Ocean this resolves to `/data/data/SmartSPIM…` — a
doubled `/data` — and the GUI dies with an HDF5 "not an HDF5 file" error.

Root cause: the `..` count and the `data/` component are computed from the
**write-time absolute layout** (outputs under `/results/…`, inputs under
`/data/…`, common ancestor `/`), which differs from the **read-time layout**
(the datapackage is now a mounted asset under `/data/<name>/…`, so the asset mount
inserts an extra level). A relative offset between two independent Code Ocean
assets is not stable across the producing and consuming runs.

Two "obvious" fixes are both wrong for the long term:

- **Relative-to-`/data`** (`<asset>/<within>`, resolved against `/data`) forces
  the GUI to hold the `/data` anchor → bakes Code Ocean into the GUI. Rejected:
  the GUI must be platform-agnostic (people run it on the desktop too).
- **Absolute `/data/<asset>/…`** stamped by the producer fixes CO→CO today, but
  it still encodes a *location* as the reference, so it does not travel to
  other platforms or survive re-hosting. Acceptable as a stop-gap, not as the
  durable shape.

The category error in all of these: **encoding a *location* as if it were a
*reference*.** Locations are the least durable property of an asset.

## 2. Principles

Separate three things the current design conflates:

- **Identity** — *what* the asset is, immutably (data-asset id, canonical URI,
  content hash). Survives platform, mount, re-hosting, time.
- **Location** — *where* a copy currently sits (`/data/<name>`, a desktop dir,
  an `s3://` key). Platform- and time-specific.
- **Resolution** — *how* a given consumer turns identity → a local path it can
  open. A consumer/platform policy, not a property of the asset.

Durable rule: **the datapackage references assets by identity; the consumer
resolves identity → location via an injected, platform-specific policy.** No
platform's location is privileged in the datapackage, and no platform logic lives
in the GUI core.

## 3. Datapackage shape

### 3.1 Every file reference is `(asset, path_within_asset)`

```json
"transforms": {
  "image_to_template_affine": {
    "asset": "smartspim",
    "path": "image_atlas_alignment/Ex_639_Em_680/ls_to_template_SyN_0GenericAffine.mat"
  },
  "template_to_ccf_affine": {
    "asset": "spim_template_to_ccf",
    "path": "syn_0GenericAffine.mat"
  }
}
```

`path` is the asset's **internal** layout — stable regardless of where the
asset lives. No `..`, no mount, ever.

### 3.2 An asset registry keyed by a logical key

```json
"platform": "code_ocean",
"external_assets": {
  "smartspim": {
    "role": "smartspim_registration",
    "name": "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04",
    "id":   "<data-asset-uuid, if available>",
    "uri":  "s3://aind-open-data/SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04",
    "checksum": null,
    "provenance": { "mounted_at": "/data/SmartSPIM_791094_..._stitched_...", "on": "code_ocean" }
  },
  "spim_template_to_ccf": {
    "role": "template_to_ccf",
    "name": "spim_template_to_ccf",
    "id":   "<uuid>",
    "uri":  "s3://aind-open-data/spim_template_to_ccf",
    "checksum": null
  }
}
```

Notes:
- **No authoritative mount/absolute path.** `provenance.mounted_at` is
  debugging/legibility only — explicitly *not* used for resolution.
- `id` + `uri` are the durable identity; `name`/`role` are human/semantic.
- `role` lets the GUI/tooling reason semantically; the reference `asset` key is
  a short logical handle into this registry.
- Do not block v1 implementation on having every identity field. `name` +
  `role` are the required, always-available reference descriptors; `id`/`uri`
  are recorded when discoverable but may be `null`. **Name as a fallback lookup
  key is explicitly acceptable** — for many deployments id/uri are not easy to
  come by, and a name the consumer resolves is still a reference, not a baked
  location. Don't let the perfect be the enemy of the good.
  - SmartSPIM: `name = asset_path.name`. If `uri` is recorded it MUST point at
    the same asset that contains the transforms (the `..._stitched_...` asset,
    i.e. `asset_path`) — NOT the raw Neuroglancer/acquisition source (that would
    reference an asset that doesn't hold `image_atlas_alignment/`). `id = null`
    if the producer cannot discover the Code Ocean data-asset id.
  - `spim_template_to_ccf`: `name = config.template_to_ccf_dir.name`;
    `uri`/`id` can be null unless/until config/discovery provides them.
  - `role` and `name` are required for every external asset.

### 3.3 Internal outputs: a null/absent `asset`

Unify internal outputs into the same reference model, but distinguish
"datapackage-local" **structurally** rather than with a reserved key: a
reference whose `asset` is **absent (or `null`)** resolves against the
datapackage's own directory.

```json
"histology": { "image_space": {
  "registration": { "asset": null, "path": "histology_img/registration.nii.gz" }
}}
```

This avoids the reserved-word pitfall — there is no magic `"self"` string that
a real external asset key could ever collide with. The distinction is simply
"does this reference name an external asset or not?" `asset: null` and an
omitted `asset` are treated identically on read; **`asset: null` is the
canonical producer output** so every file reference has the same two fields.

Decision: **unify**. One reference model (`{asset, path}` where `asset` may be
null) is easier to validate and consume; local files set `asset: null`.

If implementation pressure is high, transforms may be converted first while
internal outputs remain plain strings for one transitional version. But that is
a compatibility bridge, not the target schema.

### 3.4 Schema versioning

Changing path fields from strings to `{asset, path}` objects is a schema-shape
break, not a small extension. The clean implementation should therefore bump
the datapackage schema to **3.0.0** and have the GUI refuse older external-path
schemas with a clear "re-run preprocessing" message.

A minor bump is only acceptable if the GUI loader supports both shapes:

- legacy string paths (`2.x`)
- reference objects (`2.y+`)

Recommendation: use **3.0.0** for the unified local + external reference
schema (local = `asset: null`; omitted `asset` accepted by readers only as a
tolerance).

## 4. The resolver (consumer side)

The resolver is **code in the GUI**, plus **deployment-supplied runtime config**
(where assets live on this machine). That config may come from environment
variables, GUI settings, or a small desktop config file, but it is not copied
into the datapackage.

### 4.1 Interface (platform-agnostic; lives in the GUI core)

```python
class AssetResolver(Protocol):
    def resolve(self, asset: AssetRef | None, within: str) -> Path: ...
        # AssetRef = the registry entry (key/role/name/id/uri/checksum)
```

The GUI core only ever does `open(resolver.resolve(asset, within))`. It has no
`/data`, no mount rule, no platform branch.

### 4.2 Default implementation: locate an already-present asset

```python
class RootSearchResolver:
    """Find <root>/<asset-key>/<within> across deployment-supplied roots.
    A reference with no external `asset` resolves against the datapackage dir."""
    def __init__(
        self,
        datapackage_dir: Path,
        asset_roots: list[Path],
        asset_overrides: dict[str, Path] | None = None,
    ):
        self.datapackage_dir = datapackage_dir
        self.asset_roots = asset_roots
        self.asset_overrides = asset_overrides or {}

    def resolve(self, asset, within):
        if asset is None:                                   # datapackage-local
            return self.datapackage_dir / within
        # explicit override, keyed by logical key OR asset name
        for ovr_key in (asset.key, asset.name):
            if ovr_key and ovr_key in self.asset_overrides:
                path = self.asset_overrides[ovr_key] / within
                if path.exists():
                    return path
                raise ReferenceNotFound(asset, within, path)  # explicit: hard-fail
        for root in self.asset_roots:
            for key in (asset.name, asset.id):     # by name, then id
                if key and (root / key / within).exists():
                    return root / key / within
        raise AssetNotFound(asset, self.asset_roots)
```

Resolution order:

1. no external `asset` (absent/null) → `datapackage_dir / within`
2. explicit override, matched by **logical key or asset name** — a miss here is
   a hard error (`ReferenceNotFound`), *not* a fall-through: an explicit
   override means "it is here," so a miss is a real misconfiguration, not a hint
   to look elsewhere.
3. search roots by `asset.name`
4. search roots by `asset.id`

Override keys accept **either** the logical asset key (`"smartspim"`) **or** the
asset name. Keying by name matters for multi-subject desktop use: a single
logical-key override (`"smartspim" -> dir`) is session-scoped and becomes wrong
the moment a second subject is opened, whereas name-keyed overrides — and plain
name-based root search — generalize across subjects. The intended primary
desktop path is therefore **name-based root search** (name local asset folders
like the Code Ocean asset names); overrides are the escape hatch for
non-conventional placements.

### 4.3 Config channel (the only place platform truth enters) ⟨DECIDE⟩

The asset roots come from the **deployment**, as data:

```python
roots = [Path(p) for p in os.environ.get("IBL_ASSET_ROOTS", "").split(os.pathsep) if p]
overrides = _load_overrides()  # from IBL_ASSET_OVERRIDES / GUI settings; may be {}
resolver = RootSearchResolver(datapackage_dir, roots, overrides)
```

- **Code Ocean capsule** sets `IBL_ASSET_ROOTS=/data` — the entire
  CO-specific knowledge, living in the capsule, not the GUI.
- **Desktop** sets `IBL_ASSET_ROOTS=/home/<user>/ibl_assets` (or the app writes
  it into its own settings).

Env var vs a field in the GUI's existing settings/CLI: env is least-friction
for a capsule; a settings field is nicer for desktop. Likely support both (env
overrides settings).

Suggested config channels:

- `IBL_ASSET_CONFIG`: optional path to a JSON file with `asset_roots` and
  `asset_overrides`.
- `IBL_ASSET_ROOTS`: `os.pathsep`-separated roots such as `/data`.
- `IBL_ASSET_OVERRIDES`: optional JSON mapping from logical asset key or asset
  name to absolute directory, or the equivalent GUI settings field.

Concrete desktop/settings-file shape:

```json
{
  "asset_roots": [
    "/data",
    "/mnt/nas/aind-assets"
  ],
  "asset_overrides": {
    "SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04": "/mnt/nas/histology/SmartSPIM_791094_2025-11-12_16-36-24_stitched_2025-11-14_09-36-04",
    "spim_template_to_ccf": "/mnt/nas/shared/spim_template_to_ccf"
  }
}
```

This config is **consumer-local runtime state**. It is not copied into the
datapackage and is allowed to vary across Code Ocean, desktop, and NAS
deployments.

### 4.4 What it reproduces on Code Ocean

`root=/data`, `key=SmartSPIM_791094_…`,
`within=image_atlas_alignment/Ex_639_Em_680/ls_to_template_SyN_0GenericAffine.mat`
→ the real `/data/SmartSPIM_791094_…/…` — no `..`, no doubling. The
"producer stamps an absolute `/data` path" idea is simply *what this resolver
computes*; we moved that computation out of the (frozen, platform-bound)
datapackage into a (live, swappable) resolver.

### 4.5 Error behavior (agnostic, legible)

`AssetNotFound` prints the datapackage's own `name`/`id`/`uri`/`role`/`platform`:

> This datapackage was made on **code_ocean** and needs assets attached:
> `SmartSPIM_791094_…` (role smartspim_registration),
> `spim_template_to_ccf` (role template_to_ccf).

Reads datapackage data; encodes no platform logic. Turns the HDF5 stack trace into
an actionable message.

### 4.6 Future resolvers (same interface, no core change)

- `FetchingResolver(base, cache_root)` — on a miss, download `asset.uri`
  into `cache_root/<name>` and return it (desktop "I don't have it"). `uri`
  might be `s3://...`, `file://...`, `https://...`, or another resolvable
  scheme; S3 is only one possible transport.
- `VerifyingResolver(base)` — check `asset.checksum` after resolve.

## 5. Why reference-by-identity, not location

The principle is that the datapackage names *what* an asset is and lets the
consumer resolve *where* it is — never the reverse. The **strength** of that
reference is best-effort, and that is fine:

- **Reproducibility.** `id`/`uri`/`checksum`, when present, pin the exact
  transform/template version. When only `name` is available, `name` is still a
  lookup key the consumer resolves — weaker than an immutable id (a re-stitch
  or a template v2 could reuse a name), but far better than a baked mount path,
  and acceptable where id/uri aren't readily available.
- **Portability.** One datapackage, N platforms, each with its own resolver.
- **Survives re-hosting/re-mounting** — the exact failure mode we hit; a
  location does not survive it, whereas reference keys do.

## 6. Migration

Existing datapackages carry the broken `../../data/…` form. Do **not** teach the
GUI to heal them (that requires `/data` logic in the GUI). Instead regenerate
them with the corrected producer, or run a one-off migration script where
`/data` is real. Bump the datapackage schema major version for the unified
reference-object shape so the GUI can tell old from new and refuse old ones
with a clear "re-run preprocessing" message.

If a transitional minor schema is chosen instead, the GUI must implement both
legacy string resolution and new reference-object resolution. That is more code
and should be treated as a migration bridge only.

## 7. Producer-side validation

The producer still needs to validate that every referenced file exists at write
time. But after this change, serialized datapackage paths are references, not
local filesystem locations.

Therefore validation must use the same model as the consumer:

```python
resolver = RootSearchResolver(
    datapackage_dir=output_root,
    asset_roots=[config.data_root],
    asset_overrides={
        "smartspim": asset_info.asset_path,
        "spim_template_to_ccf": config.template_to_ccf_dir,
    },
)
```

Then `DataPackage.missing_paths(...)` resolves each `{asset, path}` through the
resolver instead of doing `root / rel`.

This preserves the useful "fail before writing a dangling datapackage" behavior
without leaking producer-local mount paths into the datapackage.

## 8. Work items (once agreed)

- **Producer** (`manifest.py`): emit `external_assets` registry
  (role/name/id/uri[/checksum]); emit transforms (and any external histology) as
  `{asset, path}`; add `platform`; convert internal paths to
  `{asset: null, path}`; bump schema to 3.0.0. Drop `_rel_up`. Add tests
  asserting no `..`, `/data`, or absolute paths in serialized references.
- **Producer validation**: replace direct `root / rel` path checks with
  resolver-backed checks. Validate both external assets and local datapackage
  outputs (`asset: null`).
- **Datapackage-only regeneration**: support Code Ocean's immutable-results
  workflow by accepting a prior results asset mounted under `/data`, copying
  its `<mouse_id>/` output tree into the current `/results`, then rewriting
  only `datapackage.json` with schema 3 references.
- **GUI** (`datapackage_loader.py`): introduce `AssetResolver` +
  `RootSearchResolver`; resolve all `(asset, within)` through it; wire roots
  and overrides from env/settings; `AssetNotFound` surfaces registry reference
  descriptors. Tests for local (`asset: null`), external hit, override hit (by
  key and by name), and missing-asset error. No `/data`, no legacy heal in GUI
  core.
- **Error messages**: replace downstream ANTs/HDF5 errors with
  `AssetNotFound`/`ReferenceNotFound` messages that include asset key, role,
  name, id, uri, requested path, and searched roots/overrides.

## 9. Open decisions

1. ⟨DECIDE⟩ **Identity primary key** — `id` vs `uri` vs content hash; store which
   / how many? Lean: store `id` + `uri` + `name`/`role`; `uri` is the most
   platform-neutral and can become a fetch source; `checksum` is
   optional-but-ideal for repro. For v1, missing `id`/`uri` may be null, but
   `role` and `name` should be present.
2. ⟨DECIDE⟩ **Config channel** — env var, GUI settings field, or both.
3. ⟨DECIDE⟩ **Integrity now or later** — compute/store checksums in the producer
   from day one?
4. ⟨DECIDE⟩ **v1 resolver scope** — CO (`/data` root) + desktop local-root,
   deferring S3 fetch/verify?
5. ⟨DECIDE⟩ **Override format** — env JSON, GUI settings, CLI flag, or all of
   them? Recommendation: env + GUI settings; env wins.

## 10. Non-goals / rejected

- **Self-contained datapackage** (bundle everything): impractical — the
  SmartSPIM histology volumes are enormous and `spim_template_to_ccf` is shared,
  so bundling duplicates large data per subject. The datapackage is a reference
  index, not a container.
- **Platform detection in the GUI** (sniffing whether `/data` exists, etc.):
  forbidden — platform binding is deployment-supplied config, not inference.
