# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.4.0 (2026-08-18)

### Feat

- **discovery**: let the manifest pin a separate registration asset
- **manifest**: declare the CSV contract in one place
- drop the CCF-space registration pass-through
- **discover**: record and enforce manifest coverage
- datapackage 4.0.0 drops the pipeline volume for its sidecar
- emit the pipeline image's geometry as a sidecar (schema 3.2.0)
- resample histology onto one 30 um isotropic grid
- pass content-resolved raw session to extract_spikes
- add ephys launch stage
- resolve session dirs by name-walk, dropping the raw/sorted sibling assumption
- forward producer_record from stage_discover for fan-out provenance
- export the datapackage contract as a bundled JSON Schema
- add pack fan-in merge for pipeline outputs
- add find_ephys_stream_config for single-read ephys dispatch
- add --stage fan-out entry points for the CO pipeline
- skip failed-sorting probes instead of aborting datapackage
- gate GUI-unused QC outputs behind emit_qc flag (default off)
- add Code Ocean asset resolver (DISCOVER) for the preprocessing trigger
- **datapackage**: expose optional channel contact IDs
- **datapackage**: emit durable asset references
- **datapackage**: encode ephys geometry mapping

### Fix

- **discover**: refuse an empty manifest, and stop resolving records by accident
- **ephys**: match sorting output on the parsed collection, not a substring
- **discover**: match sortings by input recording, not by pinned name
- **resolver**: resolve raw ecephys from sorting provenance
- bring the sync histology path back in line with the async one
- resample the base/pipeline pair together, not independently
- rehydrate geometry as a 1x1x1 stub, not a full allocation
- the resample target was in um but image spacing is mm -- 1000x
- make the geometry sidecar say which space and which axes it means
- timing spans overlap under asyncio -- say so, and emit a timeline
- fall back to the mount name when asset name and metadata disagree
- NWB internals are not raw ecephys sessions
- resolve a pinned sorting to its published capture, and report gaps
- datapackage listed no probes on a qc=0 run
- resolve ephys fan-out sessions by content, not asset name
- namespace ephys stage output by fan-out unit
- point datapackage ephys at the probe folder + validate paths before writing

### Refactor

- make discover the up-front probe-viability gate

### Perf

- build the ANTs warp domains from headers, not written volumes
- time the histology steps, and stop shipping bytes nobody needs
- target 35 um voxels so histology selects the coarser level

## v0.3.0 (2026-07-08)

### BREAKING CHANGE

- datapackage schema bumped to 1.1.0. Consumers must
resolve transform paths against the datapackage directory.

### Feat

- nest probes by recording_id (schema 2.0.0)
- record transforms as paths relative to the manifest root

## v0.2.0 (2026-04-26)

### Feat

- **validation**: probe ecephys_clipped/structure.oebin pre-flight

## v0.1.3 (2026-03-24)

### Fix

- use aind-ephys-ibl-gui-conversion>=0.2.2

## v0.1.2 (2026-03-24)

### Fix

- limit OOM by setting number of ephys extract_continuous workers

## v0.1.1 (2026-03-24)

### Fix

- publish to pypi, use published aind-ephys-ibl-gui-conversion

## v0.1.0 (2026-03-23)

### Fix

- Add Manifest (#1)
