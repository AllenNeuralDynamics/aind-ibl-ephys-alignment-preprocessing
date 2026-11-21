# Build a manifest from a mapping

Most of a manifest repeats or follows from two choices: which sorted recording each session uses,
and which stream each annotated track was recorded on. `aind-ibl-build-manifest` takes only those,
in a TOML file, and derives the rest. Track names are never parsed, because annotators follow no
naming convention.

```toml
mouseid = "837772"
annotation = "837772-probe-annotations-CK-2026-09-29.json"  # relative to this file

[[session]]
sorted_recording = "ecephys_837772_2026-06-27_14-08-21_sorted_2026-08-27_18-59-49"
# surface_finding = "..."                                     # optional

[session.tracks]
"732_MD_1" = "46101-1"   # one shank of a split quad-base probe
"733_BLA" = "46810"      # a single-stream probe
"740_PL" = { stream = "45883", shank = 1 }  # one shank of a single-stream multi-shank probe; the first shank is 1
```

```bash
aind-ibl-build-manifest 837772.toml -o manifest_837772.csv [--data-root /data]
```

A stream is the collection parsed from the sort's `postprocessed/` analyzer names. They are listed
from a copy mounted under `--data-root`, or else unsigned from that sort's prefix in
`aind-open-data`. How the columns are derived:

- `probe_file` is the annotation file's name without its extension.
- **Split quad-base probe:** a stream `<probe>-<N>`, when the sort also holds another
  `<probe>-<M>`, gets `logical_probe = <probe>`, `histology_shank = N - 1` and `ephys_shank = 0`.
- **Any other stream** is its own `logical_probe`. Its shanks stay blank, or both take the given
  shank's index.

The build writes nothing if any of these fail, and it reports every failure:
- every stream in each sort has a track, and every assigned stream is in its sort;
- every track is a layer of the annotation file with at least two points;
- no track is assigned twice;
- a split stream has exactly one track, and the tracks sharing any other stream have distinct
  shanks.

Annotated layers left out of the mapping are printed as notes, so a missed track gets noticed.
