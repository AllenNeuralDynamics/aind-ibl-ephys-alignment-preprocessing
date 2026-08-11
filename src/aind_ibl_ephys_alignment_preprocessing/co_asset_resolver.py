"""Resolve which Code Ocean data assets to attach for a preprocessing run.

This is the *pre-pipeline* (trigger-side) counterpart to :mod:`discovery`, which
resolves the file structure *inside* an already-mounted SmartSPIM asset. Given a
mouse id, the user-pinned sortings, and the Neuroglancer JSON, this module
decides which Code Ocean data assets to mount: the raw ecephys recordings, the
SmartSPIM acquisition, and the sorted-recording assets.

Design (verified live against 791094 and 781370 on 2026-07-09):

* **Raw recordings** -- taken from the chosen sorting's ``provenance.data_assets``,
  which names the recording it was produced from. Falls back to a
  ``name:ecephys_<mouseid>`` search filtered to ``tag:raw`` when a sorting
  carries no usable provenance, and warns when it does.
* **SmartSPIM** -- taken from the Neuroglancer image-layer source, which *is* the
  acquisition the tracks were drawn on (unambiguous even when a mouse was imaged
  twice, or two Code Ocean assets share a stitched name).
* **Sortings** -- **pinned by the user** (never auto-picked "newest": on real data
  the human pinned an *older* ks4 sorting over a newer ``derived`` one). The
  pinned name is resolved to an asset id by fuzzy match against the
  ``tag:<mouseid>`` pool, because asset names drift (the recording time is
  sometimes dropped before ``_sorted_``). Duplicate captures of one sorting are
  collapsed by ``provenance.computation``, preferring the external/published one.

The pure functions here operate on :class:`CandidateAsset` lists; the Code Ocean
SDK adapter (the three ``search_data_assets`` queries) lives in the trigger
capsule, which depends on ``codeocean``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

# NOTE: models live here (not in types.py) on purpose. types.py imports ants/pandas
# eagerly; keeping these local lets a lightweight consumer (e.g. the trigger capsule)
# import the resolver without pulling heavy deps. If this grows, extract the module
# into its own small package that both the converter and the trigger depend on.


@dataclass(frozen=True)
class CandidateAsset:
    """A Code Ocean data asset returned by an asset search, used by the resolver.

    Parameters
    ----------
    id : str
        Code Ocean data asset id.
    name : str
        Asset name.
    tags : tuple[str, ...]
        Asset tags.
    computation : str | None
        ``provenance.computation`` -- the id of the computation that produced
        the asset. This is the true identity of a sorting: the same run is often
        captured more than once under different names/tags.
    external : bool
        ``source_bucket.external`` -- True for the published (e.g.
        ``aind-open-data``) capture, False for an internal Code Ocean result.
    created : int
        Asset creation time (epoch seconds), used only to order competing
        computations. 0 when unknown, which degrades the ordering to the id
        tiebreak rather than making it arbitrary.
    source_assets : tuple[str, ...]
        ``provenance.data_assets`` -- the ids this asset was produced from. For
        a sorting this names its raw recording outright, which is the only
        authoritative way to find it; see :func:`_resolve_raw`.
    """

    id: str
    name: str
    tags: tuple[str, ...] = ()
    computation: str | None = None
    external: bool = False
    created: int = 0
    source_assets: tuple[str, ...] = ()


@dataclass(frozen=True)
class AssetResolution:
    """Resolved Code Ocean assets to attach for one mouse's preprocessing run.

    Parameters
    ----------
    mouseid : str
        Mouse identifier.
    smartspim : CandidateAsset | None
        SmartSPIM acquisition asset (resolved from the Neuroglancer source).
    raw : dict[str, CandidateAsset]
        Raw ecephys recording assets, keyed by recording key
        (``<mouseid>_<date>_<time>``).
    sortings : dict[str, CandidateAsset]
        Chosen sorting asset per recording key.
    warnings : tuple[str, ...]
        Human-readable warnings/notes accumulated during resolution.
    unresolved : tuple[str, ...]
        Pinned ``sorted_recording`` names that produced no sorting -- malformed,
        or matching nothing. Structured rather than left to warning-string
        matching, because the caller has to *decide* on them: a pin that resolves
        to nothing silently removes a whole session from the run, and the run
        still succeeds. See :meth:`shrunk_by`.
    """

    mouseid: str
    smartspim: CandidateAsset | None
    raw: dict[str, CandidateAsset]
    sortings: dict[str, CandidateAsset]
    warnings: tuple[str, ...] = ()
    unresolved: tuple[str, ...] = ()

    def shrunk_by(self) -> list[str]:
        """Return every reason this resolution covers less than was asked for.

        Empty means the run will process exactly what the manifest pinned.

        Returns
        -------
        list[str]
            One human-readable line per dropped pin or missing raw recording.
        """
        reasons = [f"pinned sorting {pin!r} did not resolve" for pin in self.unresolved]
        reasons += [
            f"recording {key} has a sorting but no raw asset" for key in sorted(set(self.sortings) - set(self.raw))
        ]
        return reasons

    def data_assets(self) -> list[tuple[str, str]]:
        """Return ``(asset_id, mount)`` pairs for ``RunParams.data_assets``.

        Ordered SmartSPIM, then raw recordings, then sortings; mount defaults to
        the asset name.

        Returns
        -------
        list[tuple[str, str]]
            One ``(id, mount)`` pair per asset to attach.
        """
        pairs: list[tuple[str, str]] = []
        if self.smartspim is not None:
            pairs.append((self.smartspim.id, self.smartspim.name))
        for key in sorted(self.raw):
            pairs.append((self.raw[key].id, self.raw[key].name))
        for key in sorted(self.sortings):
            pairs.append((self.sortings[key].id, self.sortings[key].name))
        return pairs


# recording key: "<mouseid>_<date>_<time>" (the identity after the "ecephys_" prefix).
# The suffix is optional because a re-uploaded raw keeps the key and adds to it
# (`..._corrected`); anchoring the pattern hid those completely.
_RAW_RE = re.compile(r"^ecephys_(?P<key>\d+_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})(?P<suffix>_.+)?$")
# Suffixes that mark a *derived* product rather than another upload of the raw.
# Without this the relaxed pattern above would also accept `..._sorted_<ts>`.
_DERIVED_SUFFIX_RE = re.compile(r"(?:^|_)(sorted|preprocessed|curated|nwb|oversplit|ccg)(?:[_-]|$)", re.IGNORECASE)
# a well-formed pinned sorted_recording: ecephys_<key>_sorted_<sort_ts>
_PINNED_RE = re.compile(
    r"^ecephys_(?P<key>\d+_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_sorted_"
    r"(?P<ts>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$"
)
# extract the SmartSPIM acquisition name from an s3 image-source uri
_SMARTSPIM_RE = re.compile(r"(SmartSPIM_[^/]+?)/")


def _acquisition_from_uri(uri: str | None) -> str | None:
    """Return the ``SmartSPIM_...`` acquisition name from an image-source URI.

    Parameters
    ----------
    uri : str | None
        An s3 (or path-like) URI pointing inside a SmartSPIM acquisition.

    Returns
    -------
    str | None
        The acquisition name, or None if the URI is missing or has no match.
    """
    if not uri:
        return None
    match = _SMARTSPIM_RE.search(uri)
    return match.group(1) if match else None


def smartspim_acquisition_from_ng(ng_data: dict[str, Any]) -> str | None:
    """Extract the SmartSPIM acquisition name from Neuroglancer image layers.

    Uses :func:`aind_zarr_utils.neuroglancer.get_image_sources` (the same parser
    the in-capsule :mod:`discovery` module uses) and returns the single acquisition
    referenced by the image layers.

    Parameters
    ----------
    ng_data : dict
        Parsed Neuroglancer JSON.

    Returns
    -------
    str | None
        The acquisition name (e.g. ``"SmartSPIM_791094_..._stitched_..."``), or
        None if no image source is present.
    """
    # Lazy import: keeps the module import-light for consumers that only call resolve().
    from aind_zarr_utils.neuroglancer import get_image_sources

    sources = get_image_sources(ng_data, remove_zarr_protocol=True)
    acquisitions = {a for a in (_acquisition_from_uri(u) for u in sources.values()) if a}
    if not acquisitions:
        return None
    return min(acquisitions)


def parse_pinned(sorted_recording: str) -> tuple[str, str] | None:
    """Parse a pinned ``sorted_recording`` into ``(recording_key, sort_ts)``.

    Parameters
    ----------
    sorted_recording : str
        A pinned sorting name, e.g.
        ``ecephys_791094_2025-10-08_16-48-57_sorted_2026-04-24_13-43-00``.

    Returns
    -------
    tuple[str, str] | None
        ``(recording_key, sort_timestamp)`` or None if malformed.
    """
    match = _PINNED_RE.match(sorted_recording)
    return (match.group("key"), match.group("ts")) if match else None


def _prefer_external(assets: list[CandidateAsset]) -> CandidateAsset:
    """Pick the published (external) capture; deterministic tiebreak on id."""
    return min(assets, key=lambda a: (not a.external, a.id))


def _with_sibling_captures(hits: list[CandidateAsset], pool: list[CandidateAsset]) -> list[CandidateAsset]:
    """Widen *hits* to every capture of the same computation found in *pool*.

    A pin names one capture, but Code Ocean often captures a single computation
    twice under *different* timestamps, and the published copy can carry the
    **older** name -- sometimes with the recording time dropped from it entirely.
    Matching on the pinned name alone can therefore only ever find the
    unpublished copy, and "prefer external" then has nothing external to prefer.

    Widening first makes the published sibling a candidate. Captures of one
    computation are the same spikes -- a published copy differs only by the AIND
    ingestion pass (``metadata.nd.json`` + ``original_metadata/``).

    Parameters
    ----------
    hits : list[CandidateAsset]
        Assets matching the pinned name.
    pool : list[CandidateAsset]
        Assets to scan for siblings (the ``tag:<mouseid>`` pool).

    Returns
    -------
    list[CandidateAsset]
        *hits* plus every same-computation sibling, de-duplicated by id.
    """
    computations = {a.computation for a in hits if a.computation}
    if not computations:
        return list(hits)
    widened = {a.id: a for a in hits}
    for asset in pool:
        if asset.computation in computations:
            widened.setdefault(asset.id, asset)
    return list(widened.values())


def _newest_computation_group(assets: list[CandidateAsset]) -> list[CandidateAsset]:
    """Return the candidates belonging to the most recently created computation.

    Used only to break a genuinely ambiguous match (one pin, several distinct
    computations). Ordering by asset creation time rather than by the timestamp
    in the name, because that timestamp is exactly what proved unreliable.

    Parameters
    ----------
    assets : list[CandidateAsset]
        Candidates spanning one or more computations.

    Returns
    -------
    list[CandidateAsset]
        The group sharing the newest computation.
    """
    groups: dict[str | None, list[CandidateAsset]] = {}
    for asset in assets:
        groups.setdefault(asset.computation, []).append(asset)
    return max(
        groups.values(),
        key=lambda g: (max(a.created for a in g), min(a.id for a in g)),
    )


def sibling_captures(asset: CandidateAsset, pool: list[CandidateAsset]) -> list[CandidateAsset]:
    """Return other assets that are captures of the same computation.

    Parameters
    ----------
    asset : CandidateAsset
        The chosen asset.
    pool : list[CandidateAsset]
        Assets to scan for siblings.

    Returns
    -------
    list[CandidateAsset]
        Assets (other than *asset*) sharing its ``provenance.computation``.
    """
    if not asset.computation:
        return []
    return [a for a in pool if a.id != asset.id and a.computation == asset.computation]


def raw_key_of(name: str) -> str | None:
    """Return the recording key an asset name denotes, if it names a raw upload.

    ``None`` for anything that is not an ``ecephys_<key>`` asset, and for the
    derived products that share the prefix (``..._sorted_<ts>`` and friends),
    which the relaxed :data:`_RAW_RE` would otherwise accept.

    Parameters
    ----------
    name : str
        Data asset name.

    Returns
    -------
    str | None
        ``<mouseid>_<date>_<time>``, or None.
    """
    match = _RAW_RE.match(name)
    if match is None:
        return None
    suffix = match.group("suffix")
    if suffix and _DERIVED_SUFFIX_RE.search(suffix):
        return None
    return match.group("key")


def _resolve_raw(
    mouseid: str,
    recording_keys: set[str],
    raw_candidates: list[CandidateAsset],
    sortings: dict[str, CandidateAsset],
    pool: list[CandidateAsset],
) -> tuple[dict[str, CandidateAsset], list[str]]:
    """Match raw ecephys recording assets to the required recording keys.

    **Provenance first.** A sorting records the assets it was produced from, so
    the raw is named outright rather than guessed. Reconstructing it from the
    sorting's *name* instead was wrong in two ways that both fail quietly:

    - a re-uploaded raw carries a suffix (``..._corrected``) the old anchored
      pattern could never match, and
    - two assets can share one name, in which case the last one seen won.

    750107's 2025-01-28 session is both at once: two same-named uncorrected
    assets plus a ``_corrected`` one, and the pinned sorting was produced from
    the corrected raw. The resolver mounted an uncorrected one -- a *different
    recording of the same session* -- which happened to fail loudly in neo only
    because the stream layouts differed. Matching layouts would have extracted
    continuous data from the wrong recording with nothing to show for it.

    The corrected asset is also ``type=result``, so it is absent from the
    ``type=Dataset`` raw search entirely; provenance ids are therefore looked up
    in the *combined* pool, not just *raw_candidates*.

    Name matching stays as the fallback for sortings with no usable provenance,
    and says so in the warnings.

    Parameters
    ----------
    mouseid : str
        Mouse identifier, for warning text.
    recording_keys : set[str]
        Recording keys the run needs a raw asset for.
    raw_candidates : list[CandidateAsset]
        Result of the ``name:ecephys_<mouseid>`` search -- the fallback pool.
    sortings : dict[str, CandidateAsset]
        Resolved sortings, keyed by recording key; the provenance source.
    pool : list[CandidateAsset]
        Every candidate seen, searched by id when following provenance.

    Returns
    -------
    tuple[dict[str, CandidateAsset], list[str]]
        ``(raw_by_recording_key, warnings)``.
    """
    by_id = {asset.id: asset for asset in pool}
    by_name_key: dict[str, list[CandidateAsset]] = {}
    for asset in raw_candidates:
        key = raw_key_of(asset.name)
        if key is not None and "raw" in asset.tags:
            by_name_key.setdefault(key, []).append(asset)

    raw: dict[str, CandidateAsset] = {}
    warnings: list[str] = []
    for key in recording_keys:
        sorting = sortings.get(key)
        # Only accept a parent that actually names this recording: a sorting's
        # provenance can also list configs, and a re-sorting lists its
        # predecessor sorting (excluded by raw_key_of).
        parents = [
            by_id[aid]
            for aid in (sorting.source_assets if sorting else ())
            if aid in by_id and raw_key_of(by_id[aid].name) == key
        ]
        if len(parents) == 1:
            raw[key] = parents[0]
            continue
        if len(parents) > 1:
            chosen = _prefer_external(parents)
            warnings.append(
                f"{key}: sorting provenance names {len(parents)} raw assets "
                f"({', '.join(p.id for p in parents)}); attaching {chosen.name!r} ({chosen.id})"
            )
            raw[key] = chosen
            continue

        hits = by_name_key.get(key, [])
        if not hits:
            warnings.append(f"{key}: no raw recording asset found (name:ecephys_{mouseid})")
            continue
        chosen = _prefer_external(hits)
        if sorting is not None and sorting.source_assets:
            warnings.append(
                f"{key}: sorting {sorting.id} records provenance {list(sorting.source_assets)} but none of it "
                f"is a raw asset for this recording; FELL BACK to the name search and attached "
                f"{chosen.name!r} ({chosen.id}) -- verify this is the recording the sorting came from"
            )
        else:
            warnings.append(f"{key}: sorting has no provenance; raw resolved by name to {chosen.name!r} ({chosen.id})")
        if len(hits) > 1:
            warnings.append(
                f"{key}: AMBIGUOUS -- {len(hits)} raw assets share this recording key "
                f"({', '.join(f'{h.name!r} ({h.id})' for h in hits)}); attached {chosen.id}"
            )
        raw[key] = chosen
    return raw, warnings


def _resolve_smartspim(
    acquisition: str | None,
    smartspim_candidates: list[CandidateAsset],
) -> tuple[CandidateAsset | None, list[str]]:
    """Match the NG-referenced SmartSPIM acquisition to a Code Ocean asset."""
    if acquisition is None:
        return None, ["SmartSPIM: no image source in the Neuroglancer file"]
    exact = [a for a in smartspim_candidates if a.name == acquisition]
    if not exact:
        return None, [f"SmartSPIM: NG references {acquisition!r} but no matching Code Ocean asset found"]
    warnings: list[str] = []
    if len(exact) > 1:
        warnings.append(
            f"SmartSPIM: {len(exact)} assets match NG acquisition {acquisition!r} "
            f"(ids {', '.join(a.id for a in exact)}); picking external/first"
        )
    return _prefer_external(exact), warnings


def _resolve_pinned_sortings(
    mouseid: str,
    pinned: list[str],
    tagged_all: list[CandidateAsset],
) -> tuple[dict[str, CandidateAsset], list[str], list[str]]:
    """Resolve each pinned ``sorted_recording`` to an asset via fuzzy match.

    Returns
    -------
    tuple[dict[str, CandidateAsset], list[str], list[str]]
        ``(sortings_by_recording_key, warnings, unresolved_pins)``. A pin in
        *unresolved_pins* removed a whole session from the run; the caller
        decides whether that is fatal.
    """
    sortings: dict[str, CandidateAsset] = {}
    warnings: list[str] = []
    unresolved: list[str] = []
    for pin in dict.fromkeys(pinned):  # dedupe, preserve order
        parsed = parse_pinned(pin)
        if parsed is None:
            warnings.append(f"pinned sorted_recording {pin!r} is not a well-formed name; skipped")
            unresolved.append(pin)
            continue
        key, sort_ts = parsed
        date = key.split("_", 1)[1].rsplit("_", 1)[0]  # <mouseid>_<date>_<time> -> <date>
        # recording time is optional because sorting asset names sometimes drop it
        fuzzy = re.compile(
            rf"^ecephys_{re.escape(mouseid)}_{re.escape(date)}"
            rf"(?:_\d{{2}}-\d{{2}}-\d{{2}})?_sorted_{re.escape(sort_ts)}$"
        )
        hits = [a for a in tagged_all if fuzzy.match(a.name)]
        if not hits:
            warnings.append(f"{key}: pinned sorting (sort_ts {sort_ts}) not found among tag:{mouseid} assets")
            unresolved.append(pin)
            continue
        # The pinned name may only match the unpublished capture, so widen to
        # every capture of the same computation before choosing.
        candidates = _with_sibling_captures(hits, tagged_all)
        computations = {a.computation for a in candidates if a.computation}
        if len(computations) > 1:
            candidates = _newest_computation_group(candidates)
            warnings.append(
                f"{key}: AMBIGUOUS -- {len(hits)} matches from {len(computations)} different "
                f"computations {sorted(computations)}; taking the newest "
                f"({candidates[0].computation}), disambiguate the pin"
            )
        elif len(candidates) > 1:
            warnings.append(
                f"{key}: {len(candidates)} captures of one computation "
                f"({next(iter(computations), '?')}); picking the external/published capture (same spikes)"
            )
        chosen = _prefer_external(candidates)
        sortings[key] = chosen
        if all(chosen.id != a.id for a in hits):
            warnings.append(
                f"{key}: pinned name matched {hits[0].name!r} ({hits[0].id}) but attaching "
                f"the published capture of the same computation: {chosen.name!r} ({chosen.id})"
            )
        elif chosen.name != pin:
            warnings.append(f"{key}: pinned name {pin!r} matched asset {chosen.name!r} (fuzzy: name inconsistency)")
    return sortings, warnings, unresolved


def resolve(
    mouseid: str,
    pinned_sortings: list[str],
    raw_candidates: list[CandidateAsset],
    smartspim_candidates: list[CandidateAsset],
    tagged_all: list[CandidateAsset],
    ng_smartspim_acquisition: str | None,
) -> AssetResolution:
    """Resolve the Code Ocean assets to attach for one mouse's preprocessing run.

    Parameters
    ----------
    mouseid : str
        Mouse identifier.
    pinned_sortings : list[str]
        User-pinned ``sorted_recording`` names (one per recording, from the
        manifest / bridge config).
    raw_candidates : list[CandidateAsset]
        Result of ``name:ecephys_<mouseid>`` (state-ready).
    smartspim_candidates : list[CandidateAsset]
        Result of ``name:SmartSPIM_<mouseid>`` (NOT ``tag:<mouseid>`` -- processed
        SmartSPIM assets are sometimes tagged only ``[smartspim, processed]``).
    tagged_all : list[CandidateAsset]
        Result of ``tag:<mouseid>`` -- the pool for fuzzy sorting resolution.
    ng_smartspim_acquisition : str | None
        SmartSPIM acquisition name from :func:`smartspim_acquisition_from_ng`.

    Returns
    -------
    AssetResolution
        Resolved SmartSPIM, raw, and sorting assets plus warnings/notes.
    """
    sortings, warn_sort, unresolved = _resolve_pinned_sortings(mouseid, pinned_sortings, tagged_all)
    recording_keys = {parsed[0] for pin in pinned_sortings if (parsed := parse_pinned(pin)) is not None}
    # Provenance ids are followed against every candidate seen: a corrected raw
    # can be type=result and so absent from the type=Dataset raw search.
    pool = list({a.id: a for a in [*raw_candidates, *tagged_all, *smartspim_candidates]}.values())
    raw, warn_raw = _resolve_raw(mouseid, recording_keys, raw_candidates, sortings, pool)
    smartspim, warn_spim = _resolve_smartspim(ng_smartspim_acquisition, smartspim_candidates)

    warnings = [*warn_sort, *warn_raw, *warn_spim]
    for key, asset in sortings.items():
        siblings = sibling_captures(asset, tagged_all)
        if siblings:
            warnings.append(
                f"{key}: note -- {len(siblings)} other capture(s) of the same computation "
                f"({', '.join(s.id for s in siblings)}); same spikes, not attached"
            )

    return AssetResolution(
        mouseid=mouseid,
        smartspim=smartspim,
        raw=raw,
        sortings=sortings,
        warnings=tuple(warnings),
        unresolved=tuple(unresolved),
    )
