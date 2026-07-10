"""Resolve which Code Ocean data assets to attach for a preprocessing run.

This is the *pre-pipeline* (trigger-side) counterpart to :mod:`discovery`, which
resolves the file structure *inside* an already-mounted SmartSPIM asset. Given a
mouse id, the user-pinned sortings, and the Neuroglancer JSON, this module
decides which Code Ocean data assets to mount: the raw ecephys recordings, the
SmartSPIM acquisition, and the sorted-recording assets.

Design (verified live against 791094 and 781370 on 2026-07-09):

* **Raw recordings** -- derived from each pinned sorting's recording name; the
  raw asset is found with a ``name:ecephys_<mouseid>`` search filtered to
  ``tag:raw``.
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
    """

    id: str
    name: str
    tags: tuple[str, ...] = ()
    computation: str | None = None
    external: bool = False


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
    """

    mouseid: str
    smartspim: CandidateAsset | None
    raw: dict[str, CandidateAsset]
    sortings: dict[str, CandidateAsset]
    warnings: tuple[str, ...] = ()

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


# recording key: "<mouseid>_<date>_<time>" (the identity after the "ecephys_" prefix)
_RAW_RE = re.compile(r"^ecephys_(?P<key>\d+_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$")
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


def smartspim_acquisition_from_ng(ng_data: dict) -> str | None:  # type: ignore[type-arg]
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
    return sorted(acquisitions)[0]


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
    return sorted(assets, key=lambda a: (not a.external, a.id))[0]


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


def _resolve_raw(
    mouseid: str,
    recording_keys: set[str],
    raw_candidates: list[CandidateAsset],
) -> tuple[dict[str, CandidateAsset], list[str]]:
    """Match raw ecephys recording assets to the required recording keys."""
    by_key: dict[str, CandidateAsset] = {}
    for asset in raw_candidates:
        match = _RAW_RE.match(asset.name)
        if match and "raw" in asset.tags:
            by_key[match.group("key")] = asset
    raw: dict[str, CandidateAsset] = {}
    warnings: list[str] = []
    for key in recording_keys:
        if key in by_key:
            raw[key] = by_key[key]
        else:
            warnings.append(f"{key}: no raw recording asset found (name:ecephys_{mouseid})")
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
) -> tuple[dict[str, CandidateAsset], list[str]]:
    """Resolve each pinned ``sorted_recording`` to an asset via fuzzy match."""
    sortings: dict[str, CandidateAsset] = {}
    warnings: list[str] = []
    for pin in dict.fromkeys(pinned):  # dedupe, preserve order
        parsed = parse_pinned(pin)
        if parsed is None:
            warnings.append(f"pinned sorted_recording {pin!r} is not a well-formed name; skipped")
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
            continue
        computations = {a.computation for a in hits if a.computation}
        if len(computations) > 1:
            warnings.append(
                f"{key}: AMBIGUOUS -- {len(hits)} matches from {len(computations)} different "
                f"computations {sorted(computations)}; picking external/first, disambiguate the pin"
            )
        elif len(hits) > 1:
            warnings.append(
                f"{key}: {len(hits)} captures of one computation "
                f"({next(iter(computations), '?')}); picking the external/published capture (same spikes)"
            )
        chosen = _prefer_external(hits)
        sortings[key] = chosen
        if chosen.name != pin:
            warnings.append(f"{key}: pinned name {pin!r} matched asset {chosen.name!r} (fuzzy: name inconsistency)")
    return sortings, warnings


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
    sortings, warn_sort = _resolve_pinned_sortings(mouseid, pinned_sortings, tagged_all)
    recording_keys = {parsed[0] for pin in pinned_sortings if (parsed := parse_pinned(pin)) is not None}
    raw, warn_raw = _resolve_raw(mouseid, recording_keys, raw_candidates)
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
    )
