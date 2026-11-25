"""Discovery integration and deployment-specific inventory metadata access."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import boto3
import pandas as pd
from aind_ephys_ibl_gui_conversion.inventory import (
    INVENTORY_FILENAME,
    MAX_METADATA_BYTES,
    LocalMetadataStore,
    MetadataStore,
    build_inventory,
    inventory_slice,
    relative_path,
    write_inventory,
)
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError

from aind_ibl_ephys_alignment_preprocessing.ephys import (
    find_raw_session_dir,
    find_session_dir,
    find_sorted_session_dirs,
    read_sorted_input_recording,
)
from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, PipelineConfig

logger = logging.getLogger(__name__)


class S3MetadataStore:
    """Read an explicitly resolved S3 asset without listing signal chunks."""

    def __init__(self, uri: str, *, anonymous: bool = False) -> None:
        """Bind storage location and credential policy outside the inventory."""
        parsed = urlparse(uri)
        if parsed.scheme != "s3" or not parsed.netloc or parsed.query or parsed.fragment:
            raise ValueError("Expected an S3 asset root")
        self.bucket = parsed.netloc
        self.prefix = parsed.path.strip("/")
        self.client = boto3.client("s3", config=Config(signature_version=UNSIGNED) if anonymous else None)
        self._json: dict[str, Any] = {}

    def _key(self, path: str) -> str:
        relative = relative_path(path)
        return "/".join(p for p in (self.prefix, "" if relative == "." else relative) if p)

    def directories(self, path: str) -> list[str]:
        """List immediate child prefixes, with a bound on metadata traversal."""
        prefix = self._key(path).rstrip("/") + "/"
        names: list[str] = []
        for page in self.client.get_paginator("list_objects_v2").paginate(
            Bucket=self.bucket, Prefix=prefix, Delimiter="/"
        ):
            names.extend(p["Prefix"][len(prefix) :].rstrip("/") for p in page.get("CommonPrefixes", []))
            if len(names) > 10000:
                raise ValueError("Metadata directory limit exceeded")
        return sorted(names)

    def read_json(self, path: str) -> Any:
        """Fetch only bounded JSON metadata and cache repeated reads."""
        key = self._key(path)
        if key not in self._json:
            try:
                response = self.client.get_object(Bucket=self.bucket, Key=key)
            except ClientError as exc:
                if exc.response["Error"]["Code"] in {"NoSuchKey", "404"}:
                    raise FileNotFoundError(path) from None
                raise
            with response["Body"] as body:
                content = body.read(MAX_METADATA_BYTES + 1)
            if len(content) > MAX_METADATA_BYTES:
                raise ValueError("Metadata exceeds size limit")
            self._json[key] = json.loads(content)
        return self._json[key]

    def read_bytes(self, path: str, *, offset: int, size: int) -> bytes:
        """Read a bounded timestamp scalar through an explicit S3 byte range."""
        if offset < 0 or not 0 < size <= MAX_METADATA_BYTES:
            raise ValueError("Invalid bounded metadata read")
        response = self.client.get_object(
            Bucket=self.bucket, Key=self._key(path), Range=f"bytes={offset}-{offset + size - 1}"
        )
        with response["Body"] as body:
            content: bytes = body.read(size + 1)
        if len(content) != size:
            raise ValueError("Timestamp response did not match requested byte range")
        return content

    def read_chunk(self, path: str) -> bytes:
        """Read a size-limited compressed timestamp chunk."""
        response = self.client.get_object(Bucket=self.bucket, Key=self._key(path))
        with response["Body"] as body:
            content: bytes = body.read(MAX_METADATA_BYTES + 1)
        if len(content) > MAX_METADATA_BYTES:
            raise ValueError("Timestamp chunk exceeds size limit")
        return content


def _local_sources(config: PipelineConfig, rows: list[ManifestRow]) -> list[dict[str, Any]]:
    """Describe only manifest-related mounts when no launcher supplies identities."""
    sources: dict[str, dict[str, Any]] = {}
    sorted_dirs = find_sorted_session_dirs(config.data_root)
    for row in rows:
        sort_key = str(row.sorted_recording)
        folder = find_session_dir(config.data_root, sort_key) or sorted_dirs.get(row.recording_id)
        raw_name = (read_sorted_input_recording(folder) if folder is not None else None) or row.recording_id
        raw_folder = find_raw_session_dir(config.data_root, recording_id=raw_name)
        sources[sort_key] = dict(key=sort_key, name=sort_key, role="sorted", source_keys=[raw_name])
        sources[raw_name] = dict(key=raw_name, name=raw_name, role="raw")
        if folder is not None:
            sources[sort_key]["local_mount"] = folder.relative_to(config.data_root).as_posix()
        if raw_folder is not None:
            sources[raw_name]["local_mount"] = raw_folder.relative_to(config.data_root).as_posix()
        if row.surface_finding is not None:
            surface = str(row.surface_finding)
            sources.setdefault(surface, dict(key=surface, name=surface, role="raw", local_mount=surface))
    return list(sources.values())


def discover_inventory(
    config: PipelineConfig,
    rows: list[ManifestRow],
    sources: list[dict[str, Any]] | None = None,
    *,
    manifest_frame: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Observe and plan before viability filtering, retaining inaccessible sources."""
    stores: list[tuple[dict[str, Any], MetadataStore | None]] = []
    for source in sources if sources is not None else _local_sources(config, rows):
        store: MetadataStore | None = None
        try:
            mount = source.get("local_mount")
            if mount:
                folder = config.data_root / relative_path(mount)
                if folder.is_dir():
                    store = LocalMetadataStore(folder)
            if store is None and source.get("uri"):
                store = S3MetadataStore(source["uri"], anonymous=bool(source.get("anonymous", False)))
        except (ValueError, OSError):
            logger.warning("[discover] inventory source %s unavailable", source["key"])
        stores.append((source, store))
    inventory: dict[str, Any] = build_inventory(stores)
    write_inventory(inventory, config.results_root / INVENTORY_FILENAME)
    from aind_ephys_ibl_gui_conversion.conversion_plan import PLAN_FILENAME

    from aind_ibl_ephys_alignment_preprocessing.conversion_planning import build_discovery_plans

    plans = build_discovery_plans(
        inventory, manifest_frame if manifest_frame is not None else pd.read_csv(config.manifest_csv)
    )
    write_inventory(plans, config.results_root / PLAN_FILENAME)
    logger.info(
        "[discover] inventory: %d recorded streams, %d compressed recordings, "
        "%d analyzers, %d findings; metadata complete=%s",
        len(inventory["recorded_streams"]),
        len(inventory["recordings"]),
        len(inventory["sortings"]),
        len(inventory["findings"]),
        inventory["metadata_complete"],
    )
    return inventory


def write_unit_inventory(inventory: dict[str, Any], config_path: Path, plans: dict[str, Any] | None = None) -> None:
    """Place portable observations and saved decisions beside each unit config."""
    cfg = json.loads(config_path.read_text())
    write_inventory(
        inventory_slice(
            inventory, cfg.get("ephys_collection"), cfg.get("sorted_recording"), cfg.get("surface_finding")
        ),
        config_path.parent / INVENTORY_FILENAME,
    )
    cfg["inventory"] = INVENTORY_FILENAME
    from aind_ibl_ephys_alignment_preprocessing.conversion_planning import read_plan_registry, write_unit_plan

    if plans is None:
        plans = read_plan_registry(config_path.parent.parent, inventory)
    write_unit_plan(inventory, cfg, config_path.parent, plans)
    config_path.write_text(json.dumps(cfg, indent=2) + "\n")


def write_unit_inventories(inventory: dict[str, Any] | None, paths: list[Path]) -> None:
    """Write slices only when discovery performed electrophysiology inventory."""
    if inventory is not None:
        for path in paths:
            write_unit_inventory(inventory, path)


def carry_inventory(config: PipelineConfig) -> None:
    """Retain canonical observations and plans beside the final datapackage."""
    source = config.manifest_csv.parent / INVENTORY_FILENAME
    if not source.is_file():
        logger.info("[pack] no discovery inventory (legacy or standalone input)")
        return
    inventory = json.loads(source.read_text())
    frame = pd.read_csv(config.manifest_csv, dtype={"mouseid": "string"})
    mouse_id = str(frame["mouseid"].iat[0])
    write_inventory(inventory, config.results_root / mouse_id / INVENTORY_FILENAME)
    from aind_ephys_ibl_gui_conversion.conversion_plan import PLAN_FILENAME

    from aind_ibl_ephys_alignment_preprocessing.conversion_planning import read_plan_registry

    plans = read_plan_registry(source.parent, inventory)
    if plans is not None:
        write_inventory(plans, config.results_root / mouse_id / PLAN_FILENAME)
