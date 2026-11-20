"""Audit the final package after the normal discovery/fan-out/pack workflow."""

from __future__ import annotations

import argparse
import json
import os
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from aind_ephys_ibl_gui_conversion.conversion_plan import inventory_digest
from aind_ephys_ibl_gui_conversion.output_validation import audit_collection

OUTPUT_MANIFEST = "ephys_output_manifest.json"


def record_output_manifests(results_root: Path) -> int:
    """Record every collection file before capture, including optional metrics."""
    receipts = set(results_root.glob("*/*/*/*/ephys_extraction.json"))
    receipts.update(results_root.glob("*/*/*/ephys_extraction.json"))
    for path in receipts:
        receipt = json.loads(path.read_text())
        if receipt["status"] != "complete":
            raise ValueError(f"Cannot record incomplete extraction: {path}")
        root = path.parent
        files = {
            p.relative_to(root).as_posix(): p.stat().st_size
            for p in root.rglob("*")
            if p.is_file() and p != root / OUTPUT_MANIFEST
        }
        manifest = dict(inventory_digest=receipt["plan"]["inventory_digest"], files=files)
        (root / OUTPUT_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    return len(receipts)


def _check_output_manifest(folder: Path, plan: dict[str, Any]) -> int:
    """Reject file loss or truncation across capture, staging, and pack."""
    manifest = json.loads((folder / OUTPUT_MANIFEST).read_text())
    if manifest["inventory_digest"] != plan["inventory_digest"]:
        raise ValueError(f"Output manifest belongs to another discovery plan: {folder}")
    for name, size in manifest["files"].items():
        path = _local_path(folder, dict(asset=None, path=name))
        if not path.is_file() or path.stat().st_size != size:
            raise ValueError(f"Packaged file changed size: {path}")
    return len(manifest["files"])


def _local_path(root: Path, reference: dict[str, Any]) -> Path:
    """Resolve a package-local reference while enforcing the package boundary."""
    if reference.get("asset") is not None:
        raise ValueError("Expected a package-local reference")
    relative = Path(reference["path"])
    path = (root / relative).resolve()
    if relative.is_absolute() or not path.is_relative_to(root.resolve()) or not path.exists():
        raise ValueError(f"Invalid or missing package reference: {reference}")
    return path


def _check_references(value: Any, root: Path, external_assets: dict[str, Any]) -> int:
    """Check local references and retain explicit external-asset membership."""
    if isinstance(value, list):
        return sum(_check_references(v, root, external_assets) for v in value)
    if not isinstance(value, dict):
        return 0
    if "asset" in value and "path" in value:
        if value["asset"] is None:
            _local_path(root, value)
        elif value["asset"] not in external_assets:
            raise ValueError(f"Unknown external asset reference: {value}")
        return 1
    return sum(_check_references(v, root, external_assets) for v in value.values())


def _check_qc(folder: Path, rows: int, blocks: int) -> None:
    """Require convention-loaded QC products and matching channel/time axes."""
    for band in ("AP", "LF"):
        tags = ("", "Main") if blocks > 1 else ("",)
        for tag in tags:
            stem = f"_iblqc_ephysTimeRms{band}{tag}"
            rms = np.load(folder / f"{stem}.rms.npy", mmap_mode="r")
            times = np.load(folder / f"{stem}.timestamps.npy", mmap_mode="r")
            assert rms.ndim == 2 and rms.shape[1] == rows, stem
            assert times.shape == (rms.shape[0],), stem
    power = np.load(folder / "_iblqc_ephysSpectralDensityLF.power.npy", mmap_mode="r")
    freqs = np.load(folder / "_iblqc_ephysSpectralDensityLF.freqs.npy", mmap_mode="r")
    assert power.shape == (len(freqs), rows)
    for band in ("delta", "theta", "alpha", "beta", "gamma"):
        for suffix in ("mean_corr", "coherency"):
            path = folder / "band_corr" / f"{band}_{suffix}.npy"
            assert np.load(path, mmap_mode="r").shape == (rows, rows), path


def audit_package(root: Path) -> dict[str, Any]:
    """Verify registry handoff, packaged arrays, QC, and histology references."""
    package = json.loads((root / "datapackage.json").read_text())
    inventory = json.loads((root / "ephys_inventory.json").read_text())
    registry = json.loads((root / "ephys_conversion_plan.json").read_text())
    if registry["inventory_digest"] != inventory_digest(inventory):
        raise ValueError("Packaged plan registry does not match discovery inventory")
    references = _check_references(package, root, package.get("external_assets", {}))
    collections = {}
    for recording, entries in package["probes"].items():
        for collection, entry in entries.items():
            unit = f"{recording}__{collection}"
            folder = _local_path(root, entry["ephys"])
            receipt = json.loads((folder / "ephys_extraction.json").read_text())
            plan = registry["units"][unit]
            if receipt["plan"] != plan:
                raise ValueError(f"Packaged execution used a different discovery plan: {unit}")
            if not receipt["execution_enabled"] or receipt["status"] != "complete":
                raise ValueError(f"Incomplete planned extraction: {unit}")
            if not entry["xyz_picks"]:
                raise ValueError(f"No histology picks for {unit}")
            files = _check_output_manifest(folder, plan)
            result = audit_collection(folder, plan)
            _check_qc(folder, result["channel_rows"], len(receipt["blocks"]))
            result["histology_pick_count"] = len(entry["xyz_picks"])
            result["collection_files_verified"] = files
            collections[unit] = result
    if not collections or collections.keys() != registry["units"].keys():
        raise ValueError("Packaged collections do not match the requested discovery units")
    return dict(
        status="passed",
        schema_version=package["schema_version"],
        global_inventory_digest_verified=True,
        discovery_units_verified=True,
        references_checked=references,
        external_assets=list(package.get("external_assets", {})),
        collections=collections,
        provenance={key: os.environ.get(key) for key in ("IBL_CONVERSION_COMMIT", "IBL_PREPROCESSING_COMMIT")},
    )


def main() -> None:
    """Audit one packed mouse tree and retain a report on success or failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=Path("/results"))
    parser.add_argument("--record-outputs", action="store_true")
    args = parser.parse_args()
    if args.record_outputs:
        count = record_output_manifests(args.results_root)
        print(f"Recorded output manifests for {count} completed ephys collections", flush=True)
        return
    packages = list(args.results_root.glob("*/datapackage.json"))
    if len(packages) != 1:
        raise ValueError(f"Expected one final datapackage, found {packages}")
    root = packages[0].parent
    report_path = root / "ephys_workflow_validation.json"
    try:
        report = audit_package(root)
    except Exception:
        report_path.write_text(json.dumps(dict(status="failed", traceback=traceback.format_exc()), indent=2))
        raise
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
