"""Pre-flight validation for the preprocessing pipeline.

This module provides comprehensive checks to detect missing files, incorrect
configurations, and resource constraints *before* the pipeline starts
processing.

Example
-------
>>> from aind_ibl_ephys_alignment_preprocessing.validation import PipelineValidator
>>> validator = PipelineValidator(config)
>>> results = validator.validate_all()
>>> if validator.has_errors(results):
...     validator.print_summary(results)
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import pandas as pd

from aind_ibl_ephys_alignment_preprocessing.types import ManifestRow, ReferencePaths

if TYPE_CHECKING:
    from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig


def _norm_shank_series(s: pd.Series[Any]) -> pd.Series[Any]:
    """Treat NaN/None as a sentinel so duplicates behave predictably."""
    return s.fillna(-1).astype(int)


@dataclass(frozen=True)
class ValidationResult:
    """Result of a single validation check.

    Attributes
    ----------
    passed : bool
        Whether the validation check passed.
    category : str
        Category of the check (e.g., ``"CLI Args"``, ``"Manifest CSV"``).
    item : str
        Specific item being checked.
    message : str
        Human-readable message about the result.
    severity : str
        ``"error"`` (must fix), ``"warning"`` (should fix), or ``"info"``.
    """

    passed: bool
    category: str
    item: str
    message: str
    severity: str = "error"


class PipelineValidator:
    """Comprehensive pre-flight validator.

    Parameters
    ----------
    config : PipelineConfig
        Fully-resolved pipeline configuration.
    skip_resource_checks : bool
        If *True*, skip resource-intensive checks (disk space, RAM).
    """

    def __init__(
        self,
        config: PipelineConfig,
        skip_resource_checks: bool = False,
    ) -> None:
        self.config = config
        self.skip_resource_checks = skip_resource_checks
        self.ref_paths = ReferencePaths.from_config(config)
        self.results: list[ValidationResult] = []

    def validate_all(self) -> list[ValidationResult]:
        """Run all validation checks in sequence.

        Returns
        -------
        list[ValidationResult]
            All validation results.
        """
        self.results = []
        self.validate_inputs()
        self.validate_manifest_structure()
        self.validate_reference_data()
        self.validate_neuroglancer_and_asset()
        self.validate_per_probe_files()
        self.validate_output_access()
        if not self.skip_resource_checks:
            self.validate_resources()
        return self.results

    # -- helpers ---------------------------------------------------------------

    def _add_result(
        self,
        passed: bool,
        category: str,
        item: str,
        message: str,
        severity: str = "error",
    ) -> None:
        self.results.append(
            ValidationResult(passed=passed, category=category, item=item, message=message, severity=severity)
        )

    def _add_unique_violation_results(
        self,
        df: pd.DataFrame,
        subset: list[str],
        label: str,
        category: str,
    ) -> None:
        dups_mask = df.duplicated(subset=subset, keep=False)
        if not dups_mask.any():
            self._add_result(True, category, f"{label}_uniqueness", f"Unique by {tuple(subset)}", severity="info")
            return
        offenders = (
            df.loc[dups_mask, subset].value_counts().reset_index(name="count").sort_values("count", ascending=False)
        )
        top_preview = offenders.head(10).to_string(index=False)
        self._add_result(
            False,
            category,
            f"{label}_uniqueness",
            f"Found non-unique keys for {label} (expected unique by {tuple(subset)}).\nTop offenders:\n{top_preview}",
        )

    # -- Category 1: Input Paths -----------------------------------------------

    def validate_inputs(self) -> None:
        """Validate that input file paths are non-empty and sane."""
        category = "Inputs"

        ng = self.config.neuroglancer_file
        if not str(ng):
            self._add_result(False, category, "neuroglancer", "Neuroglancer file path is empty")
        else:
            self._add_result(True, category, "neuroglancer", f"Neuroglancer path provided: {ng}", severity="info")

        manifest = self.config.manifest_csv
        if not str(manifest):
            self._add_result(False, category, "manifest_csv", "Manifest CSV path is empty")
        else:
            self._add_result(True, category, "manifest_csv", f"Manifest path provided: {manifest}", severity="info")

        if self.config.skip_ephys:
            self._add_result(True, category, "skip_ephys", "Ephys processing will be skipped", severity="info")

    # -- Category 2: Manifest CSV -----------------------------------------------

    def _validate_mouseid_consistency(self, df: pd.DataFrame, category: str) -> None:
        """Check that only a single mouse ID is present in the manifest."""
        if "mouseid" not in df.columns:
            return
        unique_mice = df["mouseid"].unique()
        if len(unique_mice) > 1:
            self._add_result(
                False,
                category,
                "mouseid_consistency",
                f"Multiple mouse IDs found: {', '.join(map(str, unique_mice))}. "
                "Pipeline currently supports single mouse per run.",
            )
        elif len(unique_mice) == 0 or pd.isna(unique_mice[0]):
            self._add_result(False, category, "mouseid_consistency", "Mouse ID is missing or NA in all rows")
        else:
            self._add_result(
                True, category, "mouseid_consistency", f"Single mouse ID: {unique_mice[0]}", severity="info"
            )

    def _validate_null_columns(self, df: pd.DataFrame, required_cols: list[str], category: str) -> None:
        """Warn about null values in required columns."""
        if len(df) == 0:
            return
        for col in required_cols:
            if col in df.columns:
                null_count = df[col].isna().sum()
                if null_count > 0:
                    self._add_result(
                        False,
                        category,
                        f"{col}_nulls",
                        f"Column '{col}' has {null_count} null/empty values",
                        severity="warning",
                    )

    def _validate_uniqueness_constraints(self, df: pd.DataFrame, category: str) -> None:
        """Check uniqueness constraints and recording-id multiplicity."""
        if "recording_id" not in df.columns and "sorted_recording" in df.columns:
            df = df.assign(
                recording_id=df["sorted_recording"].astype(str).str.split("_sorted", n=1, expand=True)[0],
            )

        shank_col = "probe_shank"
        if shank_col not in df.columns:
            df[shank_col] = None
        df["_probe_shank_norm"] = _norm_shank_series(df[shank_col])

        if {"mouseid", "probe_id"}.issubset(df.columns):
            self._add_unique_violation_results(
                df,
                subset=["mouseid", "probe_id", "_probe_shank_norm"],
                label="bregma_xyz (mouseid, probe_id, probe_shank)",
                category=category,
            )

        if {"recording_id", "probe_name"}.issubset(df.columns):
            self._add_unique_violation_results(
                df,
                subset=["recording_id", "probe_name", "_probe_shank_norm"],
                label="GUI (recording_id, probe_name, probe_shank)",
                category=category,
            )

        if "recording_id" in df.columns:
            dmask = df.duplicated(subset=["recording_id"], keep=False)
            if dmask.any():
                counts = (
                    df.loc[dmask, ["recording_id"]]
                    .value_counts()
                    .reset_index(name="count")
                    .sort_values("count", ascending=False)
                    .head(10)
                    .to_string(index=False)
                )
                self._add_result(
                    True,
                    category,
                    "recording_id_multiplicity",
                    "Multiple rows reference the same recording_id "
                    "(this is allowed; ephys is deduped at run time).\n"
                    f"Top counts:\n{counts}",
                    severity="warning",
                )

        if "_probe_shank_norm" in df.columns:
            df.drop(columns=["_probe_shank_norm"], inplace=True)

    def validate_manifest_structure(self) -> None:
        """Validate manifest CSV file structure and required columns."""
        category = "Manifest CSV"

        if not self.config.manifest_csv.exists():
            self._add_result(False, category, "file_exists", f"Manifest CSV not found: {self.config.manifest_csv}")
            return
        self._add_result(
            True, category, "file_exists", f"Manifest CSV exists: {self.config.manifest_csv}", severity="info"
        )

        try:
            df = pd.read_csv(self.config.manifest_csv)
        except Exception as e:
            self._add_result(False, category, "readable", f"Failed to read manifest CSV: {e}")
            return
        self._add_result(True, category, "readable", f"Manifest CSV readable ({len(df)} rows)", severity="info")

        required_cols = ["mouseid", "sorted_recording", "probe_file", "probe_id", "probe_name"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            self._add_result(
                False, category, "required_columns", f"Missing required columns: {', '.join(missing_cols)}"
            )
        else:
            self._add_result(True, category, "required_columns", "All required columns present", severity="info")

        self._validate_mouseid_consistency(df, category)
        self._validate_null_columns(df, required_cols, category)
        self._validate_uniqueness_constraints(df, category)

    # -- Category 3: Reference Data --------------------------------------------

    def validate_reference_data(self) -> None:
        """Validate that all required reference data files exist."""
        category = "Reference Data"
        checks = [
            ("template_25", self.ref_paths.template_25, "Template"),
            ("ccf_25", self.ref_paths.ccf_25, "CCF template"),
            ("ccf_labels", self.ref_paths.ccf_labels_lateralized_25, "CCF labels"),
            ("ibl_atlas", self.ref_paths.ibl_atlas_histology_path, "IBL atlas directory"),
        ]
        for item, path, label in checks:
            if path.exists():
                self._add_result(True, category, item, f"{label} found: {path}", severity="info")
            else:
                self._add_result(False, category, item, f"{label} not found: {path}")

    # -- Category 4: SmartSPIM Asset Discovery ----------------------------------

    def validate_neuroglancer_and_asset(self) -> None:
        """Validate Neuroglancer file and discover SmartSPIM asset structure."""
        category = "Asset Discovery"

        if not self.config.neuroglancer_file.exists():
            self._add_result(
                False, category, "neuroglancer_file", f"Neuroglancer file not found: {self.config.neuroglancer_file}"
            )
            return
        self._add_result(
            True,
            category,
            "neuroglancer_file",
            f"Neuroglancer file exists: {self.config.neuroglancer_file}",
            severity="info",
        )

        try:
            from aind_s3_cache.json_utils import get_json
            from aind_zarr_utils.neuroglancer import get_image_sources

            ng_data = get_json(str(self.config.neuroglancer_file))
            sources = get_image_sources(ng_data, remove_zarr_protocol=True)

            if not sources:
                self._add_result(False, category, "image_sources", "No image sources found in Neuroglancer file")
                return
            self._add_result(
                True,
                category,
                "image_sources",
                f"Found {len(sources)} image source(s) in Neuroglancer file",
                severity="info",
            )

            from aind_s3_cache.uri_utils import as_pathlike
            from aind_zarr_utils.pipeline_transformed import _asset_from_zarr_pathlike

            a_zarr_uri = next(iter(sources.values()), None)
            if a_zarr_uri:
                _, _, a_zarr_pathlike = as_pathlike(a_zarr_uri)
                asset_pathlike = _asset_from_zarr_pathlike(a_zarr_pathlike)
                asset_path = self.config.data_root / asset_pathlike

                if asset_path.exists():
                    self._add_result(True, category, "asset_path", f"Asset path exists: {asset_path}", severity="info")

                    zarr_path = asset_path / "image_tile_fusing" / "OMEZarr"
                    if zarr_path.exists():
                        self._add_result(
                            True, category, "omezarr_dir", f"OME-Zarr directory exists: {zarr_path}", severity="info"
                        )
                        zarr_files = list(zarr_path.glob("*.zarr"))
                        if zarr_files:
                            self._add_result(
                                True,
                                category,
                                "zarr_channels",
                                f"Found {len(zarr_files)} .zarr channel(s)",
                                severity="info",
                            )
                        else:
                            self._add_result(False, category, "zarr_channels", f"No .zarr files found in {zarr_path}")
                    else:
                        self._add_result(False, category, "omezarr_dir", f"OME-Zarr directory not found: {zarr_path}")

                    reg_dir = asset_path / "image_atlas_alignment"
                    if reg_dir.exists():
                        self._add_result(
                            True,
                            category,
                            "registration_dir",
                            f"Registration directory exists: {reg_dir}",
                            severity="info",
                        )
                        ccf_files = list(reg_dir.glob("*/moved_ls_to_ccf.nii.gz"))
                        if ccf_files:
                            self._add_result(
                                True,
                                category,
                                "precomputed_registration",
                                f"Found precomputed CCF registration: {ccf_files[0]}",
                                severity="info",
                            )
                        else:
                            self._add_result(
                                False,
                                category,
                                "precomputed_registration",
                                f"Precomputed registration (moved_ls_to_ccf.nii.gz) not found in {reg_dir}",
                            )

                        affine_files = list(reg_dir.glob("*/*_0GenericAffine.mat"))
                        warp_files = list(reg_dir.glob("*/*_1InverseWarp.nii.gz"))
                        if affine_files and warp_files:
                            self._add_result(
                                True,
                                category,
                                "transform_files",
                                f"Found transform files (affine: {len(affine_files)}, warp: {len(warp_files)})",
                                severity="info",
                            )
                        else:
                            missing = []
                            if not affine_files:
                                missing.append("affine (*_0GenericAffine.mat)")
                            if not warp_files:
                                missing.append("warp (*_1InverseWarp.nii.gz)")
                            self._add_result(
                                False, category, "transform_files", f"Missing transform files: {', '.join(missing)}"
                            )
                    else:
                        self._add_result(
                            False, category, "registration_dir", f"Registration directory not found: {reg_dir}"
                        )
                else:
                    self._add_result(False, category, "asset_path", f"Asset path not found: {asset_path}")
        except Exception as e:
            self._add_result(
                False, category, "asset_discovery", f"Failed to discover asset structure: {e}", severity="warning"
            )

    # -- Category 5: Per-probe Files -------------------------------------------

    def validate_per_probe_files(self) -> None:
        """Validate that per-probe annotation files and ephys data exist."""
        category = "Per-Probe Files"

        if not self.config.manifest_csv.exists():
            return

        try:
            df = pd.read_csv(self.config.manifest_csv)
        except Exception:
            return

        for idx, row in df.iterrows():
            try:
                mr = ManifestRow.from_series(row)
            except Exception as e:
                self._add_result(
                    False, category, f"row_{idx}", f"Failed to parse manifest row {idx}: {e}", severity="warning"
                )
                continue

            ext = "json" if mr.annotation_format == "json" else None
            if ext:
                pattern = f"*/{mr.probe_file}.{ext}"
                matches = list(self.config.data_root.glob(pattern))
                if matches:
                    self._add_result(
                        True,
                        category,
                        f"{mr.probe_id}_annotation",
                        f"Annotation found for {mr.probe_id}: {matches[0].relative_to(self.config.data_root)}",
                        severity="info",
                    )
                else:
                    self._add_result(
                        False,
                        category,
                        f"{mr.probe_id}_annotation",
                        f"Annotation not found for {mr.probe_id} (pattern: {pattern})",
                    )

            if not self.config.skip_ephys:
                recording_folder = self.config.data_root / mr.sorted_recording
                if recording_folder.exists():
                    self._add_result(
                        True,
                        category,
                        f"{mr.probe_id}_ephys",
                        f"Ephys folder exists for {mr.sorted_recording}",
                        severity="info",
                    )
                else:
                    self._add_result(
                        False, category, f"{mr.probe_id}_ephys", f"Ephys folder not found: {recording_folder}"
                    )

            if mr.surface_finding is not None:
                surface_path = self.config.data_root / mr.surface_finding
                if surface_path.exists():
                    self._add_result(
                        True,
                        category,
                        f"{mr.probe_id}_surface_finding",
                        f"Surface finding file exists: {surface_path}",
                        severity="info",
                    )
                else:
                    self._add_result(
                        False,
                        category,
                        f"{mr.probe_id}_surface_finding",
                        f"Surface finding file not found: {surface_path}",
                        severity="warning",
                    )

    # -- Category 6: Output Directory Access ------------------------------------

    def validate_output_access(self) -> None:
        """Validate that output directory is writable and has sufficient space."""
        category = "Output Access"
        results_dir = self.config.results_root

        if not results_dir.exists():
            self._add_result(False, category, "results_dir_exists", f"Results directory does not exist: {results_dir}")
            return

        if not os.access(results_dir, os.W_OK):
            self._add_result(False, category, "results_writable", f"Results directory is not writable: {results_dir}")
        else:
            self._add_result(
                True, category, "results_writable", f"Results directory is writable: {results_dir}", severity="info"
            )

        try:
            stat = shutil.disk_usage(results_dir)
            free_gb = stat.free / (1024**3)
            if free_gb < 10:
                self._add_result(
                    False,
                    category,
                    "disk_space",
                    f"Insufficient disk space: {free_gb:.1f} GB available (need at least 10 GB)",
                )
            elif free_gb < 50:
                self._add_result(
                    True,
                    category,
                    "disk_space",
                    f"Low disk space warning: {free_gb:.1f} GB available (recommend at least 50 GB)",
                    severity="warning",
                )
            else:
                self._add_result(
                    True, category, "disk_space", f"Sufficient disk space: {free_gb:.1f} GB available", severity="info"
                )
        except Exception as e:
            self._add_result(False, category, "disk_space", f"Failed to check disk space: {e}", severity="warning")

    # -- Category 7: Resource Checks -------------------------------------------

    def validate_resources(self) -> None:
        """Validate system resources (RAM, CPU, network mounts)."""
        category = "System Resources"

        try:
            with open("/proc/meminfo") as f:
                meminfo = f.read()
            for line in meminfo.split("\n"):
                if line.startswith("MemAvailable:"):
                    mem_kb = int(line.split()[1])
                    mem_gb = mem_kb / (1024**2)
                    if mem_gb < 8:
                        self._add_result(
                            False, category, "ram", f"Insufficient RAM: {mem_gb:.1f} GB available (need at least 8 GB)"
                        )
                    elif mem_gb < 16:
                        self._add_result(
                            True,
                            category,
                            "ram",
                            f"Low RAM warning: {mem_gb:.1f} GB available (recommend at least 16 GB)",
                            severity="warning",
                        )
                    else:
                        self._add_result(
                            True, category, "ram", f"Sufficient RAM: {mem_gb:.1f} GB available", severity="info"
                        )
                    break
        except Exception as e:
            self._add_result(False, category, "ram", f"Failed to check RAM: {e}", severity="warning")

        try:
            result = subprocess.run(
                ["df", "-T", str(self.config.data_root)],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                if len(lines) > 1:
                    fs_type = lines[1].split()[1]
                    if fs_type in ("nfs", "nfs4", "cifs", "fuse", "fuse.s3fs"):
                        self._add_result(
                            True,
                            category,
                            "network_mount",
                            f"Data directory is on network storage ({fs_type}). "
                            "Consider using local caching for better performance.",
                            severity="warning",
                        )
                    else:
                        self._add_result(
                            True,
                            category,
                            "network_mount",
                            f"Data directory is on local storage ({fs_type})",
                            severity="info",
                        )
        except Exception as e:
            self._add_result(
                False, category, "network_mount", f"Failed to check filesystem type: {e}", severity="warning"
            )

    # -- Output Methods --------------------------------------------------------

    def print_summary(self, results: list[ValidationResult] | None = None) -> None:
        """Print a formatted summary of validation results.

        Parameters
        ----------
        results : list[ValidationResult] | None
            Results to print.  Falls back to ``self.results``.
        """
        if results is None:
            results = self.results

        by_category: dict[str, list[ValidationResult]] = {}
        for r in results:
            by_category.setdefault(r.category, []).append(r)

        errors = [r for r in results if not r.passed and r.severity == "error"]
        warnings = [r for r in results if r.severity == "warning"]

        print("\n" + "=" * 80)
        print("VALIDATION SUMMARY")
        print("=" * 80)

        for category, cat_results in by_category.items():
            print(f"\n{category}:")
            print("-" * 80)
            for r in cat_results:
                if not r.passed and r.severity == "error":
                    icon = "FAIL"
                elif r.severity == "warning":
                    icon = "WARN"
                elif not r.passed:
                    icon = "FAIL"
                else:
                    icon = "OK  "
                message_lines = r.message.split("\n")
                first_line = f"{icon} {r.item}: {message_lines[0]}"
                print(f"  {first_line}")
                for line in message_lines[1:]:
                    print(f"    {line}")

        print("\n" + "=" * 80)
        if errors:
            print(f"RESULT: FAILED with {len(errors)} error(s)")
            if warnings:
                print(f"        Also {len(warnings)} warning(s)")
            print("\nPlease fix the errors above before running the pipeline.")
        elif warnings:
            print(f"RESULT: PASSED with {len(warnings)} warning(s)")
            print("\nYou can proceed, but consider addressing the warnings above.")
        else:
            print("RESULT: PASSED")
            print("\nAll checks passed! Ready to run the pipeline.")
        print("=" * 80 + "\n")

    def has_errors(self, results: list[ValidationResult] | None = None) -> bool:
        """Check if there are any error-level failures.

        Parameters
        ----------
        results : list[ValidationResult] | None
            Results to check.  Falls back to ``self.results``.

        Returns
        -------
        bool
            *True* if there are any error-level failures.
        """
        if results is None:
            results = self.results
        return any(not r.passed and r.severity == "error" for r in results)
