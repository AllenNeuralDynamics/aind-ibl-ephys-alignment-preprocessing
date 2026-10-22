"""CLI entry point for ``aind-ibl-build-manifest``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> None:
    """Build a manifest CSV from a track-to-stream mapping file."""
    parser = argparse.ArgumentParser(
        prog="aind-ibl-build-manifest",
        description="Build an IBL preprocessing manifest from a TOML mapping of tracks to streams.",
    )
    parser.add_argument("mapping", type=Path, help="Mapping TOML file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Manifest CSV to write.")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Where the sorted recordings are mounted; a sort not found here is listed from aind-open-data.",
    )
    args = parser.parse_args()

    from aind_ibl_ephys_alignment_preprocessing.manifest_builder import ManifestBuildError, build_from_files

    try:
        manifest, unmapped = build_from_files(args.mapping, args.data_root)
    except ManifestBuildError as err:
        print("Mapping is incomplete or inconsistent:", file=sys.stderr)
        for problem in err.problems:
            print(f"  - {problem}", file=sys.stderr)
        sys.exit(1)
    for track in unmapped:
        print(f"note: annotated layer {track} is not in the manifest", file=sys.stderr)
    manifest.to_csv(args.output, index=False)
    print(f"Wrote {len(manifest)} rows to {args.output}")


if __name__ == "__main__":
    main()
