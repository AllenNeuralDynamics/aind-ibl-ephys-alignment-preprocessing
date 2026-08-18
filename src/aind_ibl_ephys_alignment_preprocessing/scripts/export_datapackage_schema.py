"""Export the generated datapackage JSON Schema artifact."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from aind_ibl_ephys_alignment_preprocessing.datapackage_schema import (
    datapackage_schema_json,
    source_tree_datapackage_schema_path,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Export or check the committed datapackage JSON Schema artifact."""
    parser = argparse.ArgumentParser(description="Export the datapackage JSON Schema.")
    parser.add_argument(
        "--output",
        type=Path,
        default=source_tree_datapackage_schema_path(),
        help="Schema file to write or check. Defaults to the package's bundled schema artifact.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the output file does not match the generated schema.",
    )
    args = parser.parse_args(argv)

    schema_text = datapackage_schema_json()
    if args.check:
        try:
            current = args.output.read_text()
        except FileNotFoundError:
            print(f"schema artifact is missing: {args.output}", file=sys.stderr)
            return 1
        if current != schema_text:
            print(f"schema artifact is stale: {args.output}", file=sys.stderr)
            return 1
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(schema_text)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
