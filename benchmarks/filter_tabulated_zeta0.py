#!/usr/bin/env python3
"""Create a zeta=0-only copy of the four ChiSurfMD reference .dat files.

The input files are never changed. The output keeps the exact input columns,
row ordering, numeric formatting, and filenames expected by TabularDFT.
Only rows with |zeta| <= 1e-8 are retained.

Usage:
    python3 filter_tabulated_zeta0.py

By default the script determines which of the five columns contains zeta
using the known number of zeta=0 rows in the EP reference. Pass
--zeta-column N (1-based) to select the column explicitly instead.
"""

import argparse
import os
from pathlib import Path
from tempfile import TemporaryDirectory


DEFAULT_SOURCE = Path("/home/hadis/chimorse/data_internal/PA")
DEFAULT_OUTPUT = Path("/home/hadis/chimorse/data_internal/PA_zeta0")
TOLERANCE = 1e-8
NUMBER_OF_COLUMNS = 5

# Values previously checked by benchmark_storage.py for zeta=0.
EXPECTED_ROWS = {"EP": 303210, "EA": 305607, "OP": 321900, "OA": 323963}


def numeric_fields(raw_line: bytes, path: Path, line_number: int):
    """Return data columns, or None for blank/comment lines."""
    stripped = raw_line.strip()
    if not stripped or stripped.startswith(b"#"):
        return None
    fields = stripped.split()
    if len(fields) != NUMBER_OF_COLUMNS:
        raise ValueError(
            f"{path}:{line_number}: expected {NUMBER_OF_COLUMNS} columns, "
            f"found {len(fields)}"
        )
    return fields


def detect_zeta_column(ep_path: Path) -> int:
    """Find the unique column with the expected number of zero-zeta EP rows."""
    near_zero_counts = [0] * NUMBER_OF_COLUMNS
    with ep_path.open("rb") as source:
        for line_number, line in enumerate(source, start=1):
            fields = numeric_fields(line, ep_path, line_number)
            if fields is None:
                continue
            for index, value in enumerate(fields):
                if abs(float(value)) <= TOLERANCE:
                    near_zero_counts[index] += 1

    matches = [
        index for index, count in enumerate(near_zero_counts)
        if count == EXPECTED_ROWS["EP"]
    ]
    if len(matches) != 1:
        raise ValueError(
            "Unable to identify the zeta column uniquely from EP. "
            f"Zero-value counts per column: {near_zero_counts}. "
            "Inspect the column order and rerun with --zeta-column N "
            "(column numbers start at 1)."
        )
    return matches[0]


def filter_file(source_path: Path, output_path: Path, zeta_index: int,
                expected_rows: int) -> tuple[int, int]:
    if output_path.exists():
        raise FileExistsError(
            f"Output already exists and will not be overwritten: {output_path}"
        )

    retained_rows = 0
    with TemporaryDirectory(prefix="zeta0_", dir=output_path.parent) as temp_dir:
        temporary_path = Path(temp_dir) / output_path.name
        with source_path.open("rb") as source, temporary_path.open("wb") as output:
            for line_number, line in enumerate(source, start=1):
                fields = numeric_fields(line, source_path, line_number)
                if fields is None:
                    # Preserve comments and blank lines, if any.
                    output.write(line)
                elif abs(float(fields[zeta_index])) <= TOLERANCE:
                    output.write(line)
                    retained_rows += 1

        if retained_rows != expected_rows:
            raise ValueError(
                f"{source_path.name}: retained {retained_rows:,} rows; "
                f"expected {expected_rows:,}. No output was installed. "
                "Check the zeta column and tolerance."
            )
        output_size = temporary_path.stat().st_size
        os.replace(temporary_path, output_path)

    return retained_rows, output_size


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--zeta-column", type=int, choices=range(1, NUMBER_OF_COLUMNS + 1),
        help="1-based zeta column; automatically inferred when omitted",
    )
    args = parser.parse_args()

    source_dir = args.source.resolve()
    output_dir = args.output.resolve()
    if source_dir == output_dir:
        parser.error("The output directory must differ from the source directory.")

    source_files = {
        cls: source_dir / f"E_all_{cls}.dat" for cls in EXPECTED_ROWS
    }
    for source_file in source_files.values():
        if not source_file.is_file():
            parser.error(f"Missing source file: {source_file}")

    zeta_index = (
        args.zeta_column - 1 if args.zeta_column is not None
        else detect_zeta_column(source_files["EP"])
    )
    print(f"Zeta column: {zeta_index + 1} (1-based); tolerance: {TOLERANCE:g}")

    # Validate all destination filenames before modifying the output directory.
    for cls in EXPECTED_ROWS:
        existing_path = output_dir / f"E_all_{cls}.dat"
        if existing_path.exists():
            parser.error(f"Refusing to overwrite existing output: {existing_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    total_rows = 0
    total_bytes = 0
    for cls, expected_rows in EXPECTED_ROWS.items():
        target = output_dir / f"E_all_{cls}.dat"
        rows, size = filter_file(source_files[cls], target, zeta_index, expected_rows)
        total_rows += rows
        total_bytes += size
        print(f"{cls}: {rows:,} rows, {size / 2**20:.2f} MiB -> {target}")

    print(f"TOTAL: {total_rows:,} rows, {total_bytes / 2**20:.2f} MiB")
    print("Original files were left unchanged.")


if __name__ == "__main__":
    main()