#!/usr/bin/env python3
"""Compare ChiMorse JSON models with reference .dat data at a selected zeta.

Read-only for the input data. Requires Python 3.9+; no third-party packages.

The tabulated C++ code retains only rows with abs(zeta-target) <= tolerance.
This script reports *both* whole-source and selected-slice disk sizes.
Memory numbers are approximate numeric/container payloads, NOT measured RSS.
"""

import argparse
import csv
import json
from pathlib import Path

CLASSES = ("EP", "EA", "OP", "OA")
DEFAULT_DATA = Path("/home/hadis/chimorse/data_internal/PA")
DEFAULT_MODELS = Path("/home/hadis/chimorse/examples/models")
BYTES_PER_DOUBLE = 8          # mainMD.cpp uses double
BYTES_PER_RADIAL_POINT = 16  # C++ RadialPoint<double>: r and energy
BYTES_PER_FOURIER_TERM = 24  # typical C++ FourierTerm<double>: 4 ints + double
BYTES_PER_VECTOR_HEADER = 24  # typical 64-bit libstdc++ std::vector object
ANGULAR_GRID_NODES = 360 * 360


def kib(n):
    return f"{n / 1024:,.1f}"


def mib(n):
    return f"{n / (1024 * 1024):,.2f}"


def inspect_table(path, target_zeta, tolerance):
    """Count valid rows, rows in the selected zeta slice, and their raw text bytes."""
    if not path.is_file():
        raise FileNotFoundError(path)

    total_rows = 0
    selected_rows = 0
    selected_bytes = 0
    zeta_values = set()

    # Stream the original file; do not duplicate the potentially large table.
    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.split()
            if not fields or fields[0].startswith(b"#"):
                continue
            if len(fields) != 5:
                raise ValueError(
                    f"{path}:{line_number}: expected 5 columns: "
                    "phi1 phi2 zeta R energy"
                )
            zeta = float(fields[2])
            zeta_values.add(zeta)
            total_rows += 1
            if abs(zeta - target_zeta) <= tolerance:
                selected_rows += 1
                selected_bytes += len(raw_line)  # as if selected lines were saved unchanged

    return {
        "class": path.stem.removeprefix("E_all_"),
        "source_bytes": path.stat().st_size,
        "all_rows": total_rows,
        "zeta_slices": len(zeta_values),
        "selected_rows": selected_rows,
        "selected_text_bytes": selected_bytes,
        "radial_payload_bytes": BYTES_PER_RADIAL_POINT * selected_rows,
        "angular_lookup_headers_bytes": ANGULAR_GRID_NODES * BYTES_PER_VECTOR_HEADER,
        "estimated_cpp_floor_bytes": (BYTES_PER_RADIAL_POINT * selected_rows
                                      + ANGULAR_GRID_NODES * BYTES_PER_VECTOR_HEADER),
    }


def inspect_model(path):
    """Count actually STORED coefficients and NONZERO coefficients."""
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = []

    for interaction in CLASSES:
        model = data["interactions"][interaction]
        basis_size = len(model["basis_terms"])
        stored = 0
        nonzero = 0
        active = [False] * basis_size
        parameter_info = {}

        for parameter in ("D", "re", "alpha"):
            obj = model["parameters"][parameter]
            if obj["type"] == "constant":
                parameter_info[parameter] = "constant"
                continue
            if obj["type"] != "fourier":
                raise ValueError(f"Unexpected parameter type in {path}: {obj['type']}")

            coeffs = obj["coefficients"]
            if len(coeffs) != basis_size:
                raise ValueError(f"{path}: {interaction}/{parameter} basis length mismatch")
            stored += len(coeffs)
            count = sum(c != 0 for c in coeffs)
            nonzero += count
            parameter_info[parameter] = str(count)
            for index, coeff in enumerate(coeffs):
                if coeff != 0:
                    active[index] = True

        rows.append({
            "class": interaction,
            "basis_positions": basis_size,
            "stored_coefficients": stored,
            "nonzero_coefficients": nonzero,
            "active_basis_positions": sum(active),
            "D_nonzero": parameter_info["D"],
            "re_nonzero": parameter_info["re"],
            "alpha_nonzero": parameter_info["alpha"],
            "coefficient_only_bytes": stored * BYTES_PER_DOUBLE,
            "estimated_cpp_term_arrays_bytes": stored * BYTES_PER_FOURIER_TERM,
        })

    summary = {
        "model": str(path.parent.name) if path.parent != path.parent.parent else str(path),
        "model_path": str(path),
        "json_bytes": path.stat().st_size,
        "basis_positions": sum(r["basis_positions"] for r in rows),
        "stored_coefficients": sum(r["stored_coefficients"] for r in rows),
        "nonzero_coefficients": sum(r["nonzero_coefficients"] for r in rows),
        "active_basis_positions": sum(r["active_basis_positions"] for r in rows),
        "coefficient_only_bytes": sum(r["coefficient_only_bytes"] for r in rows),
        "estimated_cpp_term_arrays_bytes": sum(r["estimated_cpp_term_arrays_bytes"] for r in rows),
    }
    for r in rows:
        r["model"] = summary["model"]
    return summary, rows


def write_csv(path, records):
    if not records:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA,
                        help="directory containing E_all_EP/EA/OP/OA.dat")
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS,
                        help="directory to search recursively for chimorse_all.json")
    parser.add_argument("--zeta", type=float, default=0.0,
                        help="selected zeta in Angstrom (C++ default: 0)")
    parser.add_argument("--tolerance", type=float, default=1e-8,
                        help="selection tolerance in Angstrom (C++ default: 1e-8)")
    parser.add_argument("--out", type=Path, default=Path("storage_benchmark_results"),
                        help="directory for CSV reports")
    args = parser.parse_args()

    tables = [inspect_table(args.data / f"E_all_{cls}.dat", args.zeta, args.tolerance)
              for cls in CLASSES]

    print(f"\nTABULATED REFERENCE DATA (zeta = {args.zeta:g} +/- {args.tolerance:g} Angstrom)")
    print(f"{'Class':<6} {'Whole file MiB':>15} {'All rows':>12} {'Zeta values':>12} "
          f"{'Selected rows':>15} {'Slice text MiB':>16} {'C++ floor MiB*':>16}")
    for r in tables:
        print(f"{r['class']:<6} {mib(r['source_bytes']):>15} {r['all_rows']:>12,} "
              f"{r['zeta_slices']:>12,} {r['selected_rows']:>15,} "
              f"{mib(r['selected_text_bytes']):>16} "
              f"{mib(r['estimated_cpp_floor_bytes']):>16}")

    total = {k: sum(r[k] for r in tables) for k in
             ("source_bytes", "all_rows", "selected_rows", "selected_text_bytes", "radial_payload_bytes", "estimated_cpp_floor_bytes")}
    print(f"{'TOTAL':<6} {mib(total['source_bytes']):>15} {total['all_rows']:>12,} "
          f"{'—':>12} {total['selected_rows']:>15,} "
          f"{mib(total['selected_text_bytes']):>16} "
          f"{mib(total['estimated_cpp_floor_bytes']):>16}")

    model_paths = sorted(args.models.rglob("chimorse_all.json"))
    if not model_paths:
        raise FileNotFoundError(f"No chimorse_all.json found under {args.models}")
    summaries, details = [], []
    for path in model_paths:
        summary, model_rows = inspect_model(path)
        # Display a distinct name for the root-level combined model too.
        summary["model"] = str(path.parent.relative_to(args.models)) or "."
        for row in model_rows:
            row["model"] = summary["model"]
        summaries.append(summary)
        details.extend(model_rows)

    print("\nCHIMORSE MODELS (combined JSON file actually loaded by C++)")
    print(f"{'Model':<23} {'JSON KiB':>10} {'Stored':>10} {'Nonzero':>10} "
          f"{'Active basis':>14} {'Coeff KiB':>11} {'C++ terms KiB*':>16}")
    for s in summaries:
        print(f"{s['model'][:23]:<23} {kib(s['json_bytes']):>10} "
              f"{s['stored_coefficients']:>10,} {s['nonzero_coefficients']:>10,} "
              f"{s['active_basis_positions']:>14,} "
              f"{kib(s['coefficient_only_bytes']):>11} "
              f"{kib(s['estimated_cpp_term_arrays_bytes']):>16}")

    print("\nDETAILS BY CLASS (D/re/alpha = nonzero coefficient counts)")
    print(f"{'Model':<23} {'Class':<6} {'Basis':>6} {'Stored':>8} {'Nonzero':>9} "
          f"{'Active':>7} {'D':>9} {'re':>9} {'alpha':>10}")
    for r in details:
        print(f"{r['model'][:23]:<23} {r['class']:<6} {r['basis_positions']:>6} "
              f"{r['stored_coefficients']:>8} {r['nonzero_coefficients']:>9} "
              f"{r['active_basis_positions']:>7} {r['D_nonzero']:>9} "
              f"{r['re_nonzero']:>9} {r['alpha_nonzero']:>10}")

    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "tabulated_data.csv", tables)
    write_csv(args.out / "chimorse_models.csv", summaries)
    write_csv(args.out / "chimorse_by_class.csv", details)
    print(f"\nCSV tables saved in: {args.out.resolve()}")
    print("\nInterpretation:")
    print("- Whole file: original .dat including EVERY zeta value (disk size).")
    print("- Slice text: bytes of only matching zeta rows, still with 5 columns; not a created file.")
    print("- *C++ floor: approximate radial payload (r,energy doubles) + 360x360 vector headers;")
    print("  assumes 24-byte std::vector headers and excludes map nodes, look-up keys, allocations.")
    print("  NOT measured total C++ RAM; see CSV for radial and lookup parts separately.")
    print("- Stored: ALL coefficients in JSON arrays, including zeros after pruning.")
    print("- Active basis: indices nonzero in at least one Fourier parameter.")
    print("- *C++ term arrays: illustrative 24 bytes per stored FourierTerm<double>.")
    print("  Excludes vector capacities/headers and other objects; NOT measured total RAM.")
    print("- During MD, the current mainMD.cpp constructs separate force and potential objects;")
    print("  both load their own model copy. Do not compare these numbers to whole-process RAM.")


if __name__ == "__main__":
    main()
