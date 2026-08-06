import argparse
from pathlib import Path

import numpy as np
from astropy.table import Table


# Normalize truthy string-like values into real booleans for a specified column.
# Useful because CSV round-trips may store booleans as strings.
def _coerce_bool_column(lc_table, col_name):
    if col_name not in lc_table.colnames:
        return lc_table

    values = np.array(lc_table[col_name]).astype(str)
    truthy = {"true", "1", "t", "yes"}
    coerced = np.array([v.strip().lower() in truthy for v in values], dtype=bool)
    lc_table[col_name] = coerced
    return lc_table


# Print a compact schema/summary report for every column.
# Optionally includes full array dumps for key light-curve columns.
def inspect_table(lc_table, show_arrays=True):
    print("Columns:", lc_table.colnames)

    for col in lc_table.colnames:
        arr = np.array(lc_table[col])
        if np.issubdtype(arr.dtype, np.number):
            print(
                f"{col:20s} dtype={arr.dtype} shape={arr.shape} min={arr.min():.4f} max={arr.max():.4f}"
            )
        else:
            unique_vals = np.unique(arr)
            preview = ", ".join(str(v) for v in unique_vals[:5])
            print(
                f"{col:20s} dtype={arr.dtype} shape={arr.shape} unique={len(unique_vals)} [{preview}]"
            )

    if show_arrays:
        for col in (
            "T_flux",
            "uncertainty",
            "x_centroid",
            "y_centroid",
            "norm_t",
            "norm_flux",
            "norm_flux_error",
        ):
            if col in lc_table.colnames:
                print(f"\n{col}")
                print(lc_table[col])

    print("\nFirst row")
    print(lc_table[0])
    print("\nLast row")
    print(lc_table[-1])


# Read CSV output into an Astropy Table, clean types, then print inspection output.
def inspect_csv(csv_path, show_arrays=True):
    lc_table = Table.read(csv_path, format="ascii.csv")
    lc_table = _coerce_bool_column(lc_table, "oot_mask")
    inspect_table(lc_table, show_arrays=show_arrays)
    return lc_table


# Run the target-star photometry phase and return its table output.
def run_target_pipeline(output_csv):
    import target_star_photometry

    return target_star_photometry.run_pipeline(output_csv=output_csv)


# Load the generated CSV into q/kdb+ using the q script provided.
def maybe_load_to_q(output_csv, q_script):
    import target_star_photometry

    print(f"Loading into q using {q_script}...")
    target_star_photometry.load_into_q(output_csv=output_csv, q_script=q_script)


# Run the bycatch phase, which performs field-source (non-target) photometry.
def run_bycatch_phase():
    import bycatch_photometry

    print("Running bycatch photometry phase...")
    bycatch_photometry.execute_bycatch_photometry()


# CLI entrypoint that supports three modes:
# - inspect: inspect existing CSV only
# - target: run target pipeline
# - full: run target pipeline and optional bycatch phase
def main(argv=None):
    parser = argparse.ArgumentParser(description="Photometry pipeline entrypoint")
    parser.add_argument(
        "mode",
        choices=["inspect", "target", "full"],
        nargs="?",
        default="inspect",
        help="inspect: inspect CSV only, target: run target photometry, full: target + bycatch",
    )
    parser.add_argument(
        "--csv",
        default="target_star_photometry.csv",
        help="Path to photometry CSV output",
    )
    parser.add_argument(
        "--q-script",
        default="targetStarOnly.q",
        help="q loader script path used with --load-q",
    )
    parser.add_argument(
        "--load-q",
        action="store_true",
        help="Load generated CSV into q after target pipeline runs",
    )
    parser.add_argument(
        "--skip-inspect",
        action="store_true",
        help="Skip table inspection printout",
    )
    parser.add_argument(
        "--no-arrays",
        action="store_true",
        help="Hide full array dumps and print summaries only",
    )
    parser.add_argument(
        "--skip-bycatch",
        action="store_true",
        help="Only for full mode: skip bycatch phase",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would run without executing the pipeline",
    )

    args = parser.parse_args(argv)
    csv_path = Path(args.csv)

    # Dry-run mode prints planned actions and exits without processing data.
    if args.dry_run:
        print(f"Mode: {args.mode}")
        print(f"CSV path: {csv_path}")
        if args.mode == "inspect":
            print(f"CSV exists: {csv_path.exists()}")
            if not csv_path.exists():
                print("Warning: inspect mode will fail because CSV file does not exist.")
        print(f"Load to q: {args.load_q}")
        if args.load_q:
            print(f"q script: {args.q_script}")
        print(f"Inspect output: {not args.skip_inspect}")
        print(f"Show arrays: {not args.no_arrays}")
        if args.mode == "full":
            print(f"Run bycatch: {not args.skip_bycatch}")
        return

    # Inspect mode validates file existence, then prints table diagnostics.
    if args.mode == "inspect":
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV not found: {csv_path}")
        inspect_csv(str(csv_path), show_arrays=not args.no_arrays)
        return

    # For target/full modes, run target photometry to generate/refresh the CSV.
    print("Running target star photometry pipeline...")
    lc_table = run_target_pipeline(str(csv_path))

    # Optional q load step for database ingestion.
    if args.load_q:
        maybe_load_to_q(str(csv_path), args.q_script)

    # Full mode optionally runs the bycatch stage after target processing.
    if args.mode == "full" and not args.skip_bycatch:
        run_bycatch_phase()

    # Final inspection gives quick QA feedback on output values/types.
    if not args.skip_inspect:
        lc_table = _coerce_bool_column(lc_table, "oot_mask")
        inspect_table(lc_table, show_arrays=not args.no_arrays)


# Script execution guard.
if __name__ == "__main__":
    main()
