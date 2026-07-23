import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def infer_target_from_path(input_path: Path) -> str:
    stem = input_path.stem.upper()
    return stem.replace("_", "-").replace(" ", "")


# Alias groups mirror the q loader column resolver in targetStarOnly.q.
CSV_ALIAS_GROUPS = {
    "time": ["time", "BJD_TDB"],
    "flux": ["T_flux", "Flux", "flux"],
    "uncertainty": ["uncertainty", "Uncertainty"],
    "norm_t": ["norm_t"],
    "oot_mask": ["oot_mask"],
    "norm_flux": ["norm_flux"],
    "norm_flux_error": ["norm_flux_error"],
    "orbital_phase": ["orbital_phase", "Orbital Phase"],
    "model": ["model", "Model"],
    "airmass": ["airmass", "Airmass"],
    "amplitude": ["amplitude"],
    "offset": ["offset"],
    "x_centroid": ["x_centroid"],
    "y_centroid": ["y_centroid"],
    "sigma_x": ["sigma_x"],
    "sigma_y": ["sigma_y"],
    "rotation": ["rotation"],
}

# Minimum columns needed for useful ingestion.
REQUIRED_GROUPS = ["time", "flux", "uncertainty"]


def _read_csv_headers(csv_path: Path) -> list[str]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"CSV is empty: {csv_path}") from exc

    headers = [h.strip() for h in header if h is not None]
    if not headers:
        raise ValueError(f"CSV header row is empty: {csv_path}")
    return headers


def validate_csv_schema(csv_path: Path) -> None:
    headers = set(_read_csv_headers(csv_path))

    missing_required = []
    for group in REQUIRED_GROUPS:
        aliases = CSV_ALIAS_GROUPS[group]
        if not any(alias in headers for alias in aliases):
            missing_required.append((group, aliases))

    if missing_required:
        accepted = "; ".join(
            f"{group}: {', '.join(aliases)}" for group, aliases in missing_required
        )
        raise ValueError(
            "CSV schema not compatible with targetStarOnly loader. "
            f"Missing required column groups -> {accepted}."
        )

    # Optional visibility: tell user which optional groups are present.
    optional_groups = [g for g in CSV_ALIAS_GROUPS if g not in REQUIRED_GROUPS]
    matched_optional = [
        group
        for group in optional_groups
        if any(alias in headers for alias in CSV_ALIAS_GROUPS[group])
    ]
    print(
        "CSV schema check passed. "
        f"Required groups found: {', '.join(REQUIRED_GROUPS)}. "
        f"Optional groups found: {', '.join(matched_optional) if matched_optional else 'none'}."
    )


def build_q_program(input_path: Path, target: str | None, q_script: Path) -> str:
    ext = input_path.suffix.lower()
    path_posix = input_path.resolve().as_posix()
    q_script_posix = q_script.resolve().as_posix()

    lines = [f"\\l {q_script_posix}"]

    if ext == ".csv":
        if target:
            # Explicit target is best for generic CSV names.
            lines.append(f'n:sortTargetInput["{path_posix}";"{target}"]')
        else:
            # Fallback to inferred target from CSV file name.
            lines.append(f'n:loadTargetCSV["{path_posix}"]')
    elif ext in {".fits", ".fit", ".fts"}:
        resolved_target = target or infer_target_from_path(input_path)
        lines.append(f'n:sortTargetInput["{path_posix}";"{resolved_target}"]')
    else:
        raise ValueError(
            f"Unsupported file type: {ext or '<no extension>'}. Use CSV or FITS/FIT/FTS."
        )

    lines.append("show n")
    lines.append("showStats[]")
    lines.append("\\\\")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Route any CSV/FITS file into targetStarOnly.q without hardcoded target names"
    )
    parser.add_argument("--input", required=True, help="Path to CSV/FITS input file")
    parser.add_argument(
        "--target",
        default=None,
        help="Optional target name (for example HATP-32-B). If omitted for CSV, inferred from file name.",
    )
    parser.add_argument(
        "--q-script",
        default="targetStarOnly.q",
        help="Path to q loader script",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print q commands without executing",
    )

    args = parser.parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists() or not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    q_script = Path(args.q_script)
    if not q_script.exists() or not q_script.is_file():
        raise FileNotFoundError(f"q script not found: {q_script}")

    q_binary = shutil.which("q")
    if not q_binary:
        raise RuntimeError("q executable not found on PATH")

    if input_path.suffix.lower() == ".csv":
        validate_csv_schema(input_path)

    q_program = build_q_program(input_path=input_path, target=args.target, q_script=q_script)

    if args.dry_run:
        print("=== q program ===")
        print(q_program)
        return

    result = subprocess.run(
        [q_binary],
        input=q_program,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())

    if result.returncode != 0:
        raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
