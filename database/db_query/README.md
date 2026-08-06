# Benchmark Data Pipeline (kdb+/q + CSV Ingestion)

This repository provides tools for:

1. Inspecting light-curve style CSV files from Python.
2. Routing CSV/FITS inputs into kdb+/q.
3. Running reusable q queries for target-level analysis.

The documentation in this repository is intentionally scoped to the maintained ingestion/query workflow.

## Repository structure

```text
.
|-- db_pipeline.py
|-- input_router.py
|-- targetStarOnly.q
|-- config/
|   |-- __init__.py
|   |-- config_hatp32.py
|   |-- config_hatp36.py
|   |-- config_wasp12.py
|-- data files/
|-- images/
|-- requirements.txt
```

## Core scripts

1. `db_pipeline.py`
CLI entrypoint used here for table/CSV inspection.

2. `input_router.py`
Schema-aware CSV router/loader into q.

3. `targetStarOnly.q`
Table schema, ingestion helpers, and query utilities in q.

## Installation

### 1. Create and activate virtual environment

```in powershell terminal:
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

### 2. Install Python dependencies

```in powershell terminal:
pip install -r requirements.txt
```

### 3. Verify q is available

```in powershell terminal:
q -q
```

If this fails, install kdb+ and/or add q to PATH.

## Configuration

The owned workflow is configured through CLI options and q script paths.

## Data availability

Light-curve style CSV data files are included in this repository.

Raw and aligned FITS image datasets are maintained outside this repository and are not versioned on GitHub.

This repository focuses on code, ingestion/query tooling, and light-curve style CSV workflow. To run full photometry from FITS, provide access to the external image dataset.

## Quick start

### Inspect existing CSV

```in powershell terminal:
.\.venv\Scripts\python db_pipeline.py inspect --csv target_star_photometry.csv --no-arrays
```

### Preview ingestion commands (dry run)

```in powershell terminal:
.\.venv\Scripts\python input_router.py --input target_star_photometry.csv --target HAT-P-36b --dry-run
```

### Load CSV into q

```in powershell terminal:
.\.venv\Scripts\python input_router.py --input target_star_photometry.csv --target HAT-P-36b --q-script targetStarOnly.q
```

### Load FITS path via router

```in powershell terminal:
.\.venv\Scripts\python input_router.py --input images/raw_images_hatp36/sample.fits --target HAT-P-36b --q-script targetStarOnly.q
```

## q usage

In a q session:

```in q terminal:
\l targetStarOnly.q
n:loadTargetCSVAs["target_star_photometry.csv";"HAT-P-36b"]
show count getObsByTarget["HAT-P-36b"]
show 5#getObsByTarget["HAT-P-36b"]
```

If q is launched from a different folder, use an absolute path:

```in q terminal:
n:loadTargetCSVAs["C:/kdb/db/db query/target_star_photometry.csv";"HAT-P-36b"]
```

Useful functions:

1. `showStats[]`
2. `getObsByTarget[target]`
3. `getObsByTimeRange[start;end]`
4. `getObsByFlux[min;max]`
5. `getObsByPhase[min;max]`
6. `saveDB["benchmark.qdb"]`
7. `loadDB["benchmark.qdb"]`

## Outputs

Typical artifacts:

1. q loader command output and validation logs.
2. q in-memory tables and optional persisted qdb file.

## Automation status

### Automated

1. CSV schema validation before ingestion is automated in `input_router.py`.
2. CSV/FITS routing into q loader commands is automated in `input_router.py`.
3. Dataset inspection summaries are automated in `db_pipeline.py inspect`.
4. Target-level query helpers are automated in `targetStarOnly.q` once data is loaded.

### Manual

1. Choosing the correct input file and target label still requires user input.
2. Opening q and running ad hoc analysis queries is manual.
3. Persisting and restoring q databases (`saveDB`/`loadDB`) is manual.
4. End-to-end scientific photometry generation from raw FITS to CSV is outside this maintained workflow and must be handled separately.

## Troubleshooting

### CSV schema errors in `input_router.py`

Ensure at least one alias from each required group exists:

1. Time: `time` or `BJD_TDB`
2. Flux: `T_flux` or `Flux` or `flux`
3. Uncertainty: `uncertainty` or `Uncertainty`

### q commands in PowerShell

Commands such as `\l` and `show` must run in q (`q)`), not PowerShell (`PS>`).

## Reproducibility checklist

1. Pin dependency versions in `requirements.txt`.
2. Keep raw data immutable; write generated artifacts separately.
3. Record the exact ingestion command used for each dataset.
4. Save qdb snapshots for published results.

## Contributing

Contributions are welcome via pull requests with a clear summary of changes and validation steps.

## Credits

Repository owner and maintainer: current project author.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
