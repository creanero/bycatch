 `THE PRIMARY CONTRIBUTOR OF THIS WORK IS *D.LAWLEY* WITH ADAPTATIONS MADE FROM *DEAN WINTERS*`

`- Dean Winters - - dean.winters4@mail.dcu.ie - - May 2026`


# Astronomical Bycatch Pipeline

Photometry pipeline for transit detection from archived FITS imgs (bycatch method).

## File structure

```
bycatch-pipeline/
config.py                             # user-specific paths and params.
target_star_photometry.py             # pahse 1: alignment ; target/comp star photometry
bycatch-photometry.py                 # phase 2: source detection ; bycatch photometry
run_pipeline.py                       # run phases in order
photutils-photometry-commented.ipynb  # reference notebook
requirements.txt                      # the dependencies
```

## How the files relate

The notebook (`photutils-photometry-commented.ipynb`) is the main reference and highly annotated version of the pipeline. The `.py` scripts are developed derivative scripts:

- `target_star_photometry.py`, from the notebook, executes: imports ; alignment loop ; target/comp photometry loop ; normalisation ; binned plot. It also adds error calculations, SNR, and an error bar plot (not found in notebook).
- `bycatch-photometry.py` covers: source detection (DAOStarFinder) ; bycatch photometry loop ; flux matrix ; plotings.

## Setup

Install dependencies:
```
pip install -r requirements.txt
```

Edit `config.py` with user-specific filepaths and params.

## Running

Execute pipeline end to end:
```
python run_pipeline.py
```

Or can run indivdiually:
```
python target_star_photometry.py
python bycatch-photometry.py
```

From my understanding of the work, `target_star_photometry.py` has to be ran `bycatch-photometry.py` since it makes the aligned imgs which the bycatch script uses.

## Notes for further work (from DL)

- EXOTIC outputs (particularly the `temp` folder) has useful info on optimised params for aperture and annulus sizes, and comp star suggestions. can use these in `config.py`.
- AstroImageJ (AIJ) may be useful for identifying target and comp pixel positions, for previewing image stacks, and to estimate FWHMs.
- Notebook uses rawe MJD time values in its OOT mask; `target_star_photometry.py` shifts to time=0 start for its mask. Both domain sets are separately denoted in the `config.py`.
