"""
TODO: verify this works for all EXOTIC data available in database
Combines the two EXOTIC output files into a single DataFrame to match the column format of targetStarOnly.q.

EXOTIC produces two separate files per target-star run:

    FinalLightCurve CSV:
        time-series photometry
        (BJD_TDB, Orbital Phase, Flux, Uncertainty, Model, Airmass)

    PSF data CSV:
        per-frame data
        (x_centroid, y_centroid, amplitude, sigma_x, sigma_y, rotation, offset)

Both files have one row per FITS frame and no shared key column; they're just aligned by row pos. (frame idx).

-----------------

PIPELINE:

Currently have Python to disk, Q read from disk.

PyKX/otherwise path (Python->Q directly):
    Photometry pipeline produces DataFrames in memory
    ->merge_datasets(lc_df, psf_df) called directly (no disk read)
    ->combined DataFrame pushed to Q with pykx

load_and_combine() and merge_datasets() are split such that when PyKX is implemented load_and_combine can be discarded (is 
a temporary disk-reading helper);
and thus will just work off merge_datasets() with in-memory data.

--------------------

PSF.csv HEADER:

PSF.cv has an empty # which pandas will read as the first col. header=None and reading by position corrects this.
"""



import pandas as pd
import numpy as np
from pathlib import Path



# -- DATA FORMAT --



# targetStarOnly.q column order (13 columns, 1 row per FITS frame)

# observations: BJD_TDB, orbital_phase, flux, uncertainty, model, airmass, amplitude, offset
# centroids: x_centroid, y_centroid
# psf_params: sigma_x, sigma_y, rotation

# TODO: alter loadObservations() in Q to match this nomenclature

COMBINED_COLUMNS = [
    "BJD_TDB",
    "orbital_phase",
    "flux",
    "uncertainty",
    "model",
    "airmass",
    "amplitude",
    "offset",
    "x_centroid",
    "y_centroid",
    "sigma_x",
    "sigma_y",
    "rotation"
]

LC_COLUMNS = ["BJD_TDB", "orbital_phase", "flux", "uncertainty", "model", "airmass"]

 # Note that the EXOTIC csv outputs have a trailing comma, which pandas reads as an empty col
 # Will ignore later by specifying usecols=range(len(PSF_COLUMNS))
PSF_COLUMNS =[
    "x_centroid",
    "y_centroid",
    "amplitude",
    "sigma_x",
    "sigma_y",
    "rotation",
    "offset"
]



# -- CSV LOADERS --



def import_lightcurves(lc_path:str) -> pd.DataFrame:
    """
    Read an EXOTIC FinalLightCurve CSV.

    EXOTIC format:
    Skip Line 1
    Line 2: # BJD_TDB,Orbital Phase, etc.
    Line 3+: data rows (float)

    Returns: 
        padas.DataFrame with cols that match LC_COLUMNS
    """

    path = Path(lc_path)

    # Column formatting: EXOTIC stores col names on the secodn line
    with open(path) as f:
        raw_lines = f.readlines()
    header = raw_lines[1].lstrip("# ").strip()
    col_names = [c.strip() for c in header.split(",")]

    # Read
    df = pd.read_csv(path, skiprows=2, header=None, names=col_names)

    # Match cols to targetStarOnly.q names (no spaces & lowercase)
    df = df.rename(columns={
        "Orbital Phase":"orbital_phase",
        "Flux":"flux",
        "Uncertainty":"uncertainty",
        "Model":"model",
        "Airmass":"airmass"
    })

    return df[LC_COLUMNS].reset_index(drop=True)


def import_PSF_data(psf_path:str) -> pd.DataFrame:
    """
    Read EXOTIC PSF/centroid CSV by col index.

    Skip first line "#" header

    Trailing empty column (pos 7) is dropped via usecols.

    Returns: 
        pandas.DataFrame with cols matching PSF_COLUMNS.
    """
    path = Path(psf_path)

    df = pd.read_csv(
        path,
        comment="#", # drop first line
        header=None, # read by position
        usecols=range(len(PSF_COLUMNS)),
        names=PSF_COLUMNS,
    )

    return df.reset_index(drop=True)



# -- MERGE DATA --



def merge_datasets(lc_df:pd.DataFrame, psf_df:pd.DataFrame) -> pd.DataFrame:
    """
    Merge light-curve and PSF DataFrames by row pos.

    There's no index/key column for alignment here, so the dataframes must be the same length (x1 row/FITS frame)

    TODO:
    When the pipeline moves to direct PyKX/otherwise the load_* execs will be discarded and merge_datasets()
    will work on stored in-memory DataFrames from the data pipeline instead (after photometry processing).

    Return:
        pandas.DataFrame with cols in COMBINED_COLUMNS order.
    """
    if len(lc_df) != len(psf_df):
        raise ValueError("data mismatch!")

    merged_data = pd.concat([lc_df.reset_index(drop=True), psf_df.reset_index(drop=True)], axis=1) #

    return merged_data[COMBINED_COLUMNS].reset_index(drop=True)



def load_and_combine(lc_path:str, psf_path:str, verbose:bool=True) -> pd.DataFrame:
    """
    Load both EXOTIC CSV files from disk and return a combined DataFrame.
    
    Disk-based approach. When using PyKX/otherwise will replace any call to this func with merged_data = merge_datasets(lc_df, psf_df)
    where lc_df and psf_df are stored in-memroy

    Note on ARG: verbose prints a summary of the loaded data

    Return:
        Combined DataFrame, 13 columns, 1 row per FITS frame.
    """
    lc_df  = import_lightcurves(lc_path)
    psf_df = import_PSF_data(psf_path)
    merged_data = merge_datasets(lc_df, psf_df)

    if verbose:
        print_debug(merged_data)

    return merged_data



# -- DEBUG --



def print_debug(df:pd.DataFrame) -> None:
    print("="*15)
    print(f"Combined obserations: {len(df)} frames, {len(df.columns)} columns")
    print(f"BJD range: {df['BJD_TDB'].min():.6f} ... {df['BJD_TDB'].max():.6f}")
    print(f"Flux range: {df['flux'].min():.6f} ... {df['flux'].max():.6f}")
    print(f"Centroid X: {df['x_centroid'].min():.2f} ... {df['x_centroid'].max():.2f} px")
    print(f"Centroid Y: {df['y_centroid'].min():.2f} ... {df['y_centroid'].max():.2f} px")
    print(f"Offset range: {df['offset'].min():.2f} ... {df['offset'].max():.2f}")
    null_cols = df.columns[df.isna().any()].tolist()
    if null_cols:
        print(f"Columns with nulls: {null_cols}")
    else:
        print("No null values.")
    print("="*15)



# debug/test


def save_combination(df:pd.DataFrame, out_path:str):
    df.to_csv(out_path, index=False)

if __name__ == "__main__":
    import sys

    LC_PATH = r"C:/Users/dwint/Downloads/FinalLightCurve_HAT-P-32 b_2017-12-19.csv"
    PSF_PATH = r"C:/Users/dwint/Downloads/HAT-P-32b_psf_data_target.csv"
    OUT_PATH = r"C:/Users/dwint/Downloads/data_restructure_test.csv"
    df = load_and_combine(LC_PATH,PSF_PATH)
    save_combination(df,OUT_PATH)

    print("\nFirst 3 rows:")
    print(df.head(3).to_string())
    print("\nData types:")
    print(df.dtypes.to_string())
