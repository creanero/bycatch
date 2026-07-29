# THE PRIMARY CONTRIBUTOR OF THIS WORK IS *D.LAWLEY* WITH ADAPTATIONS MADE FROM *DEAN WINTERS*
# - Dean Winters - - dean.winters4@mail.dcu.ie - - May 2026

# This is the configuration file in which necessary paths and parameters used across the pipeline are hardcoded. 
# This is an amalgamation of work of DL without much further inference/contribution. Some notes in quotations "" about choice of parameter are given and are
# mostly derived from what is mentioned in GoogleDrive and the main commented notebook; these are not direct quotes from any one person and are currently there to note for
# future use cases



# Edit the values in this file before running the pipeline.
# Each script imports from here so settings here are global.
# The notebook has its own small config cell that mirrors this.




# ============== Paths ==============

# Path to the directory containing the raw FITS img folder

# in target_star_photometry.py this was os.chdir("/Users/davidlawley/...")
# in bycatch-photometry.py this was embedded in the glob.glob() path directly
# in the notebook this was the placeholder string "working directory filepath"

# Process currently reads from disk
WORKING_DIR = "C:/Users/dwint/Downloads"



# Name of the subfolder inside WORKING_DIR which contains the raw FITS images

# in original target_star_photometry.py work this was hard coded as "data" inside glob.glob("data/*.FITS")
# in the master notebook this was the placehoder "image folder filepath"

RAW_DATA_FOLDER = "aligned-images" # TODO ... this is already aligned in source code


# Name the folder the aligned imgs will be saved into (need not previously exist, will be created automatically)

# in target_star_photometry.py this was hard coded as "aligned-images"
# in the notebook this was "name of new directory for aligned images"

ALIGNED_FOLDER = "aligned-images-out"




# ============== TARGET AND COMP STAR POSITIONS ==============

# PX positions for target and comp stars, targ first.

# ... "AAVSO + AIJ can help find good comp stars" ...

TARGET_POSITIONS = (
  (424.4, 286.8), 
  (348, 215.5), 
  (465, 182.6),
)





# ============== APERTURE PARAMS ==============

# FWHM, aperture radius, and annulus radii are derived per-run based on the reference frame
# using fwhm_calibration.py instead of being hardcoded via AIJ/EXOTIC.

# Bootstrap initial guess
FWHM_INIT = 3.0

# Aperture radius = APERTURE_COEFF * aggregated FWHM
# This is the dominant driver of SNR (=1.5 vs =2.6 on HATP32b doubled the SNR) 
# TODO: a proper automated way to derive this from data would need to avoid 
# re-reading/re-processing the full aligned stack per candidate value
APERTURE_COEFF = 2.6

# Annulus radii are defined relative to the winning aperture radius, not as
# independent ratios of FWHM: much less sensitive than APERTURE_COEFF - they just
# need to (a) clear the target's PSF wings so the background estimate isn't
# biased by stellar flux, and (b) span enough pixels for a stable median.
#   annulus_r_in  = aperture_r + ANNULUS_GAP_COEFF * measured_FWHM
#   annulus_r_out = annulus_r_in + ANNULUS_WIDTH_COEFF * measured_FWHM
ANNULUS_GAP_COEFF = 1.55
ANNULUS_WIDTH_COEFF = 2.9




# ============== OUT-OF-TRANSIT (OOT) MASK(S) ==============

# The notebook uses raw MJD time values for its OOT mask.
# original target_star_photometry.py work uses normalised time starting from 0.
# Both are seperately kept here since they are/may be intentionally different.



# Notebook OOT mask 

# raw MJD values — ingress and egress times
# EXOTIC output files have mid transit time... approximated ingress/egress 

OOT_START_MJD = 0.34
OOT_END_MJD = 0.46



# target_star_photometry.py 

# normalised time, shifted to start from 0

OOT_START_NORM = 0.06
OOT_END_NORM = 0.2




# ============== BINS ==============

# Bin sizes for the light curves.

# Notebook used 7, original .py work used 6. Not overly important to distinguish, perhaps down to what the user prefers possible by visual inspection

BIN_SIZE_NOTEBOOK = 7
BIN_SIZE_SCRIPT = 6




# ============== ERROR PROPAGATION ==============

# e- read noise
SIGMA_READ = 5


# =============== STELLAR INDEX =================

# Chosen star in original bycatch_photometry.py work

STAR_INDEX = 177 #DEAN
