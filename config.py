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

#WORKING_DIR = "/Users/davidlawley/Downloads/Research_Project/photmetry_data/Hat_P_32_Dec202017.FITS"
# C:\Users\dwint\Downloads\images
WORKING_DIR = "C:/Users/dwint/Downloads"



# Name of the subfolder inside WORKING_DIR which contains the raw FITS images

# in target_star_photometry.py this was hard coded as "data" inside glob.glob("data/*.FITS")
# in the master notebook this was the placehoder "image folder filepath"

RAW_DATA_FOLDER = "aligned-images" # TODO ... why is this already aligned in source code


# Name the folder the aligned imgs will be saved into (need not previously exist, will be created automatically)

# in target_star_photometry.py this was hard coded as "aligned-images"
# in the notebook this was "name of new directory for aligned images"
# C:/Users/dwint/Downloads/aligned-images
ALIGNED_FOLDER = "aligned-images-out"




# ============== TARGET AND COMP STAR POSITIONS ==============

# PX positions for target and comp stars, targ first.

# ... "AAVSO + AIJ can help find good comp stars" ...

TARGET_POSITIONS = (
  (424.4, 286.8), 
  (348, 215.5), 
  (465, 182.6)
)




# ============== APERTURE PARAMS ==============

# Aperture radius. 

#EXOTIC gives optimised aperture radius in the output file.

APERTURE_R = 4.44



# Annulus inner radius.

#EXOTIC provides optimised again.

ANNULUS_R_IN = 7.15



# Annulus outer radius. 

#"r_out" is not provided by EXOTIC. "usually ~< 2x r_in is good".
# notebook target cell used 12, notebook bycatch cell used 12.15,
# bycatch-photometry.py used 12.15 — using 12.15 to match bycatch convention

ANNULUS_R_OUT = 12.15




# ============== SOURCE DETECTION ==============

# FWHM for DAOStarFinder source detection.

# "FWHM remains fairly const for all sources" ... "good approx using AIJ".

FWHM = 2.6




# ============== OUT-OF-TRANSIT (OOT) MASK(S) ==============

# The notebook uses raw MJD time values for its OOT mask.
# target_star_photometry.py uses normalised time starting from 0.
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

# Notebook used 7, .py used 6. Not overly important to distinguish, perhaps down to what the user prefers possible by visual inspection

BIN_SIZE_NOTEBOOK = 7
BIN_SIZE_SCRIPT = 6




# ============== ERROR PROPAGATION ==============

# e- read noise
SIGMA_READ = 5


# =============== STELLAR INDEX =================

# Chosen star in bycatch_photometry.py

STAR_INDEX = 177 #DEAN
