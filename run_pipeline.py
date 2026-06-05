# THE PRIMARY CONTRIBUTOR OF THIS WORK IS *D.LAWLEY* WITH ADAPTATIONS MADE FROM *DEAN WINTERS*
# - Dean Winters - - dean.winters4@mail.dcu.ie - - May 2026


# Run the pipeline end to end

# Order of execution:


# 1. target_star_photometry.py
#    aligns the FITS imgs and saves to ALIGNED_FOLDER
#    perform aperture photometry on the target and comp stars
#    plot raw flux, norm. flux, and error plots (RMS/SNR)
#    binned curves

# 2. bycatch-photometry.py
#    reads ALIGNED_FOLDER
#    detect sources in the field using DOAStarFinder
#    perform aperture photometry on detections
#    plots flux of a particular detected source


# The notebook has the steps of these scripts annotated and describes the pipeline in more depth

import time
import target_star_photometry
import bycatch_photometry #DEAN: need to rename bycatch-photometry to remove the hyphen otherwise this cannot import

print("Phase 1: Target star photometry and alignment:")
target_star_photometry.execute_target_star_photometry()

print("\n Phase 2: Bycatch photometry:")
bycatch_photometry.execute_bycatch_photometry()

print("\n=== Pipeline complete ===")
