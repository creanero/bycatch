# Measure FWHM from actual image data
# Uses IRAFStarFinder to detect sources (within 5 sigmas),
# then fit a Gaussian to the N brightest cases.
# Take FWHM_final as the median of each source fitted FWHM
# Uses some outlier rejection (sigma_clip; sigma=3)

# Aperture/annuli vals are then derived as ratios of FWHM (still lives in config.py)

import numpy as np
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils.detection import IRAFStarFinder


def determine_fwhm(imgs, fwhm_init, threshold_sigma=5.0, n_stars=20):
    # TODO: this function currently fits the Gaussian to every source in the image data and only does
    # the top N and clipping rejections after the fact. Might be more data-efficient to do a cheap
    # source detection pass then fit N number of Gaussians only.

    # Find background stats to pass to finder

    mean, median, std = sigma_clipped_stats(imgs, sigma=3.0)

    # IRAFStarFinder is used (unlike eg: DAOStarFinder) since it fits a
    # Gaussian to each detection and report fwhm per source. 
    # Before, DAO was used and yielded one FWHM for the entire image.

    # threshold_sigma is the source threshold above background.
    # threshold_sigma=5.0 is just an estimate good starter for this. 
    # TODO: verify if this needs to be exported to configpy or stay hardcoded here

    finder = IRAFStarFinder(threshold = threshold_sigma*std, fwhm=fwhm_init, exclude_border=True)

    # background subtraction
    sources = finder(imgs - median)

    # For FWHM calc, will use the top N brightest source cases since dimmer sources (lower SNR)
    # will probably give a more noisy FWHM fit that would also be high variance.

    # another inescapble param here is n_stars, the top N=n amount to use
    # TODO: verify if this needs to be exported to configpy or stay hardcoded here

    sources.sort("flux")
    sample = sources[-n_stars:] if len(sources) > n_stars else sources
    raw_fwhms = np.array(sample["fwhm"])

    # This might not actually be necessary but the FWHM vals are then themselves sigma clipped
    # in case of non-stellar detection or some other kind of outlier/error
    fwhms = sigma_clip(raw_fwhms, sigma=3.0)

    # Debug printout of the number of sources, per-source FWHM and how many outliers were rejected.
    print(f"determine_fwhm: {len(sources)} sources detected, using brightest {n_stars}")
    print(f"per-star fwhm (px): {np.array2string(raw_fwhms, precision=2)}") # string formatting error. use array2string
    print(f"clipped {np.ma.count_masked(fwhms)} outlier(s) before taking median")

    # Median considered a fairer average for single peak distrubtions
    return float(np.ma.median(fwhms))


def derive_apertures(fwhm, aperture_coeff, annulus_gap_coeff, annulus_width_coeff):

    # Aperture radius is the dominant lever for SNR - aperture_coeff is
    # currently a manually-chosen constant in config.py (matched against
    # EXOTIC's own optimised radius for this dataset via compare_calibration.py)
    aperture_r = aperture_coeff * fwhm

    # Define annulus relative to fwhm-scaled aperture, not from fwhm directly
    # annulus_gap_coeff needs to go past PSF tails so theres no residual flux bias
    # annulus_width_coeff just needs enough pixels for a stable median
    annulus_r_in = aperture_r + annulus_gap_coeff*fwhm
    annulus_r_out = annulus_r_in + annulus_width_coeff*fwhm

    return {
        "aperture_r": aperture_r,
        "annulus_r_in": annulus_r_in,
        "annulus_r_out": annulus_r_out
    }