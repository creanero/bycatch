# THE PRIMARY CONTRIBUTOR OF THIS WORK IS *D.LAWLEY* WITH ADAPTATIONS MADE FROM *DEAN WINTERS*
# - Dean Winters - - dean.winters4@mail.dcu.ie - - May 2026


# "-The vast majority of this code was built using the documentation from the photutils user guide
# -There is some very useful information in the Exotic output files (particularly the "temp" folder) that include optimised parameters,
# such as aperture size, annulus size, comp stars and so on, raw flux plots for comparison etc. "

# import libraries

import astroalign as aa
from astropy.io import fits
from pathlib import Path
import glob
#DEAN: from astropy.io import fits ...duplicate
from photutils.aperture import CircularAnnulus, CircularAperture
#DEAN: import matplotlib as plt ... duplicate and error-laden
import matplotlib.pyplot as plt
from photutils.aperture import ApertureStats
from photutils.aperture import aperture_photometry
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.table import Table
#DEAN: from photutils.aperture import aperture_photometry ... duplicate
from astropy.nddata import Cutout2D
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)
#DEAN: import matplotlib.pyplot as plt ... duplicate
from photutils.detection import DAOStarFinder
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats
import os
import time

#DEAN: pull settings from config.py instead of defining locally
from config import (RAW_DATA_FOLDER, TARGET_POSITIONS, ALIGNED_FOLDER,
                    APERTURE_R, ANNULUS_R_IN, ANNULUS_R_OUT,
                    SIGMA_READ,
                    OOT_START_NORM, OOT_END_NORM,
                    BIN_SIZE_SCRIPT,
                    WORKING_DIR)

#DEAN: no longer run this at top-level
#start_time = time.time()

#DEAN: no longer run this at top-level
# change working directory to directory containing fits folder
#os.chdir(WORKING_DIR) #DEAN: was os.chdir("/Users/xxx/Downloads/Research_Project/photmetry_data/Hat_P_32_Dec202017.FITS")

def align_images(RAW_DATA_FOLDER):

    # load in the raw fits files folder
    # this cell aligns images to a chosen reference image. all coordinates used should be based off this reference image
    # The raw image files can be viewed in AIJ as a stack to chose a clean reference file

    # Be sure to leave the /*.FITS at the end of the file path as this is used as an identifier when iterating through the folder
    # It is also case sensitive so you might need to check if your .fits is caps or not

    # iterates through image files, loads all images ending in .FITS to memory
    fits_files = sorted(glob.glob(RAW_DATA_FOLDER + "/*.FITS")) #DEAN: was sorted(glob.glob("data/*.FITS"))

    # gets the flux values and header data for the chosen fits file, this is your reference image
    reference = fits.getdata(fits_files[0])

    # makes a new fodler in your working directory for the aligned images to be saved to
    os.makedirs(ALIGNED_FOLDER, exist_ok=True) #DEAN: was os.makedirs("aligned-images", exist_ok=True)

    # iterate through the fits file folder
    for f in fits_files:

        # pulls out the flux values and header info for each fits file and stores them in memory
        data, header = fits.getdata(f, header=True)
        # data (flux values) is turned into an array, this is not essential but easier to work with
        data = np.asarray(data, dtype=np.float64)

        # this is an astropy module that aligns each frame in the loop based off the reference frame that we chose above
        aligned_data, footprint = aa.register(data, reference)

        # writes the each of the aligned files to the new folder we created above.
        base = os.path.basename(f)
        new_name = os.path.join(ALIGNED_FOLDER, base) #DEAN: was os.path.join("aligned-images", base)

        fits.writeto(new_name, aligned_data, header, overwrite=True)

        # after you have run this cell if you check your working directory you should see an additional folder with the new name you have chosen, containing the aligned images

        # this cell performs the photometry for just the target and comp stars
        # I have commented out the background gradient reduction in the loop, as it seems to add a slight bit of noise if anything. however you can uncomment it and play around with it.
        #

    # now we are iterating through the new aligned fits images we have jsut created, copy the complete filepath into the inverted commas, being sure to include the /*.FITS still.
    aligned_fits = sorted(glob.glob(ALIGNED_FOLDER + "/*.FITS")) #DEAN: was sorted(glob.glob("aligned-images/*.FITS"))

    return aligned_fits

#DEAN: debug align_images() func
#aligned_fits = align_images(RAW_DATA_FOLDER)

def perform_photometry(aligned_fits, TARGET_POSITIONS):

    # choose positions for your target and comp stars, putting the target star first. AAVSO chart finder used in conjunction with AIJ is extremely helpful for finding good comp stars
    positions = TARGET_POSITIONS #DEAN: was positions = ((424.4, 286.8), (348, 215.5), (465, 182.6))

    # define your aperture, positions are defined above so just chose a radius value. Exotic provides the optimised aperture radius in their output files, found in the "final parameters" file
    # or you can make a rough guess from observing the target star in AIJ
    aperture = CircularAperture(positions, r=APERTURE_R) #DEAN: was CircularAperture(positions, r=4.44)

    # define annulus aperture, this is required for the backgorund reduction. Exotic provides the  optimised "r_in" parameter but not "r_out", usually something a little less than double r_in is good.
    annulus_aperture = CircularAnnulus(positions, r_in=ANNULUS_R_IN, r_out=ANNULUS_R_OUT) #DEAN: was CircularAnnulus(positions, r_in=7.15, r_out=12)

    # the variable rows is created to append the desired values later in the loop.
    rows = []

    # iterate through aligned fits files
    for i, f in enumerate(aligned_fits, start=1):

        # grabs the data (flux values) and header data
        data, header = fits.getdata(f, header=True)
        # time is defined as the "MJD-OBS" value from the header, there are a few different astronomical time systems present in the headers that you can chose from
        time_MJD = header["MJD-OBS"]
        # sigma_clip = SigmaClip(sigma=3.0)
        # bkg_estimator = MedianBackground()
        # bkg = Background2D(data, (50, 50), filter_size=(3, 3),
        #                sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        # reduced_data = data - bkg.background

        # "aperture_photometry" is a photutils module that is used to perform our aperture photometry
        phot = aperture_photometry(data, aperture)
        # we define the fluxes as the value from the "aperture_sum" column from our photometry table
        fluxes = phot[("aperture_sum")]  # single aperture
        # aperstats is another photutils module, that provides us with statistics from our apertures such as mean, median ,std etc. The full list of stats that can be called can be found in the user guide
        aperstats = ApertureStats(data, annulus_aperture)

        # we now have multiple flux values for each iteration, due to having multiple positions defined. I suggest defining the target star position first for this reason
        # Target flux defined as the first flux value
        T_flux = fluxes[0]
    
        # comp fluxes are defined as the rest of the flux values that are present in the table
        comp_fluxes = fluxes[1:]

        # using the aperstats module we find the median background value from the annulus.
        bkg_median = aperstats.median #DEAN: notebook uses bkg_median; target_star_photometry.py used annulus_median ... kept notebook naming
        # Multiply the median value by our aperture area to get the total background count within our aperture. this background method works well because it is very local to the Target star
        bkg_aperture = bkg_median * aperture.area
        # Each background is deinfed accordingly like it was done for the fluxes
        T_bkg = bkg_aperture[0]
        comp_bkg = bkg_aperture[1:]

        # subtract the background from the flux counts
        T_reduced_flux = T_flux - T_bkg
        comp_reduced_fluxes = comp_fluxes - comp_bkg

        # this is where the comp stars are defined, it saves a few lines of extra code by doing this at the end
        # if you have defined more than two comp stars you will of course need to define these here
        C1_reduced_flux = comp_reduced_fluxes[0]
        C2_reduced_flux = comp_reduced_fluxes[1]

        # this is a useful check to see the quality of your chosen comparison stars. the less noise in this plot the better, ideally you want the comp stars plot to roughly match eachother aswell as the target star.
        comp_star_check = C1_reduced_flux / C2_reduced_flux
        # the reference flux is used for the differential photometry. it is the sum of the comp stars flux
        reference_flux = np.sum(comp_reduced_fluxes)
        # the differential flux is the target star flux divided by the reference flux. This helps enormously with reducing any atmosherpic/seeing effects that are being produced across the image
        differential_T_flux = T_reduced_flux / reference_flux

        #DEAN: error calculations below are from target_star_photometry.py — not in the notebook, added here as an extension
        aperture_sum = fluxes - bkg_aperture  # reduced fluxes for all positions
        n_pix = aperture.area
        sigma_read = SIGMA_READ #DEAN: was sigma_read = 5
        std_ann = aperstats.std
        n_ann = annulus_aperture.area

        # error calculations
        var_ap = aperture_sum + n_pix * sigma_read**2
        var_bkg = (n_pix**2) * (std_ann**2 / n_ann)
        var_total = var_ap + var_bkg

        sigma_e = np.sqrt(var_total)

        sigma_e_T = sigma_e[0]
        sigma_e_C1 = sigma_e[1]
        sigma_e_C2 = sigma_e[2]
        comp_errors = sigma_e[:1] #DEAN: this line only captures the first comp star error ... sigma_e[1:] would capture all comp stars ... perhaps this was intended?
        #comp_errors = sigma_e[1:] #DEAN: try this TODO

        sigma_comp = np.sqrt(np.sum(np.array(comp_errors)**2, axis=0))

        rel_err = np.sqrt(
            (sigma_e_T / T_reduced_flux)**2 +
            (sigma_comp / reference_flux)**2
        )

        sigma_diff = differential_T_flux * rel_err

        # we now append the deisred values to rows. this creates a list of rows, where each row contains a value that is defined below
        # you can add or remove variables below
        #DEAN: COL set follows the notebook (C1, C2, comp_star_check, reference_flux kept in)
        # error columns from target_star_photometry.py added on top
        rows.append((i, time_MJD, T_reduced_flux, differential_T_flux, C1_reduced_flux, C2_reduced_flux, comp_star_check, reference_flux, sigma_e_T, sigma_diff))

    # convert rows into a table, makes it easier to query
    lc_table = Table(
        rows=rows,
        names=("id", "time", "T_flux", "diff_T_flux", "C1", "C2", "comp_star_check", "reference_flux", "t_error", "diff_error")
    )


    # plot desired variables from the table created
    plt.plot(lc_table["time"], lc_table["diff_T_flux"], marker='.', color='gray')
    plt.xlabel("MJD")
    plt.ylabel("flux")
    plt.show()

    return lc_table

#DEAN: debug perform_photometry()
#lc_table = perform_photometry(TARGET_POSITIONS)

def normalise_light_curves(lc_table):
        
    # this cell normailises the differential flux, so that the baseline is roughly one.
    # the time parameters are the ingress and egress times, Exotic output files contain the mid transit time so the ingress and egress are rough approximations based off this.

    #DEAN: t now defined to begin at time=0 and count up to make it easier when creating time masks
    # this is from target_star_photometry.py ... the notebook used raw MJD for its mask, the .py shift-started at time=0
    t = lc_table["time"] - np.min(lc_table["time"])
    lc_table["norm_t"] = t

    # define ingress and egress times respectively
    oot_mask = (lc_table["norm_t"] < OOT_START_NORM) | (lc_table["norm_t"] > OOT_END_NORM)  #DEAN: was (lc_table["norm_t"] < 0.06) | (lc_table["norm_t"] > 0.2)
    # uses the mask to normailise the baseline time so that it is now roughly 1
    baseline = np.median(lc_table["diff_T_flux"][oot_mask])
    # defines the normailised baseline as "norm_flux"
    norm_flux = lc_table["diff_T_flux"] / baseline
    # adds these values to the table
    lc_table["norm_flux"] = norm_flux

    #DEAN: error normalisation from target_star_photometry.py
    norm_flux_error = lc_table["diff_error"] / baseline
    lc_table["norm_flux_error"] = norm_flux_error


    return lc_table



def plot_target_star_photometry(lc_table):

    #DEAN: initialise the mask
    oot_mask = (lc_table["norm_t"] < OOT_START_NORM) | (lc_table["norm_t"] > OOT_END_NORM)

    #DEAN: plot of the normalised flux
    plt.plot(lc_table["norm_t"], lc_table["norm_flux"], marker='.', color='gray')
    plt.xlabel("MJD")
    plt.ylabel("normailised flux")
    plt.show()

    #DEAN: error bar plot with RMS and SNR annotation from target_star_photometry.py
    fig1, ax = plt.subplots()

    ax.errorbar(lc_table["norm_t"],
                lc_table["norm_flux"],
                yerr=lc_table["norm_flux_error"],
                marker='o',
                color='blue')

    ax.set_xlabel("MJD")
    ax.set_ylabel("normalised flux")

    oot_flux = lc_table["norm_flux"][oot_mask]
    rms = np.std(oot_flux)

    it_mask = (lc_table["norm_t"] > OOT_START_NORM) & (lc_table["norm_t"] < OOT_END_NORM)  #DEAN: was (lc_table["norm_t"] > 0.06) & (lc_table["norm_t"] < 0.2)
    transit_depth = 1 - np.median(lc_table["norm_flux"][it_mask])
    scatter = np.std(lc_table["norm_flux"][oot_mask])
    N = np.sum(it_mask)
    SNR = transit_depth / scatter * np.sqrt(N)

    info_text = f"""RMS = {rms:.5f}
    SNR = {SNR:.2f}"""

    ax.text(0.02, 0.98, info_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round',
                    facecolor='white',
                    alpha=0.8))

    plt.show()

    # bins the data and plots average of binned points, useful to see the overall shape of the curve

    # I have chosen a bin size of 14 so that i have 10 binned data points. the exact bin size is not important, but the best bin size would be the one that shows the shape of the trasnit the best
    bin_size = BIN_SIZE_SCRIPT  #DEAN: was bin_size = 6 in target_star_photometry.py, 7 in the notebook ... using script value here, see config.py

    lc_time = lc_table["norm_t"]
    flux = lc_table["norm_flux"]

    # divides the number of time values (frames) by the bin size, to give the number of bins
    n_bins = len(lc_time) // bin_size

    # create a list for the binned time and flux values
    binned_time = []
    binned_flux = []

    for i in range(n_bins):
        start = i * bin_size
        end = start + bin_size

        # avergaes the time values contained in each bin
        binned_time.append(np.mean(lc_time[start:end]))
        # averages the flux values conatined in each bin
        binned_flux.append(np.mean(flux[start:end]))

    binned_time = np.array(binned_time)
    binned_flux = np.array(binned_flux)

    plt.plot(binned_time, binned_flux, marker='s')
    plt.scatter(lc_table["norm_t"], lc_table["norm_flux"], marker='.', color='gray')

    plt.xlabel("time")
    plt.ylabel("T_flux")
    plt.show()

    oot_flux = lc_table["norm_flux"][oot_mask]
    rms = np.std(oot_flux)
    print(rms)

def execute_target_star_photometry():

    # change working directory to directory containing fits folder
    os.chdir(WORKING_DIR) #DEAN: was os.chdir("/Users/xxx/Downloads/Research_Project/photmetry_data/Hat_P_32_Dec202017.FITS")

    aligned_fits = align_images(RAW_DATA_FOLDER)
    lc_table     = perform_photometry(aligned_fits, TARGET_POSITIONS)
    lc_table     = normalise_light_curves(lc_table)
    plot_target_star_photometry(lc_table)


if __name__ == "__main__":
    start_time = time.time()
    execute_target_star_photometry()
    print("Runtime:", time.time() - start_time, "seconds")


#end_time = time.time()
#print("Runtime:", end_time - start_time, "seconds")
