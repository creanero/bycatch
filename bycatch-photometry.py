# THE PRIMARY CONTRIBUTOR OF THIS WORK IS *D.LAWLEY* WITH ADAPTATIONS MADE FROM *DEAN WINTERS*
# - Dean Winters - - dean.winters4@mail.dcu.ie - - May 2026



# import libraries

import time
start_time = time.time()
import astroalign as aa
from astropy.io import fits
from pathlib import Path
import glob
#DEAN: from astropy.io import fits ... duplicate import commented out
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
#DEAN: from photutils.aperture import aperture_photometry  ... duplicate
from astropy.nddata import Cutout2D
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)
#DEAN: import matplotlib.pyplot as plt ... duplicate
from photutils.detection import DAOStarFinder
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats
import os

#DEAN: hyperparams sourced from config.py instead of being hardcoded in the individual files
from config import WORKING_DIR, ALIGNED_FOLDER, APERTURE_R, ANNULUS_R_IN, ANNULUS_R_OUT, FWHM

#DEAN: added os.chdir so ALIGNED_FOLDER resolves correctly as a relative path, like in target_star_photometry.py
os.chdir(WORKING_DIR)

# Load in the aligned fits images, if you do not have the images aligned they can be created using the more general "photutils-photometry" script or notebook
# Be sure to leave the /*.FITS at the end of the file path as this is used as an identifier when iterating through the folder
# It is also case sensitive so you might need to check if your .fits is caps or not
aligned_fits = sorted(glob.glob(ALIGNED_FOLDER + "/*.FITS")) #DEAN: was sorted(glob.glob("/Users/davidlawley/.../aligned-images/*.FITS"))

# detects sources present, using a reference image again to find sources
# the fwhm and std parameters are important here, the fwhm remains fairly constant for all the sources in the image and a good approximations can be made using tools in AIJ.
# std 5 is usually a safe bet for most fits images

# defines the reference image that will be used to detect sources present in the image stack
source_image = fits.getdata(aligned_fits[0])

# define the mean, median, and std as the values obtained from "sigma clipped stats", which is a photutils module
mean, median, std = sigma_clipped_stats(source_image, sigma=3)

# daofind is the module from photutils that is used to detect sources present in the images
# define the fwhm and threshold you would like to use.
daofind = DAOStarFinder(fwhm=FWHM, threshold=5.*std) #DEAN: was DAOStarFinder(fwhm=2.6, threshold=5.*std)
# subtract the mediaan from the "source_image" data. the vast majority of the data points in the image is taken up by background, so it is apropriate in this case to just use the median value here to background reduce.
sources = daofind(source_image - median)
# the loop below was copied from the phoutils user guide, it neatly formats the output table
for col in sources.colnames:
    if col not in ('id', 'npix'):
        sources[col].info.format = '%.2f'
sources.pprint(max_width=76)

# if you want to view the table of sources along with their positions and so on uncomment the line below
# sources.pprint(max_width=76)  #DEAN: bycatch-photometry.py had this commented out; notebook had it active ... keeping notebook version active above

# to perform aperture photometry on the sources positions need to be in column order as opposed to row.
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
# creates an aperture around all detected sources
apertures = CircularAperture(positions, r=4.32)
norm = ImageNormalize(stretch=SqrtStretch())

# plt.imshow(data, cmap='Greys', origin='lower', norm=norm,
#            interpolation='nearest')
# apertures.plot(color='blue', lw=1.5, alpha=0.5)


# uses the positions of the sources present to perform aperture photometry on all sources (bycatch photometry) the same way we did earlier for target and comp stars.

aperture = CircularAperture(positions, r=APERTURE_R) #DEAN: was CircularAperture(positions, r=4.44)
annulus_aperture = CircularAnnulus(positions, r_in=ANNULUS_R_IN, r_out=ANNULUS_R_OUT) #DEAN: was CircularAnnulus(positions, r_in=7.15, r_out=12.15)

# the Loop below works exactly the same way as it does when we are concerned with only the target and comp stars, except here we just iterate through all the sources present.
# for information on how this loop works, see the commented photutils photometry python file / notebook.
rows = []

for i, f in enumerate(aligned_fits, start=1):
    data, header = fits.getdata(f, header=True)
    time = header["MJD-OBS"] #DEAN: notebook uses "time" as variable name here; bycatch-photometry.py used "time_MJD" ... use notebook

    phot = aperture_photometry(data, aperture)
    fluxes = phot[("aperture_sum")] # single aperture
    aperstats = ApertureStats(data, annulus_aperture)

    bkg_median = aperstats.median
    bkg = bkg_median * aperture.area

    reduced_fluxes = fluxes - bkg

    rows.append((i, time, reduced_fluxes))

bycatch_table = Table(
    rows=rows,
    names=("id", "time", "fluxes") #DEAN: notebook uses "time" as variable name here; bycatch-photometry.py used "time_MJD" ... use notebook
)

print(bycatch_table)

flux_matrix = np.array(bycatch_table['fluxes'])

end_time = time.time()
print("Runtime:", end_time - start_time, "seconds")

star_flux = flux_matrix[:, 177] #DEAN: notebook used index 177; bycatch-photometry.py used 68 ... use notebook
plt.plot(bycatch_table["time"], star_flux, '.-') #DEAN: updated col name to match table above
plt.xlabel("Time")
plt.ylabel("Flux")
plt.show()
