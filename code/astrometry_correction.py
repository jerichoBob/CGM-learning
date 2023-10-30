descr = """
Path: code/astrometry_correction.py
This code will confirm/reject that the files that I am using 
for the 7 separate observations of J1429+1202 
have been correctly astrometry corrected :D.
"""
import bobutils.utils as bu

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from photutils.centroids import centroid_sources, centroid_com
from astropy import units as u
from astropy.coordinates import SkyCoord

cmap = 'gnuplot'


def get_hst_data():
    fh = fits.open("/Users/robertseaton/School/github_repos/CGM-learning/data/data_idpw02010_WFC3_IR_F160W_drz.fits")
    hdr = fh[0].header    # Reference image header
    data = fh[0].data     # Reference image 2D data
    wcs = WCS(hdr)        # WCS of the reference image
    return hdr, data, wcs


def show_whitelight_image(ax, wl_image):
    ax.imshow(wl_image, origin='lower', interpolation='nearest', cmap=cmap, vmin=0, vmax=0.3*wl_image.max())
    ax.set_title("White Light from the data cube")

    # %matplotlib qt
    # %matplotlib inline

    fig = plt.figure(1, figsize=(5,5))
    ax1 = fig.add_subplot(111)
    show_whitelight_image(ax1, wl_k)
    fig.tight_layout()
    plt.show()

def main():
    hdr_hst, data_hst, wcs_hst = get_hst_data()
    observations = bu.get_corrected_kcwi_data()
    for ii,o in enumerate(observations):
        if ii > 0: break
        print(f"=======================  Observation {ii+1}  =============================")
        hdr_f, hdr_v, flux, var, wcs_f, wcs_v, wave, wl_k = o.get_all()
        # print("="*40)
        # for k in hdr_f.keys():
        #     print(f"   {k} = {hdr_f[k]}")

if __name__ == "__main__":
    main()