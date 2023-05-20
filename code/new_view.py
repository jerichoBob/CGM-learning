# read in the flux and var data
# show the flux image 
# plot them with common wavelength axis so you can see them both together

# my crap
import utils as u;
import layout_utils as lu;

import kcwitools.io as kio
import kcwitools.utils as ku
from kcwitools.image import build_whitelight
from astropy.io import fits, ascii
from astropy.nddata import Cutout2D
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, TextBox, RangeSlider, RectangleSelector

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429"
flux_fn = base_path+"/J1429_rb_flux.fits"
var_fn = base_path+"/J1429_rb_var.fits"
# flux_file = fits.open(base_path+"/SGAS1429+1202.fits") # J1429_flux
# fluxfile = fits.open(flux_fn) # J1429_flux
# varfile = fits.open(var_fn) # J1429_rb_var.fits

#normally a (55,64) x,y data cube
# ku.trim_kcwi_cube_with_header(flux_fn, var_fn, size = (66, 25), position = (16.5, 47.5)) # works
# this rb cube is whacked, so trim out the bad stuff
newflux, newvar, newfhdr, newvhdr = ku.trim_kcwi_cube_with_header(flux_fn, var_fn, size = (64, 54), position = (35.5, 31.5))


# flux_hdr=fluxfile[0].header
# for i in newfhdr:
#     print(i,": ", newfhdr[i])
# flux_data=fluxfile[0].data
# flux_wave=ku.build_wave(flux_hdr)

# wavedim, ydim, xdim = newflux.shape
# print("wavedim: ", wavedim)
# print("ydim: ", ydim)
# print("xdim: ", xdim)
# image = flux_data[1149]

# position = (5.5, 5.5)
# size = (10, 10) # ny,nx
# cutout = Cutout2D(image, size = (65, 54), position = (35.5, 31.5))
ax = lu.make_image_flux_var_layout(plt)
wl_image=build_whitelight(newfhdr, newflux)
wave=ku.build_wave(newfhdr) # our trust, air2vac adjusted wavelength array

ax[1].imshow(wl_image,origin="lower",interpolation="nearest",cmap="afmhot",vmin=0)
c = '#0f0f0f60' # one color, for now -- or colors[j];
flux = []
for i in range(wave.size):
    flux.append(newflux[i][38][27]); # center

    # flux.append(newflux[i])
    # flux.append(u.acc_circle(newflux, x_cen, y_cen, i, 1))
    # flux.append(data[i][y_cen[j]][x_cen[j]]) # need to use kcwi tools to get the circular 
ax[3].plot(wave, flux, '-', color=c)

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95, wspace=0.129, hspace=0.17)
plt.show()