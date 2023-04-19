import utils;
import sys

import kcwitools.io as kio
import kcwitools.utils as ku
from kcwitools.image import build_whitelight
from astropy.io import fits, ascii
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, TextBox, RangeSlider, RectangleSelector

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS data/J1429"

# flux_filename = base_path+"/SGAS1429+1202.fits"
# flux_filename = base_path+"/J1429_flux.fits"
flux_filename = base_path+"/J1429_rb_flux.fits"

f1 = fits.open(flux_filename)
hdr=f1[0].header

# print(flux_filename)
# def check_key(hdr, key):
#     if key in hdr:
#         print(key,"=",hdr[key])
#         return(True)
#     else:
#         print(key," not present")
#         return(False)

# cd_cdelta3 = -1; # this takes the place of CD3_3
# if check_key(hdr, 'CD3_3'):
#     cd_cdelta3 = hdr['CD3_3']
# elif check_key(hdr, 'CDELT3'):
#     cd_cdelta3 = hdr['CDELT3']
# else:
#     print("Can't read ", flux_filename)
#     print("aborting...")
#     sys.exit()

# print("cd_cdelta3: ", cd_cdelta3)
# crval3 = hdr['CRVAL3']
# crpix3 = hdr['CRPIX3']
# wavedim = hdr['NAXIS3']
# # Do it
# #wave = crval3 + (crpix3 + np.arange(0, wavedim, 1.0)) * cd3_3
# wave = crval3 + cd_cdelta3 * (np.arange(wavedim) + 1. - crpix3)
# print(cd_cdelta3)
# print(hdr.size())
data=f1[0].data
wave=ku.build_wave(hdr)
wl_image=build_whitelight(hdr, data)

# our list of x,y locations on the source image to extract the wavelength
x_cen = [30, 32, 34, 37, 35, 51];
y_cen = [37, 38, 39, 39, 24, 32];

labels = ["1","2","3","4","5","6"];
alignment = [ # [ha, va] :: ha={'left', 'center', 'right'} va={'bottom', 'baseline', 'center', 'center_baseline', 'top'}
    ['center', 'top'],
    ['center', 'top'],
    ['center', 'top'],
    ['center', 'top'],
    ['center', 'top'],        
    ['center', 'top']
]; 


fig = plt.figure()
grid_size = (12,7)
fig = plt.figure(figsize=(14,8))

ax0 = plt.subplot2grid(grid_size, (1, 0), colspan=3, rowspan=4)
ax = [];
ax.append(plt.subplot2grid(grid_size, (0, 3), colspan=4, rowspan=2))
ax.append(plt.subplot2grid(grid_size, (2, 3), colspan=4, rowspan=2))
ax.append(plt.subplot2grid(grid_size, (4, 3), colspan=4, rowspan=2))
ax.append(plt.subplot2grid(grid_size, (6, 3), colspan=4, rowspan=2))
ax.append(plt.subplot2grid(grid_size, (8, 3), colspan=4, rowspan=2))
ax.append(plt.subplot2grid(grid_size, (10, 3), colspan=4, rowspan=2))


ax[0].tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
ax[1].tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
ax[2].tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)


axnext = fig.add_axes([0.81, 0.05, 0.10, 0.03]) # left, bottom, width, height
# bnext = Button(axnext, 'Extract')
# bnext.on_clicked(extract_data)
# sslider.on_changed(update_sourceRange)
# cslider.on_changed(update_continuumRange)


ax0.imshow(wl_image,origin="lower",interpolation="nearest",cmap="gnuplot",vmin=0)
utils.plotbox(ax0, x_cen, y_cen, labels, alignment, 2, 'c')

colors = ['r','g','b','y'];
j = 0
for j in range(6):
    flux = []
    for i in range(wave.size):
        # print("i: ",i)
        # print("wavelength: ",wave[i], "flux: ",data[i][y][x])
        flux.append(data[i][y_cen[j]][x_cen[j]])
    ax[j].plot(wave, flux, '-', color='black')

plt.subplots_adjust(left=0, bottom=0.0, right=0.979, top=0.993, wspace=0.0, hspace=0.05)
plt.show()
