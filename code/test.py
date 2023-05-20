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

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429"

# flux_filename = base_path+"/SGAS1429+1202.fits"
# flux_filename = base_path+"/J1429_flux.fits"
flux_filename = base_path+"/J1429_rb_flux.fits"

f1 = fits.open(flux_filename)
hdr=f1[0].header
data=f1[0].data
wave=ku.build_wave(hdr)  # updated kcwitools.utils.build_wave to handle a missing CD3_3 from the header
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
grid_size = (16,8)
fig = plt.figure(figsize=grid_size)

ax0 = plt.subplot2grid(grid_size, (1, 0), colspan=8, rowspan=4)
ax = [];
graph_count = 6 
graph_x = 4
graph_rows = 2
graph_cols = 4
for i in range(graph_count):
    print(i)
    ax.append(plt.subplot2grid(grid_size, (graph_rows*i, graph_x), colspan=graph_cols, rowspan=graph_rows))
    if i != graph_count-1: 
        ax[i].tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)

ax0.imshow(wl_image,origin="lower",interpolation="nearest",cmap="gnuplot",vmin=0)

# utils.plotbox(ax0, x_cen, y_cen, labels, alignment, 2, 'c') # identify/label the locations on the plot
utils.plotcircle(pax0, x_cen, y_cen, labels, 1, 'c')

colors = ['r','g','b','y'];
j = 0
for j in range(6):
    flux = []
    for i in range(wave.size):
        flux.append(data[i][y_cen[j]][x_cen[j]])
    ax[j].plot(wave, flux, '-', color='black')

plt.subplots_adjust(left=0, bottom=0.0, right=0.979, top=0.993, wspace=0.0, hspace=0.05)
plt.show()
