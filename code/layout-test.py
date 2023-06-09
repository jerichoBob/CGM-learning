# my crap
import utils;
import layout_utils as lu;

import kcwitools.io as kio
import kcwitools.utils as ku
from kcwitools.image import build_whitelight
from astropy.io import fits
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, TextBox, RangeSlider, RectangleSelector

# mpl.rc('image', cmap='gray')
mpl.rc('image', cmap='RdYlBu_r')


base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429"
# flux_filename = base_path+"/J1429_rb_flux.fits"
flux_filename = base_path+"/SGAS1429+1202.fits"

f1 = fits.open(flux_filename)
hdr=f1[0].header
data=f1[0].data
wave=ku.build_wave(hdr)  # updated kcwitools.utils.build_wave to handle a missing CD3_3 from the header
wl_image=build_whitelight(hdr, data)

# our list of x,y locations on the source image to extract the wavelength
# points for J1429_rb_flux.fits
# x_cen = [30, 32, 34, 37, 35, 51];
# y_cen = [37, 38, 39, 39, 24, 32];
# points for SGAS1429+1202.fits
x_cen = [20, 22, 24, 27, 25, 40];
y_cen = [35, 36, 37, 37, 22, 30];

labels = ["1","2","3","4","5","6","7"];
alignment = [ # [ha, va] :: ha={'left', 'center', 'right'} va={'bottom', 'baseline', 'center', 'center_baseline', 'top'}
    ['left', 'top'],
    ['left', 'top'],
    ['center', 'top'],
    ['right', 'top'],
    ['center', 'top'],        
    ['center', 'top'],        
    ['center', 'top']
]; 

ax = lu.make_6speclayout(plt)

ax[1].imshow(wl_image,origin="lower",interpolation="nearest",cmap="afmhot",vmin=0);
colors = ['b','g','r','c','m','y','w'];
utils.plotcircle(ax[1], x_cen, y_cen, labels, 2, colors) # identify/label the locations on the plot

# accumulate the flux at a point
def acc_point():
    loc_flux = [];

def acc_circle(data, x_cen, y_cen, z_index, rad):
    # assuming rad = 1 for now (KISS)
    loc_flux = [];
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]]); # center
    loc_flux.append(data[i][y_cen[z_index]+1][x_cen[z_index]]); # north
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]+1]); # east
    loc_flux.append(data[i][y_cen[z_index]-1][x_cen[z_index]]); # south
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]-1]); # west
    return loc_flux


j = 0
for j in range(6):
    c = '#0f0f0f60' # one color, for now -- or colors[j];
    flux = []
    for i in range(wave.size):
        flux.append(acc_circle(data, x_cen, y_cen, j, 1));
        # flux.append(data[i][y_cen[j]][x_cen[j]]) # need to use kcwi tools to get the circular 
    ax[3+j].plot(wave, flux, '-', color=c)


plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95, wspace=0.129, hspace=0.17)
plt.show()


