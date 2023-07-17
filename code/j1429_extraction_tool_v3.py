import bobutils.layout_utils as lu
import bobutils.utils as bu

from kcwitools import io as kcwi_io
from kcwitools import utils as kcwi_u
from kcwitools import spec as kcwi_s
from kcwitools import plot as kp
from kcwitools import image as im
from kcwitools import extract_weighted_spectrum as ke

import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
import warnings


# our new sightlines identified using QFitsView using square 3x3 aperture
# SNR was waaay too low for a 3x3 (mostly less than 1) when treating 3700-5500 as the continuum
sz = 3  # length of one side of the square box
points = [
    (29,36,3), #0
    (29,39,3),    
    (32,37,3),
    (32,40,3),
    (35,39,3),
    # (35,41,3),
    (38,39,3), 
    # (38,41,3), 

    # (49,31,3),
    # (49,34,3),
    # (52,31,3),

    (50,33,7),
    
    (35,24,7),
]

x_coords, y_coords, sizes = zip(*points)
# convert to lists
x_coords = list(x_coords)
y_coords = list(y_coords)
sizes    = list(sizes)

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
flux_filename = base_path+"/J1429_rb_flux.fits"
var_filename = base_path+"/J1429_rb_var.fits"

# Load flux and variance data cubes
_, var = kcwi_io.open_kcwi_cube(var_filename)
hdr, flux = kcwi_io.open_kcwi_cube(flux_filename)

# first do a little data cleanup
var[np.isnan(flux)]=1.
flux[np.isnan(flux)]=0.0000
header =kp.tweak_header(deepcopy(hdr)) #why do this if it is never used?

wave = kcwi_u.build_wave(hdr) # should I be using header instead of hdr??

# First create a narrow-band whitelight image to plot
narrowband_center = 4686
del_lambda  = 10
whitelight= im.build_whitelight(hdr, flux, minwave=narrowband_center - del_lambda, maxwave=narrowband_center + del_lambda)



fig_size=8

fig = plt.figure(figsize=(18, 10))
# plt.rcParams['text.usetex'] = True
# plt.rcParams['font.family'] = 'serif'
# ax = fig.add_subplot(111)
ax_info, ax_image, ax_spectra = lu.make_image_and_12_plots(plt)

# the continuum, place in a rlatively stable region of the spectra
co_begin=4864
co_end=4914

info_str = f'J1429 Analysis Metadata\n\nContinuum: {co_begin} - {co_end}\nAperture Size: Position Dependent\nNarrowband Center: {narrowband_center} +/- {del_lambda} Angstroms\n\n'

ax_info.text(0.05, 0.9, 
             info_str,
             color='k', 
             fontsize = 12, 
             ha='left', va='top',
             transform=ax_info.transAxes)

ax_image.imshow(whitelight/1e4, interpolation='nearest',cmap=plt.get_cmap('gray'),origin="lower")#,vmin=1., vmax=8.)



font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
color_sightline = 'r'
color_error = 'gray'
color_snr = 'k'

output_path = "./analysis/j1429/"
for i in range(len(points)):
    x = x_coords[i]
    y = y_coords[i]
    sz = sizes[i]
    # print("center: (", x,",",y,")")
    xx_A, yy_A, ff, vv = bu.corrected_corner_define(x, y, flux, var, deltax=sz, deltay=sz)

    # the bounding box and text label on the image
    ax_image.plot(xx_A, yy_A, color=color_sightline, linewidth=0.75)
    ax_image.text(x-1, y-1, str(i), color=color_sightline, 
                  fontsize = fig_size * 1.25, 
                  ha='center', va='center')
    filename = f'{i}_1d_spectra_{x}.{y}-{sz}x{sz}-{co_begin}-{co_end}.fits'
    filepath = output_path + filename
    print("file: ",filename)

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        sp=ke.extract_weighted_spectrum(ff, vv, wave, weights='Data', verbose=False)
        # sp=kcwi_s.extract_square(x, y, wave, flux, var, squaresize=sz, outfile=None)


        #plot the extracted spectrum before writing it out to file.
    ax_spectra[i].plot(sp.wavelength, sp.flux, '-', color=color_sightline)
    ax_spectra[i].plot(sp.wavelength, sp.sig, '-', color=color_error, linewidth=0.5)
    snr = bu.sig_figs(bu.signal_to_noise(sp.wavelength, sp.flux, co_begin=co_begin, co_end=co_end), 5)
    print("SNR: ",snr)

    ax_spectra[i].text(0.05, 0.9, 
               "Sightline: "+str(i),
               color=color_sightline, 
               fontsize = 12, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)
    ax_spectra[i].text(0.05, 0.8, 
               "SNR: "+str(snr),
               color=color_snr, 
               fontsize = 12, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)    
    ax_spectra[i].text(0.05, 0.7, 
               f'Coord: ({x},{y})',
               color=color_snr, 
               fontsize = 12, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)   
    ax_spectra[i].text(0.05, 0.6, 
               f'Aperture: {sz} x {sz}',
               color=color_snr, 
               fontsize = 12, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)  
    sp.write_to_fits(filepath)


plt.subplots_adjust(left=0.038, bottom=0.057, right=0.96, top=0.929, wspace=0.174, hspace=0.041)
plt.show()
