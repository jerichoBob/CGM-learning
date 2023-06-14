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
    (29,36), #0
    (29,39),    
    (32,37),
    (32,40),
    (35,38),
    (35,41),
    (38,38), 
    (38,41), 

    (49,31),
    (49,34),
    (52,31),

    (35,24),
]

x_coords, y_coords = zip(*points)
# convert to lists
x_coords = list(x_coords)
y_coords = list(y_coords)

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
flux_filename = base_path+"/J1429_rb_flux.fits"
var_filename = base_path+"/J1429_rb_var.fits"

# Load flux and variance data cubes
_, var = kcwi_io.open_kcwi_cube(var_filename)
hdr, flux = kcwi_io.open_kcwi_cube(flux_filename)

# first do a little data cleanup
var[np.isnan(flux)]=1.
flux[np.isnan(flux)]=0.0000
header =kp.tweak_header(deepcopy(hdr))

wave = kcwi_u.build_wave(hdr) # should I be using header instead of hdr??

# First create a white light image to plot
wl_center = 4686
wl_halfwidth  = 10
whitelight= im.build_whitelight(hdr, flux, minwave=wl_center - wl_halfwidth, maxwave=wl_center + wl_halfwidth)



fig_size=10

fig = plt.figure(figsize=(18, 12))
# ax = fig.add_subplot(111)
ax_info, ax_image, ax_spectra = lu.make_image_and_12_plots(plt)

co_begin=4750
co_end=5500

info_str = f'J1429 analysis metadata\n\nContinuum Start: {co_begin}\nContinuum Stop: {co_end}\nAperture Size: {sz}x{sz}'

ax_info.text(0.05, 0.9, 
             info_str,
             color='k', 
             fontsize = 12, 
             ha='left', va='top',
             transform=ax_info.transAxes)

ax_image.imshow(whitelight/1e4, interpolation='nearest',cmap=plt.get_cmap('gray'),origin="lower")#,vmin=1., vmax=8.)



spec_path = base_path+"/"

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
color_sightline = 'r'
color_error = 'gray'
color_snr = 'k'

for i in range(len(points)):
    x = x_coords[i]
    y = y_coords[i]
    # print("center: (", x,",",y,")")
    xx_A, yy_A, ff, vv = bu.corrected_corner_define(x, y, flux, var, deltax=sz, deltay=sz)

    # the bounding box and text label on the image
    ax_image.plot(xx_A, yy_A, color=color_sightline, linewidth=0.75)
    ax_image.text(x-1, y-1, str(i), color=color_sightline, 
                  fontsize = fig_size * 1.25, 
                  ha='center', va='center')
    filename = f'1d_spectra_{x}.{y}-{sz}x{sz}-{co_begin}-{co_end}.fits'
    filepath = spec_path + filename
    print("file: ",filename)
    # print(filepath)

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        sp=ke.extract_weighted_spectrum(ff, vv, wave, weights='Data', verbose=False)

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
    ax_spectra[i].text(0.05, 0.75, 
               "SNR: "+str(snr),
               color=color_snr, 
               fontsize = 12, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)    
    ax_spectra[i].text(0.05, 0.6, 
               f'({x},{y})',
               color=color_snr, 
               fontsize = 10, 
               ha='left', va='top',
               transform=ax_spectra[i].transAxes)   

    # fspec.plot()
    sp.write_to_fits(filepath)

    # kcwi_s.extract_square(x, y, wave, flux, var, radius, outfile=outfile_name)
    # spectra.append(kcwi_s.extract_circle(x, y, wave, flux, var, radius, outfile=outfile_name))
    # spectra.append(kcwi_s.extract_circle(x, y, wave, flux, var, radius))


plt.subplots_adjust(left=0.038, bottom=0.057, right=0.96, top=0.929, wspace=0.174, hspace=0.057)
# plt.show()
