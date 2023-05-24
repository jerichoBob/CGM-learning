from kcwitools import io as kcwi_io
from kcwitools import spec as kcwi_s
from kcwitools import utils as kcwi_u
import numpy as np

points = [
    (30,37),
    (32,38),
    (34,39),
    (37,39),
    (51,32),
    (35,24)
]

x_coords, y_coords = zip(*points)
# convert to lists
x_coords = list(x_coords)
y_coords = list(y_coords)

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429"
flux_filename = base_path+"/J1429_rb_flux.fits"
var_filename = base_path+"/J1429_rb_var.fits"

# Load
hdr, flux = kcwi_io.open_kcwi_cube(flux_filename)
wave = kcwi_u.build_wave(hdr)

# our var file is busted right now (2023.05.15) - 
# trying to figure out a good substitute so we don't wind up with a whacky spectrum 
# _, var = kcwi_io.open_kcwi_cube(var_filename)
# var = 0.01 * flux
var = np.copy(flux)

spec_basename = base_path+"/"+"spectrum_"
# Spectrum
# print("len: ", len(points))
radius = 1
spectra = []
for i in range(len(points)):
    x = x_coords[i]
    y = y_coords[i]

    outfile_name = spec_basename + str(x)+"."+str(y)+".fits"
    # print("outfile_name: ",outfile_name)
    # spectra.append(kcwi_s.extract_circle(x, y, wave, flux, var, radius, outfile=outfile_name))
    spectra.append(kcwi_s.extract_circle(x, y, wave, flux, var, radius))

