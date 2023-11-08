description = """
File I/O utilities - especially when dealing with FITS files.
"""
import os, sys
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools import image as im
from kcwitools import utils as kcwi_u
from kcwitools import io as kcwi_io

# modify the python path to include the current directory
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from observations import Observation

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
reference_flux_filename = f"{base_path}/J1429+1202_KCWI_corrected_flux.fits"

def find_fits(basedir, endswith):
    return (os.path.join(root, file)
        for root, dirs, files in os.walk(basedir)
            for file in files if root == basedir and file.lower().endswith(endswith))

def find_files_ending_with(dir, endswith):
    # List all files in the specified directory
    files = os.listdir(dir)
    # Filter the list to only include files that end with the specified suffix
    return [os.path.join(dir, file) for file in files if file.endswith(endswith)]

def read_and_prep_flux_var_data(flux_file, var_file, minwave, maxwave, ybot, ytop):
    """ This method reads the flux and var cubes for a specific observation, and cleans them up before returning a new Observation object"""
    hdr_f, flux = open_kcwi_cube(flux_file)
    hdr_v, var = open_kcwi_cube(var_file)
    wave = build_wave(hdr_f)
    # print(f"flux file: {flux_file}")
    # print(f"flux.shape={flux.shape}") # (lambda, y, x)
    # print(f"var file: {var_file}")
    # print(f"var.shape={var.shape}") # (lambda, y, x)

    # first do a little data cleanup
    var[np.isnan(var)]=1.
    flux[np.isnan(flux)] = 0.0000

    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]
    wave = wave[slices]
    flux = flux[slices,:,:] 
    var  =  var[slices,:,:] 
    
    var_is_stddev = False # set this to true if we need to interpret the variance data as standard deviation data
    if var_is_stddev:
        stddev = var
        var = stddev**2
    else:
        stddev = np.sqrt(var)
    
    # print(f"hdr_v: {hdr_v}")

    # if we are going to crop the data, we need to modify the header so that everything tracks correctly
    # flux = flux[slices,int(ybot):int(ytop),:] # NOTE: Cropping is busted.  Need to fix this.
    # var  =  var[slices,int(ybot):int(ytop),:] # NOTE: Cropping is busted.  Need to fix this.
    # print(f"cropped flux.shape={flux.shape}") # (slices, y, x)
    
    return Observation(hdr_f=hdr_f, flux=flux, hdr_v=hdr_v, var=var, wave=wave, flux_file=flux_file, stddev=stddev)

def load_observations():
    """returns a set of ready-to-use (hdr, flux, var, wave) observation data"""
    o_dir = base_path+"/observations"

    flux_file_ends_with = "_icubes_corrected_flux.fits"
    var_file_ends_with = "_icubes_corrected_var.fits"
    files = find_files_ending_with(o_dir, flux_file_ends_with)
    # make sure that the flux and var files are in the same order
    # get the observation name from the flux files
    obs_names = []
    for i in range(len(files)):
        obs_names.append(files[i].split("/")[-1].split(flux_file_ends_with)[0])

    obs_names.sort()
    o_cnt = len(obs_names)

    ybot = 15.5
    ytop = 80.5

    observations = []
    for i in range(o_cnt):
        ffile = o_dir+"/"+obs_names[i]+flux_file_ends_with
        vfile = o_dir+"/"+obs_names[i]+var_file_ends_with
    
        minwave = 3500.
        maxwave = 5500.

        observations.append(read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop))

    return observations

def load_original_cube(nb_min, nb_max):
    """ 
    The assumption for this function is that we're only using this as a reference for the sightlines.
    The actual extraction is performed on the 7 observational cubes.
    Therefore we only need the astrometry-corrected flux cube.
    """

    hdr, flux = kcwi_io.open_kcwi_cube(reference_flux_filename)
    wcs = WCS(hdr).celestial  # the World Coordinate System (WCS) of the flux cube

    # Removing all the NaNs from the flux cube before making white-light or narrow-band image
    q_nan = np.isnan(flux) 
    flux[q_nan] = 0.0

    wl_image = im.build_whitelight(hdr, flux, minwave=nb_min, maxwave=nb_max)
    return wl_image, wcs
