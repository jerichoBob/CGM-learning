""" Utility functions for working with Jupyter notebooks. """

import os, sys, re
import pickle
import numpy as np
from typing import Any, Union
from pydantic import BaseModel

base = '/Users/robertseaton/School/github_repos/CGM-learning/code'
if base not in sys.path:
    sys.path.insert(0, os.path.abspath(base))
from bobutils import data_analysis as da

def find_pickle_files(basedir):
    return (os.path.join(root, file)
        for root, dirs, files in os.walk(basedir)
            for file in files if file.lower().endswith('.p'))
    
def get_sightline_from_path(path: str) -> Union[str, None]:
    """Extract the sightline from the path."""
    sightline_match = re.search(r'[\w\d\_\+]+/(\d+)_', path)
    if sightline_match:
        sightline = sightline_match.group(1)
    else:
        sightline = None
    return sightline

def get_redshift_from_path(path: str) -> Union[float, None]:
    """Extract the redshift from the path."""
    redshift_match = re.search(r'z_([\d.]+)', path)
    if redshift_match:
        redshift = float(redshift_match.group(1))
    else:
        redshift = None
    return redshift

def get_ion_index_from_path(path: str) -> Union[int, None]:
    """Extract the ion index from the path."""
    ion_index_match = re.search(r'Ions(\d+)', path)
    if ion_index_match:
        ion_index = int(ion_index_match.group(1))
    else:
        ion_index = None
    return ion_index

def get_tuple_from_path(path: str) -> Union[tuple, None]:
    """Extract the sightline, redshift, and ion index from the path."""
    sightline = get_sightline_from_path(path)
    redshift = get_redshift_from_path(path)
    if sightline and redshift:
        return (sightline, redshift)
    else:
        return None
    
def load_pickled_ion_data() -> dict:
    """Load the pickled ion data from files."""
    basedir = "/Users/robertseaton/School/github_repos/CGM-learning/code/analysis/J1429+1202_240131/combined_spectra"
    files = list(sorted(find_pickle_files(basedir)))
    
    ion_data = {} # ion_data[z][sl] = ions - and ions[]
    
    for pfile in files:
        # if i > 1: break
        (sl, z) = get_tuple_from_path(pfile)
        if z not in ion_data: 
            ion_data[z] = {}
        with open(pfile,'rb') as pickle_file:
            ions = pickle.load(pickle_file)
            pickle_file.close()
            ion_data[z][sl] = ions
    return ion_data

def count_sightlines(z, ion_data):
    """ how many sightlines do we have for a specific redshift """
    return len(ion_data[z])

class fitDataClass(BaseModel):
    xdata: Any = None
    nflux: Any = None
    nerror: Any = None
    gxdata_rest: Any = None
    gnflux: Any = None
    gnerror: Any = None
    ewlims_rest: list = []
    wlims: list = []
    wingspan: float = 0.0
    # vasic fit data
    g: Any = None
    fit_g: Any = None
    ew: float = 0.0
    ew_sig: float = 0.0
    good_fit: bool = False
    # broken out data for downstream processing
    cont: float = 0.0 # continuum
    abs1_amp: float = 0.0 # amplitude of the first gaussian model component
    abs2_amp: float = 0.0 # amplitude of the second gaussian model component
    abs1_mean: float = 0.0 # mean of the first gaussian model component
    abs2_mean: float = 0.0  # mean of the second gaussian model component
    abs1_stdev: float = 0.0  # standard deviation of the first gaussian model component
    abs2_stdev: float = 0.0  # standard deviation of the second gaussian model component
    
    def __init__(self, **data):
        """ initialize the simulation based on input values (things are implicitly set via the super class __init__)"""
        super().__init__(**data)  # Call the super class __init__
    
    
DEBUG=False
def create_metadata_for_gaussian_fitting(ion, use_vel=False):
    """prepare the ion ready for the gaussian fitting."""
    """Prepare the data: Choose how you want to do the analysis, and set the analysis range/limits"""
    vel = ion['vel']
    x_lim = ion['window_lim']
    wave  = ion['wave']
    ew_lims = ion['EWlims']
    flux = ion['flux']
    y_error = ion['error']
    cont = ion['cont']
    name = ion['name']
    zabs = ion['z']    
    
    x1_index = np.where(vel >= ew_lims[0])[0][0] # the x index where the velocity is greater than the lower limit
    x2_index = np.where(vel >= ew_lims[1])[0][0] # the x index where the velocity is greater than the upper limit
    if DEBUG: print(f"x1_index = {x1_index}   x2_index = {x2_index}")
    if use_vel:
        if DEBUG: print('Using velocity')
        adj = 0
        xdata = vel
        ewlims = [ew_lims[0]-adj, ew_lims[1]+adj]
        wlims = [-2000,2000]
    else: 
        if DEBUG: print('Using wavelength')
        xdata = wave
        ewlims = [wave[x1_index-1],wave[x2_index+1]]
        right_bound = x2_index+min(20, len(wave)-x2_index-1)
        if DEBUG: print(f"Right Bound: {right_bound}")
        wlims = [wave[x1_index-20],wave[right_bound]]

    # normalize the flux and error to the continuum
    nflux = flux/cont
    nerror = y_error/cont

    # print(f"ewlims = {ewlims}")
    ewlims_rest = [ewlims[0]/(1. + zabs),ewlims[1]/(1. + zabs)]
    # print(f"ewlims-obs: {ewlims}")
    # print(f"ewlims_rest: {ewlims_rest}")

    delwv       = np.double(ewlims[1])-np.double(ewlims[0])
    delwv_rest  = np.double(ewlims_rest[1])-np.double(ewlims_rest[0])
    # print(f"delwv = {delwv}   delvw_rest = {delwv_rest}")  

    # to limit the influence of points away from the signal of interest, we limit the analysis range
    wingspan = delwv / 2.0

    new_limits = (xdata >= ewlims[0]-wingspan) & (xdata <= ewlims[1]+wingspan)
    gxdata = xdata[new_limits]
    gxdata_rest = gxdata/(1. + zabs)
    # print(f"gxdata-obs: {gxdata}")
    # print(f"gxdata-rest: {gxdata_rest}")

    gnflux = nflux[new_limits]
    gnerror = nerror[new_limits]
    if DEBUG: print(f"limiting fit analysis range to {ewlims[0]-delwv}, {ewlims[1]+delwv}")  

    return xdata, nflux, nerror, gxdata_rest, gnflux, gnerror, ewlims_rest, wlims, wingspan

def fit_and_compute_ew(gauss_metadata):
    """Compute the equivalent width of an ion."""
    """wave is either the rest-frame wavelength or the velocity."""
    """flux is the normalized flux, ewlims_rest is the rest-frame equivalent width limits in Angstroms."""
    g, fit_g, good_fit = da.create_gaussian_fitter(gauss_metadata.gxdata_rest, 
                                                   gauss_metadata.gnflux, 
                                                   gauss_metadata.ewlims_rest)
    if good_fit:
        ew, ew_sig = da.calculate_ew_from_gaussian(g, fit_g, gauss_metadata.ewlims_rest)
        ew = 1000 * ew # convert to mAngstroms
    else:
        ew = 0
        ew_sig = 0
    return g, fit_g, ew, ew_sig, good_fit



DEBUG=False
def fit_ion(ion, use_vel=False):
    """create the necessary metadata, and then fit the ion."""

    vel = ion['vel']
    # x_lim = ion['window_lim']
    wave  = ion['wave']
    ew_lims = ion['EWlims']
    flux = ion['flux']
    y_error = ion['error']
    cont = ion['cont']
    # name = ion['name']
    zabs = ion['z']    

    fitData = fitDataClass()
    
    x1_index = np.where(vel >= ew_lims[0])[0][0] # the x index where the velocity is greater than the lower limit
    x2_index = np.where(vel >= ew_lims[1])[0][0] # the x index where the velocity is greater than the upper limit
    if DEBUG: print(f"x1_index = {x1_index}   x2_index = {x2_index}")
    if use_vel:
        if DEBUG: print('Using velocity')
        adj = 0
        xdata = vel
        ewlims = [ew_lims[0]-adj, ew_lims[1]+adj]
        wlims = [-2000,2000]
    else: 
        if DEBUG: print('Using wavelength')
        xdata = wave
        ewlims = [wave[x1_index-1],wave[x2_index+1]]
        right_bound = x2_index+min(20, len(wave)-x2_index-1)
        if DEBUG: print(f"Right Bound: {right_bound}")
        wlims = [wave[x1_index-20],wave[right_bound]]

    fitData.xdata = xdata
    fitData.wlims = wlims

    # normalize the flux and error to the continuum
    nflux = flux/cont
    nerror = y_error/cont
    
    fitData.nflux = nflux
    fitData.nerror = nerror

    # print(f"ewlims = {ewlims}")
    ewlims_rest = [ewlims[0]/(1. + zabs),ewlims[1]/(1. + zabs)]
    # print(f"ewlims-obs: {ewlims}")
    # print(f"ewlims_rest: {ewlims_rest}")
    fitData.ewlims_rest = ewlims_rest
    
    delwv       = np.double(ewlims[1])-np.double(ewlims[0])
    delwv_rest  = np.double(ewlims_rest[1])-np.double(ewlims_rest[0])
    # print(f"delwv = {delwv}   delvw_rest = {delwv_rest}")  

    # to limit the influence of points away from the signal of interest, we limit the analysis range
    wingspan = delwv / 2.0

    new_limits = (xdata >= ewlims[0]-wingspan) & (xdata <= ewlims[1]+wingspan)
    gxdata = xdata[new_limits]
    gxdata_rest = gxdata/(1. + zabs)
    # print(f"gxdata-obs: {gxdata}")
    # print(f"gxdata-rest: {gxdata_rest}")

    gnflux = nflux[new_limits]
    gnerror = nerror[new_limits]
    if DEBUG: print(f"limiting fit analysis range to {ewlims[0]-delwv}, {ewlims[1]+delwv}")  



    g, fit_g, good_fit = da.create_gaussian_fitter(gxdata_rest, 
                                                   gnflux, 
                                                   ewlims_rest)
    if good_fit:
        ew, ew_sig, amp1, amp2, mean, stdev = da.calculate_ew_from_gaussian(g, fit_g, ewlims_rest)
        ew = 1000 * ew # convert to mAngstroms
    else:
        ew, ew_sig, amp1, amp2, mean, stdev = 0, 0, 0, 0, 0, 0
        
    fitData.g = g
    fitData.fit_g = fit_g
    fitData.good_fit = good_fit
    fitData.ew = ew
    fitData.ew_sig = ew_sig
    fitData.cont
        
    # return xdata, nflux, nerror, gxdata_rest, gnflux, gnerror, ewlims_rest, wlims, wingspan, 
    return fitData