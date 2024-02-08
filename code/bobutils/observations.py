description = """
Observation utilities - especially when dealing with FITS files.
"""
import os, sys
from pydantic import BaseModel, conlist
from typing import List, Union, Callable, Any, Optional
import numpy as np

# modify the python path to include the current directory
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from bobutils import sightlines as bus

from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning
import warnings

from kcwitools import extract_weighted_spectrum as ke

warnings.simplefilter('ignore')
warnings.filterwarnings('ignore', category=FITSFixedWarning)

class Observation(BaseModel):
    """ A container for the data associated with an observation """
    id: Any = None          # the observation id
    hdr_f: Any = None       # FITS header for the flux image
    hdr_v: Any =None        # FITS header for the variance image
    flux: Any = None        # the flux image
    var: Any = None         # the variance image
    wcs_flux: Any = None    # the WCS for the flux image
    wcs_var: Any = None     # the WCS for the variance image
    wave: Any = None        # the wavelength array
    wl: Any = None          # the whitelight image
    flux_file: Any = None   # the full pathname of the flux file
    var_file: Any = None    # the full pathname of the variance file
    stddev: Any = None      # the calculated stddev


    def __init__(self, **data):
        """ initialize the class based on input values (things are implicitly set via the super class __init__)"""
        super().__init__(**data)  # Call the super class __init__

    def __repr__(self):
        return f"""Observation hdr_f: {self.hdr_f}: """

    # --- get_all -----------------------------------------------
    def get_all(self):
        return self._hdr_f, self._hdr_v, self._flux, self._var, self.wcs_flux, self.wcs_var, self._wave, self._wl, self.flux_file, self.var_file

def extract_spectra_from_obs(sl_radec, observations):
    """
    Extracts (but doesn't combine) spectra from the collection of observations, contained within the sl_radec box.   
    ie. observations = [class Observation, class Observation, ...]
    params:
        * observations: an list of Observation objects
        * sl_radec[0]: the ras of the extraction box
        * sl_radec[1]: the decs of the extraction box
    
    returns: 
        * an array of XSpectrum1D flux and variance objects, one for each observation
    """
    specs = []
    ras = sl_radec[0]
    decs = sl_radec[1]
    # print("=-"*40)
    print(f"=============== BEGINNING EXTRACTION FOR RA:{ras[0]} DEC:{decs[0]} ==================")
    for ondx,ob in enumerate(observations):    
        x_min, y_min = bus.world_to_pixel_int(ob.wcs_flux, ras[0], decs[0])
        x_max, y_max = bus.world_to_pixel_int(ob.wcs_flux, ras[1], decs[1])
        # print(f"ras={ras}, decs={decs}  -->  x_min={x_min}, x_max={x_max}   y_min={y_min}, y_max={y_max}")
        flux_cut = ob.flux[:, y_min:y_max, x_min:x_max]
        var_cut  =  ob.var[:, y_min:y_max, x_min:x_max]
        var_cut[var_cut < 0] = 0.0
        var_cut[np.isnan(var_cut)] = 0.0
        # if we are in observation #6, sqrt the variance
        if ondx == 6:
            print(f"SQRT the variance for observation {ondx}")
            var_cut = np.sqrt(var_cut)
        var_cut
        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore', AstropyWarning)
            # warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, ob.wave, weights='Data')
            
        specs.append(sp) # add the spectrum to our list
        
    return specs