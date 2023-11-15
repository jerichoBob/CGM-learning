from math import floor, log10
import numpy as np
import warnings
from math import log10, floor, isnan

from astropy import units as u
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyWarning

from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools.spec import extract_square, extract_rectangle
from kcwitools import extract_weighted_spectrum as ke
from kcwitools import image as im

from linetools.spectra.xspectrum1d import XSpectrum1D

from shapely.geometry import Polygon

import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import bobutils.fileio as bio


def box_corners(pt_x, pt_y, deltax=5, deltay=5, coord_corr=None):
    """
    Given a center point (pt_x,pt_y) (possibly defined in QFitsView) and a box of size (deltax,deltay), 
    return x1,x2,y1,y2 of the extraction box
    - if defined in QFitsView, coord_corr should be -1.5
    - if not, coord_corr should be -0.5
    # - - - - - 
    # NOTE:  to get the right plot alignment when ploting, you need to subtract both x and y by 1.5. WHY?
    # 1. QFitsView starts in the lower-left corner as 1,1, whereas a python/np array is base 0.
    # 2. The 0,0 point is at the center of the pixel, not at the lower left hand corner, which adds 0.5 to the coordinates.
    """
    # assumes deltax & deltay are odd
    x_halfbox = (deltax-1)//2
    y_halfbox= (deltay-1)//2

    x1=pt_x-x_halfbox  
    x2=pt_x+x_halfbox+1
    y1=pt_y-y_halfbox  
    y2=pt_y+y_halfbox+1
    
    if coord_corr is not None:
        x1+=coord_corr
        x2+=coord_corr
        y1+=coord_corr
        y2+=coord_corr
        
    xs = [x1,x2]
    ys = [y1,y2]

    return xs, ys

def flux_var_cutout(flux, var, xs, ys):
    x0 = int(xs[0])
    x1 = int(xs[1])
    # if x1-x0==1: x1+=1
    y0 = int(ys[0])
    y1 = int(ys[1])
    # if y1-y0==1: y1+=1
    print(f"flux_var_cutout x0:{x0}, x1:{x1}, y0:{y0}, y1:{y1}")
    return flux[:, y0:y1, x0:x1], var[:, y0:y1, x0:x1] # axis0 = lambda, axis1 = y, axis2 = x



def calc_overlap_area(box_corners, pixel_corners):
    """
    Calculate the area of overlap between the extraction box and a pixel.
    
    :param box_corners: Corners of the extraction box
    :param pixel_corners: Corners of the pixel
    :return: The overlap area
    """
    # Create shapely polygons for the box and the pixel
    box_poly = Polygon(box_corners)
    pixel_poly = Polygon(pixel_corners)
    
    # Calculate the area of overlap
    overlap_area = box_poly.intersection(pixel_poly).area
    return overlap_area


# Let's redefine the create_fractional_mask_rotated function to correctly calculate the weights
def create_fractional_mask(box_corners, shape):
    """
    Create a fractional mask for a given bounding box with subpixel accuracy, handling rotation.
    
    :param box_corners: The corners of the bounding box (possibly rotated)
    :param shape: The shape of the 2D array to apply the mask to 
       - assuming this is the shape of the whitelight image (but could be subset of that)
    :return: A 2D numpy array representing the mask
    """
    mask = np.zeros(shape)
    
    # Iterate over the coordinates of the pixels
    for y in range(shape[0]):
        for x in range(shape[1]):
            # Define the corners of the current pixel
            pixel_corners = np.array([
                [x, y],
                [x+1, y],
                [x+1, y+1],
                [x, y+1]
            ])
            # Calculate the area of overlap and set it as the weight for the pixel
            mask[y, x] = calc_overlap_area(box_corners, pixel_corners)
    
    return mask

def fractional_flux_var_spectra(wave, flux, var, xs, ys):
    """
    Extracts the flux and variance data from the observation cubes using a fractional mask.
    """
    # Calculate the 2D mask
    shape = flux.shape[1:]  # Assuming flux shape is (wavelength, y, x)
    print(f"flux shape={shape}")
    box_corners = np.array([
        [xs[0], ys[0]],
        [xs[1], ys[0]],
        [xs[1], ys[1]],
        [xs[0], ys[1]]
    ])
    mask = create_fractional_mask(box_corners, shape)
    print(f"mask shape={mask.shape}")
    experiment = False
    if experiment:
        # get the s and y bounds of the non-zero part of the mask
        ys = np.nonzero(mask)[0]
        xs = np.nonzero(mask)[1]
        #NOTE: shrink the mask down to it's non-zero values
        mask_cut = mask[np.nonzero(mask)]
        print(f"mask shape={mask.shape}")
        # get the x and y bounds of the mask
        flux_cut = flux[:, ys[0]:ys[1], xs[0]:xs[1]]
        var_cut = var[:, ys[0]:ys[1], xs[0]:xs[1]]
        weighted_flux = np.nansum(flux_cut * mask_cut, axis=(1, 2))
        weighted_var = np.nansum(var_cut * mask_cut, axis=(1, 2))
    else:
        # Apply the mask to the flux and variance
        weighted_flux = np.nansum(flux * mask, axis=(1, 2))
        # weighted_var = np.nansum(var * mask**2, axis=(1, 2))  # Variance scales with the square of the mask
        weighted_var = np.nansum(var * mask, axis=(1, 2)) 


    # Replace NAN - no warnings
    bad = np.isnan(weighted_var)
    weighted_var[bad] = 100. # make it obvious that we have a problem

    print(f"weighted_flux.shape={weighted_flux.shape}")
    print(f"weighted_var.shape={weighted_var.shape}")
    print("="*40)

    # Create and return XSpectrum1D
    return XSpectrum1D.from_tuple((wave, weighted_flux, weighted_var))
    

def corrected_corner_define(pt_x, pt_y, flux, var, deltax=5, deltay=5):
    """
    convenience function which returns the flux and var cutouts for a box around the specified point
    """

    xs, ys = box_corners(pt_x, pt_y, deltax=deltax, deltay=deltay)
    return flux_var_cutout(flux, var, xs, ys)

def combine_spectra_ivw(specs):
    """Combine a collection of 1D spectra into a single spectrum using inverse variance weighting"""
    """See https://en.wikipedia.org/wiki/Inverse-variance_weighting"""
    print(f"# of spectra={len(specs)}")
    fluxlen = len(specs[0].flux)
    print(f"fluxlen={fluxlen}")
    flux_tot = np.zeros(fluxlen)
    var_tot = np.zeros(fluxlen)
    wave_tot = specs[0].wavelength
    for lndx in range(fluxlen): # for each lambda, apply inverse variance weighting
        sum_ysigma = 0.0
        sum_1sigma = 0.0
        for sp in specs: # for each observation (that is, each spectrum)...
            sum_ysigma += sp.flux[lndx] / (sp.sig[lndx] ** 2)
            sum_1sigma += 1.0 / (sp.sig[lndx] ** 2)
        flux_tot[lndx] = sum_ysigma / sum_1sigma
        var_tot[lndx] = 1.0 / sum_1sigma

    return flux_tot, var_tot, wave_tot


def combine_spectra_ivw2(specs):
    """ refined implementation of above -- UNTESTED"""
    """Combine a collection of 1D spectra into a single spectrum using inverse variance weighting"""
    """assumes all spectra have the same dimensionality, and have .wavelength, .flux and .sig attributes"""
    n_spectra = len(specs)
    print(f"# of spectra={n_spectra}")
    
    # Assuming all spectra have the same wavelength grid
    fluxlen = len(specs[0].flux)
    print(f"fluxlen={fluxlen}")
    
    # Extract all fluxes and variances into a 2D array (spectra x wavelength)
    all_fluxes = np.array([sp.flux for sp in specs])
    all_variances = np.array([sp.sig ** 2 for sp in specs])  # Squaring the sigmas to get variances
    
    # Calculate weights for each wavelength across all spectra
    weights = 1 / all_variances
    weighted_fluxes = all_fluxes * weights
    
    # Sum along the spectra axis to get weighted sums
    sum_ysigma = np.sum(weighted_fluxes, axis=0)
    sum_1sigma = np.sum(weights, axis=0)
    
    # Calculate the total flux and variance for each wavelength
    flux_tot = sum_ysigma / sum_1sigma
    var_tot = 1 / sum_1sigma

    wave_tot = specs[0].wavelength
    
    return flux_tot, var_tot, wave_tot


def simple_extract_rectangle(wave, sub_flux, sub_var):
    """
    This is a drastically simplified version of kcwitools.spec.extract_rectangle().
    Assumes that flux and var are already cut (hence sub_flux, sub_var).
    """

    flux_spec = np.nansum(sub_flux, axis=(1,2))
    err_spec = np.sqrt(np.sum(sub_var, axis=(1,2))) # assumes that var is really variance data

    # Replace NAN - no warnings
    bad = np.isnan(err_spec)
    err_spec[bad] = 0.

    # Create and return XSpectrum1D
    return XSpectrum1D.from_tuple((wave, flux_spec, err_spec))
    
def extract_spectra(flux_files, var_files, ra, dec, box_size):
    """
    Extracts a combined spectra from the collection of flux and var files, centered at the specified RA and Dec, with the specified box size.   
    params:
        * flux_files: a list of flux files to extract spectra from
        * var_files: a list of var files to extract spectra from
        * ra: the RA of the sightline
        * dec: the Dec of the sightline
        * box_size: the size of the aperature to extract the spectra from
    
    returns: 
        * XSpectrum1D flux and variance objects
    """
    file_cnt = len(flux_files)
 
    ybot = 15.5
    ytop = 80.5
    specs = []
    for i in range(file_cnt):
        ffile = flux_files[i]
        vfile = var_files[i]

        minwave = 3500.
        maxwave = 5500.
        
        o = bio.read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop)
        wcs_cur = WCS(o.hdr_f).celestial
        o.wcs_f = wcs_cur
        xc, yc = wcs_cur.world_to_pixel_values(ra, dec)
        xc = np.int32(xc)
        yc = np.int32(yc)
        hb = box_size // 2

        flux_cut = o.flux[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        var_cut  =  o.var[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore', AstropyWarning)
            # warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, o.wave, weights='Data')
        specs.append(sp) # add the spectrum to our list

    spec, var, wave = combine_spectra_ivw(specs)            
    return spec, var, wave


def signal_to_noise2(wave, flux, co_begin=3500, co_end=5500):
    '''
    Calculates the signal to noise ratio (SNR) against a continuum region.
    Continuum range is inclusive of the endpoints
    this generally assumes that we are working with a subset or window of the flux produced by 'corrected_corner_define()'
    '''
    
    # Select the data within this range
    wave_min = co_begin * u.AA  # Adjust unit as per your data
    wave_max = co_end * u.AA  # Adjust unit as per your data
    mask = (wave >= wave_min) & (wave <= wave_max)    
    selected_flux = flux[mask]

    # Compute the mean flux and its standard deviation within this range
    mean_flux = np.mean(selected_flux, axis=0)
    stddev_flux = np.std(selected_flux, axis=0)

    # Compute the SNR
    snr = mean_flux / stddev_flux
    #..................................................
    # print(f"mean flux between {wave_min}-{wave_max}: ", mean_flux)
    # print(f"stddev flux between {wave_min}-{wave_max}: ", stddev_flux)
    # print(f"snr between {wave_min}-{wave_max}: ", snr)
    return snr

def signal_to_noise3(wave, flux_spec, var_spec, co_begin=3500, co_end=5500):
    """
    Calculate the signal-to-noise ratio for a given wavelength range in 1D spectra.
    
    :param wave: The 1D numpy array representing the wavelength axis
    :param flux_spec: The 1D numpy array representing the flux spectrum
    :param var_spec: The 1D numpy array representing the variance spectrum
    :param co_begin: The start of the wavelength range
    :param end_wavelco_endength: The end of the wavelength range
    :return: The SNR for the specified wavelength range
    """
    # Find the indices for the specified wavelength range
    wave_min = co_begin * u.AA  # Adjust unit as per your data
    wave_max = co_end * u.AA  # Adjust unit as per your data
    wavelength_indices = np.where((wave >= wave_min) & (wave <= wave_max))[0]

    # Extract the flux and variance for the wavelength range
    # create a local copy of the flux and variance arrays
    flux_copy = np.copy(flux_spec)
    var_copy = np.copy(var_spec)

    flux_range = flux_copy[wavelength_indices]
    variance_range = var_copy[wavelength_indices]


    # Calculate SNR
    snr = np.sum(flux_range) / np.sqrt(np.sum(variance_range))
    return snr

def sig_figs(x: float, precision: int):
    """
    Rounds a number to number of significant figures
    Parameters:
    - x - the number to be rounded
    - precision (integer) - the number of significant figures
    Returns:
    - float
    """

    if isnan(x) or x == 0:
        # Return a default value or raise an error
        return x  # or use `raise ValueError("Input cannot be NaN or zero")`
    precision = int(precision)

    return round(x, -int(floor(log10(abs(x)))) + (precision - 1))


