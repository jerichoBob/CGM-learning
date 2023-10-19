import numpy as np
from math import floor, log10
from astropy import units as u
import warnings
from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools import extract_weighted_spectrum as ke

from astropy.wcs import WCS, FITSFixedWarning

# This just draws a single box centered at (x,y) and sz from that center point in the n/e/s/w directions 
def plotbox(plt, x, y, sz, c):
    ax = x - 0.5
    ay = y - 0.5
    plt.plot(
        [ax, ax,    ax-sz, ax-sz, ax],
        [ay, ay-sz, ay-sz, ay,    ay],
        '-', color = c)

def plotcircle(plt, x, y, labels, sz, c):
  for i in range(len(x)):
        colour = 'c' #fix it at cyan for now c[i];
        x_ = x[i]
        y_ = y[i]
        plotbox(plt, x_, y_,  1, colour) # center
        plotbox(plt, x_, y_+1, 1, colour) # north
        plotbox(plt, x_+1, y_, 1, colour) # east
        plotbox(plt, x_, y_-1, 1, colour) # south
        plotbox(plt, x_-1, y_, 1, colour) # west

        plt.text(x[i]-1.5*sz, y[i], labels[i], color=colour)

def read_and_prep_flux_var_data(flux_file, var_file, minwave, maxwave, ybot, ytop):
    """ This method reads the flux and var cubes for a specific observation, and cleans them up before returning """
    hdr, flux = open_kcwi_cube(flux_file)
    _, var = open_kcwi_cube(var_file)
    wave = build_wave(hdr)

    # first do a little data cleanup
    var[np.isnan(flux)]=1.
    flux[np.isnan(flux)] = 0.0000

    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]
    wave = wave[slices]
    flux = flux[slices,int(ybot):int(ytop),:]
    var = var[slices,int(ybot):int(ytop),:]
    print(f"flux.shape={flux.shape}") # (slices, y, x)

    return hdr, flux, var, wave

def corrected_corner_define(xx, yy, flux, var, deltax=5, deltay=5):
    """corrected variation of "corner_define" function found in Extract_1D_Spectra_J2222.ipynb"""
# 
    x_halfbox = (deltax-1)//2
    y_halfbox= (deltay-1)//2

    # extract flux and variance for the region specified
    sub_flux = flux[:, yy-y_halfbox:yy+y_halfbox+1, xx-x_halfbox:xx+x_halfbox+1]
    sub_var  =  var[:, yy-y_halfbox:yy+y_halfbox+1, xx-x_halfbox:xx+x_halfbox+1]

    # - - - - - 
    # NOTE:  to get the right plot alignment when ploting, you need to subtract both x and y by 1.5. WHY?
    # 1. QFitsView starts in the lower-left corner as 1,1, whereas a python/np array is base 0.
    # 2. The 0,0 point is at the center of the pixel, not at the lower left hand corner, 
    coordinate_correction = -1.5

    x1=yy-y_halfbox   + coordinate_correction
    x2=yy+y_halfbox+1 + coordinate_correction
    
    y1=xx-x_halfbox   + coordinate_correction
    y2=xx+x_halfbox+1 + coordinate_correction
    x=[y1,y1,y2,y2,y1]
    y=[x1,x2,x2,x1,x1]
    return x, y, sub_flux, sub_var

def combine_spectra_ivw(specs):
    """Combine the collection of 1D spectra into a single spectrum using inverse variance weighting"""
    """See https://en.wikipedia.org/wiki/Inverse-variance_weighting"""
    print(f"# of spectra={len(specs)}")
    fluxlen = len(specs[0].flux)
    print(f"fluxlen={fluxlen}")
    flux_tot = np.zeros(fluxlen)
    var_tot = np.zeros(fluxlen)
    for lndx in range(fluxlen): # for each lambda, apply inverse variance weighting
        sum_ysigma = 0.0
        sum_1sigma = 0.0
        for sp in specs: # for each observation (that is, each spectrum)...
            sum_ysigma += sp.flux[lndx] / (sp.sig[lndx] ** 2)
            sum_1sigma += 1.0 / (sp.sig[lndx] ** 2)
        flux_tot[lndx] = sum_ysigma / sum_1sigma
        var_tot[lndx] = 1.0 / sum_1sigma

    return flux_tot, var_tot


def extract_spectra(flux_files, var_files, ra, dec, box_size):
    """
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
        hdr, flux, var, wave = read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop)
        wcs_cur = WCS(hdr).celestial
        xc, yc = wcs_cur.world_to_pixel_values(ra, dec)
        xc = np.int32(xc)
        yc = np.int32(yc)
        hb = box_size // 2

        flux_cut = flux[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        var_cut  =  var[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, wave, weights='Data')
        specs.append(sp) # add the spectrum to our list

    
    spec, var = combine_spectra_ivw(specs)            
    return spec, var

def signal_to_noise(wave, flux, co_begin=3500, co_end=5500):
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
def signal_to_noise(wave, flux, var, co_begin=3500, co_end=5500):
    """"""
    pass


def sig_figs(x: float, precision: int):
    """
    Rounds a number to number of significant figures
    Parameters:
    - x - the number to be rounded
    - precision (integer) - the number of significant figures
    Returns:
    - float
    """

    x = float(x)
    precision = int(precision)

    return round(x, -int(floor(log10(abs(x)))) + (precision - 1))
