import os
from math import floor, log10
import numpy as np
import warnings

from astropy import units as u
from astropy.wcs import WCS, FITSFixedWarning
from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools import extract_weighted_spectrum as ke

class Observation:
    """ A container for the info associated with a single observation """
    
    # def __init__(self):
    #     self.hdr_f = None
    #     self.hdr_v = None
    #     self.flux = None
    #     self.var = None
    #     self.wcs_f = None
    #     self.wcs_v = None
    #     self.wave = None
    #     self.wl = None
    #     self.f_file = None
    #     self.v_file = None

    # def __init__(self, hdr_f, flux, hdr_v, var, wave, flux_file):
    #     self.hdr_f = hdr_f
    #     self.flux = flux
    #     self.hdr_v = hdr_v
    #     self.var = var
    #     self.wave = wave
    #     self.flux_file = flux_file

    def __init__(self, hdr_f=None, hdr_v=None, flux=None, var=None, wcs_flux=None, wcs_var=None, wave=None, wl=None, flux_file=None, var_file=None):
        self._hdr_f = hdr_f
        self._hdr_v = hdr_v
        self._flux = flux
        self._var = var
        self._wcs_f = wcs_flux
        self._wcs_v = wcs_var
        self._wave = wave
        self._wl = wl          # the whitelight image\
        self._f_file = flux_file
        self._v_file = var_file

    
    # --- hdr_f -----------------------------------------------
    @property
    def hdr_f(self):
        """Get the 'hdr_f' property."""
        return self._hdr_f

    @hdr_f.setter
    def hdr_f(self, value):
        self._hdr_f = value

    @hdr_f.deleter
    def hdr_f(self):
        del self._hdr_f

    # --- hdr_v -----------------------------------------------
    @property
    def hdr_v(self):
        """Get the 'hdr_v' property."""
        return self._hdr_v

    @hdr_v.setter
    def hdr_v(self, value):
        self._hdr_v = value

    @hdr_v.deleter
    def hdr_v(self):
        del self._hdr_v

    # --- flux -----------------------------------------------
    @property
    def flux(self):
        """Get the 'flux' property."""
        return self._flux

    @flux.setter
    def flux(self, value):
        self._flux = value

    @flux.deleter
    def flux(self):
        del self._flux

    # --- var -----------------------------------------------
    @property
    def var(self):
        """Get the 'var' property."""
        return self._var

    @var.setter
    def var(self, value):
        self._var = value

    @var.deleter
    def var(self):
        del self._var
        
    # --- wcs_f -----------------------------------------------
    @property
    def wcs_f(self):
        """Get the 'wcs_f' property."""
        return self._wcs_f

    @wcs_f.setter
    def wcs_f(self, value):
        self._wcs_f = value

    @wcs_f.deleter
    def wcs_f(self):
        del self._wcs_f
        
    # --- wcs_v -----------------------------------------------
    @property
    def wcs_v(self):
        """Get the 'wcs_v' property."""
        return self._wcs_v

    @wcs_v.setter
    def wcs_v(self, value):
        self._wcs_v = value

    @wcs_v.deleter
    def wcs_v(self):
        del self._wcs_v
        
    # --- wave -----------------------------------------------
    @property
    def wave(self):
        """Get the 'wave' property."""
        return self._wave

    @wave.setter
    def wave(self, value):
        self._wave = value

    @wave.deleter
    def wave(self):
        del self._wave
        
    # --- wl -----------------------------------------------
    @property
    def wl(self):
        """Get the 'wl' property."""
        return self._wl

    @wl.setter
    def wl(self, value):
        self._wl = value

    @wl.deleter
    def wl(self):
        del self._wl
        
    # --- f_file -----------------------------------------------
    @property
    def f_file(self):
        """Get the 'f_file' property."""
        return self._f_file

    @f_file.setter
    def f_file(self, value):
        self._f_file = value

    @f_file.deleter
    def f_file(self):
        del self._f_file
        
    # --- v_file -----------------------------------------------
    @property
    def v_file(self):
        """Get the 'v_file' property."""
        return self._v_file

    @v_file.setter
    def v_file(self, value):
        self._v_file = value

    @v_file.deleter
    def v_file(self):
        del self._v_file
        
    def get_all(self):
        return self._hdr_f, self._hdr_v, self._flux, self._var, self._wcs_f, self._wcs_v, self._wave, self._wl, self._f_file, self._v_file

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

    # flux = flux[slices,int(ybot):int(ytop),:] # NOTE: Cropping is busted.  Need to fix this.
    # var  =  var[slices,int(ybot):int(ytop),:] # NOTE: Cropping is busted.  Need to fix this.
    print(f"cropped flux.shape={flux.shape}") # (slices, y, x)
    
    return Observation(hdr_f=hdr_f, flux=flux, hdr_v=hdr_v, var=var, wave=wave, flux_file=flux_file)

def box_corners(pt_x, pt_y, deltax=5, deltay=5, qfitsview_correction=-1.5):
    """
    Given a center point (pt_x,pt_y) defined in QFitsView and a box of size (deltax,deltay), 
    return x1,x2,y1,y2 of the extraction box
    """
    # assumes deltax & deltay are odd
    x_halfbox = (deltax-1)//2
    y_halfbox= (deltay-1)//2

    # - - - - - 
    # NOTE:  to get the right plot alignment when ploting, you need to subtract both x and y by 1.5. WHY?
    # 1. QFitsView starts in the lower-left corner as 1,1, whereas a python/np array is base 0.
    # 2. The 0,0 point is at the center of the pixel, not at the lower left hand corner, 
    # qfitsview_correction = -1.5

    x1=pt_x-x_halfbox   + qfitsview_correction
    x2=pt_x+x_halfbox+1 + qfitsview_correction
    xs = [x1,x2]

    y1=pt_y-y_halfbox   + qfitsview_correction
    y2=pt_y+y_halfbox+1 + qfitsview_correction
    ys =[y1,y2]

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

def corrected_corner_define(pt_x, pt_y, flux, var, deltax=5, deltay=5):
    """
    convenience function which returns the flux and var cutouts for a box around the specified point
    """

    xs, ys = box_corners(pt_x, pt_y, deltax=deltax, deltay=deltay)
    return flux_var_cutout(flux, var, xs, ys)

def combine_spectra_ivw(specs):
    """Combine the collection of 1D spectra into a single spectrum using inverse variance weighting"""
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
        o = read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop)
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
            warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, o.wave, weights='Data')
        specs.append(sp) # add the spectrum to our list

    
    spec, var, wave = combine_spectra_ivw(specs)            
    return spec, var, wave

def extract_spectra_from_observations(observations, sl_radec):
    """
    Extracts a combined spectra from the collection of observations, contained within the sl_radec box.   
    ie. observations = [class Observation, class Observation, ...]
    params:
        * observations: an list of Observation objects
        * sl_radec[0]: the ras of the extraction box
        * sl_radec[1]: the decs of the extraction box
    
    returns: 
        * XSpectrum1D flux and variance objects
    """
    specs = []
    ras = sl_radec[0]
    decs = sl_radec[1]
    print("=-"*40)
    print("=============== BEGINNING EXTRACTION ==================")
    for ob in observations:
        # hdr_f, flux, hdr_v, var, wave, flux_file = o
        wcs_cur = WCS(ob.hdr_f).celestial
        xs, ys = wcs_cur.world_to_pixel_values(ras, decs)
        print(f"before ras={ras}, decs={decs}")
        print(f"before xs={xs}, ys={ys}")
        xs = np.round(xs)
        ys = np.round(ys)
        print(f"after xs={xs}, ys={ys}")        
        # _, _, flux_cut, var_cut = corrected_corner_define(xc, yc, flux, var, deltax=box_size, deltay=box_size)
        flux_cut, var_cut = flux_var_cutout(ob.flux, ob.var, xs, ys)
        print(f"flux_cut.shape={flux_cut.shape}")
        print(f"var_cut.shape={var_cut.shape}")

        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, ob.wave, weights='Data')
        specs.append(sp)
    
    spec, var, wave = combine_spectra_ivw(specs)            
    return spec, var, wave

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

base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
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

def get_corrected_kcwi_data():
    """ returns an array of Observation objects"""
    ob_array = load_observations()
    """ now chocked full with more goodness... """
    for ob in ob_array:
        ob.wcs_f = WCS(ob.hdr_f).celestial

        # Removing all the NaNs from the flux cube before making white-light or narrow-band image
        q_nan = np.isnan(ob.flux) 
        ob.flux[q_nan] = 0.0

        # making a white-light image from the flux data-cube by summing 
        # the flux along the wavelength axis
        ob.wl_k = np.sum(ob.flux, axis=0)
        # ob.wl_k = im.build_whitelight(hdr_k, flux, minwave=global_nb_min, maxwave=global_nb_max)

        # Loading the variance cube from index 2 of the KCWI fits file
        ob.wcs_v = WCS(ob.hdr_v).celestial # The WCS of the variance cube
        # ob = Observation(hdr_f, hdr_v, flux, var, wcs_flux, wcs_var, wave, wl_k)
        # results.append(ob)

    return ob_array
