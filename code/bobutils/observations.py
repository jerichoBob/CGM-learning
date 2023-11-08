description = """
Observation utilities - especially when dealing with FITS files.
"""

class Observation:
    """ A container for the data associated with an observation """

    def __init__(self, hdr_f=None, hdr_v=None, flux=None, var=None, wcs_flux=None, wcs_var=None, wave=None, wl=None, flux_file=None, var_file=None, stddev=None):
        self._hdr_f = hdr_f
        self._hdr_v = hdr_v
        self._flux = flux
        self._var = var
        self._wcs_f = wcs_flux
        self._wcs_v = wcs_var
        self._wave = wave
        self._wl_k = wl          # the whitelight image\
        self._f_file = flux_file
        self._v_file = var_file
        self._stddev = stddev

    
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
    def wl_k(self):
        """Get the 'wl' property."""
        return self._wl_k

    @wl_k.setter
    def wl_k(self, value):
        self._wl_k = value

    @wl_k.deleter
    def wl_k(self):
        del self._wl_k
        
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

    # --- stddev -----------------------------------------------
    @property
    def stddev(self):
        """Get the 'stddev' property."""
        return self._stddev
    
    @stddev.setter
    def stddev(self, value):
        self._stddev = value

    @stddev.deleter
    def stddev(self):
        del self._stddev

    # --- get_all -----------------------------------------------
    def get_all(self):
        return self._hdr_f, self._hdr_v, self._flux, self._var, self._wcs_f, self._wcs_v, self._wave, self._wl, self._f_file, self._v_file

