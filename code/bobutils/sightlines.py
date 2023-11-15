description = """
Sightline convenience class.
"""
import utils as bu
import numpy as np

class Sightline:
    """Class for managing Sightlines"""
    def __init__(self, x=None, y=None, w=None, h=None, color=None, radecs=[],
                 label_alignment=None, snr=None, label=None):
        self._x = x
        self._y = y
        self._w = w
        self._h = h
        self._color = color
        self._label = label
        self._label_alignment = label_alignment
        self._radecs = radecs

        self._snr = snr
        
    def __repr__(self):
        return f"""Sightline {self.label}: x={self.x}, y={self.y}, w={self.w}, h={self.h}, radecs={self.radecs} """

    # --- x -----------------------------------------------
    @property
    def x(self):
        """Get the 'x' property."""
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @x.deleter
    def x(self):
        del self._x    

    # --- y -----------------------------------------------
    @property
    def y(self):
        """Get the 'y' property."""
        return self._y
    
    @y.setter
    def y(self, value):
        self._y = value
        
    @y.deleter
    def y(self):
        del self._y
        
    # --- disp_x -----------------------------------------------
    @property
    def disp_x(self):
        """Get the 'disp_x' property."""
        return self._disp_x
    
    @disp_x.setter
    def disp_x(self, value):
        self._disp_x = value
        
    @disp_x.deleter
    def disp_x(self):
        del self._disp_x
        
    # --- disp_y -----------------------------------------------
    @property
    def disp_y(self):
        """Get the 'disp_y' property."""
        return self._disp_y
    
    @disp_y.setter
    def disp_y(self, value):
        self._disp_y = value
        
    @disp_y.deleter
    def disp_y(self):
        del self._disp_y
        
    # --- radius -----------------------------------------------
    @property
    def radius(self):
        """Get the 'radius' property."""
        return self._radius
    
    @radius.setter
    def radius(self, value):
        self._radius = value
        
    @radius.deleter
    def radius(self):
        del self._radius
        
    # --- w -----------------------------------------------
    @property
    def w(self):
        """Get the 'w' property."""
        return self._w
    
    @w.setter
    def w(self, value):
        self._w = value
        
    @w.deleter
    def w(self):
        del self._w
        
    # --- h -----------------------------------------------
    @property
    def h(self):
        """Get the 'h' property."""
        return self._h

    @h.setter
    def h(self, value):
        self._h = value
        
    @h.deleter
    def h(self):
        del self._h
        
    # --- color -----------------------------------------------
    @property
    def color(self):
        """Get the 'color' property."""
        return self._color
    
    @color.setter
    def color(self, value):
        self._color = value
        
    @color.deleter
    def color(self):
        del self._color
        
    # --- label_alignment -----------------------------------------------
    @property
    def label_alignment(self):
        """Get the 'label_alignment' property."""
        return self._label_alignment
    
    @label_alignment.setter
    def label_alignment(self, value):
        self._label_alignment = value
        
    @label_alignment.deleter
    def label_alignment(self):
        del self._label_alignment
        
    # --- snr -----------------------------------------------
    @property
    def snr(self):
        """Get the 'snr' property."""
        return self._snr
    
    @snr.setter
    def snr(self, value):
        self._snr = value
        
    @snr.deleter
    def snr(self):
        del self._snr
        
    # --- label -----------------------------------------------
    @property
    def label(self):
        """Get the 'label' property."""
        return self._label
    
    @label.setter
    def label(self, value):
        self._label = value
        
    @label.deleter
    def label(self):
        del self._label
    # --- radecs -----------------------------------------------
    @property
    def radecs(self):
        """Get the 'radecs' property."""
        return self._radecs
    
    @radecs.setter
    def radecs(self, value):
        self._radecs = value
        
    @radecs.deleter
    def radecs(self):
        del self._radecs


def add_radec_to_sightlines(wcs_ref, sightlines):
    """ Converts the sightline box corners to an array of ra,dec tuples"""
    radecs = []
    for sl in sightlines: # loop through the sightlines

        # print(f"=== sightline #{sl_ndx}   x: {x} y: {y} w: {w} h: {h}")
        xs, ys = bu.box_corners(sl.x, sl.y, deltax=sl.w, deltay=sl.h,coord_corr=None)
        # print(f"=== extraction box corners   xs: {xs} ys: {ys}")
        radecs = wcs_ref.pixel_to_world_values(xs, ys)
        # print(f"radec: {radec}")
        sl.radecs = radecs # [[ra_min, ra_max], [dec_min, dec_max]]
    return sightlines
    
def load_OLD_sightlines():
    """ we'll figure out a way to make this more generic, but trying to consolidate the set of sightlines here"""
    points = [
        (29,36, 3), # just doing 1 sightline for now
        # (29,39,3),    
        # (32,37,3),
        # (32,40,3),
        # (35,39,3),
        # (38,39,3), 
        # (50,33,7),
        # (35,24,7),
    ]
    x_coords, y_coords, sizes = zip(*points)
    # convert to lists
    x_coords = list(x_coords)
    y_coords = list(y_coords)
    sizes    = list(sizes)
    return x_coords, y_coords, sizes

def load_231113_sightlines():
    sightlines = []
    
    sightlines.append(Sightline(x=27, y=36, w=3, h=9, label="1"))
    sightlines.append(Sightline(x=30, y=37, w=3, h=11, label="2"))
    sightlines.append(Sightline(x=33, y=38, w=3, h=11, label="3"))
    sightlines.append(Sightline(x=36, y=39, w=3, h=11, label="4"))
    sightlines.append(Sightline(x=40, y=39, w=5, h=7, label="5"))
    
    sightlines.append(Sightline(x=49, y=31, w=11, h=7, label="6"))
    sightlines.append(Sightline(x=34, y=23, w=11, h=7, label="7"))
    
    return sightlines

def get_xywhs(sightlines):
    """ returns the x,y,w,h of the sightlines"""
    points = []
    for sl in sightlines:
        points.append((sl.x, sl.y, sl.w, sl.h))
    xs, ys, ws, hs = zip(*points)
    # convert to lists
    xs = list(xs)
    ys = list(ys)
    ws = list(ws)
    hs = list(hs)
    return xs, ys, ws, hs

def world_to_pixel_int(wcs, ra, dec):
    x, y = wcs.world_to_pixel_values(ra, dec)  
    x=int(abs(np.rint(x)))
    y=int(abs(np.rint(y)))
    return x, y

def sightline_cuts(sl, wcs):
    """for the sightline passed in, calculate the x and y cuts to be used when cutting the flux and var cubes"""
    ras = sl.radecs[0]
    decs = sl.radecs[1]
    
    x_min, y_min = world_to_pixel_int(wcs, ras[0], decs[0])
    x_max, y_max = world_to_pixel_int(wcs, ras[1], decs[1])
    
    return x_min, x_max, y_min, y_max