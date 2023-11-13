description = """
Sightline convenience class.
"""
import utils as bu

class Sightline:
    """Class for managing Sightlines"""
    def __init__(self, x=None, y=None, disp_x=None, disp_y=None, radius=None, box_width=None, box_height=None, 
                 color=None, label_alignment=None, snr=None):
        self._x = x
        self._y = y
        self._disp_x = disp_x
        self._disp_y = disp_y
        self._radius = radius
        self._box_width = box_width
        self._box_height = box_height
        self._color = color
        self._label_alignment = label_alignment
        self._snr = snr


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
        
    # --- box_width -----------------------------------------------
    @property
    def box_width(self):
        """Get the 'box_width' property."""
        return self._box_width
    
    @box_width.setter
    def box_width(self, value):
        self._box_width = value
        
    @box_width.deleter
    def box_width(self):
        del self._box_width
        
    # --- box_height -----------------------------------------------
    @property
    def box_height(self):
        """Get the 'box_height' property."""
        return self._box_height

    @box_height.setter
    def box_height(self, value):
        self._box_height = value
        
    @box_height.deleter
    def box_height(self):
        del self._box_height
        
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
        
    # --- methods -----------------------------------------------

def radecs_from_sightline_boxes(wcs_ref, pt_xs, pt_ys, szs):
    """ Converts the sightline box corners to an array of ra,dec tuples"""
    radecs = []
    sightline_count = len(pt_xs)
    for sl_ndx in range(sightline_count): # loop through the sightlines
        print(f"=== sightline   pt_xs: {pt_xs[sl_ndx]} pt_ys: {pt_ys[sl_ndx]} szs: {szs[sl_ndx]}")
        xs, ys = bu.box_corners(pt_xs[sl_ndx], pt_ys[sl_ndx], deltax=szs[sl_ndx], deltay=szs[sl_ndx]) #NOTE: still a question whether we should use these for extraction, or just drawing the box against the wl image
        print(f"=== sightline box corners   xs: {xs} ys: {ys}")
        radec = wcs_ref.pixel_to_world_values(xs, ys)
        print(f"radec: {radec}")
        radecs.append(radec)
    return radecs
    
def load_sightlines():
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