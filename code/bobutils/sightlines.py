description = """
Sightline convenience class.
"""

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
    