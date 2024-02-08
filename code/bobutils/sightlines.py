description = """
Sightline convenience class.
"""
import utils as bu
import numpy as np
from pydantic import BaseModel, conlist
from typing import List, Union, Callable, Any, Optional

class RaDec(BaseModel):
    """ a simple RA/DEC class"""
    ra: float = None
    dec: float = None    
    def __str__(self):
        return f"RA: {self.ra}, DEC: {self.dec}"
    def __repr__(self):
        return f"RA: {self.ra}, DEC: {self.dec}"    
class Sightline(BaseModel):
    """Class for managing Sightlines"""
    x: Any = None
    y: Any = None
    w: Any = None
    h: Any = None
    color: Any = None
    radecs: List[RaDec] = []
    label: Any = None    
    label_alignment: Any = None    
    srn: Any = None
            
    def __init__(self, **data):
        """ initialize the class based on input values (things are implicitly set via the super class __init__)"""
        super().__init__(**data)  # Call the super class __init__
                
    def __repr__(self):
        return f"""Sightline {self.label}: x={self.x}, y={self.y}, w={self.w}, h={self.h}, radecs={self.radecs} """


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

def radecs_from_sightline_array(wcs_ref, sightlines: List[Sightline]):
    """ Converts the sightline box corners to an array of ra,dec tuples"""
    radecs = []
    sightline_count = len(sightlines)
    for sl_ndx in range(sightline_count): # loop through the sightlines
        # print(f"=== sightline   pt_xs: {sightlines[sl_ndx].x} pt_ys: {sightlines[sl_ndx].y} szs: {sightlines[sl_ndx].w}")
        xs, ys = bu.box_corners(pt_x=sightlines[sl_ndx].x, 
                                pt_y=sightlines[sl_ndx].y, 
                                deltax=sightlines[sl_ndx].w, 
                                deltay=sightlines[sl_ndx].h) 
        # print(f"=== sightline box corners   xs: {xs} ys: {ys}")
        radec = wcs_ref.pixel_to_world_values(xs, ys)
        # print(f"radec: {radec}")
        radecs.append(radec)
    return radecs

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

def sightline_cuts(sl_radec, wcs):
    """ for the sightline passed in, calculate the x and y cuts to be used when cutting the flux and var cubes"""
    ras = sl_radec[0]
    decs = sl_radec[1]
    
    x_min, y_min = world_to_pixel_int(wcs, ras[0], decs[0])
    x_max, y_max = world_to_pixel_int(wcs, ras[1], decs[1])
    
    return x_min, x_max, y_min, y_max
# radec: (array([217.4790609 , 217.47881255]), array([12.04378713, 12.04403001]))