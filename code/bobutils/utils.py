# print("<hi>")
import numpy as np
from math import floor, log10
from astropy import units as u

# This just draws a single box centered at (x,y) and sz from that center point in the n/e/s/w directions 
def plotbox(plt, x, y, sz, c):
    ax = x - 0.5;  # I don't know why this is yet.... :()
    ay = y - 0.5;
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


def corrected_corner_define(xx, yy, flux, var, deltax=5, deltay=5):
# corrected variation of Rongmon's "corner_define" function in Extract_1D_Spectra_J2222.ipynb
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
