#Read in packages 
from astropy.io import fits
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   
import numpy as np

# Read in the 1D spectrum to be analyzed
# filename='1d_spectra_29.36-3x3-4750-5500.fits'
# filename = 'analysis/j1429/4_1d_spectra_35.39-3x3-4864-4914.fits'
filename = '../4_1d_spectra_35.39-3x3-4864-4914.fits'
a=fits.open(filename)

flux=a[0].data
error=a[1].data
wave=a[2].data

#------------------------------
#Specify redshift at which to perform analysis
# z = 2.34708 (value Rongmon suggested)
z=2.3482  # seems a better fit for my data

# Give approximate absorption line rest frame wavelengths to be analyzed
N_I   = np.array([ 1200.7098 ]) # OK
D_I   = np.array([ 1215.3394 ]) # PERFECT
H_I   = np.array([ 1215.6701 ]) # OK

Si_II_1 = np.array([ 1264.7377 ]) # PERFECT

Si_II_2 = np.array([ 1304.3702 ]) # PERFECT 



lines = np.concatenate((N_I, D_I, H_I, Si_II_1, Si_II_2)).tolist()

# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines, window_lim=[-2000,2000])   
Abs=absys.ions

#Run the main GUI
M.Transitions(Abs)
