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
# Rongmon's value: z = 2.22245 
# z=2.2238
z= 2.2265

# Give approximate absorption line rest frame wavelengths to be analyzed
D_I   = np.array([ 1215.3394 ]) # PERFECT
H_I   = np.array([ 1215.6701 ]) # OK

Si_II_1 = np.array([ 1264.7377 ]) # PERFECT -- actually it's Si_II*

N_V     = np.array([ 1242.804 ]) # OK, but shouldn't this be a doublet?
O_I     = np.array([ 1302.1685 ]) # OK? 
Si_II_2 = np.array([ 1304.3702 ]) # PERFECT 

Si_II_3 = np.array([ 1526.70698 ]) # PERFECT
C_IV    = np.array([ 1548.2049, 1550.77845 ]) # PERFECT, but weak


lines = np.concatenate((D_I, H_I, Si_II_1, N_V, O_I, Si_II_2, Si_II_3, C_IV)).tolist()

# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines, window_lim=[-2000,2000])   
Abs=absys.ions

#Run the main GUI
M.Transitions(Abs)
