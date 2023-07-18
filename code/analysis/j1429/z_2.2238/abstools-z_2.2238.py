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
# z=2.18025 # Rongmon's value
z=2.18110 # Seems a better fit for Si_III, Si_II_2 & C_IV

# Give approximate absorption line rest frame wavelengths to be analyzed
Si_II_1 = np.array([ 1190.4158, 1193.2897 ]) # ???
Si_III  = np.array([ 1206.500 ])
N_V     = np.array([ 1238.821, 1242.804 ])
Si_II_2 = np.array([ 1526.70698 ])
C_IV    = np.array([ 1548.2049, 1550.77845 ])
Fe_II   = np.array([ 1608.45085 ]) # ??
Al_II   = np.array([ 1670.7886 ]) # ??

lines = np.concatenate((Si_II_1, Si_III, N_V, Si_II_2, C_IV, Fe_II, Al_II)).tolist()

# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines, window_lim=[-2000,2000])   
Abs=absys.ions

#Run the main GUI
M.Transitions(Abs)
