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
z=2.1812 # Seems a better fit for Si_III, Si_II_2 & C_IV

# Give approximate absorption line rest frame wavelengths to be analyzed
Fe_II   = np.array([ 1144.9379 ]) # ok @ z=2.1812
Si_II_1 = np.array([ 1190.4158, 1193.2897 ]) # 1190 ok and 1193 perfect @ z=2.1812
Si_III  = np.array([ 1206.500 ]) # perfect @ z=2.1812
N_V     = np.array([ 1238.821, 1242.804 ]) # 1238 perfect, 1242 ok @ z=2.1812

Si_II_4 = np.array([ 1260.4221 ]) # ok @ z=2.1812

O_I     = np.array([ 1302.1685 ]) # perfect @ z=2.1812
Si_II_3 = np.array([ 1304.3702 ]) # perfect @ z=2.1812
C_II    = np.array([ 1334.5323 ]) # perfect @ z=2.1812
Si_IV   = np.array([ 1393.76018 ]) # ok @ z=2.1812

Si_II_2 = np.array([ 1526.70698 ]) # perfect @ z=2.1812
C_IV    = np.array([ 1548.2049, 1550.77845 ]) # perfect @ z=2.1812
Fe_II   = np.array([ 1608.45085 ]) # perfect @ z=2.1812
Al_II   = np.array([ 1670.7886 ]) # ok @ z=2.1812

# lines = np.concatenate((Fe_II, Si_II_1, Si_III, N_V, Si_II_2, C_IV, Fe_II, Al_II)).tolist()
lines = np.concatenate((Fe_II, Si_II_1, Si_III, N_V,   Si_II_4, O_I, Si_II_3, C_II, Si_IV,    Si_II_2, C_IV, Fe_II, Al_II)).tolist()

# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines, window_lim=[-2000,2000])   
Abs=absys.ions

#Run the main GUI
M.Transitions(Abs)
