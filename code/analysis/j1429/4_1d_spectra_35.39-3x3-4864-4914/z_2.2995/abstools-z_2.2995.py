#Read in packages 
from astropy.io import fits
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   
import numpy as np
import sys, getopt

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

def load
def main(argv):
    z = 0
    try:
        opts, args = getopt.getopt(argv,"hp:",["p="])
    except getopt.GetoptError:
        print('abstools.py -p <pickle_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('abstools.py -p <pickle_file>')
            sys.exit()
        elif opt in ("-p", "--p"):
            inputfile = str(arg)
    
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
    # z = 2.29849 # Rongmon's value
    z = 2.2995 # Seems a better fit for Si_III, Si_II_2 & C_IV

    # Give approximate absorption line rest frame wavelengths to be analyzed

    # Si_II_1 = np.array([ 1193.2897 ]) # Si_II 1190 wasn't good
    Si_III  = np.array([ 1206.500 ])

    D_I   = np.array([ 1215.3394 ]) # OK
    H_I   = np.array([ 1215.6701 ]) # OK

    N_V     = np.array([ 1238.821, 1242.804 ]) # PERFECT

    Si_II_2 = np.array([ 1260.4221 ]) # ok @ z=2.2995, but missing doublet  1264

    C_IV    = np.array([ 1548.2049, 1550.77845 ]) # PERFECT, but 1548 shares the bed with another absorption line


    lines = np.concatenate(( Si_III, D_I, H_I, N_V, Si_II_2, C_IV)).tolist()

    # Create an absorber class to feed into the main GUI
    absys=A.Absorber(z,wave,flux,error,lines=lines, window_lim=[-2000,2000])   
    Abs=absys.ions

    #Run the main GUI
    M.Transitions(Abs)

if __name__ == "__main__":
   main(sys.argv[1:])