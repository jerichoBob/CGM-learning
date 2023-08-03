# A "generic" version of abstools that can read in a line list (absorbers at specific redshifts) 
# from rbspecgui and display the stack plots for the absorbers.

from astropy.io import fits
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   
import numpy as np
import csv
import pandas as pd
import bobutils.utils as bu
import sys, os
import argparse

analysisdir = './analysis/j1429'

def read_linelist(specdir):
    linelist_file = specdir + '/SpecPlot_Projects/Identified_LineList.txt'
    reader = csv.reader(open(linelist_file), delimiter=" ")
    redshifts = {}
    next(reader) # skip header
    for row in sorted(reader):
        index, ion, observed, redshift = row
        # print(f"index: {index}, ion: {ion}, observed: {observed}, redshift: {redshift}")
        if redshift not in redshifts: redshifts[redshift] = []
        redshifts[redshift].append({ "ion": ion, "observed": bu.sig_figs(float(observed),8) })
    return redshifts

def read_spec(specdir):
    spec_filename = specdir + '/7_1d_spectra_35.24-7x7-4864-4914.fits'
    a=fits.open(spec_filename)

    flux=a[0].data
    error=a[1].data
    wave=a[2].data
    return wave, flux, error

def get_rest_wavelength(ion, wave_obs, redshift):
    if type(wave_obs) == str:  wave_obs = float(wave_obs)
    if type(redshift) == str:  redshift = float(redshift)

    lambda_rest = wave_obs / (1 + redshift)
    # print(f"lamda_rest: {lambda_rest} ion: {ion} wave_obs: {wave_obs} redshift: {redshift}")
    return lambda_rest

# print("-"*60)
def main():
    for redshift, line_list in redshifts.items():
        lines = []
        for el in line_list:
            ion = el["ion"]
            obs = el["observed"]
            line = get_rest_wavelength(ion, obs, redshift)
            lines.append(bu.sig_figs(float(line),8))

        z = float(redshift)
        lines = sorted(lines)
        last_12_lines = lines[-12:] # have to do this because Metal_Plot.Transitions() can only take up to 12 lines

        # Create an absorber class to feed into the main GUI
        print("."*80)
        print(f"z: {z}  lines: {last_12_lines}")
        print(f"""
    A.Absorber(z={z}, 
            wave={wave}, 
            flux={flux}, 
            error={error}, 
            lines={lines}, 
            window_lim={[-2000,2000]})
    """)
        try:
            absys=A.Absorber(z=z,wave=wave,flux=flux,error=error,lines=last_12_lines, window_lim=[-2000,2000])  
            ions=absys.ions
            M.Transitions(ions)
        except Exception as e:
            print(f"Exception: {e}")
        print(". "*60)

parser = argparse.ArgumentParser(
                    prog='abscli.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This wrapper around abstools takes the output of rbspecgui.py (LineList_Identified.txt) and 
creates a stack plot of the absorbers at specific redshifts.

Example usage:
    python abscli.py -s <specdir> -l - to get a list of the available redshifts
    python abscli.py -s <specdir> -z <redshift> - to create a stack plot of the absorbers at the specified redshift

""", 
                    epilog='----')
parser.add_argument('-s', '--specdir', required=True, help='the directory containing the 1d spectra & LineList_Identified.txt')
parser.add_argument('-l', '--list',     action='store_true', help='list the available redshifts')
parser.add_argument('-z', '--redshift', help='the redshift to analyze')

if __name__ == "__main__":

    args = parser.parse_args()
    if args.specdir is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.specdir:
        specdir = args.specdir
        redshifts = read_linelist(specdir)

    if args.redshift:
        # assumes that the specdir has already been set
        wave, flux, error = read_spec(specdir)
        z = float(args.redshift)
        line_list = redshifts[str(z)]
        lines = []
        for el in line_list:
            ion = el["ion"]
            obs = el["observed"]
            line = get_rest_wavelength(ion, obs, z)
            lines.append(bu.sig_figs(float(line),8))

        z = float(z)
        lines = sorted(lines)
        last_12_lines = lines[-12:]
        print(f"z: {z}  lines: {last_12_lines}")
        try:
            absys=A.Absorber(z=z,wave=wave,flux=flux,error=error,lines=last_12_lines, window_lim=[-2000,2000])  
            ions=absys.ions
            basedir = os.path.join(specdir, 'z_' + str(z))
            print(f"basedir: {basedir}")
            M.Transitions(ions, basedir=basedir)
        except Exception as e:
            print(f"Exception: {e}")

    if args.list:
        # assumes that the specdir has already been set
        print(f"Available redshifts: {list(redshifts.keys())}") 

    # main()