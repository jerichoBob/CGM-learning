# A "generic" version of abstools that can read in a line list (absorbers at specific redshifts) 
# from rbspecgui and display the stack plots for the absorbers.
# make sure you first establish your own venv like this:
#   python3 -m venv astroresearch
#   source astroresearch/bin/activate

import sys, os

from astropy.io import fits
# PYTHONPATH=/Users/robertseaton/School/github_repos/rbcodes
rbcodes = '/Users/robertseaton/School/github_repos/rbcodes'
rbcodes_ui = '/Users/robertseaton/School/github_repos/rbcodes-ui2-2'
for path in sys.path:
    if rbcodes in path:
        sys.path.remove(rbcodes)

if rbcodes_ui not in sys.path:
    sys.path.append(rbcodes_ui)

# for path in sys.path:
#     print(f"sys.path: {path}")

from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   
import pickle
import csv
import bobutils.utils as bu
import argparse

analysisdir = './analysis/j1429'

def get_ion_redshifts(specdir):
    """ Creates a redshifts dictionary, indexed by a specific redshift, and contains the ions (and their observed wavelength) at that redshift"""
    cwd = os.getcwd()
    # print(f"cwd: {cwd}")
    # print(f"specdir: {specdir}")
    full_specdir = os.path.join(cwd, specdir)

    identified_linelist_file = os.path.join(full_specdir,'Identified_LineList.txt')
    linelist_identified_file = os.path.join(full_specdir,'LineList_Identified.txt')
    if os.path.exists(identified_linelist_file): 
        filepath=identified_linelist_file
    elif os.path.exists(linelist_identified_file):
        filepath=linelist_identified_file
    else:
        print(f"""\nERROR: Can\'t find:
{identified_linelist_file} 
              or 
{linelist_identified_file}!
""")
        sys.exit(-1)
    reader = csv.reader(open(filepath), delimiter=" ")
    redshifts = {}
    next(reader) # skip header
    for row in sorted(reader):
        index, ion, observed, redshift = row
        # print(f"index: {index}, ion: {ion}, observed: {observed}, redshift: {redshift}")
        if redshift not in redshifts: redshifts[redshift] = []
        redshifts[redshift].append({ "ion": ion, "observed": bu.sig_figs(float(observed),8) })
    return redshifts

def read_spec(specdir):
    """grabs the 1d spectrum from `specdir` and returns the wave, flux, and error arrays"""
    # print(f"specdir: {specdir}")
    elements = os.path.split(specdir)
    fits_filename = elements[-1] + '.fits'
    # print(f"last element: {elements[-1]}")

    spec_filename = os.path.join(specdir, fits_filename)
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

def load_transitions_from_pickle(basedir, z):
    pfile = basedir + '/' + "Spectrum_Analysis_z_" + str(z) + ".p"
    # print(f"pfile: {pfile}")
    if not os.path.exists(pfile): 
        print(f"\nERROR: pickle file {pfile} does not exist!\n")
        sys.exit(-1)

    with open(pfile,'rb') as pklfile: absys=pickle.load(pklfile)
    M.Transitions(absys, basedir=basedir)


def load_transitions_from_scratch(specdir, basedir, z):
    wave, flux, error = read_spec(specdir)
    line_list = redshifts[str(z)]

    lines = []
    for el in line_list:
        ion = el["ion"]
        obs = el["observed"]
        line = get_rest_wavelength(ion, obs, z)
        lines.append(bu.sig_figs(float(line),8))

    lines = sorted(lines)
    print(f"z: {z}  lines: {lines}")
    try:
        absys=A.Absorber(z=z,wave=wave,flux=flux,error=error,lines=lines, window_lim=[-2000,2000], order_init=0, mask_init=[-200,200])  
        # absys=A.Absorber(z=z,wave=wave,flux=flux,error=error,lines=lines, window_lim=[-2000,2000],  mask_init=[-200,200])  
        ions=absys.ions

        M.Transitions(ions, basedir=basedir)
        # M.Transitions(ions)
    except Exception as e:
        print(f"Exception: {e}")

def process_transitions(specdir, z, use_pickle = False):

    basedir = os.path.join(specdir, 'z_' + str(z))
    if not os.path.exists(basedir): os.makedirs(basedir)

    if use_pickle:
        load_transitions_from_pickle(basedir, z)
    else:
        load_transitions_from_scratch(specdir, basedir, z)


parser = argparse.ArgumentParser(
                    prog='abscli.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This wrapper around abstools takes the output of rbspecgui.py (LineList_Identified.txt) and 
creates a stack plot of the absorbers at specific redshifts.

Example usage:
    python abscli.py -s <specdir> -l               ::: print a list of the available redshifts
    >>> Available redshifts: ['2.1812', '2.2265', '2.2995', '2.3842', '2.3852']

    python abscli.py -s <specdir> -z <redshift>    ::: create a stack plot of the absorbers at the specified redshift
    >>> z: 2.2265  lines: [1206.5, 1215.6701, 1242.804, 1302.1685, 1304.3702, 1526.7066] & then presents the GUI

    python abscli.py -s <specdir> -z <redshift> -p ::: display a previously saved stack plot of the absorbers at the specified redshift
    >>> 

""", 
                    # epilog='----'
                    )
parser.add_argument('-s', '--specdir', required=True, help='the directory containing the 1d spectra & LineList_Identified.txt')
parser.add_argument('-l', '--list', action='store_true', help='list the available redshifts')
parser.add_argument('-z', '--redshift', help='the redshift to analyze')
parser.add_argument('-a', '--all', action='store_true', help='EXPERIMENTAL: loop over all redshifts and create stack plots for each one')
parser.add_argument('-p', '--pickle', action='store_true', help='used with --specdir and --redshift; load pickle file from redshift directory')
# parser.add_argument('-n', '--next', action='store_true', help='used with --specdir and --redshift; load up the next redshift that has not yet been analyzed (no pickle file yet)')

if __name__ == "__main__":

    args = parser.parse_args()
    # print(f"args.specdir: {args.specdir} args.list: {args.list}  args.pickle: {args.pickle}  args.all: {args.all}")

    if args.specdir is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        specdir = args.specdir
        redshifts = get_ion_redshifts(specdir)

    if args.specdir and args.redshift:
        z = float(args.redshift)
        process_transitions(specdir, z, args.pickle)

    if args.specdir and args.all:
        for redshift in redshifts.keys():
            # print(f"z: {z}")
            # basedir = os.path.join(specdir, 'z_' + str(redshift))
            # print(f"basedir: {basedir}")
            z = float(redshift)
            print(f"process_transitions(specdir={specdir}, z={z}, use_pickle={args.all})")

    if args.specdir and args.list:
        # assumes that the specdir has already been set
        # print(f"Available redshifts: {list(redshifts.keys())}") 
        for z in redshifts.keys():
            print(f"z: {z}  lines: {len(redshifts[z])}")

    if args.redshift == None and args.list == None and args.pickle == None:
        parser.print_help(sys.stderr)
        sys.exit(1)


    # main()