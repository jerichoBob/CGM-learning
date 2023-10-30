doc = """
This is an updated extraction tool to handle extraction of sightlines across multiple (7) different observations.
"""
import bobutils.layout_utils as lu
import bobutils.utils as bu

from astropy.wcs import WCS, FITSFixedWarning

from kcwitools import io as kcwi_io
from kcwitools import utils as kcwi_u
from kcwitools import spec as kcwi_s
from kcwitools import image as im
from kcwitools import extract_weighted_spectrum as ke

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse

from copy import deepcopy
import warnings
import os
# sightlines were identified using QFitsView using square 3x3 aperture

# global defines
base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"

def old_stuff():
    # SNR was waaay too low for a 3x3 (mostly less than 1) when treating 3700-5500 as the continuum
    sz = 3  # length of one side of the square box


    # First create a narrow-band whitelight image to plot
    narrowband_center = 4686
    del_lambda  = 10
    whitelight= im.build_whitelight(hdr, flux, minwave=narrowband_center - del_lambda, maxwave=narrowband_center + del_lambda)



    fig_size=8

    fig = plt.figure(figsize=(18, 10))
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams['font.family'] = 'serif'
    # ax = fig.add_subplot(111)
    ax_info, ax_image, ax_spectra = lu.make_image_and_12_plots(plt)

    # the continuum, place in a rlatively stable region of the spectra
    co_begin=4864
    co_end=4914

    info_str = f"""J1429+1202 Analysis Metadata
    Continuum: {co_begin} - {co_end}
    Aperture Size: Position Dependent
    Narrowband Center: {narrowband_center} +/- {del_lambda} Angstroms
    """

    ax_info.text(0.05, 0.9, 
                info_str,
                color='k', 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_info.transAxes)

    ax_image.imshow(whitelight/1e4, interpolation='nearest',cmap=plt.get_cmap('gray'),origin="lower")#,vmin=1., vmax=8.)



    font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 14,
            }
    color_sightline = 'yellow'
    color_error = 'red'
    color_snr = 'k'

    output_path = "./analysis/j1429/"
    for i in range(len(points)): # cycle through each sightline
        x = x_coords[i]
        y = y_coords[i]
        sz = sizes[i]
        # print("center: (", x,",",y,")")

        xx_A, yy_A, ff, vv = bu.corrected_corner_define(x, y, flux, var, deltax=sz, deltay=sz)

        # the bounding box and text label on the image
        ax_image.plot(xx_A, yy_A, color=color_sightline, linewidth=1)
        ax_image.text(x-1, y-1, str(i), color=color_sightline, 
                    fontsize = fig_size * 1.25, 
                    ha='center', va='center')
        filename = f'{i}_1d_spectra_{x}.{y}-{sz}x{sz}-{co_begin}-{co_end}.fits'
        filepath = output_path + filename
        print("file: ",filename)

        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp=ke.extract_weighted_spectrum(ff, vv, wave, weights='Data', verbose=False)
            # sp=kcwi_s.extract_square(x, y, wave, flux, var, squaresize=sz, outfile=None)


            #plot the extracted spectrum before writing it out to file.
        ax_spectra[i].plot(sp.wavelength, sp.flux, '-', color=color_sightline)
        ax_spectra[i].plot(sp.wavelength, sp.sig, '-', color=color_error, linewidth=0.5)
        snr = bu.sig_figs(bu.signal_to_noise(sp.wavelength, sp.flux, co_begin=co_begin, co_end=co_end), 5)
        print("SNR: ",snr)

        ax_spectra[i].text(0.05, 0.9, 
                "Sightline: "+str(i),
                color=color_sightline, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_spectra[i].transAxes)
        ax_spectra[i].text(0.05, 0.8, 
                "SNR: "+str(snr),
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_spectra[i].transAxes)    
        ax_spectra[i].text(0.05, 0.7, 
                f'Coord: ({x},{y})',
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_spectra[i].transAxes)   
        ax_spectra[i].text(0.05, 0.6, 
                f'Aperture: {sz} x {sz}',
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_spectra[i].transAxes)  
        ax_spectra[i].set_facecolor('darkgrey')
        # sp.write_to_fits(filepath)


    plt.subplots_adjust(left=0.038, bottom=0.057, right=0.96, top=0.929, wspace=0.174, hspace=0.041)
    plt.show()

def load_ref_cube():
    """Our reference data cube is where we locate the sightlines and construct reference WCS."""
    flux_filename = base_pat+"/J1429_rb_flux.fits" # NOTE: you want the corrected flux and variance cubes here
    var_filename = base_path+"/J1429_rb_var.fits"

    # Load flux and variance data cubes
    _, var = kcwi_io.open_kcwi_cube(var_filename)
    hdr, flux = kcwi_io.open_kcwi_cube(flux_filename)

    # first do a little data cleanup
    var[np.isnan(flux)]=1.
    flux[np.isnan(flux)]=0.0000
    wave = kcwi_u.build_wave(hdr) 
    wcs_ref = WCS(hdr).celestial

    nb_min = 4676.
    nb_max = 4696.
    wl_image = im.build_whitelight(hdr, flux, minwave=nb_min, maxwave=nb_max)

    return flux, var, hdr, wave, wcs_ref, wl_image
    

def load_sightlines():
    """"""
    points = [
        (29,36, 3), # 
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
    # return points


def extract_sightline_spectra(observations, sl_radec):
    """for a given sightline specified by ra,dec, extract the spectra from each observation and combine them"""
    """observations = [(hdr, flux, var, wave)]"""
    spec, var, wave = bu.extract_spectra_from_observations(observations, sl_radec)
    return spec, var, wave

def label_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def create_mpl_layout(plt, count=8):
    """returns an array of {count} axes for images to be placed within"""
    fig = plt.figure(figsize=(10,7), layout="constrained", dpi=150)
    # fig.suptitle(title)
    axs_spec = []
    rows = 17
    cols = 13
    gs = gridspec.GridSpec(rows, cols, figure=fig)
    ax_image = fig.add_subplot(gs[9:11, :])
    for i in range(count):
        print(f"index:{i}")
        axs_spec.append(fig.add_subplot(gs[0:9, i]))
        label_axes(fig)
    # return fig, ax_spectrum, axs_image
    return fig, ax_image, axs_spec

def show_ref_image(axs_image, wl_image):
    """"""
    pass
def show_sightline_spectra(ax_specs):
    """"""
    pass

parser = argparse.ArgumentParser(
                    prog='j1429_extraction_tool_v4.py',
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description="""
This tool extract the spectra for the sightlines identified in the J1429+1202 system.
For each sightline identified, it extracts the spectra at the same ra,dec for each of the 7 observations, 
and combines the spectra using inverse variance weighting.

""", 
                    # epilog='----'
                    )
# parser.add_argument('-s', '--specdir', required=True, help='the directory containing the 1d spectra & LineList_Identified.txt')
# parser.add_argument('-l', '--list', action='store_true', help='list the available redshifts')
# parser.add_argument('-z', '--redshift', help='the redshift to analyze')
# parser.add_argument('-a', '--all', action='store_true', help='EXPERIMENTAL: loop over all redshifts and create stack plots for each one')
# parser.add_argument('-p', '--pickle', action='store_true', help='used with --specdir and --redshift; load pickle file from redshift directory')
# parser.add_argument('-n', '--next', action='store_true', help='used with --specdir and --redshift; load up the next redshift that has not yet been analyzed (no pickle file yet)')

def handle_commandline_args():
    """ """
    pass
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

def main(): 
    """ """
    flux, var, hdr, wave, wcs_ref, wl_image = load_ref_cube()
    xs, ys, szs = load_sightlines()
    sl_radecs = radecs_from_sightline_boxes(wcs_ref, xs, ys, szs)

    observations = bu.load_observations()
    sightline_count = len(xs)
    for sl_ndx in range(sightline_count): # loop through the sightlines
        spec, var, wave = extract_sightline_spectra(observations, sl_radecs[sl_ndx])

    # definte graphics environment
    _, ax_image, ax_specs = create_mpl_layout(plt, count=sightline_count )

    show_ref_image(ax_image, wl_image)
    show_sightline_spectra(ax_specs)

if __name__ == "__main__":
    main()

