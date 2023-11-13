doc = """
This is an updated-updated extraction tool to handle extraction of spectra 
using flux/variance cubes from 7 different exposures (observations). 
The extraction boxes are no longer square (independent xw and yh sizes).
This extraction also takes into account that when the extraction boxes are 
WCS-converted to RA/DEC and then to xs,ys, we MUST CONFIRM that:
 - working left to right, confirm that if a box overlaps a pixel, 
   it is included in the sum for that sightline, 
   EXCEPT when it is already in another sightline 
   (probably a good idea to colorize the pixels with alpha=0.5 to visually confirm)
"""
import bobutils.layout_utils as lu
import bobutils.utils as bu
import bobutils.fileio as bio

from kcwitools import image as im
from kcwitools import extract_weighted_spectrum as ke

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse

import warnings
# sightlines were identified using QFitsView using square 3x3 aperture

# global defines
base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
global_nb_min = 4676. 
global_nb_max = 4696. 
global_cmap = 'gnuplot'
global_lw = 0.5

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
    ax_info, ax_image, ax_spectra = lu.make_image_and_8_plots(plt)

    ax_image.imshow(whitelight/1e4, interpolation='nearest',cmap=plt.get_cmap('gray'),origin="lower")#,vmin=1., vmax=8.)

    font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 14,
            }
    color_sightline = 'blue'
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
  
        filename = f'{i}_1d_spectra_{x}.{y}-{sz}x{sz}-{co_begin}-{co_end}.fits'
        filepath = output_path + filename
        print("file: ",filename)

        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp=ke.extract_weighted_spectrum(ff, vv, wave, weights='Data', verbose=False)
            # sp=kcwi_s.extract_square(x, y, wave, flux, var, squaresize=sz, outfile=None)


            #plot the extracted spectrum before writing it out to file.
        # sp.write_to_fits(filepath)


def load_sightlines():
    """"""
    points = [
        (29,36, 3), # 
        (29,39,3),    
        (32,37,3),
        (32,40,3),
        (35,39,3),
        (38,39,3), 
        (50,33,7),
        (35,24,7),
    ]
    x_coords, y_coords, sizes = zip(*points)
    # convert to lists
    x_coords = list(x_coords)
    y_coords = list(y_coords)
    sizes    = list(sizes)
    return x_coords, y_coords, sizes
    # return points


# def extract_sightline_spectra(observations, sl_radec):
#     """----------------------- UNUSED -----------------------"""
#     """for a given sightline specified by ra,dec, extract the spectra from each observation and combine them"""
#     """observations = [(hdr, flux, var, wave)]"""
#     spec, var, wave = bu.extract_spectra_from_observations(observations, sl_radec)
#     return spec, var, wave

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

def show_ref_image(ax, image, title="White Light", xh_lim=None, yh_lim=None):
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=0)
    ax.set_title(title)
    if xh_lim is not None: ax.set_xlim(xh_lim)
    if yh_lim is not None: ax.set_ylim(yh_lim)

    ax.grid()

def show_sightline_spectra(ax_specs, specs, color_sightline='blue', color_error='red', color_snr='k', co_begin=4864, co_end=4914, xs=[], ys=[], szs=[]):
    """ show all of the sightline spectra on their own axes """

    for i, sp in enumerate(specs):
        spec, var, wave = sp
        ax_specs[i].plot(wave, spec, '-', color=color_sightline, linewidth=global_lw)
        twin = ax_specs[i].twinx()
        twin.plot(wave, var, '-', color=color_error, linewidth=global_lw)
        # snr = bu.sig_figs(bu.signal_to_noise2(wave, spec, co_begin=co_begin, co_end=co_end), 5)
        snr = bu.sig_figs(bu.signal_to_noise3(wave, spec, var, co_begin=co_begin, co_end=co_end), 5)
        print(f"#{i}  SNR: {snr}")

        ax_specs[i].text(0.05, 0.9, 
                "Sightline: "+str(i),
                color=color_sightline, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_specs[i].transAxes)
        ax_specs[i].text(0.05, 0.8, 
                "SNR: "+str(snr),
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_specs[i].transAxes)    
        ax_specs[i].text(0.05, 0.7, 
                f'Coord: ({xs[i]},{ys[i]})',
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_specs[i].transAxes)   
        ax_specs[i].text(0.05, 0.6, 
                f'Aperture: {szs[i]} x {szs[i]}',
                color=color_snr, 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_specs[i].transAxes)  
        ax_specs[i].set_facecolor('darkgrey')

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

def show_info(ax_info, nb_min, nb_max):
    """ """
    # the continuum, place in a rlatively stable region of the spectra
    co_begin=4864
    co_end=4914

    info_str = f"""J1429+1202 Analysis Metadata
    Continuum: {co_begin} - {co_end}
    Aperture Size: Position Dependent
    Narrowband Min: {nb_min} to {nb_max} Angstroms
    """

    ax_info.text(0.05, 0.9, 
                info_str,
                color='k', 
                fontsize = 12, 
                ha='left', va='top',
                transform=ax_info.transAxes)

def main(): 
    """ Do all the things """
    print("""
    
    ================= Loading Sightlines =================
    """)
    xs, ys, szs = load_sightlines()
    cnt = len(xs)
    print("""
    
    ================= Loading Reference Image =================
    """)    
    wl_image, wcs_ref = bio.load_original_cube(global_nb_min, global_nb_max)
    print("""
    
    ================= Computing Extraction Box RA/DECs =================
    """)
    sl_radecs = radecs_from_sightline_boxes(wcs_ref, xs, ys, szs)

    fig = plt.figure(figsize=(18, 10))
    print("""
    
    ================= Creating Matplotlib Layout =================
    """)
    ax_info, ax_image, ax_spectra = lu.make_image_and_8_plots(plt, wcs=wcs_ref)
    fig.suptitle(f"J1429+1202 Extracted Spectra for {cnt} sightlines", fontsize=16)

    show_info(ax_info, global_nb_min, global_nb_max)

    print("""
    
    ================= Loading Observations =================
    """)

    observations = bio.load_observations()
    specs = []
    print("""
    
    ================= Spectra Extraction for Each Sightline =================
    """)    
    for sl_ndx in range(cnt):  # loop through the sightlines
        sl_radec = sl_radecs[sl_ndx] # let's plot for the first sightline only
        # box_ras, box_decs = sl_radec

        # spec, var, wave = bu.extract_spectra_from_observations(sl_radec)
        spec = bu.extract_spectra_from_observations(sl_radec, observations)
        specs.append(spec)

    show_ref_image(ax_image, wl_image)
    bu.plot_sightlines_wcs(ax_image, wcs_ref, sl_radecs, color='w-', lw=0.5)

    show_sightline_spectra(ax_spectra, specs, 
                           color_sightline='blue', 
                           color_error='red', 
                           color_snr='k', 
                           co_begin=4864, co_end=4914, 
                           xs=xs, ys=ys, szs=szs)

    # definte graphics environment
    mng = plt.get_current_fig_manager()
    if hasattr(mng, 'window'):
        mng.window.state('zoomed') #works fine on Windows!    
    plt.subplots_adjust(left=0.038, bottom=0.057, right=0.917, top=0.929, wspace=0.165, hspace=0.244)
    plt.show()


if __name__ == "__main__":
    main()

