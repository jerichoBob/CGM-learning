description = """
This tool tests out a few things:
1. if I am combining the spectra from the 7 observations correctly (using inverse variance weighting)
2. if there is noise in a specific observation I need to be aware of
3. it will help me debug the fractional-pixel weighting algorithm for spectra extraction
4. error(stddev) or variance coming from the variance cube? remember when you combine the spectra, you add the variance, not the stddev.

"""

import os, sys
import numpy as np

import bobutils.layout_utils as lu
import bobutils.utils as bu
import bobutils.fileio as bio


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors


import argparse

from astropy.wcs import WCS, FITSFixedWarning
import warnings

warnings.simplefilter('ignore')
warnings.filterwarnings('ignore', category=FITSFixedWarning)

global_nb_min = 4676. 
global_nb_max = 4696. 
global_cmap = 'gnuplot'
global_lw = 0.5


def label_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def create_subplot(fig, gs, row, col, obs, oindex):
    gs_sub = gs[row, col].subgridspec(1, 3)
    ax_container = fig.add_subplot(gs[row, col])
    # ax_container.set_facecolor('lightblue')
    ax_container.set_xticks([])  # Remove x ticks
    ax_container.set_yticks([])  # Remove y ticks
    ax_container.set_xticklabels([])  # Remove x tick labels
    ax_container.set_yticklabels([])  # Remove y tick labels
    ax_container.spines.right.set_visible(False)
    ax_container.spines.bottom.set_visible(False)
    ax_container.spines.left.set_visible(False)

    ax_image = fig.add_subplot(gs_sub[0], projection=obs[oindex].wcs_f)
    ax_image.coords[0].set_ticks_visible(False)  # Hide RA ticks
    ax_image.coords[1].set_ticks_visible(False)  # Hide Dec ticks
    ax_image.coords[0].set_ticklabel_visible(False)  # Hide RA tick labels
    ax_image.coords[1].set_ticklabel_visible(False)  # Hide Dec tick labels
    ax_image.set_xlabel('')
    ax_image.set_ylabel('')
    # Optional: Hide the spines
    for spine in ax_image.spines.values():
        spine.set_visible(False)

    ax_flux = fig.add_subplot(gs_sub[1:])
    # ax_flux.set_xlabel('wavelength (Angstroms)')
    ax_flux.set_ylabel('Flux', color='k', rotation=90)  # Primary y-axis label, rotated
    ax_flux.tick_params(axis='y', labelcolor='k')

    ax_error = ax_flux.twinx()

    # Plot the secondary data and customize the secondary y-axis
    ax_error.set_ylabel('Error', color='r', rotation=90)  # Secondary y-axis label, rotated
    ax_error.tick_params(axis='y', labelcolor='r')  # Set the tick color to red

    return [ax_image, ax_container, ax_flux, ax_error]


def create_mpl_layout(plt, wcs_ref, obs):
    """ top level gridspec is 3 colums by 4 rows"""
    fig = plt.figure(figsize=(10,7), dpi=150)
    rows = 4
    cols = 3
    gs = gridspec.GridSpec(rows, cols, figure=fig)
    ax_ref_image = fig.add_subplot(gs[0:2, 0:1], projection=wcs_ref)
    obs_count = 7 # 7 observations
    ax_obs = []

    ax_obs.append(create_subplot(fig, gs, row=0, col=1, obs=obs, oindex=0))
    ax_obs.append(create_subplot(fig, gs, row=0, col=2, obs=obs, oindex=1))
    ax_obs.append(create_subplot(fig, gs, row=1, col=1, obs=obs, oindex=2))
    ax_obs.append(create_subplot(fig, gs, row=1, col=2, obs=obs, oindex=3))
    ax_obs.append(create_subplot(fig, gs, row=2, col=0, obs=obs, oindex=4))
    ax_obs.append(create_subplot(fig, gs, row=2, col=1, obs=obs, oindex=5))
    ax_obs.append(create_subplot(fig, gs, row=2, col=2, obs=obs, oindex=6))


    ax_combined_spec = fig.add_subplot(gs[3, :])

    # label_axes(fig)
    return fig, ax_ref_image, ax_obs, ax_combined_spec

def load_sightlines():
    """"""
    points = [
        (29,36, 3), # just doing 1 sightline for now
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

def show_ref_image(ax, image, title="White Light", xh_lim=None, yh_lim=None):
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=0)
    ax.coords[0].set_ticks_visible(False)  # Hide RA ticks
    ax.coords[1].set_ticks_visible(False)  # Hide Dec ticks
    ax.coords[0].set_ticklabel_visible(False)  # Hide RA tick labels
    ax.coords[1].set_ticklabel_visible(False)  # Hide Dec tick labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    # Optional: Hide the spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title)
    if xh_lim is not None: ax.set_xlim(xh_lim)
    if yh_lim is not None: ax.set_ylim(yh_lim)
    ax.grid()


def extract_spectra_from_obs(sl_radec, observations):
    """
    Extracts (but doesn't combine) spectra from the collection of observations, contained within the sl_radec box.   
    ie. observations = [class Observation, class Observation, ...]
    params:
        * observations: an list of Observation objects
        * sl_radec[0]: the ras of the extraction box
        * sl_radec[1]: the decs of the extraction box
    
    returns: 
        * XSpectrum1D flux and variance objects
    """
    specs = []
    ras = sl_radec[0]
    decs = sl_radec[1]
    print("=-"*40)
    print(f"=============== BEGINNING EXTRACTION FOR RA:{ras[0]} DEC:{decs[0]} ==================")
    for ob in observations:
        # hdr_f, flux, hdr_v, var, wave, flux_file = ob
        wcs_cur = WCS(ob.hdr_f).celestial
        xs, ys = wcs_cur.world_to_pixel_values(ras, decs)
        print(f"ras={ras}, decs={decs}  -->  xs={xs}, ys={ys}")
        sp = bu.fractional_flux_var_spectra(ob.wave, ob.flux, ob.var, xs, ys)
        specs.append(sp)
    return specs

def show_combined_spectra(ax_spec, combined, color_sightline='blue', color_error='red', color_snr='k', co_begin=4864, co_end=4914, xs=[], ys=[], szs=[]):
    """ show all of the sightline spectra on their own axes """

    spec, var, wave = combined

    ax_spec.plot(wave, spec, '-', color=color_sightline, linewidth=global_lw)
    ax_spec.set_xlabel('wavelength (Angstroms)')
    ax_spec.set_ylabel('Flux', color='k', rotation=90)  # Primary y-axis label, rotated
    ax_spec.tick_params(axis='y', labelcolor='k')
    ax_error = ax_spec.twinx()
    ax_error.plot(wave, var, '-', color=color_error, linewidth=global_lw)
    ax_error.set_ylabel('Error', color='r', rotation=90)  # Secondary y-axis label, rotated
    ax_error.tick_params(axis='y', labelcolor='r')  # Set the tick color to red
    snr2 = bu.sig_figs(bu.signal_to_noise2(wave, spec, co_begin=co_begin, co_end=co_end), 5)
    # snr3 = bu.sig_figs(bu.signal_to_noise3(wave, spec, var, co_begin=co_begin, co_end=co_end), 5)
    # print(f"#0  SNR2: {snr2} SNR3: {snr3}")
    print(f"#0  SNR2: {snr2}")

    ax_spec.text(0.05, 0.9, 
            "Sightline: "+str(0),
            color=color_sightline, 
            fontsize = 12, 
            ha='left', va='top',
            transform=ax_spec.transAxes)
    ax_spec.text(0.05, 0.8, 
            # f"#0  SNR2: {snr2} SNR3: {snr3}",
            f"#0  SNR2: {snr2}",
            color=color_snr, 
            fontsize = 12, 
            ha='left', va='top',
            transform=ax_spec.transAxes)    
    ax_spec.text(0.05, 0.7, 
            f'Coord: ({xs[0]},{ys[0]})',
            color=color_snr, 
            fontsize = 12, 
            ha='left', va='top',
            transform=ax_spec.transAxes)   
    ax_spec.text(0.05, 0.6, 
            f'Aperture: {szs[0]} x {szs[0]}',
            color=color_snr, 
            fontsize = 12, 
            ha='left', va='top',
            transform=ax_spec.transAxes)  
    ax_spec.set_facecolor('darkgrey')

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
parser.add_argument('-l', '--list', action='store_true', help='list the available redshifts')

def handle_commandline_args():
    """ """
    pass

def main(): 
    args = parser.parse_args()

    if args.list:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    xs, ys, szs = load_sightlines()
    cnt = len(xs)
    wl_image, wcs_ref = bio.load_original_cube(global_nb_min, global_nb_max)
    sl_radecs = radecs_from_sightline_boxes(wcs_ref, xs, ys, szs)
    obs = bu.get_corrected_kcwi_data(global_nb_min, global_nb_max)
    ocnt = len(obs)

    specs = []

    fig, ax_ref_image, ax_obs, ax_combined_spec = create_mpl_layout(plt, wcs_ref, obs)
    fig.suptitle(f"J1429+1202 Single Sightline Extraction", fontsize=16)

    sl_radec = sl_radecs[0] # let's plot for the first sightline only
    specs = extract_spectra_from_obs(sl_radec, obs)
    combined = bu.combine_spectra_ivw(specs)            

    show_ref_image(ax_ref_image, wl_image)
    bu.plot_sightlines_wcs(ax_ref_image, wcs_ref, sl_radecs, color='w-', lw=0.5)
    
    # plot the flux and error/variance/stddev/wtf for each observation
    for i in range(ocnt):
        ax_image, ax_container, ax_flux, ax_error = ax_obs[i]
        ax_container.set_title(f"Obs {i+1}")
        ax_image.imshow(obs[i].wl_k, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=0)
        ax_image.grid()

        sp = specs[i]
        ax_flux.plot(sp.wavelength, sp.flux, '-', color='k', lw=0.5)
        ax_error.plot(sp.wavelength, sp.sig, '-', color='r', lw=0.5)

        bu.plot_sightlines_wcs(ax_image, obs[i].wcs_f, sl_radecs, color='w-', lw=0.5, show_label=False)

    show_combined_spectra(ax_combined_spec, combined, 
                           color_sightline='blue', 
                           color_error='red', 
                           color_snr='k', 
                           co_begin=4864, co_end=4914, 
                           xs=xs, ys=ys, szs=szs)    
    mng = plt.get_current_fig_manager()
    if hasattr(mng, 'window'):
        mng.window.state('zoomed') 
    plt.subplots_adjust(left=0.038, bottom=0.057, right=0.917, top=0.929, wspace=0.152, hspace=0.244)
    plt.show()

def show_max_against_obs():
    xs, ys, szs = load_sightlines()
    cnt = len(xs)
    wl_image, wcs_ref = bio.load_original_cube(global_nb_min, global_nb_max)
    sl_radecs = radecs_from_sightline_boxes(wcs_ref, xs, ys, szs)
    obs = bu.get_corrected_kcwi_data(global_nb_min, global_nb_max)
    ocnt = len(obs)
    specs = []

    fig = plt.figure(figsize=(10,7), dpi=150)
    rows = 1
    cols = 2
    obs_index = 0
    ob = obs[obs_index]
    gs = gridspec.GridSpec(rows, cols, figure=fig)
    ax_image = fig.add_subplot(gs[0], projection=ob.wcs_f)
    ax_image.imshow(ob.wl_k, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=0)
    ax_image.grid()
    bu.plot_sightlines_wcs(ax_image, ob.wcs_f, sl_radecs, color='r-', lw=0.5)



    sl_radec = sl_radecs[0] # the RA/DECs bounds for the first sightline
    ras = sl_radec[0]
    decs = sl_radec[1]
    wcs_cur = WCS(ob.hdr_f).celestial
    xs, ys = wcs_cur.world_to_pixel_values(ras, decs)
    # ----- from bu.fractional_flux_var_spectra -----------------------
    shape = ob.flux.shape[1:]  # Assuming flux shape is (wavelength, y, x)
    print(f"flux shape={shape}")
    box_corners = np.array([
        [xs[0], ys[0]],
        [xs[1], ys[0]],
        [xs[1], ys[1]],
        [xs[0], ys[1]]
    ])
    mask = bu.create_fractional_mask(box_corners, shape)
    print(f"mask shape={mask.shape}")
    
    ax_image = fig.add_subplot(gs[1], projection=ob.wcs_f)
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_gray', ['white', 'black'])
    x_corners = box_corners[:, 0]
    y_corners = box_corners[:, 1]
    ax_image.scatter(x_corners, y_corners)


    ax_image.imshow(mask, origin='lower', interpolation='nearest', cmap=cmap, vmin=0)
    ax_image.grid()
    bu.plot_sightlines_wcs(ax_image, ob.wcs_f, sl_radecs, color='r-', lw=0.5)


    mng = plt.get_current_fig_manager()
    if hasattr(mng, 'window'):
        mng.window.state('zoomed') 
    plt.subplots_adjust(left=0.038, bottom=0.057, right=0.917, top=0.929, wspace=0.152, hspace=0.244)
    plt.show()


if __name__ == "__main__":
    main()
    # show_max_against_obs()
