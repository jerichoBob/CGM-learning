import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy import units as u
from astropy.utils.exceptions import AstropyWarning

from kcwitools import io as kcwi_io
from kcwitools import image as im
from kcwitools import utils as kcwi_u

from linetools.utils import radec_to_coord
import bobutils.utils as bu
import bobutils.fileio as bfi

import warnings
# global variables
base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/Analysis/J1429"
# cmap = 'gnuplot'
# cmap = plt.get_cmap('gray')
# cmap = plt.cm.viridis
# valid values for cmap; supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'
# 'magma', 'inferno', 'plasma', 'viridis', 'cividis', 'twilight', 'twilight_shifted', 'turbo',
global_cmap = 'gnuplot'
global_vmin = 0
global_vmax = 80
# for narrow band whitelight images
global_nb_min = 4676. 
global_nb_max = 4696. 

def ahmeds_mpl_params():
    #### Modifying Matplotlib parameters
    mpl.rcParams.update(mpl.rcParamsDefault)
    # mpl.rc('text',usetex=True) # don't want to install latex for this environment (yet)
    XSMALL_SIZE = 10
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    LARGE_SIZE = 16
    XLARGE_SIZE = 18
    XXLARGE_SIZE = 24

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=XLARGE_SIZE)
    plt.rc('axes',labelsize=XLARGE_SIZE)
    plt.rc('axes',labelweight=700)
    plt.rc('axes',titleweight=700)
    plt.rc('xtick',labelsize=XLARGE_SIZE)
    plt.rc('ytick',labelsize=XLARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.rc('figure',titlesize=XXLARGE_SIZE)

    tdir = 'out'

    major = 5.0
    minor = 3.0
    plt.rcParams["xtick.minor.visible"] =  True
    plt.rcParams["ytick.minor.visible"] =  True
    plt.rcParams["xtick.major.size"] = major
    plt.rcParams["ytick.major.size"] = major
    plt.rcParams["xtick.minor.size"] = minor
    plt.rcParams["ytick.minor.size"] = minor

    plt.rcParams['xtick.direction'] = tdir
    plt.rcParams['ytick.direction'] = tdir
    mpl.rcParams['axes.linewidth'] = 1.0
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams["savefig.facecolor"] = 'white'

def bobs_mpl_params():
    #### Modifying Matplotlib parameters
    mpl.rcParams.update(mpl.rcParamsDefault)
    # mpl.rc('text',usetex=True) # don't want to install latex for this environment (yet)
    XSMALL_SIZE = 10
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    LARGE_SIZE = 16
    XLARGE_SIZE = 18
    XXLARGE_SIZE = 24

    plt.rc('font',size=MEDIUM_SIZE)
    plt.rc('axes',titlesize=SMALL_SIZE)
    plt.rc('axes',labelsize=XSMALL_SIZE)
    plt.rc('axes',labelweight=300)
    plt.rc('axes',titleweight=300)
    plt.rc('xtick',labelsize=XSMALL_SIZE)
    plt.rc('ytick',labelsize=XSMALL_SIZE)
    plt.rc('legend',fontsize=SMALL_SIZE)
    plt.rc('figure',titlesize=MEDIUM_SIZE)

    tdir = 'out'

    major = 4.0
    minor = 2.0
    plt.rcParams["xtick.minor.visible"] =  True
    plt.rcParams["ytick.minor.visible"] =  True
    plt.rcParams["xtick.major.size"] = major
    plt.rcParams["ytick.major.size"] = major
    plt.rcParams["xtick.minor.size"] = minor
    plt.rcParams["ytick.minor.size"] = minor

    plt.rcParams['xtick.direction'] = tdir
    plt.rcParams['ytick.direction'] = tdir
    mpl.rcParams['axes.linewidth'] = 1.0
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams["savefig.facecolor"] = 'white'

def transform_pixels(x1, y1, wcs1, wcs2):
    ra, dec = wcs1.pixel_to_world_values(x1, y1) # from this
    x2, y2 = wcs2.world_to_pixel_values(ra, dec) # to this
    return x2, y2

def find_index_exceeding_value(lst, val):
    indices = [index for index, value in enumerate(lst) if value >= val]
    return indices[0] if indices else None

def label_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def create_images_layout(plt, count=6):
    """returns an array of {count} axes for images to be placed within"""
    fig = plt.figure(figsize=(14,10), layout="constrained", dpi=150)
    # fig.suptitle(title)
    axs_image = []
    gs = gridspec.GridSpec(12, count, figure=fig)
    for i in range(count):
        print(f"index:{i}")
        axs_image.append(fig.add_subplot(gs[0:9, i]))
        label_axes(fig)
    ax_spectrum = fig.add_subplot(gs[9:11, :])
    # return fig, ax_spectrum, axs_image
    return fig, axs_image

def plot_sightlines(ax, xs, ys, box_size, color='w-', lw=0.5):
    # ax.plot(xs, ys, 'gs')
    hb = box_size // 2
    for i in range(len(xs)):
        ax.plot([xs[i] - hb, xs[i] + hb, xs[i] + hb, xs[i] - hb, xs[i] - hb], 
                [ys[i] - hb, ys[i] - hb, ys[i] + hb, ys[i] + hb, ys[i] - hb], color, lw=lw)

def plot_sightlines_wcs(mpl_ax, target_wcs, radecs, color='w-', lw=0.5):
    """ this isn't quite right yet. need to figure out how ras and decs are packed """
    for i in range(len(radecs)):
        ras, decs = radecs[i] # ras and decs define the edges of the extraction box
        ax, ay = target_wcs.world_to_pixel_values(ras[1], decs[1])  
        bx, by = target_wcs.world_to_pixel_values(ras[0], decs[1])
        cx, cy = target_wcs.world_to_pixel_values(ras[0], decs[0])
        dx, dy = target_wcs.world_to_pixel_values(ras[1], decs[0]) 
        mpl_ax.plot([ax, bx, cx, dx, ax], 
                    [ay, by, cy, dy, ay], color, lw=lw)

    
def display_wl_image(ax, title, wl_image, xs, ys, xh_lim, yh_lim, hdr_f, draw_axis=False):
    """displays the whitelight image with the sightline(s) overlaid"""
    # box_size_arcsecond = 2
    # box_size = int(box_size_arcsecond / (hdr_f["CD2_2"] * 3600))
    ax.imshow(wl_image, origin='lower', interpolation="nearest", cmap=global_cmap, vmin=0)
    # plt.colorbar(im)

    box_size = 3
    # plot_sightlines(ax, xs, ys, box_size, color='w-', lw=1.0)

    ax.set_title(title, fontstyle='oblique', fontsize='small')
    ax.set_xlim(xh_lim)
    ax.set_ylim(yh_lim)
    if draw_axis:
        ax.set_xlabel("RA")
        ax.set_ylabel("Decl.")

    ax.grid()


configs = {
    # "base_dir":"/Users/robertseaton/rwseaton@ncsu.edu - Google Drive/My Drive/Astro Research/KCWI_1429_Data/Fluxed",
    "base_dir":"/Users/robertseaton/Library/CloudStorage/GoogleDrive-rwseaton@ncsu.edu/My Drive/Astro Research/KCWI_1429_Data/Fluxed",
    "observations":[
        {
            "flux":"KB.20180708.23055_icubes_corrected_flux.fits",
            "var":"KB.20180708.23055_icubes_corrected_var.fits",
            "pt_x": 17, 
            "pt_y": 39,
            "nb_wave_lo": 4640.0,
            "nb_wave_hi": 4660.0,
            "nb_index_lo": 1418,
            "nb_index_hi": 1432,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20180708.24349_icubes_corrected_flux.fits",
            "var":"KB.20180708.24349_icubes_corrected_var.fits",
            "pt_x": 16, 
            "pt_y": 43,
            "nb_index_lo": 1420,
            "nb_index_hi": 1434,
            "nb_wave_lo": 4646.02,
            "nb_wave_hi`": 4660.08,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20190601.32100_icubes_corrected_flux.fits",
            "var":"KB.20190601.32100_icubes_corrected_var.fits",
            "pt_x": 18,
            "pt_y": 36,
            "nb_index_lo": 1420,
            "nb_index_hi": 1431,
            "nb_wave_lo": 4646.26,
            "nb_wave_hi`": 4657.01,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20190601.33358_icubes_corrected_flux.fits",
            "var":"KB.20190601.33358_icubes_corrected_var.fits",
            "pt_x": 18, 
            "pt_y": 36,
            "nb_index_lo": 1419,
            "nb_index_hi": 1435,
            "nb_wave_lo": 4645.07,
            "nb_wave_hi`": 4661.1,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20200520.37578_icubes_corrected_flux.fits",
            "var":"KB.20200520.37578_icubes_corrected_var.fits",
            "pt_x": 17, 
            "pt_y": 39,
            "nb_index_lo": 1418,
            "nb_index_hi": 1433,
            "nb_wave_lo": 4644.62,
            "nb_wave_hi`": 4660.42,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20200520.38847_icubes_corrected_flux.fits",
            "var":"KB.20200520.38847_icubes_corrected_var.fits",
            "pt_x": 17, 
            "pt_y": 39,
            "nb_index_lo": 1418,
            "nb_index_hi": 1434,
            "nb_wave_lo": 4645.45,
            "nb_wave_hi`": 4660.39,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        },
        {
            "flux":"KB.20200520.40119_icubes_corrected_flux.fits",
            "var":"KB.20200520.40119_icubes_corrected_var.fits",
            "pt_x": 17, 
            "pt_y": 39,
            "nb_index_lo": 1416,
            "nb_index_hi": 1432,
            "nb_wave_lo": 4643.4,
            "nb_wave_hi`": 4658.92,
            "window_size_x": 0,
            "window_size_y": 0,
            "image_shift_x": 0,
            "image_shift_y": 0
        }
    ]
}

def display_test():
    count = 7
    fig, _ = create_images_layout(plt, count=count)
    fig.suptitle(f"{count} observations of 1429")
    plt.show()




def show_whitelight_image(ax, image, title="White Light", xh_lim=None, yh_lim=None):
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=0)
    ax.set_title(title)
    if xh_lim is not None: ax.set_xlim(xh_lim)
    if yh_lim is not None: ax.set_ylim(yh_lim)

    ax.grid()


def load_hst_image():
    fh = fits.open("/Users/robertseaton/School/github_repos/CGM-learning/data/data_idpw02010_WFC3_IR_F160W_drz.fits")
    hdr = fh[0].header    # Reference image header
    data = fh[0].data     # Reference image 2D data
    wcs = WCS(hdr)        # WCS of the HST image
    return data, wcs

def show_hst_image(ax, image):
    # ax.imshow(image, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=1, vmax=8)
    # print(f"HST image.max(): {image.max()}")
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=global_cmap, vmin=1, vmax=4)
    ax.set_title("HST Image")
    ax.set_xlim([470, 594])
    ax.set_ylim([606, 720])
    ax.grid()

def load_sightlines():
    """"""
    points = [ # these were defined against the original kcwi datacube
        (29,36, 3), # 0
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

def load_KB20180708_23055_sightlines(wcs):
    """ defines a new set of boxes and returns the radecs"""
    boxes = [ # these are defined against the new kcwi observations, specifically KB.20180708.23055_icubes_corrected_flux.fits
        (13,50,3), #0
        (13,53,3), #1
        (16,52,3), #2
        (16,55,3), #3
        (19,53,3), #4
        (19,56,3), #5
        (23,47,5), #6
        (16,39,5), #7
    ]
    radecs = []
    box_count = len(boxes)
    for ndx in range(box_count): # loop through the sightlines
        xval,yval,szval = boxes[ndx]
        print(f"=== sightline #{ndx} pt_xs: {xval} pt_ys: {yval} szs: {szval}")
        xs, ys = bu.box_corners(xval, yval, deltax=szval, deltay=szval) #NOTE: still a question whether we should use these for extraction, or just drawing the box against the wl image
        print(f"=== sightline box corners   xs: {xs} ys: {ys}")
        radec = wcs.pixel_to_world_values(xs, ys)
        print(f"radec: {radec}")
        radecs.append(radec)
    return radecs

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

def stack_spectra():
    bobs_mpl_params()
    # get the 7 observation flux and var files
    observations = bu.get_corrected_kcwi_data(global_nb_min, global_nb_max)
    file_cnt = len(observations)


    # only one sightline for now -- just checking our numbers
    # xc0 = [17]
    # yc0 = [39]
    global_lw = 1.0
    fig = plt.figure(figsize=(14,10))
    fig.suptitle(f"{file_cnt} observations of 1429", fontsize=16)
    gs = gridspec.GridSpec(21, 2*file_cnt, figure=fig)
    wl_image, wcs_kcwi = bfi.load_narrowband_reference_image(global_nb_min, global_nb_max)
    ax_kcwi = fig.add_subplot(gs[0:6, 3:6], projection=wcs_kcwi)
    show_whitelight_image(ax_kcwi, wl_image, title="KCWI White Light")

    original_method = True
    if original_method:
        xs, ys, szs = load_sightlines()
        sl_radecs = radecs_from_sightline_boxes(wcs_kcwi, xs, ys, szs) # these are the reference radecs for the sightlines
    else:
        wcs_ref = WCS(observations[0].hdr_f).celestial # using the first observation as the reference WCS
        sl_radecs = load_KB20180708_23055_sightlines(wcs_ref)

    plot_sightlines_wcs(ax_kcwi, wcs_kcwi, sl_radecs, color='w-', lw=global_lw)

    hst_image, wcs_hst = load_hst_image()
    ax_hst = fig.add_subplot(gs[0:6, 8:11], projection=wcs_hst)
    show_hst_image(ax_hst, hst_image)
    plot_sightlines_wcs(ax_hst, wcs_hst, sl_radecs, color='w-', lw=global_lw)
    
    wcs_ref = None # our first observiation will go here, and we will use it to transform the other WCS's
    wcs_cur = None # this will be the current WCS we are working with
    
    ybot = 15.5 # this doesn't work as it causes the result to not match with what's in the header - need to fix if we want to crop
    ytop = 80.5

    for ii,ob in enumerate(observations):
        # if ii > 0: break
        print(f"=======================  Observation {ii+1}  =============================")
    
        # about the points we are interested in
        with warnings.catch_warnings():
            # Ignore a warning on using DATE-OBS in place of MJD-OBS
            warnings.simplefilter('ignore', AstropyWarning)

            # warnings.filterwarnings('ignore', 
            #                         message="""'datfix' made the change 'Set MJD-OBS to 58989.000000 from DATE-OBS.""",
            #                         category=FITSFixedWarning)
            if wcs_ref is None:
                wcs_ref = WCS(ob.hdr_f).celestial
                # ra, dec = wcs_ref.pixel_to_world_values(xc0, yc0)
                # print(f"ra={ra}, dec={dec} xc0={xc0}, yc0={yc0}")
                # sky = radec_to_coord((ra[0], dec[0]))
                # ra_str = sky.ra.to_string(u.hour, precision=2, alwayssign=True, pad=True)
                # dec_str = sky.dec.to_string(u.degree, precision=2, alwayssign=True, pad=True)
                # print(f"sky ra={ra_str}, dec={dec_str}")
                wcs_cur = wcs_ref

            else:
                wcs_cur = WCS(ob.hdr_f).celestial

        sh = ob.wl_k.shape
        # print(f"sh={sh}")
        xh_lim, yh_lim = transform_pixels([0, sh[1]], [0, sh[0]], wcs_ref, wcs_cur)
        # xc, yc = wcs_cur.world_to_pixel_values(ra, dec)
        # xc = np.int32(xc)
        # yc = np.int32(yc)

        ax2 = fig.add_subplot(gs[7:16, 2*ii:2*(ii+1)], projection= wcs_cur)
        # title = ob.f_file.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        title = ob.f_file.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        show_whitelight_image(ax2, ob.wl_k, title=title, xh_lim=xh_lim, yh_lim=yh_lim)
        # display_wl_image(ax2, title, wl_image, xc, yc, xh_lim, yh_lim, ob.hdr_f)
        plot_sightlines_wcs(ax2, wcs_cur, sl_radecs, color='w-', lw=global_lw)

    # outside the loop
    plot_spectrum = True
    if plot_spectrum:
        
        sl_radec = sl_radecs[1] # let's plot for the first sightline only
        box_ras, box_decs = sl_radec
        spec, var, wave = bu.extract_spectra_from_observations(observations, sl_radec)

        lw = 0.5
        ax_spec = fig.add_subplot(gs[17:21, :])    
        # ax_err = ax_spec.twinx()
        ax_err = ax_spec
        ax_err.plot(wave, var, 'r-', lw=lw, label="Error")
        ax_spec.plot(wave, spec, 'b-', lw=lw, label="Flux")
        sky = radec_to_coord((box_ras[0], box_decs[0]))
        ra_str = sky.ra.to_string(u.hour, precision=2, alwayssign=True, pad=True)
        dec_str = sky.dec.to_string(u.degree, precision=2, alwayssign=True, pad=True)

        title = f"""Sightline {ra_str} {dec_str}
    Inverse Variance-weighted Spectrum summed over {file_cnt} Observations"""
        ax_spec.set_title(title)
        ax_spec.set_xlabel("Wavelength")
        ax_spec.set_ylabel("Flux & Error")
        # ax_err.set_ylabel("Error")
        ax_spec.legend()
    # ax_spec.legend(loc='upper right', bbox_to_anchor=(1, 1))
    # ax_err.legend(loc='upper right', bbox_to_anchor=(1, .9))

    # display the spectrum and variance
    plt.subplots_adjust(left=0.079, bottom=0.043, right=0.967, top=0.952, wspace=0.638, hspace=0.35)
    
    # mng = plt.get_current_fig_manager()
    mng = plt.get_current_fig_manager()
    if hasattr(mng, 'window'):
        mng.window.state('zoomed') #works fine on Windows!

        # fig_manager.window.showMaximized()
    plt.show()
    # mng.full_screen_toggle()


def read_headers():
    basedir = configs["base_dir"]
    observations = configs["observations"]
    flux_files = [os.path.join(basedir,obs['flux']) for obs in observations if 'flux' in obs]
    var_files = [os.path.join(basedir,obs['var']) for obs in observations if 'var' in obs]
    count = len(flux_files)
    for i in range(count):
        ffile = flux_files[i]
        vfile = var_files[i]
    
        # Create a narrow-band whitelight image to plot
        hdr_f, _ = kcwi_io.open_kcwi_cube(ffile)
        f_name = ffile.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        v_name = vfile.split("/")[-1].split("_icubes_corrected_var.fits")[0]

        hdr_v, _ = kcwi_io.open_kcwi_cube(vfile)
        print("-="*40)
        print( f"f_name={f_name} hdr_f[DATE-BEG]:{hdr_f['DATE-BEG']} hdr_f[DATE-END]:{hdr_f['DATE-END']}")
        print( f"v_name={v_name} hdr_v[DATE-BEG]:{hdr_v['DATE-BEG']} hdr_v[DATE-END]:{hdr_v['DATE-END']}")


if __name__ == "__main__":
    stack_spectra()
