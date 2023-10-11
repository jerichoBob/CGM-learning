import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from astropy.wcs import WCS
from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools import extract_weighted_spectrum as ke
from kcwitools import image as im
from kcwitools import utils as kcwi_u

import warnings

# def ahmeds_mpl_params():
    #### Modifying Matplotlib parameters
    # mpl.rcParams.update(mpl.rcParamsDefault)
    # # mpl.rc('text',usetex=True) # don't want to install latex for this environment (yet)
    # XSMALL_SIZE = 10
    # SMALL_SIZE = 12
    # MEDIUM_SIZE = 14
    # LARGE_SIZE = 16
    # XLARGE_SIZE = 18
    # XXLARGE_SIZE = 24

    # plt.rc('font',size=LARGE_SIZE)
    # plt.rc('axes',titlesize=XLARGE_SIZE)
    # plt.rc('axes',labelsize=XLARGE_SIZE)
    # plt.rc('axes',labelweight=700)
    # plt.rc('axes',titleweight=700)
    # plt.rc('xtick',labelsize=XLARGE_SIZE)
    # plt.rc('ytick',labelsize=XLARGE_SIZE)
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.rc('figure',titlesize=XXLARGE_SIZE)

    # tdir = 'out'

    # major = 5.0
    # minor = 3.0
    # plt.rcParams["xtick.minor.visible"] =  True
    # plt.rcParams["ytick.minor.visible"] =  True
    # plt.rcParams["xtick.major.size"] = major
    # plt.rcParams["ytick.major.size"] = major
    # plt.rcParams["xtick.minor.size"] = minor
    # plt.rcParams["ytick.minor.size"] = minor

    # plt.rcParams['xtick.direction'] = tdir
    # plt.rcParams['ytick.direction'] = tdir
    # mpl.rcParams['axes.linewidth'] = 1.0
    # mpl.rcParams['font.family'] = 'serif'
    # mpl.rcParams["savefig.facecolor"] = 'white'

def transform_pixels(x1, y1, wcs1, wcs2):
    r, d = wcs1.pixel_to_world_values(x1, y1)
    x2, y2 = wcs2.world_to_pixel_values(r, d)
    return x2, y2
def find_fits(basedir, endswith):
    return (os.path.join(root, file)
        for root, dirs, files in os.walk(basedir)
            for file in files if root == basedir and file.lower().endswith(endswith))
def find_files_ending_with(directory, endswith):
    # List all files in the specified directory
    files = os.listdir(directory)
    # Filter the list to only include files that end with the specified suffix
    matching_files = [file for file in files if file.endswith(endswith)]
    return matching_files

def find_index_exceeding_value(lst, val):
    indices = [index for index, value in enumerate(lst) if value >= val]
    return indices[0] if indices else None

def read_and_prep_flux_var_data(flux_file, var_file, minwave, maxwave):
    """ This method reads the flux and var cubes for a specific observation, and cleans them up before returning """
    hdr, flux = open_kcwi_cube(flux_file)
    _, var = open_kcwi_cube(var_file)
    wave = build_wave(hdr)

    # first do a little data cleanup
    var[np.isnan(flux)]=1.
    flux[np.isnan(flux)]=0.0000

    # print("="*50)
    low_index = find_index_exceeding_value(wave, minwave)
    high_index = find_index_exceeding_value(wave, maxwave)
    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]

    print(f"wave[{low_index}]:{wave[low_index]}  wave[{high_index}]:{wave[high_index]}")
    print(f"shape of wave: {wave.shape}")
    print(f"shape of flux: {flux.shape}")
    print(f"shape of var: {var.shape}")
    print("-"*50)
    whitelight= im.build_whitelight(hdr, flux, minwave=minwave, maxwave=maxwave)
    wave = wave[slices]
    flux = flux[slices,:,:]
    var = var[slices,:,:]
    print(f"shape of wave: {wave.shape}")
    print(f"shape of flux: {flux.shape}")
    print(f"shape of var: {var.shape}")    
    print("-="*25)

    # print(f"wavelength range:{wave[0]}-{wave[-1]}")
    # print("="*50)

    # Create a narrow-band whitelight image to plot
    return flux, var, wave, whitelight, hdr

def label_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def create_images_layout(plt, count=6):
    """returns an arround of {count} axes for images to be placed within"""
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
    ax.plot(xs, ys, 'gs')
    hb = box_size // 2
    for i in range(len(xs)):
        ax.plot([xs[i] - hb, xs[i] + hb, xs[i] + hb, xs[i] - hb, xs[i] - hb], 
                [ys[i] - hb, ys[i] - hb, ys[i] + hb, ys[i] + hb, ys[i] - hb], color, lw=lw)
    

def display_wl_image(ax, title, wl_image, xs, ys, xh_lim, yh_lim, hdr_f, draw_axis=False):
    """displays the whitelight image with the sightline(s) overlaid"""
    box_size_arcsecond = 2
    box_size = int(box_size_arcsecond / (hdr_f["CD2_2"] * 3600))

    cmap = 'gnuplot'
    # valid values for cmap; supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'
    # playing around with the color map
    # print(f"wl_image.shape={wl_image.shape}")
    ax.imshow(wl_image, origin='lower', interpolation="nearest", cmap=cmap, vmin=0, vmax=50)
    # plt.colorbar(im)

    plot_sightlines(ax, xs, ys, box_size, color='w-', lw=1.5)

    ax.set_title(title, fontstyle='oblique', fontsize='small')
    ax.set_xlim(xh_lim)
    ax.set_ylim(yh_lim)
    if draw_axis:
        ax.set_xlabel("RA")
        ax.set_ylabel("Decl.")

    ax.grid()


configs = {
    "base_dir":"/Users/robertseaton/rwseaton@ncsu.edu - Google Drive/My Drive/Astro Research/KCWI_1429_Data/Fluxed",
    "observations":[
        {
            "flux":"KB.20180708.23055_icubes_corrected_flux.fits",
            "var":"KB.20180708.23055_icubes_corrected_var.fits",
            "pt_x": 17, 
            "pt_y": 39,
            "nb_wave_lo": 4645.0,
            "nb_wave_hi": 4659.0,
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


def main():
    basedir = configs["base_dir"]
    observations = configs["observations"]
    flux_files = [os.path.join(basedir,obs['flux']) for obs in observations if 'flux' in obs]
    var_files = [os.path.join(basedir,obs['var']) for obs in observations if 'var' in obs]
    wl_low = 3500.
    wl_high = 5500. 
 
    # only one sightline for now -- just checking our numbers
    xc0 = [17]
    yc0 = [39]
    spec_sum = None # this needs to have the same shape as the spectrum shape
    err_sum = None # this needs to have the same the spectrum shape
    count = len(flux_files)

    fig = plt.figure(figsize=(14,10))
    fig.suptitle(f"{count} observations of 1429")
    gs = gridspec.GridSpec(12, count, figure=fig)
    ax_spec = fig.add_subplot(gs[9:11, :])    
    
    wcs_ref = None # our first observiation will go here, and we will use it to transform the other WCS's
    wcs_cur = None # this will be the current WCS we are working with

    for i in range(count):
        ffile = flux_files[i]
        vfile = var_files[i]
    
        flux, var, wave, wl_image, hdr_f = read_and_prep_flux_var_data(ffile, vfile, wl_low, wl_high)
        if wcs_ref is None:
            wcs_ref = WCS(hdr_f).celestial
            ra, dec = wcs_ref.pixel_to_world_values(xc0, yc0)
            wcs_cur = wcs_ref
        else:
            wcs_cur = WCS(hdr_f).celestial

        sh = wl_image.shape
        print(f"sh={sh}")
        xh_lim, yh_lim = transform_pixels([0, sh[1]], [0, sh[0]], wcs_ref, wcs_cur)
        xc, yc = wcs_cur.world_to_pixel_values(ra, dec)
        xc = np.int32(xc)
        yc = np.int32(yc)
        print(f"adjusted pixels for wcs_cur: ({xc}, {yc}) =============================")

        ax2 = fig.add_subplot(gs[0:9, i], projection= wcs_cur)
        # ax2 = fig.add_subplot(gs[0:9, i])
        # ax2.set_axis_off()
        # display the WL image
        title = ffile.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        display_wl_image(ax2, title, wl_image, xc, yc, xh_lim, yh_lim, hdr_f)
        box_size_arcsecond = 2
        box_size = int(box_size_arcsecond / (hdr_f["CD2_2"] * 3600))
        hb = box_size // 2

        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            print(f"extracting spectrum for {xc0[0]},{yc0[0]}")
            print(f"shape of flux: {flux.shape}")
            flux_cut = flux[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
            var_cut  =  var[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
            print(f"shape of flux_cut: {flux_cut.shape}")
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, wave, weights='Data')
            flux_1d = sp.flux
            wave_1d = sp.wavelength
            err_1d = sp.sig
            if spec_sum is None:
                spec_sum = flux_1d
                err_sum = err_1d
            else:
                spec_sum = np.add(spec_sum, flux_1d) 
                err_sum  = np.add(err_sum, err_1d) 
            
        # add spectrum / variance to our stack

    # outside the loop
    lw = 0.5
    ax_spec.plot(wave_1d, spec_sum, 'b-', lw=lw, label="Flux")
    ax_spec.plot(wave_1d, err_sum, 'r-', lw=lw, label="Error")
    ax_spec.set_title("Stacked Light-Weighted Spectrum for ("+str(xc0[0])+","+str(yc0[0])+f") across {count} Observations")

    # display the spectrum and variance
    plt.subplots_adjust(left=0.079, bottom=0.02, right=0.967, top=0.952, wspace=0.071, hspace=0.317)

    plt.show()


main()
