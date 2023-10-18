import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy import units as u

from kcwitools.io import open_kcwi_cube
from kcwitools.utils import build_wave
from kcwitools import extract_weighted_spectrum as ke
from kcwitools import image as im
from linetools.utils import radec_to_coord
# from kcwitools import utils as kcwi_u
import bobutils.utils as bu

import warnings

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

def read_and_prep_flux_var_data(flux_file, var_file, minwave, maxwave, ybot, ytop):
    """ This method reads the flux and var cubes for a specific observation, and cleans them up before returning """
    hdr, flux = open_kcwi_cube(flux_file)
    _, var = open_kcwi_cube(var_file)
    wave = build_wave(hdr)

    # first do a little data cleanup
    var[np.isnan(flux)]=1.
    flux[np.isnan(flux)] = 0.0000

    slices = np.where((wave >= minwave) & (wave <= maxwave))[0]
    wave = wave[slices]
    flux = flux[slices,int(ybot):int(ytop),:]
    var = var[slices,int(ybot):int(ytop),:]
    print(f"flux.shape={flux.shape}") # (slices, y, x)

    return hdr, flux, var, wave

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
    # ax.plot(xs, ys, 'gs')
    hb = box_size // 2
    for i in range(len(xs)):
        ax.plot([xs[i] - hb, xs[i] + hb, xs[i] + hb, xs[i] - hb, xs[i] - hb], 
                [ys[i] - hb, ys[i] - hb, ys[i] + hb, ys[i] + hb, ys[i] - hb], color, lw=lw)
    

def display_wl_image(ax, title, wl_image, xs, ys, xh_lim, yh_lim, hdr_f, draw_axis=False):
    """displays the whitelight image with the sightline(s) overlaid"""
    # box_size_arcsecond = 2
    # box_size = int(box_size_arcsecond / (hdr_f["CD2_2"] * 3600))
    # cmap = 'gnuplot'
    # cmap = plt.get_cmap('gray')
    cmap = plt.cm.viridis
    # valid values for cmap; supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'
    # playing around with the color map
    # print(f"wl_image.shape={wl_image.shape}")
    # ax.imshow(wl_image, origin='lower', interpolation="nearest")
    ax.imshow(wl_image, origin='lower', interpolation="nearest", cmap=cmap, vmin=0)
    # ax.imshow(wl_image, origin='lower', interpolation="nearest", cmap=cmap, vmin=0, vmax=80)
    # plt.colorbar(im)
    # ax.set_xlabel(' ')
    # ax.set_ylabel(' ')

    # projection needs to be turned on in order for this to work
    # ax.coords.grid(True, color='white', ls='solid')
    # ax.coords[0].set_axislabel('Galactic Longitude')
    # ax.coords[1].set_axislabel('Galactic Latitude')

    # overlay = ax.get_coords_overlay('fk5')
    # overlay.grid(color='white', ls='dotted')
    # overlay[0].set_axislabel('Right Ascension (J2000)')
    # overlay[1].set_axislabel('Declination (J2000)')


    box_size = 3
    plot_sightlines(ax, xs, ys, box_size, color='w-', lw=1.0)

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


def combine_spectra_ivw(specs):
    """Combine the collection of 1D spectra into a single spectrum using inverse variance weighting"""
    """See https://en.wikipedia.org/wiki/Inverse-variance_weighting"""
    print(f"# of spectra={len(specs)}")
    fluxlen = len(specs[0].flux)
    print(f"fluxlen={fluxlen}")
    flux_tot = np.zeros(fluxlen)
    var_tot = np.zeros(fluxlen)
    for lndx in range(fluxlen): # for each lambda, apply inverse variance weighting
        sum_ysigma = 0.0
        sum_1sigma = 0.0
        for sp in specs: # for each observation (that is, each spectrum)...
            sum_ysigma += sp.flux[lndx] / (sp.sig[lndx] ** 2)
            sum_1sigma += 1.0 / (sp.sig[lndx] ** 2)
        flux_tot[lndx] = sum_ysigma / sum_1sigma
        var_tot[lndx] = 1.0 / sum_1sigma

    return flux_tot, var_tot

def deal_with_wl_image(ax, title, wl_image, xs, ys, xh_lim, yh_lim, hdr_f, draw_axis=False):
    """displays the whitelight image with the sightline(s) overlaid"""
    pass

def stack_spectra():
    bobs_mpl_params()
    # get the 7 observation flux and var files
    basedir = configs["base_dir"]
    observations = configs["observations"]
    flux_files = [os.path.join(basedir,obs['flux']) for obs in observations if 'flux' in obs]
    var_files = [os.path.join(basedir,obs['var']) for obs in observations if 'var' in obs]
    file_cnt = len(flux_files)
 
    # only one sightline for now -- just checking our numbers
    xc0 = [17]
    yc0 = [39]

    # define y-axis cropping bounds for the flux and var cubes
    ybot = 15.5
    ytop = 80.5
    
    fig = plt.figure(figsize=(14,10))
    fig.suptitle(f"{file_cnt} observations of 1429")
    gs = gridspec.GridSpec(12, file_cnt, figure=fig)
    ax_spec = fig.add_subplot(gs[9:12, :])    
    
    wcs_ref = None # our first observiation will go here, and we will use it to transform the other WCS's
    wcs_cur = None # this will be the current WCS we are working with

    specs = []

    for i in range(file_cnt):
        ffile = flux_files[i]
        vfile = var_files[i]
    
        minwave = 3500.
        maxwave = 5500.
        hdr, flux, var, wave = read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop)

        # Create a narrow-band whitelight image to plot
        # nb_min = 4630.
        # nb_max = 4670.
        nb_min = 4676. #4637. #3500. #4546.02
        nb_max = 4696. # 5500. #4760.08
        wl_image = im.build_whitelight(hdr, flux, minwave=nb_min, maxwave=nb_max)
    
        # about the points we are interested in
        with warnings.catch_warnings():
            # Ignore a warning on using DATE-OBS in place of MJD-OBS
            warnings.filterwarnings('ignore', 
                                    message="""'datfix' made the change 'Set MJD-OBS to 58989.000000 from DATE-OBS.""",
                                    category=FITSFixedWarning)
            if wcs_ref is None:
                wcs_ref = WCS(hdr).celestial
                ra, dec = wcs_ref.pixel_to_world_values(xc0, yc0)
                print(f"ra={ra}, dec={dec} xc0={xc0}, yc0={yc0}")
                sky = radec_to_coord((ra[0], dec[0]))
                ra_str = sky.ra.to_string(u.hour, precision=2, alwayssign=True, pad=True)
                dec_str = sky.dec.to_string(u.degree, precision=2, alwayssign=True, pad=True)
                print(f"sky ra={ra_str}, dec={dec_str}")
                wcs_cur = wcs_ref

            else:
                wcs_cur = WCS(hdr).celestial

        sh = wl_image.shape
        # print(f"sh={sh}")
        xh_lim, yh_lim = transform_pixels([0, sh[1]], [0, sh[0]], wcs_ref, wcs_cur)
        xc, yc = wcs_cur.world_to_pixel_values(ra, dec)
        xc = np.int32(xc)
        yc = np.int32(yc)

        ax2 = fig.add_subplot(gs[0:9, i], projection= wcs_cur)
        title = ffile.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        display_wl_image(ax2, title, wl_image, xc, yc, xh_lim, yh_lim, hdr)

        # ax2.axhline(y=ybot, linewidth=1, linestyle="--")
        # ax2.axhline(y=ytop, linewidth=1, linestyle="--")

        box_size_arcsecond = 2
        box_size = int(box_size_arcsecond / (hdr["CD2_2"] * 3600))
        hb = box_size // 2
        print(f"hb={hb} pixels")

        flux_cut = flux[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        var_cut  =  var[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        # extract the weighted spectrum and variance
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, wave, weights='Data')
        specs.append(sp) # add the spectrum to our list

    spec, var = combine_spectra_ivw(specs)            

    # outside the loop
    lw = 0.5
    # ax_err = ax_spec.twinx()
    ax_err = ax_spec
    ax_err.plot(wave, var, 'r-', lw=lw, label="Error")
    ax_spec.plot(wave, spec, 'b-', lw=lw, label="Flux")
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
    plt.subplots_adjust(left=0.079, bottom=0.043, right=0.967, top=0.952, wspace=0.064, hspace=0.0)

    plt.show()

def cropper():
    """this is me figuring out how to crop down a data cube too just the region I want"""
    bobs_mpl_params()
    basedir = configs["base_dir"]
    observations = configs["observations"]
    flux_files = [os.path.join(basedir,obs['flux']) for obs in observations if 'flux' in obs]
    var_files = [os.path.join(basedir,obs['var']) for obs in observations if 'var' in obs]
 
    # only one sightline for now -- just checking our numbers
    xc = [17]
    yc = [39]
    # cropping bounds
    ybot = 15.5
    ytop = 80.5
    spec_sum = None # this needs to have the same shape as the spectrum shape
    err_sum = None # this needs to have the same the spectrum shape
    
    ffile = flux_files[0]
    vfile = var_files[0]
    
    minwave = 3500.
    maxwave = 5500.
    hdr, flux, var, wave = read_and_prep_flux_var_data(ffile, vfile, minwave, maxwave, ybot, ytop)
    
    # Create a narrow-band whitelight image to plot
    nb_min = 4546.02
    nb_max = 4760.08
    wl_image = im.build_whitelight(hdr, flux, minwave=nb_min, maxwave=nb_max)

    sh = wl_image.shape
    xh_lim = [0, sh[1]]
    yh_lim = [0, sh[0]]

    title = ffile.split("/")[-1].split("_icubes_corrected_flux.fits")[0]

    fig = plt.figure(figsize=(14,10))
    fig.suptitle("Cropping 1429")
    gs = gridspec.GridSpec(5, 3, figure=fig)


    # ax_image = fig.add_subplot(gs[:, 0:1])
    wcs_cur = WCS(hdr).celestial
    ax_image = fig.add_subplot(gs[:, 0:1], projection=wcs_cur)
    display_wl_image(ax_image, title, wl_image, xc, yc, xh_lim, yh_lim, hdr)

    # ax_image.axhline(y=ybot, linewidth=1, linestyle="--")
    # ax_image.axhline(y=ytop, linewidth=1, linestyle="--")

    show_spectra = True
    if show_spectra:
        # sightline and spectrum
        box_size_arcsecond = 2
        box_size = int(box_size_arcsecond / (hdr["CD2_2"] * 3600))
        hb = box_size // 2

        # extract the weighted spectrum and variance
        flux_cut = flux[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        var_cut  =  var[:, yc[0]-hb:yc[0]+hb+1, xc[0]-hb:xc[0]+hb+1]
        # print(f"shape of flux_cut: {flux_cut.shape}")
        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter
            warnings.simplefilter('ignore')
            sp = ke.extract_weighted_spectrum(flux_cut, var_cut, wave, weights='Data')

        flux_1d = sp.flux
        wave_1d = sp.wavelength
        err_1d = sp.sig
        spec_sum = flux_1d
        err_sum = err_1d

        lw = 0.5
        ax_spec = fig.add_subplot(gs[2:3, 1:3]) 
        ax_spec.plot(wave_1d, spec_sum, 'b-', lw=lw, label="Flux")
        ax_spec2 = ax_spec.twinx()
        ax_spec2.plot(wave_1d, err_sum, 'r-', lw=lw, label="Error")
        ax_spec.set_title("Flux and Error Spectrum")
        ax_spec.set_xlabel("Wavelength")
        ax_spec.set_ylabel("Flux")
        ax_spec2.set_ylabel("Error")
        # ax_spec.legend(loc='upper right', bbox_to_anchor=(1, 1))
        # ax_spec2.legend(loc='upper right', bbox_to_anchor=(1, .9))
        ax_spec.legend()

    plt.subplots_adjust(left=0.079, bottom=0.043, right=0.967, top=0.952, wspace=0.064, hspace=0.0)

    plt.show()

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
        hdr_f, _ = open_kcwi_cube(ffile)
        f_name = ffile.split("/")[-1].split("_icubes_corrected_flux.fits")[0]
        v_name = vfile.split("/")[-1].split("_icubes_corrected_var.fits")[0]

        hdr_v, _ = open_kcwi_cube(vfile)
        print("-="*40)
        print( f"f_name={f_name} hdr_f[DATE-BEG]:{hdr_f['DATE-BEG']} hdr_f[DATE-END]:{hdr_f['DATE-END']}")
        print( f"v_name={v_name} hdr_v[DATE-BEG]:{hdr_v['DATE-BEG']} hdr_v[DATE-END]:{hdr_v['DATE-END']}")


if __name__ == "__main__":
    stack_spectra()
    # cropper()
