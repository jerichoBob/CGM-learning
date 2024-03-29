import sys, os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.patches as patches
import seaborn as sns
import numpy as np

# modify the python path to include the current directory
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from bobutils import sightlines as bus

# This just draws a single box centered at (x,y) and sz from that center point in the n/e/s/w directions 
def plotbox(plt, x, y, sz, c):
    ax = x - 0.5
    ay = y - 0.5
    plt.plot(
        [ax, ax,    ax-sz, ax-sz, ax],
        [ay, ay-sz, ay-sz, ay,    ay],
        '-', color = c)

def plotcircle(plt, x, y, labels, sz, c):
  for i in range(len(x)):
        colour = 'c' #fix it at cyan for now c[i];
        x_ = x[i]
        y_ = y[i]
        plotbox(plt, x_, y_,  1, colour) # center
        plotbox(plt, x_, y_+1, 1, colour) # north
        plotbox(plt, x_+1, y_, 1, colour) # east
        plotbox(plt, x_, y_-1, 1, colour) # south
        plotbox(plt, x_-1, y_, 1, colour) # west

        plt.text(x[i]-1.5*sz, y[i], labels[i], color=colour)
            
def show_wl_image(ax, image, title="White Light", xh_lim=None, yh_lim=None, cmap='gnuplot', axis_vis=True, show_labels=True):
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=cmap, vmin=0)

    if axis_vis:    
        ax_type = type(ax).__name__
        # print(f"ax_type: {ax_type}")
        if not show_labels:
            if ax_type == 'WCSAxesSubplot' or ax_type == 'WCSAxes':
                ax.coords[0].set_ticks_visible(False)  # Hide RA ticks
                ax.coords[1].set_ticks_visible(False)  # Hide Dec ticks
                ax.coords[0].set_ticklabel_visible(False)  # Hide RA tick labels
                ax.coords[1].set_ticklabel_visible(False)  # Hide Dec tick labels
            elif ax_type == 'AxesSubplot':
                ax.set_xticks([])
                ax.set_yticks([])
            ax.set_xlabel('')
            ax.set_ylabel('')
            # Optional: Hide the spines
            for spine in ax.spines.values():
                spine.set_visible(False)

    ax.set_title(title)
    if xh_lim is not None: ax.set_xlim(xh_lim)
    if yh_lim is not None: ax.set_ylim(yh_lim)
    ax.grid()

def quick_show_wl_image_and_sightlines(image, title="White Light", wcs_ref=None, sightlines=None, xh_lim=None, yh_lim=None, cmap='gnuplot', axis_vis=True, show_label=True):
    fig, ax = plt.subplots(1, 1, figsize=(8,8))
    ax = fig.add_axes((0,0,1,1), projection=wcs_ref, polar=False)

    show_wl_image(ax, image, title, xh_lim, yh_lim, cmap, axis_vis)
    plot_sightlines_wcs(ax, wcs_ref, sightlines, lw=0.5, show_label=show_label)
    
    plt.show()
    
def show_wl_image_and_sightlines(ax, image, title="White Light", wcs_ref=None, sightlines=None, xh_lim=None, yh_lim=None, cmap='gnuplot', axis_vis=True, show_label=True):
    show_wl_image(ax, image, title, xh_lim, yh_lim, cmap, axis_vis)
    plot_sightlines_wcs(ax, wcs_ref, sightlines, lw=0.5, show_label=show_label)
    
    
def shift_coords(xs, ys, dx, dy):
    """ shifts the x,y coords by dx,dy """
    xs = [x+dx for x in xs]
    ys = [y+dy for y in ys]
    return xs, ys


def plot_sightlines_wcs(mpl_ax, wcs, sightlines, lw=0.5, color='w', label="default", show_label=True):
    """ plots the wcs converted sightline extraction boxes on the given mpl axis """

    my_cmap = sns.color_palette("husl", len(sightlines))

    for sl in sightlines:
        x_min, x_max, y_min, y_max = bus.sightline_cuts(sl.radecs, wcs)
        
        xs = [x_min,  x_max, x_max, x_min,  x_min]
        ys = [y_min,  y_min, y_max, y_max,  y_min]
        xs, ys = shift_coords(xs, ys, -0.5, -0.5)
        
        mpl_ax.plot(xs, ys, color=color, lw=lw, alpha=1.0, zorder=1)
        
        if show_label:
            mpl_ax.text(xs[2], ys[2], str(sl.label), color='w', fontsize = 14, ha='center', va='center')

def plot_sightlines_wcs_old(mpl_ax, wcs, sl_radecs, lw=0.5, color='w', label="default", show_label=True):
    """ plots the wcs converted sightline extraction boxes on the given mpl axis """

    my_cmap = sns.color_palette("husl", len(sl_radecs))

    for sl_radec in sl_radecs:
        # print(f"sl_radec: {sl_radec}") 
        # print(f"type(sl_radec): {type(sl_radec)}")
        if type(sl_radec) == bus.Sightline:
            sl_radec = sl_radec.radecs
            x_min, x_max, y_min, y_max = bus.sightline_cuts(sl_radec, wcs)
        else:
            x_min, x_max, y_min, y_max = bus.sightline_cuts(sl_radec, wcs)
        
        xs = [x_min,  x_max, x_max, x_min,  x_min]
        ys = [y_min,  y_min, y_max, y_max,  y_min]
        xs, ys = shift_coords(xs, ys, -0.5, -0.5)
        
        mpl_ax.plot(xs, ys, color=color, lw=lw, alpha=1.0, zorder=1)
        
        if show_label:
            mpl_ax.text(xs[2], ys[2], str(label), color='w', fontsize = 14, ha='center', va='center')
            
def plot_pixel_wcs(mpl_ax, x,y, color):
    """ 
    a little helper that fills in a selected set of pixels
     - pixels array is a list of tuples (x,y) x,y>=0 of the pixels to fill in
    """
    """ 
    let's make a box around a pixel: so for (0,0), we have the following corners
    """

    xs = [x, x+1, x+1, x,   x]
    ys = [y, y,   y+1, y+1, y]
    mpl_ax.fill(xs, ys, alpha=0.5, 
                facecolor=color,
                edgecolor=color, 
                linewidth=1, 
                zorder=1)

def plot_spectra(sl, xspec, label=f"Sightline", figsize=(8,2), dpi=200, fontsize=10):
    """plot_spectra plots the associated XSpectrum1D with a single sightline. 
    Args:   sl (Sightline): _description_
            xspec object (spec,var,wave): Xspectrum1D object
    """
    wave = xspec.wavelength.value
    spec = xspec.flux.value
    err = xspec.sig.value
        
    plt.close()
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax_spec = fig.add_axes((0,0,1,1))

    ax_spec.plot(wave, spec, '-', color="k", linewidth=0.5, label="flux")
    ax_spec.set_xlabel('$\lambda$ ($\AA$)', )
    ax_spec.set_ylabel('Flux', color='k', rotation=90)  # Primary y-axis label, rotated
    ax_spec.tick_params(axis='y', labelcolor='k')
    ax_spec.plot(wave, err, '-', color="r", linewidth=0.5, label="error")
    ax_spec.legend(loc='upper right', fontsize=fontsize*.8)

    ax_spec.text(0.05, 0.9, 
            f"{label}",
            color="black", 
            fontsize = fontsize, 
            ha='left', va='top',
            transform=ax_spec.transAxes)
 
    ax_spec.set_facecolor('darkgrey')
    
    display(fig) # plots the figure as a static png
        
def zoom_plot():
    """ I would love to figure out how to zoom into a wcs projection plot"""
    zoom_pct = 0.1
    # Center of the image
    w2, h2 = w / 2, h / 2


    # Convert pixel coordinates to WCS coordinates
    ra_center, dec_center = ob._wcs_f.all_pix2world(w2, h2, 0)

    # Assuming you want to zoom in on a region around the center
    # Here, we need to determine how much RA and Dec should change
    # for a 10% zoom. This is highly specific to your data and WCS.
    # For example, you might need to calculate this based on the field of view, pixel scale, etc.

    # For demonstration, let's assume a small change in RA and Dec
    delta_ra, delta_dec = 0.01, 0.01  # These values should be adjusted based on your data

    ra_min = 500.0
    ra_max = -1.0
    dec_min = 500.0
    dec_max = -1.0
    for sl in sightlines:
        print(f"{sl}")
        ras, decs = sl.radecs
        if ras[0] < ra_min: ra_min = ras[0]
        if ras[1] > ra_max: ra_max = ras[1]
        if decs[0] < dec_min: dec_min = decs[0]
        if decs[1] > dec_max: dec_max = decs[1]
    # Set new limits in WCS coordinates
    print(f"ra_min, ra_max = {ra_min}, {ra_max}")   
    print(f"dec_min, dec_max = {dec_min}, {dec_max}")
    ax.set_xlim(ra_min, ra_max)
    ax.set_ylim(dec_min, dec_max)
    # ax.set_ylim(dec_center - delta_dec, dec_center + delta_dec)

