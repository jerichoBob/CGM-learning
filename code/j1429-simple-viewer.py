import layout_utils as lu

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from matplotlib.widgets import SpanSelector
from linetools.spectra.xspectrum1d import XSpectrum1D  

# get some data
base_path = "/Users/robertseaton/Desktop/Physics-NCState/---Research/FITS-data/J1429"
filename = base_path+"/spectrum_30.37.fits"

sp=XSpectrum1D.from_file(filename)
wave=sp.wavelength.value
flux=sp.flux.value
error=0.1*flux

#..................................................
# Define the range of interest
continuum_range = (4700, 4800)

# Select the data within this range
mask = (wave >= continuum_range[0]) & (wave <= continuum_range[1])
selected_flux = flux[mask]

# Compute the mean flux and its standard deviation within this range
mean_flux = np.mean(selected_flux, axis=0)
stddev_flux = np.std(selected_flux, axis=0)

# Compute the SNR
snr = mean_flux / stddev_flux
#..................................................
print("mean flux between 4700-4800: ", mean_flux)
print("stddev flux between 4700-4800: ", stddev_flux)
print("snr between 4700-4800: ", snr)


# let's draw some stuff
fig = plt.figure(figsize=(14,8))
gs = GridSpec(10, 5, figure=fig)
clear_axes = lu.clear_axes

# left panel
ax = []
for i in range(5):
    ax.append(clear_axes(fig.add_subplot(gs[2*i:2*i+2, 0])))
    ax.append(clear_axes(fig.add_subplot(gs[2*i:2*i+2, -1])))

# right panel
ax1 = fig.add_subplot(gs[0:5, 1:-1])
ax2 = fig.add_subplot(gs[5:10, 1:-1])


flux_color = '#0055ff99'
error_color = '#ff330099'
mean_color = '#00ff3399'

ax1.plot(wave, flux, '-', color=flux_color)
ax1.plot(wave, error, '-', color=error_color)
line2, = ax2.plot([], [])

def onselect(wave_min, wave_max):
    indmin, indmax = np.searchsorted(wave, (wave_min, wave_max))
    indmax = min(len(wave) - 1, indmax)

    wave_range = wave[indmin:indmax]
    flux_range = flux[indmin:indmax]
    error_range = error[indmin:indmax]

    if len(wave_range) < 2:
        ax[1].clear()
        ax[1].cla()
        print("yep")

    else:
        mean_flux = np.mean(flux_range, axis=0)
        snr = mean_flux / stddev_flux
        wmin = wave_range[0]
        wmax = wave_range[-1]

        ax[1].clear()
        clear_axes(ax[1])
        ax[1].text(0, .9, f"Selected Range ({wmin:.2f}, {wmax:.2f})", va="top", ha="left")
        ax[1].text(0, .7, f"Mean Flux: {mean_flux:.4f}", va="top", ha="left")
        ax[1].text(0, .5, f"S/N: {snr:.4f}", va="top", ha="left")

        ax2.clear()
        ax2.plot(wave_range, flux_range, '-', color=flux_color)
        ax2.plot(wave_range, error_range, '-', color=error_color)
        ax2.axhline(y=mean_flux,linewidth=1, color=mean_color)

        ax2.set_xlim(wave_range[0], wave_range[-1])
        ax2.set_ylim(min(flux_range.min(),error_range.min()),
                     max(flux_range.max(), error_range.max()))
        fig.canvas.draw_idle()
    


span = SpanSelector(
    ax1,
    onselect,
    "horizontal",
    useblit=True,
    props=dict(alpha=0.5, facecolor="tab:blue"),
    interactive=True,
    drag_from_anywhere=True
)

fig.tight_layout()

plt.subplots_adjust(left=0.02, bottom=0.075, right=0.97, top=0.95, wspace=0.258, hspace=0.614)
plt.show()
