# import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Makes a 1 row : 2 column even split (with initial size set to 10x8 (idk what units))
# left column is TextBox, Image, Slider Controls
# right column is 6 spectral lines (one for each of the 6 points on the image)
def make_6speclayout(plt):
    fig = plt.figure(constrained_layout=True,figsize=(10,8))
    gs0 = fig.add_gridspec(1, 2)

    # create an axes object that we can stuff full and send back to the caller
    ax = []
    gsleft  = gs0[0].subgridspec(10, 10)
    gsright = gs0[1].subgridspec(6, 1)

    ax_text = fig.add_subplot(gsleft[0:1, :]) # text box
    ax_text.set_xticks([])
    ax_text.set_xticklabels([])
    ax_text.set_yticks([])
    ax_text.set_yticklabels([])
    ax.append(ax_text)

    ax_imag = fig.add_subplot(gsleft[1:9, :]) # image box
    ax_imag.set_xticks([])
    ax_imag.set_xticklabels([])
    ax_imag.set_yticks([])
    ax_imag.set_yticklabels([])
    ax.append(ax_imag)

    ax_cntl = fig.add_subplot(gsleft[9:, :]) # control box
    ax_cntl.set_xticks([])
    ax_cntl.set_xticklabels([])
    ax_cntl.set_yticks([])
    ax_cntl.set_yticklabels([])
    ax.append(ax_cntl)

    ax1 = fig.add_subplot(gsright[0, :])
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax.append(ax1)
    
    ax2 = fig.add_subplot(gsright[1, :])
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    ax.append(ax2)

    ax3 = fig.add_subplot(gsright[2, :])
    ax3.set_xticks([])
    ax3.set_xticklabels([])
    ax.append(ax3)

    ax4 = fig.add_subplot(gsright[3, :])
    ax4.set_xticks([])
    ax4.set_xticklabels([])
    ax.append(ax4)

    ax5 = fig.add_subplot(gsright[4, :])
    ax5.set_xticks([])
    ax5.set_xticklabels([])
    ax.append(ax5)

    ax6 = fig.add_subplot(gsright[5, :])
    ax.append(ax6)

    fig.tight_layout()
    return ax



# Makes a 1 row : 2 column even split (with initial size set to 10x8 (idk what units))
# left column is TextBox, Image, Slider Controls
# right column is single spectral plot with 
def make_image_flux_var_layout(plt):
    # fig = plt.figure(constrained_layout=True,figsize=(14,8))
    fig = plt.figure(figsize=(14,8))
    gs0 = fig.add_gridspec(1, 2)

    # create an axes object that we can stuff full and send back to the caller
    ax = []
    gsleft  = gs0[0].subgridspec(10, 10)
    gsright = gs0[1].subgridspec(10, 10)

    ax_text = fig.add_subplot(gsleft[0:1, :]) # text box
    ax_text.set_xticks([])
    ax_text.set_xticklabels([])
    ax_text.set_yticks([])
    ax_text.set_yticklabels([])
    ax.append(ax_text)

    ax_imag = fig.add_subplot(gsleft[1:9, :]) # image box
    # ax_imag.set_xticks([])
    # ax_imag.set_xticklabels([])
    # ax_imag.set_yticks([])
    # ax_imag.set_yticklabels([])
    ax.append(ax_imag)

    ax_cntl = fig.add_subplot(gsleft[9:, :]) # control box
    ax_cntl.set_xticks([])
    ax_cntl.set_xticklabels([])
    ax_cntl.set_yticks([])
    ax_cntl.set_yticklabels([])
    ax.append(ax_cntl)


    ax5 = fig.add_subplot(gsright[4:6, :])
    ax5.set_xticks([])
    ax5.set_xticklabels([])
    ax.append(ax5)

    fig.tight_layout()
    return ax

def make_simple_image_flux_var_layout(plt, image_size=6):
    fig = plt.figure(figsize=(14,8))
    gs0 = fig.add_gridspec(1, 2)

    # create an axes object that we can stuff full and send back to the caller
    ax = []
    gsleft  = gs0[0].subgridspec(10, 10)
    gsright = gs0[1].subgridspec(10, 10)

    ax_text = fig.add_subplot(gsleft[0:1, :]) # text box
 
    clear_axes(ax_text)
    ax.append(ax_text)

    ax_imag = fig.add_subplot(gsleft[1:image_size, :]) # image box
    ax.append(ax_imag)

    ax_spec = fig.add_subplot(gsright[4:6, :])
    ax_spec.set_xticks([])
    ax_spec.set_xticklabels([])
    ax.append(ax_spec)

    fig.tight_layout()
    return ax

def make_image_and_12_plots(plt):
    gs = gridspec.GridSpec(6, 3)
    # Create a subplot for the image on the left
    ax_text = plt.subplot(gs[0:1, :1]) # text box
    ax_text = clear_axes(ax_text)
    ax_image = plt.subplot(gs[1:5, :1])

    # Create 12 subplots on the right for the spectra
    ax_spectra = [plt.subplot(gs[i // 2, i % 2 + 1]) for i in range(12)]

    plt.tight_layout()
    return ax_text, ax_image, ax_spectra

def clear_axes(ax):
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])    
    ax.set_axis_off()
    ax.tick_params(labelbottom=False, labelleft=False)
    return ax

def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)