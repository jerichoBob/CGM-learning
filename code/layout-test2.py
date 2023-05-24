import matplotlib.pyplot as plt

def clear_axes(ax):
    # ax.set_axis_off()
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    return ax

fig = plt.figure(figsize=(14,8))
gs0 = fig.add_gridspec(1, 2)

# create an axes object that we can stuff full and send back to the caller
gsleft  = gs0[0].subgridspec(20, 10)
gsright = gs0[1].subgridspec(10, 10)

ax_text1 = clear_axes(fig.add_subplot(gsleft[0:1, :])) # text box
ax_text2 = clear_axes(fig.add_subplot(gsleft[1:2, :])) # text box
ax_text3 = clear_axes(fig.add_subplot(gsleft[2:3, :])) # text box
ax_imag  = clear_axes(fig.add_subplot(gsleft[3:17, :])) # image box
ax_cntl1 = clear_axes(fig.add_subplot(gsleft[17:18:, :])) # control box
ax_cntl2 = clear_axes(fig.add_subplot(gsleft[18:19:, :])) # control box
ax_cntl3 = clear_axes(fig.add_subplot(gsleft[19:, :])) # control box

ax5 = fig.add_subplot(gsright[4:6, :])

fig.tight_layout()

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95, wspace=0.129, hspace=0.17)
plt.show()
