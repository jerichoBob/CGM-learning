import matplotlib.pyplot as plt


fig = plt.figure(constrained_layout=True,figsize=(10,8))
gs0 = fig.add_gridspec(1, 2)

gsleft  = gs0[0].subgridspec(10, 10)
gsright = gs0[1].subgridspec(6, 1)

ax_text = fig.add_subplot(gsleft[0:1, :]) # text box
ax_text.set_xticks([]);
ax_text.set_xticklabels([]);
ax_text.set_yticks([]);
ax_text.set_yticklabels([]);
ax_imag = fig.add_subplot(gsleft[1:9, :]) # image box
ax_imag.set_xticks([]);
ax_imag.set_xticklabels([]);
ax_imag.set_yticks([]);
ax_imag.set_yticklabels([]);
ax_cntl = fig.add_subplot(gsleft[9:, :]) # control box
ax_cntl.set_xticks([]);
ax_cntl.set_xticklabels([]);
ax_cntl.set_yticks([]);
ax_cntl.set_yticklabels([]);

ax1 = fig.add_subplot(gsright[0, :])
ax1.set_xticks([]);
ax1.set_xticklabels([]);
ax2 = fig.add_subplot(gsright[1, :])
ax2.set_xticks([]);
ax2.set_xticklabels([]);
ax3 = fig.add_subplot(gsright[2, :])
ax3.set_xticks([]);
ax3.set_xticklabels([]);
ax4 = fig.add_subplot(gsright[3, :])
ax4.set_xticks([]);
ax4.set_xticklabels([]);
ax5 = fig.add_subplot(gsright[4, :])
ax5.set_xticks([]);
ax5.set_xticklabels([]);
ax6 = fig.add_subplot(gsright[5, :])

fig.tight_layout()

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95, wspace=0.129, hspace=0.17)
plt.show()
