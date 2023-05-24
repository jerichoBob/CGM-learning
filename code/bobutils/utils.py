# print("<hi>")

# This just draws a single box centered at (x,y) and sz from that center point in the n/e/s/w directions 
def plotbox(plt, x, y, sz, c):
    ax = x - 0.5;  # I don't know why this is yet.... :()
    ay = y - 0.5;
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

        plt.text(x[i]-1.5*sz, y[i], labels[i], color=colour);

def acc_circle(data, x_cen, y_cen, z_index, rad):
    # assuming rad = 1 for now (KISS)
    loc_flux = [];
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]]); # center
    loc_flux.append(data[i][y_cen[z_index]+1][x_cen[z_index]]); # north
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]+1]); # east
    loc_flux.append(data[i][y_cen[z_index]-1][x_cen[z_index]]); # south
    loc_flux.append(data[i][y_cen[z_index]][x_cen[z_index]-1]); # west
    return loc_flux