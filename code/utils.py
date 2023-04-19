print("<hi>")
def test():
    print("hi")

# This just draws a box centered at (x,y) and sz from that center point in the n/e/s/w directions 
def plotbox(plt, x, y, labels, align, sz, c):
    for i in range(len(x)):
        plt.plot(
            [x[i]-sz, x[i]-sz, x[i]+sz, x[i]+sz, x[i]-sz], 
            [y[i]-sz, y[i]+sz, y[i]+sz, y[i]-sz, y[i]-sz], 
            '-', color=c)
        ha_ = align[i][0]
        va_ = align[i][1]
        plt.text(x[i], y[i]+1.5*sz, labels[i], color=c);

