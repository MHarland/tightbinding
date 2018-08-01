import numpy as np
from tightbinding.tightbinding import TightBinding
try:
    from mpl_to_latex import PRL
    mplcfg = PRL()
    ax = mplcfg.ax
    mpltolat = True
    plt = mplcfg.plt
except ImportError:
    from matplotlib import pyplot as plt
    ax = plt.gca()
    mpltolat = False

    
nk = 128
tzs = [0, -.15]
ntbs = len(tzs)
colors = [plt.cm.viridis(float(i)/ntbs) for i in range(ntbs)]
for tz, color  in zip(tzs, colors):
    mu, t, s,  = 0., -1., 0
    t_r = {(0,0,0)  :[[-mu]],
           (1,0,0)  :[[t]],
           (1,1,0)  :[[s]],
           (0,1,0)  :[[t]],
           (-1,1,0) :[[s]],
           (-1,0,0) :[[t]],
           (-1,-1,0):[[s]],
           (0,-1,0) :[[t]],
           (1,-1,0) :[[s]],
           (0,0,1)  :[[tz]],
           (0,0,-1) :[[tz]]}
    symmetrypath = [[0,0,0],[.5,0,0],[.5,.5,0],[0,0,0],[.5,.5,.5],[.5,0,0]]
    pathlabels = ["$\Gamma$","$X$","$M$","$\Gamma$","$R$","$X$"]
    latticedisp = TightBinding(t_r,symmetrypath,nk)
    latticedisp.calculate_dispersion()
    latticedisp.plot_dispersion(ax, pathlabels, color = color, label = '$'+str(tz)+'$')

ax.set_xlabel("$k$")
ax.set_ylabel("$\\epsilon(k)$")
if mpltolat:
    leg = ax.legend(labelspacing = .1, title = '$\\mathrm{lattice}$, $t_\\perp$', **mplcfg.legendkwargs)
    mplcfg.save("cubic_tz")
else:
    leg = ax.legend(labelspacing = .1, title = '$t_\\perp$')
    plt.tight_layout()
    plt.savefig('cubic_tz.pdf')
