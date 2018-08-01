import numpy as np
from tightbinding.tightbinding import TightBinding
try:
    from mpl_to_latex import PRL
    mplcfg = PRL(add_ax = False)
    mpltolat = True
    plt = mplcfg.plt
    fig = mplcfg.fig
except ImportError:
    from matplotlib import pyplot as plt
    fig = plt.gcf()
    mpltolat = False


ax = fig.add_axes([.14,.15,.6,.83])
ax2 = fig.add_axes([.74,.15,.25,.83])
col = plt.cm.viridis(0.)

t, s,  = -1., .15
mu = 4*s
t_r = {(0,0)  :[[-mu]],
       (1,0)  :[[t]],
       (1,1)  :[[s]],
       (0,1)  :[[t]],
       (-1,1) :[[s]],
       (-1,0) :[[t]],
       (-1,-1):[[s]],
       (0,-1) :[[t]],
       (1,-1) :[[s]]}
symmetrypath = [[0,0],[.5,0],[.5,.5],[0,0]]
pathlabels = ["$\Gamma$","$X$","$M$","$\Gamma$"]
latticedisp = TightBinding(t_r,symmetrypath,128)
latticedisp.calculate_dispersion()
latticedisp.plot_dispersion(ax, pathlabels, color = col)
latticedisp.calculate_dos(200,32,2**(-2))
latticedisp.plot_dos(ax2, color = col)

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xticks([])
ax2.set_ylim(ax.get_ylim())
ax.set_xlabel("$k$")
ax.set_ylabel("$\\epsilon(k)$")
if mpltolat:
    mplcfg.save("square_dos")
else:
    plt.tight_layout()
    plt.savefig('square_dos.pdf')
