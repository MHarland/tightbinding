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
    t = -1.
    mu = 0
    tnnn = s = .3
    u,v,w = tz, tz/4., -tz/2.
    t_r = {(0,0,0):[[-mu,t,t,s],[t,-mu,s,t],[t,s,-mu,t],[s,t,t,-mu]],
           (1,0,0):[[0,t,0,s],[0,0,0,0],[0,s,0,t],[0,0,0,0]],
           (1,1,0):[[0,0,0,s],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           (0,1,0):[[0,0,t,s],[0,0,s,t],[0,0,0,0],[0,0,0,0]],
           (-1,1,0):[[0,0,0,0],[0,0,s,0],[0,0,0,0],[0,0,0,0]],
           (-1,0,0):[[0,0,0,0],[t,0,s,0,],[0,0,0,0],[s,0,t,0]],
           (-1,-1,0):[[0,0,0,0],[0,0,0,0],[0,0,0,0],[s,0,0,0]],
           (0,-1,0):[[0,0,0,0],[0,0,0,0],[t,s,0,0],[s,t,0,0]],
           (1,-1,0):[[0,0,0,0],[0,0,0,0],[0,s,0,0],[0,0,0,0]],

           (0,0,1):[[u,0,0,w],[0,u,w,0],[0,w,u,0],[w,0,0,u]],
           (0,0,-1):[[u,0,0,w],[0,u,w,0],[0,w,u,0],[w,0,0,u]],

           (1,1,1):[[0,0,0,w],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           (1,1,-1):[[0,0,0,w],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           (-1,1,1):[[0,0,0,0],[0,0,w,0],[0,0,0,0],[0,0,0,0]],
           (-1,1,-1):[[0,0,0,0],[0,0,w,0],[0,0,0,0],[0,0,0,0]],
           (-1,-1,1):[[0,0,0,0],[0,0,0,0],[0,0,0,0],[w,0,0,0]],
           (-1,-1,-1):[[0,0,0,0],[0,0,0,0],[0,0,0,0],[w,0,0,0]],
           (1,-1,1):[[0,0,0,0],[0,0,0,0],[0,w,0,0],[0,0,0,0]],
           (1,-1,-1):[[0,0,0,0],[0,0,0,0],[0,w,0,0],[0,0,0,0]],

           (1,0,-1):[[v,0,0,w],[0,v,0,0],[0,w,v,0],[0,0,0,v]],
           (1,0,1):[[v,0,0,w],[0,v,0,0],[0,w,v,0],[0,0,0,v]],
           (0,1,-1):[[v,0,0,w],[0,v,w,0],[0,0,v,0],[0,0,0,v]],
           (0,1,1):[[v,0,0,w],[0,v,w,0],[0,0,v,0],[0,0,0,v]],
           (-1,0,-1):[[v,0,0,0],[0,v,w,0],[0,0,v,0],[w,0,0,v]],
           (-1,0,1):[[v,0,0,0],[0,v,w,0],[0,0,v,0],[w,0,0,v]],
           (0,-1,-1):[[v,0,0,0],[0,v,0,0],[0,w,v,0],[w,0,0,v]],
           (0,-1,1):[[v,0,0,0],[0,v,0,0],[0,w,v,0],[w,0,0,v]]}
    symmetrypath = [[0,0,0],[.5,0,0],[.5,.5,0],[0,0,0],[.5,.5,.5],[.5,0,0]]
    pathlabels = ["$\Gamma$","$X$","$M$","$\Gamma$","$R$","$X$"]
    latticedisp = TightBinding(t_r,symmetrypath,nk)
    latticedisp.calculate_dispersion()
    latticedisp.plot_dispersion(ax, pathlabels, color = color, label = '$'+str(tz)+'$')

ax.set_xlabel("$\\tilde{k}$")
ax.set_ylabel("$\\epsilon(\\tilde{k})$")
if mpltolat:
    mplcfg.legendkwargs['loc'] = 'upper right'
    leg = ax.legend(labelspacing = .1, title = '$t_\\perp$', **mplcfg.legendkwargs)
    mplcfg.save("cubic_andersen_2x2x1_tz")
else:
    leg = ax.legend(labelspacing = .1, title = '$t^\\prime$')
    plt.tight_layout()
    plt.savefig('cubic_andersen_2x2x1_tz.pdf')
