from tightbinding import TightBinding as TBD
from matplotlib import pyplot as plt, cm, patches
from hoppinggenerator import HoppingGenerator
import numpy as np


# square lattice
sl = HoppingGenerator([[1,0],[0,1]],[[0,0],[.5,0],[0,.5],[.5,.5]])
symmetrypath = [[0,0],[.5,0],[.5,.5],[0,0]]
pathlabels = ["$\Gamma$","$X$","$M$","$\Gamma$"]

"""
# pyrochlore lattice
sl = HoppingGenerator([[1, 0, 0],[.5, .5 * np.sqrt(3), 0],[.5, .5 * 1/np.sqrt(3), np.sqrt(2. / 3)]],[[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5]])
#symmetrypath = [[0,0,0],[0.0,.5,.5],[.25,.75,.5],[3./8,.75,3./8],[0,0,0],[.5,.5,.5]]
#pathlabels = ["$\Gamma$","$X$","$W$","$K$","$\Gamma$","$L$"]
symmetrypath = [[0,0,0],[0.0,.5,.5],[.25,.75,.5],[.5,.5,.5],[0,0,0],[3./8,.75,3./8]]#,[0.0,.5,.5]]
pathlabels = ["$\Gamma$","$X$","$W$","$L$","$\Gamma$","$K$"]#,"$X$"]
"""

# parameters
hw = 21*.5
mus = [-1.1]
t1s = [-1]*len(mus)
t2s = [-.4]*len(mus) # lower lower band at G
t3s = [-.35]*len(mus) # lower upper band at G
t4s =[-.25]*len(mus) # bend X-W
t5s=[0]*len(mus) # breaks X-w
t6s=[0]*len(mus)
t7s=[0]*len(mus)
hw=1
mus = [0]*4
t1s = [-1]*len(mus)
t2s = [0,-.5,-1,-1.5]#*len(mus) # lower lower band at G
t3s = [0]*len(mus) # lower upper band at G
t4s =[0]*len(mus) # bend X-W
t5s=[0]*len(mus) # breaks X-w
t6s=[0]*len(mus)
t7s=[0]*len(mus)

colors = [cm.jet(i/float(max(1,len(mus)-1))) for i in range(len(mus))]
linestyles = ["-","--",":","-."]
#linestyles = ["-"]*len(mus)

fig = plt.figure()
ax = fig.add_axes([.1,.1,.6,.85])
ax2 = fig.add_axes([.7,.1,.25,.85])
ax2.set_yticklabels([])
for i,mu,t1, t2, t3, t4,t5,t6,t7,col, ls in zip(range(len(mus)),mus,t1s,t2s,t3s,t4s,t5s,t6s,t7s,colors,linestyles):
    sl.generate([mu/hw,t1/hw,t2/hw,t3/hw,t4/hw,t5/hw,t6/hw,t7/hw])
    latticedisp = TBD(sl.get_hopping(),symmetrypath,200)
    latticedisp.calculate_dispersion()
    latticedisp.plot_dispersion(ax, pathlabels, color = col, ls = ls)
    latticedisp.calculate_dos(200,64,2**(-3))
    #latticedisp.calculate_dos(300,16,2**(-3))
    latticedisp.plot_dos(ax2, color = col, ls = ls)
    nl,nu = latticedisp.get_N_T0()
    ax2.text(.98,.03*i,str(nl)[:4]+"$<N_\sigma(T=0)<$"+str(nu)[:4], transform=ax2.transAxes, horizontalalignment="right", verticalalignment="bottom", color = col)
    #ax2.set_xticklabels([])
    #ax2.set_xticks([])
    """
    nf = 10
    ef = None
    for i in np.arange(2,-3,-.01):
        n = abs(latticedisp.get_N_T0(i)[1]-2)
        if n < nf:
            nf = n
            ef = i
    print ef
    print latticedisp.bandwidth
    """
#ax.plot([0,10**6],[0]*2,alpha=.5,color="gray",ls="--")
#ax2.set_xlim(0,3)
#ax.legend(handles = [patches.Patch(color=cm.jet(i), label=lab) for i,lab in zip([0,.5,.999],["-0.1","0","+0.1"])], bbox_to_anchor=[1.02,1],loc="upper left", title="$t_6$")
plt.savefig("dispersion.pdf", dpi=300)
print "dispersion.pdf written."



