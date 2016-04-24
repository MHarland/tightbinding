from tightbinding import TightBindingDispersion as TBD
from matplotlib import pyplot as plt, cm
from hoppinggenerator import HoppingGenerator
import numpy as np

"""
# square lattice
sl = HoppingGenerator([[1,0],[0,1]],[[0,0],[.5,0],[0,.5],[.5,.5]])
symmetrypath = [[0,0],[.5,0],[.5,.5],[0,0]]
pathlabels = ["$\Gamma$","$X$","$M$","$\Gamma$"]
"""

# pyrochlore lattice
sl = HoppingGenerator([[1, 0, 0],[.5, .5 * np.sqrt(3), 0],[.5, .5 * 1/np.sqrt(3), np.sqrt(2. / 3)]],[[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5]])
symmetrypath = [[0,0,0],[0.0,.5,.5],[.25,.75,.5],[3./8,.75,3./8],[0,0,0],[.5,.5,.5]]
pathlabels = ["$\Gamma$","$X$","$W$","$K$","$\Gamma$","$L$"]

# parameters
tnns = [-1]*4
tnnns = [0,-.05,-.1,-.5]
colors = [cm.jet(i/float(max(1,len(tnnns)-1))) for i in range(len(tnnns))]
linestyles = ["-","--",":","-."]
#tnns = [-1]
#tnnns = [0]
#colors = [cm.jet(i/float(max(1,len(tnnns)-1))) for i in range(len(tnnns))]
#linestyles = ["-"]


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for t, tp, col, ls in zip(tnns,tnnns,colors,linestyles):
    sl.generate([0,t,0,tp])
    """
    for r, h in sl.get_hopping().items():
        print r
        print h
    print
    """
    latticedisp = TBD(sl.get_hopping(),symmetrypath,100)
    latticedisp.calculate_dispersion()
    latticedisp.plot_dispersion(ax, pathlabels, color = col, ls = ls)
plt.savefig("dispersion.pdf", dpi=300)
print "dispersion.pdf written."
