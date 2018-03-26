from tight_binding.tightbinding import TightBinding as TBD
from matplotlib import pyplot as plt, cm, patches
import numpy as np


symmetrypath = [[0,0],[.5,0],[2./3.,-1./3.],[0,0]]
pathlabels = ["$\Gamma$","$M$","$K$","$\Gamma$"]
t = -1
s = -.2
m = -0
hop = {(1, 1):
       [[0,  t,  0,  0],
        [0,  0,  0,  0],
        [0,  0,  0,  t],
        [0,  0,  0,  0]],
       (0, 0):
       [[-m,  t,  s,  0],
        [t,  -m,  0,  s],
        [s,  0,  -m,  t],
        [0,  s,  t,  -m]],
       (1, 0):
       [[0,  t,  0,  0],
        [0,  0,  0,  0],
        [0,  0,  0,  t],
        [0,  0,  0,  0]],
       (-1, -1):
       [[0,  0,  0,  0],
        [t,  0,  0,  0],
        [0,  0,  0,  0],
        [0,  0,  t,  0]],
       (-1, 0):
       [[0,  0,  0,  0],
        [t,  0,  0,  0],
        [0,  0,  0,  0],
        [0,  0,  t,  0]]}

fig = plt.figure()
ax = fig.add_axes([.1,.1,.6,.85])
ax2 = fig.add_axes([.7,.1,.25,.85])
ax2.set_yticklabels([])
latticedisp = TBD(hop, symmetrypath, 500)
latticedisp.calculate_dispersion()
latticedisp.plot_dispersion(ax, pathlabels)
latticedisp.calculate_dos(100,32,2**(-4))
latticedisp.plot_dos(ax2)
ax2.set_xticklabels([])
ax2.set_xticks([])
plt.savefig("dispersion.pdf", dpi=300)
print "dispersion.pdf written."



