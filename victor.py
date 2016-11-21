from tightbinding import TightBinding
from hoppinggenerator import HoppingGenerator


hgen = HoppingGenerator([[1,0],[1.5,0]], [[0,0]])
hgen.generate([0, -1,-.5,-.25]) # hopping energies sorted by distance
tb = TightBinding(hgen.get_hopping())
mesh = tb.calculate_k_mesh(20)
for k in mesh:
    print k, tb.epsilon_at(k)
