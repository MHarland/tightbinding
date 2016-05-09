from hoppinggenerator import HoppingGenerator
import numpy as np
from itertools import product

sl = HoppingGenerator([[1, 0, 0],[.5, .5 * np.sqrt(3), 0],[.5, .5 * 1/np.sqrt(3), np.sqrt(2. / 3)]],[[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5]])
sl.generate([0,2,3,4,5])
for r, h in sl.get_hopping().items():
    is_zero = True
    for i, j in product(*[range(len(h))]*2):
        if not np.allclose(h[i,j],0):
            is_zero = False
            break
    if not is_zero:
        r = list(r)
        for k in range(len(r)):
            r[k] = int(round(r[k]))
        print tuple(r),":"
        print h,","
