from hoppinggenerator import HoppingGenerator

hg = HoppingGenerator([[1,0],[0,1]],[[0,0],[.5,0],[0,.5],[.5,.5]])
hg.generate([0,-1])
#print hg.translations()
#print
print hg.get_hopping_dict()