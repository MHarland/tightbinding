import numpy as np
from itertools import product

class HoppingGenerator:
    def __init__(self, translationvectors, basisvectors):
        """translationvectors cartesian, basisvectors in direct coordinates"""
        self.translationvectors = np.array(translationvectors)
        self.basisvectors = np.array(basisvectors)
        self.basis = range(len(basisvectors))
        self.t_cartesian = {}
        self.t = {}
        self.n_basis_sites = len(self.basis)

    def generate(self, hoppings):
        """hoppings: list, [mu,tnn,tnnn,...]"""
        distances = {}
        for r in self.translations():
            distances[tuple(r)] = np.zeros([self.n_basis_sites]*2)
            self.t_cartesian[tuple(r)] = np.zeros([self.n_basis_sites]*2)
            for i, j in product(*[self.basis]*2):
                b1 = self.basisvectors[i]
                b2 = self.basisvectors[j]
                distances[tuple(r)][j,i] = np.linalg.norm(r+b2.dot(self.translationvectors) - b1.dot(self.translationvectors))
        hoppings_distance_sorted = self.classify_distances(distances)
        for k,t in enumerate(hoppings):
            for r,d in distances.items():
                for i,j in product(*[range(len(d))]*2):
                    if np.allclose(hoppings_distance_sorted[k],d[i,j]):
                        self.t_cartesian[r][i,j] = t
        for r, t in self.t_cartesian.items():
            self.t[tuple(self.transform_to_displacement(r))] = t
                
    def get_hopping(self, *args, **kwargs):
        return self.t

    def classify_distances(self, distance_dict):
        distances = []
        for r, d in distance_dict.items():
            for i,j in product(*[range(len(d))]*2):
                distances.append(d[i,j])
        distances_unique = self.make_unique(np.array(distances))
        return np.sort(distances_unique)

    def translations(self):
        all_translations = []
        for r in self.translationvectors:
            for sign in [+1,-1]:
                r = sign*r
                all_translations.append(r)
                if len(self.translationvectors) > 1:
                    for r2 in self.translationvectors:
                        for sign2 in [+1,-1]:
                            r2 = sign2*r2
                            all_translations.append(r+r2)
                            if len(self.translationvectors) > 2:
                                for r3 in self.translationvectors:
                                    for sign3 in [+1,-1]:
                                        r3 = sign3*r3
                                        all_translations.append(r+r2+r3)
        return self.make_unique(all_translations)

    def make_unique(self, listofv):
        unique = []
        for v in listofv:
            is_unique = True
            for v2 in unique:
                if np.allclose(v,v2):
                    is_unique = False
            if is_unique:
                unique.append(v)
        return np.array(unique)

    def transform_to_displacement(self, v):
        return self.reciprocal_latticevectors().dot(v)

    def reciprocal_latticevectors(self):
        rlvs = []
        lvs = self.translationvectors
        if len(lvs) == 3:
            volume = lvs[0].dot(np.cross(lvs[1],lvs[2]))
            for i in range(3):
                rlvs.append(np.cross(lvs[(i+1)%3],lvs[(i+2)%3])/volume)
        return np.array(rlvs)
