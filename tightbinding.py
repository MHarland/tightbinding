import numpy as np, scipy as sp
from matplotlib import pyplot as plt

class TightBindingDispersion:
    def __init__(self, t_r, path, n_k):
        """
        r and k normalized to 1
        """
        self.dimension = len(path[0])
        if self.dimension == 1:
            local = (0)
        elif self.dimension == 2:
            local = (0,0)
        elif self.dimension == 3:
            local = (0,0,0)
        self.n_bands = len(t_r[local])
        self.t_r = t_r
        self.path = [np.array(p) for p in path]
        self.n_k = n_k
        self.n_kpoints_segment = None
        self.k_mesh = None
        self.epsilon_k = None
        self.t_k = None

    def calculate_dispersion(self):
        self.k_mesh, self.n_kpoints_segment = self.calculate_k_mesh()
        self.t_k = []
        self.epsilon_k = []
        for k in self.k_mesh:
            t_k = np.sum([np.exp(complex(0,2*np.pi*np.array(r).dot(k)))*np.array(t) for r, t in self.t_r.items()], axis=0, dtype=complex)
            self.t_k.append(t_k)
            self.epsilon_k.append(np.linalg.eigvalsh(t_k))
        self.epsilon_k = np.array(self.epsilon_k)
        self.t_k = np.array(self.t_k)

    def calculate_k_mesh(self):
        length_path = np.sum([np.linalg.norm(self.path[i+1]-self.path[i]) for i in range(len(self.path)-1)], axis=0)
        k_mesh = []
        n_points_per_segment = []
        for i,p in enumerate(self.path[:-1]):
            n_points_this_segment = int(np.ceil(self.n_k * np.linalg.norm(self.path[i+1]-self.path[i])/length_path))
            n_points_per_segment.append(n_points_this_segment)
            for j in range(n_points_this_segment):
                steplength = j/float(n_points_this_segment)
                k = self.path[i] + steplength * (self.path[i+1] - self.path[i])
                k_mesh.append(k)
        return np.array(k_mesh), np.array(n_points_per_segment)

    def get_data_positions_of_path_points(self):
        pos = [0]
        for interval in self.n_kpoints_segment:
            pos.append(pos[-1]+interval)
        return pos

    def plot_dispersion(self, ax, path_labels, **kwargs):
        assert len(path_labels) == len(self.path), "number of path labels wrong"
        for bandnr in range(self.n_bands):
            ax.plot(range(len(self.k_mesh)), self.epsilon_k[:,bandnr], **kwargs)
        ax.set_xticks(self.get_data_positions_of_path_points())
        ax.set_xticklabels(path_labels)
        ax.set_ylabel("$\epsilon (k)$")
        ax.set_xlabel("$k$")
        ax.set_xlim(0,len(self.k_mesh+1))

