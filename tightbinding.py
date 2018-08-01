import numpy as np, scipy as sp, itertools
from matplotlib import pyplot as plt

class TightBinding:
    def __init__(self, t_r, path = None, n_k = None):
        """
        t_r is a dict that maps translation-tuples to rank 2 arrays over the hopping's orbitals
        r and k normalized to 1
        """
        #self.dimension = len(path[0])
        self.dimension = len([r for r in t_r.keys()][0])
        if self.dimension == 1:
            local = (0)
        elif self.dimension == 2:
            local = (0,0)
        elif self.dimension == 3:
            local = (0,0,0)
        self.n_bands = len(t_r[local])
        self.t_r = t_r
        if path is not None:
            self.path = [np.array(p) for p in path]
        self.n_k = n_k
        self.n_kpoints_segment = None
        self.k_mesh_path = None
        self.epsilon_k_path = None
        self.k_mesh = None
        self.dos = None

    def epsilon_at(self, k):
        return np.sum([np.exp(complex(0,2*np.pi*np.array(r).dot(k)))*np.array(t) for r, t in self.t_r.items()], axis=0, dtype=complex).real
        
    def calculate_dispersion(self):
        self.k_mesh_path, self.n_kpoints_segment = self.calculate_k_mesh_path()
        self.epsilon_k_path = []
        for k in self.k_mesh_path:
            t_k = self.epsilon_at(k)
            self.epsilon_k_path.append(np.linalg.eigvalsh(t_k))
        self.epsilon_k_path = np.array(self.epsilon_k_path)

    def calculate_k_mesh_path(self):
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
            ax.plot(range(len(self.k_mesh_path)), self.epsilon_k_path[:,bandnr], **kwargs)
            if 'label' in kwargs.keys():
                del kwargs['label']
        ax.set_xticks(self.get_data_positions_of_path_points())
        ax.set_xticklabels(path_labels)
        ax.set_ylabel("$\epsilon (k)$")
        ax.set_xlabel("$k$")
        ax.set_xlim(0,len(self.k_mesh_path+1))

    def plot_dos(self, ax, rotated = True, **kwargs):
        if rotated:
            ax.plot(self.dos[1], self.dos[0], **kwargs)
            ax.set_xlabel("$N_\sigma(\omega)$")
        else:
            ax.plot(self.dos[0], self.dos[1], **kwargs)
            ax.set_ylabel("$N_\sigma(\omega)$")
            ax.set_xlabel("$\omega$")

    def calculate_dos(self, n_omega, n_k_per_dim, eps = 2**(-4)):
        """
        The continuous dispersion is chunked into discrete energy levels that are
        approximated with lorentzians of widths defined by parameter eps
        """
        self.dos = np.zeros([2, n_omega])
        self.k_mesh = self.calculate_k_mesh(n_k_per_dim)
        self.epsilon_k = np.empty([n_k_per_dim**self.dimension, self.n_bands])
        for i, k in enumerate(self.k_mesh):
            t_k = self.epsilon_at(k)
            self.epsilon_k[i,:] = np.linalg.eigvalsh(t_k)
        self.epsilon_k_min = self.epsilon_k.min()
        self.epsilon_k_max = self.epsilon_k.max()
        self.bandwidth = round(self.epsilon_k_max - self.epsilon_k_min)
        delta_omega = abs(self.epsilon_k_max - self.epsilon_k_min) / float(n_omega -1)
        eps_norm = eps/np.pi/len(self.k_mesh)
        for w in range(n_omega):
            omega = self.epsilon_k_min + w * delta_omega
            self.dos[0, w] = omega
            for i, n in itertools.product(range(len(self.k_mesh)),range(self.n_bands)):
                self.dos[1, w] += eps_norm / ((omega - self.epsilon_k[i, n])**2 + eps**2)

    def calculate_k_mesh(self, n_points):
        line = [-.5+i/float(n_points) for i in range(n_points)]
        k_pts = np.empty([n_points**self.dimension, self.dimension])
        i = 0
        for k in itertools.product(*[line]*self.dimension):
            k_pts[i,:] = k
            i += 1
        return k_pts

    def get_N_T0(self, eps_F = 0):
        """return lower and upper estimate (riemannian integral appr.)"""
        #not written for mu above largest epsilon_k
        eps_fermi_arg = np.argmin(abs(self.dos[0,:]-eps_F))
        integral = np.zeros([eps_fermi_arg+1])
        for i, d in enumerate(self.dos[1,:eps_fermi_arg+1]):
            integral[i] = d * abs(self.dos[0, i+1] - self.dos[0, i])
        nl = np.sum(integral[:-1], axis = 0)
        nu = nl + integral[-1]
        return nl, nu
    
class TightBindingK(TightBinding):
    """
    t_k must be a callable function that returns a rank 2 array
    dimension is that of the lattice
    Don't use factor of 2*pi in t_k!
    """
    def __init__(self, t_k, dimension, path = None, n_k = None):
        self.t_k = t_k
        self.dimension = dimension
        self.n_bands = t_k(np.array([0]*dimension)).shape[0]
        if path is not None:
            self.path = [np.array(p) for p in path]
        self.n_k = n_k
        self.n_kpoints_segment = None
        self.k_mesh_path = None
        self.epsilon_k_path = None
        self.k_mesh = None
        self.dos = None

    def epsilon_at(self, k):
        return self.t_k(2*np.pi*k)
