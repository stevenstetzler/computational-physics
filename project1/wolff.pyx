import numpy as np
import matplotlib.pyplot as plt

class Lattice():
    def __init__(self, size_x, size_y):
        self._size_x = size_x
        self._size_y = size_y
        self._lattice = np.zeros((self._size_x, self._size_y))
        for i in range(self._size_x):
            for j in range(self._size_y):
                spin = 2 * np.random.uniform() - 1
                if spin > 0:
                    spin = 1
                else:
                    spin = -1
                self._lattice[i, j] = spin
        
    def get(self, i, *args):
        if len(args) == 0:
            i, j = i
        else:
            i, j = i, args[0]
            
        new_i = i % self._size_x
        new_j = j % self._size_y
        return self._lattice[new_i, new_j]

    def flip(self, i, *args):
        if len(args) == 0:
            i, j = i
        else:
            i, j = i, args[0]
            
        new_i = i % self._size_x
        new_j = j % self._size_y
        self._lattice[new_i, new_j] = -self._lattice[new_i, new_j]
    
    def print(self):
        print(self._lattice)
        
    def get_size(self):
        return self._size_x, self._size_y
    
    def get_neighbors(self, i, *args):
        if len(args) == 0:
            i, j = i
        else:
            i, j = i, args[0]
        
        size_x, size_y = self.get_size()
        
        left_x, left_y = (i - 1) % size_x, j % size_y
        right_x, right_y = (i + 1) % size_x, j % size_y
        top_x, top_y = i % size_x, (j - 1) % size_y
        bottom_x, bottom_y = i % size_x, (j + 1) % size_y
        
        left_result = self.get(left_x, left_y)
        right_result = self.get(right_x, right_y)
        top_result = self.get(top_x, top_y)
        bottom_result = self.get(bottom_x, bottom_y)
        
        return (left_result, right_result, top_result, bottom_result), ((left_x, left_y), (right_x, right_y), (top_x, top_y), (bottom_x, bottom_y))

    def get_lattice(self):
        return self._lattice

def compute_H(lattice, J):
    H = 0
    size_x, size_y = lattice.get_size()
    for i in range(size_x):
        for j in range(size_y):
            sigma_ij = lattice.get(i, j)
            neighbor_spins, neighbor_locations = lattice.get_neighbors(i, j)
            for neighbor in neighbor_spins:
                H += -J * sigma_ij * neighbor
    return H

def metropolis(lattice, J, beta, num_trials):
    size_x, size_y = lattice.get_size()
    for n in range(num_trials):
        i, j = int(np.random.uniform(0, size_x)), int(np.random.uniform(0, size_y))

        H_mol = 0
        neighbor_spins, neighbor_locations = lattice.get_neighbors(i, j)
        for neighbor in neighbor_spins:
            H_mol += J * neighbor

        delta_E = 2 * lattice.get(i, j) * H_mol

        if np.random.uniform() < np.exp(-beta * delta_E):
            lattice.flip(i, j)

def metropolis_sim(N, J, Ts, num_mc_trials):
    energies = []
    lattice = Lattice(N, N)
    for T in Ts:
        beta = 1 / T
        metropolis(lattice, J, beta, num_mc_trials)
        H = compute_H(lattice, J)
        energies.append(H / (2 * N**2))
    return energies

def wolff(lattice, J, beta, num_mc_trials):
    size_x, size_y = lattice.get_size()
    for n in range(num_mc_trials):
        i, j = int(np.random.uniform(0, size_x)), int(np.random.uniform(0, size_y))

        C = []
        F_old = [(i, j)]

        while len(F_old) != 0:
            F_new = []
            for site in F_old:
                site_spin = lattice.get(site)
                neighbor_spins, neighbor_locations = lattice.get_neighbors(site)
                
                same_neighbors = []
                for spin, location in zip(neighbor_spins, neighbor_locations):
                    x, y = location
                    if spin == site_spin and location not in C:
                        same_neighbors.append((spin, location))
                        
        
                for spin, location in same_neighbors:
                    if np.random.uniform() < 1 - np.exp(-2 * beta * J):
                        F_new.append(location)
                        x, y = location
                        C.append(location)
            F_old = F_new

        for site in C:
            lattice.flip(site)

def wolff_fast(lattice, J, beta, num_mc_trials):
    size_x, size_y = lattice.get_size()
    for n in range(num_mc_trials):
        i, j = int(np.random.uniform(0, size_x)), int(np.random.uniform(0, size_y))
        
        # All the points in the cluster of spins to flip
        cluster = np.zeros((size_x, size_y))
        cluster[i, j] = 1

        
        # All of the sites who's neighbors to consider for flipping at each round
        F_old = np.zeros((size_x, size_y))
        F_old[i, j] = 1
        
        num_sites_to_add = -1

        spin = lattice.get(i, j)
        like_spin = lattice.get_lattice() == spin
        
        while num_sites_to_add != 0:           
            # Shift around the old edges of the cluster to find new ones
            left = np.roll(F_old, -1, axis=1)
            right = np.roll(F_old, 1, axis=1)
            top = np.roll(F_old, -1, axis=0)
            bottom = np.roll(F_old, 1, axis=0)

            edges = np.zeros((size_x, size_y))
            edges = (left + right + top + bottom)
            # Don't include new edges that are inside the cluster
            edges[np.where(cluster == 1)] = 0
            # Identify all of the edges of the cluster
            on_edge = edges >= 1
            
            # Identfy where all of the new edges are and which have spins that match the cluster's
            neighbor_sites = np.where(like_spin & on_edge)
            
            locations_to_add = np.where(np.random.uniform(size=neighbor_sites[0].size) < 1 - np.exp(-2 * beta * J))
            
            sites_to_add = (neighbor_sites[0][locations_to_add], neighbor_sites[1][locations_to_add])
            num_sites_to_add = sites_to_add[0].shape[0]
            # Add the sites to the cluster
            cluster[sites_to_add] = 1
            # Add the added edges to the list of considered edges for next time
            F_new = np.zeros((size_x, size_y))
            F_new[sites_to_add] = 1
            
            F_old = F_new
            
        # Flip the spins that are inside the cluster
        lattice.flip(np.where(cluster == 1))

def wolff_sim(N, J, Ts, num_mc_trials, fast=True):
    energies = []
    for T in Ts:
        lattice = Lattice(N, N)
        
        beta = 1 / T
        if fast:
            wolff_fast(lattice, J, beta, num_mc_trials)
        else:
            wolff(lattice, J, beta, num_mc_trials)
        H = compute_H(lattice, J)
        energies.append(H / (2*N**2))
    return energies

# metropolis simulation
N = 30
J = 1
Ts = np.linspace(0.01, 5, 200)
num_mc_trials = 500000

energies = metropolis_sim(N, J, Ts, num_mc_trials)

out_file = open("energies_metropolis.out", "w")
for t, e in zip(Ts, energies):
    out_file.write("{} {}\n".format(t, e))
out_file.close()

# plt.plot(Ts, energies, linestyle="None", marker="^", color="green")

# Wolff simulation
N = 30
J = 1
Ts = np.linspace(0.01, 5, 200)
num_mc_trials = 20000

energies = wolff_sim(N, J, Ts, num_mc_trials)

out_file = open("energies_wolff.out", "w")
for t, e in zip(Ts, energies):
    out_file.write("{} {}\n".format(t, e))
out_file.close()

# plt.plot(Ts, energies, linestyle="None", marker="^", color="green")

