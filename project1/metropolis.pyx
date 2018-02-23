import numpy as np

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
