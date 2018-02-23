import numpy as np
# Ignore math warnings
import warnings
warnings.filterwarnings('ignore')
import math
import argparse

# The following lines allow us to adjust our plots, made with `matplotlib`.
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 18})
rcParams.update({'font.sans-serif':'Arial'})
rcParams.update({'font.family':'sans-serif'})
rcParams.update({'backend':'Agg'})
import matplotlib.pyplot as plt

def get_sigma_matrices_x(N):
    mats = []
    pauli_x = np.array([[0, 1], [1, 0]])
    
    for i in range(N):
        dim_left = 2 ** (i)
        dim_right = 2 ** (N - i - 1)

        identity_left = np.eye(dim_left)
        identity_right = np.eye(dim_right)

        mat = np.kron(identity_left, pauli_x)
        mat = np.kron(mat, identity_right)
        mats.append(mat)
    
    return mats

def get_sigma_matrices_z(N):
    mats = []
    pauli_z = np.array([[1, 0], [0, -1]])

    for i in range(N):
        dim_left = 2 ** (i)
        dim_right = 2 ** (N - i - 1)

        identity_left = np.eye(dim_left)
        identity_right = np.eye(dim_right)

        mat = np.kron(identity_left, pauli_z)
        mat = np.kron(mat, identity_right)
        mats.append(mat)
        
    return mats

def get_hamiltonian_terms(N):
    # Get Sigma X Matrices for each spin
    x_mats = get_sigma_matrices_x(N)
    # Sum up those matrices to get them x-term
    x_term = np.sum(x_mats, axis=0)
    # Get Sigma Z Matrices for each spin
    z_mats = get_sigma_matrices_z(N)
    # Impose periodic boundary conditions
    z_mats.append(z_mats[0])
    # mulitply adjacent matrices and sum them all up
    z_term = np.sum([np.dot(z_mats[i], z_mats[i + 1]) for i in range(len(z_mats) - 1)], axis=0)

    return z_term, x_term

def construct_hamiltonian(z_term, x_term, h_J):
    return -z_term + h_J * x_term

def compute_Z(e_vals, T_J):
    return np.sum([math.exp(-e_val * 1. / T_J) for e_val in e_vals])

def compute_E(e_vals, T_J, N):
    Z = compute_Z(e_vals, T_J)
    return (-1. / N) * np.sum([-e_val * math.exp(-e_val * 1. / T_J) for e_val in e_vals]) / Z

def compute_C(e_vals, T_J, N):
    Z = compute_Z(e_vals, T_J)
    sum_1 = np.sum([e_val**2 * math.exp(-e_val * 1. / T_J) for e_val in e_vals]) / Z
    sum_2 = np.sum([-e_val * math.exp(-e_val * 1. / T_J) for e_val in e_vals]) / Z
    return (1. / (N * T_J**2)) * (sum_1 - sum_2**2)

def make_plots(N, h_Js, T_Js, out, titles=None):
    z_term, x_term = get_hamiltonian_terms(N)
    
    fs = 12
    f1, ax1 = plt.subplots(figsize=(fs, fs))
    f2, ax2 = plt.subplots(figsize=(fs, fs))
    
    for h_J in h_Js:
        H = construct_hamiltonian(z_term, x_term, h_J)
        e_vals, e_vecs = np.linalg.eig(H)
        
        E_Js = []
        C_Js = []

        for T_J in reversed(T_Js):            
            try:
                E_J = compute_E(e_vals, T_J, N)
            except:
                E_J = E_Js[-1]
                
            try:
                C_J = compute_C(e_vals, T_J, N)
            except:
                C_J = C_Js[-1]
                
            E_Js.append(E_J)
            C_Js.append(C_J)
        
        if h_J == h_Js[0]:
            l = "h/J = {}".format(h_J)
        else:
            l = str(h_J)
        
        E_Js = list(reversed(E_Js))
        C_Js = list(reversed(C_Js))
        
        ax1.plot(T_Js, E_Js, label=l)
        ax2.plot(T_Js, C_Js, label=l)
    
    ax1.set_xlabel("T/J")
    ax1.set_ylabel("E/J")
    
    ax2.set_xlabel("T/J")
    ax2.set_ylabel("C/J")
    
    ax1.legend()
    ax2.legend()

    ax1.set_xlim(xmin=min(T_Js), xmax=max(T_Js))
    ax2.set_xlim(xmin=min(T_Js), xmax=max(T_Js))
    
    if titles is not None:
        ax1.set_title(titles[0])
        ax2.set_title(titles[1])
    
    if out is not None:
        plt.figure(f1.number)
        out_name = "Energy_{}".format(out)
        print("Saving plot to {}".format(out_name))
        plt.savefig(out_name)

        plt.figure(f2.number)
        out_name = "Heat_Capacity_{}".format(out)
        print("Saving plot to {}".format(out_name))
        plt.savefig(out_name)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("list_of_N", type=int, nargs="+", help="The list of the number of spins to test for the Ising ring model. e.g. 3 4 5 10")

    args = parser.parse_args()
    Ns = args.list_of_N

    # Choose a large number of evenly spaced values of T/J
    T_Js = np.linspace(0, 3.5, 1000)
    # Choose a selection of values of h/J
    H_Js = [0, 0.5, 1.0, 1.5, 2.0]

    for N in Ns:
        make_plots(N, H_Js, T_Js, out="N_{}.png".format(N), titles=['N = {} Ising Chain Energy'.format(N), 'N = {} Ising Chain Heat Capacity'.format(N)])

if __name__ == '__main__':
    main()

