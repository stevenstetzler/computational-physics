#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
using namespace std;

const int size_x = 30;
const int size_y = 30;

double num_accepted = 0.;

double compute_H(int lattice[size_x][size_y], int J) {
    double H = 0;
    int i, j, left, right, top, bottom, sigma_ij, neighbor_left, neighbor_right, neighbor_top, neighbor_bottom;

    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            left = i - 1;
            right = i + 1;
            top = j - 1;
            bottom = j + 1;

            sigma_ij = lattice[i][j];
        
            if (left % size_x < 0) {
                left = left % size_x + size_x;
            } else {
                left = left % size_x;
            }


            if (right % size_x < 0) {
                right = right % size_x + size_x;
            } else {
                right = right % size_x;
            }
            
            if (top % size_y < 0) {
                top = top % size_y + size_y;
            } else {
                top = top % size_y;
            }
            
            if (bottom % size_y < 0) {
                bottom = bottom % size_y + size_y;
            } else {
                bottom = bottom % size_y;
            }
            
            neighbor_left = lattice[left][j];
            neighbor_right = lattice[right][j];
            neighbor_top = lattice[i][top];
            neighbor_bottom = lattice[i][bottom];

            H += -J * sigma_ij * ( neighbor_left + neighbor_right + neighbor_top + neighbor_bottom );
        }
    }

    return H / 2;
}

void initialize_lattice(int lattice[size_x][size_y]) {
    int i, j;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            lattice[i][j] = 2 * (rand() % 2) - 1;
        }
    }
}

void metropolis_step(int lattice[size_x][size_y], int i, int j, double J, double beta) {
    double H_mol = 0;
    int left, right, top, bottom;

    left = i - 1;
    right = i + 1;
    top = j - 1;
    bottom = j + 1;

    if (left % size_x < 0) {
        left = left % size_x + size_x;
    } else {
        left = left % size_x;
    }


    if (right % size_x < 0) {
        right = right % size_x + size_x;
    } else {
        right = right % size_x;
    }
    
    if (top % size_y < 0) {
        top = top % size_y + size_y;
    } else {
        top = top % size_y;
    }
    
    if (bottom % size_y < 0) {
        bottom = bottom % size_y + size_y;
    } else {
        bottom = bottom % size_y;
    }

    H_mol = J * ( lattice[left][j] + lattice[right][j] + lattice[i][top] + lattice[i][bottom] );

    double delta_E = 2 * lattice[i][j] * H_mol;

    if (drand48() < exp( -beta * delta_E)) {
        lattice[i][j] = -lattice[i][j];
        num_accepted += 1;
    }
}

void sweep_metropolis(int lattice[size_x][size_y], double J, double beta) {
    int i, j;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            metropolis_step(lattice, i, j, J, beta);
        }
    }
}

void thermalize(int lattice[size_x][size_y], double J, double beta) {
    int num_sweeps = 10000;
    int i;
    for (i = 0; i < num_sweeps; i++) {
        sweep_metropolis(lattice, J, beta);
    }
}

void metropolis_simulation(int lattice[size_x][size_y], double J, double min_T, double max_T, int num_T) {

    double Ts[num_T];
    double Es[num_T];
    double Cs[num_T];
    double acceptance[num_T];

    for (int i = 0; i < num_T; i++) {
        Ts[i] = min_T + (i / (double) (num_T - 1)) * (max_T - min_T);
    }

    double beta;
    int num_snapshots = 100000;
    int num_sweeps = 150;

    double avg_e, avg_e2, e;

    for (int i = 0; i < num_T; i++) {
        printf("> %0.2f%% done       \r", 100. * (double) i / num_T);
        fflush(stdout);

        beta = 1. / Ts[i];

        thermalize(lattice, J, beta);

        num_accepted = 0;
        avg_e = 0;
        avg_e2 = 0;

        for (int j = 0; j < num_snapshots; j++) {
            e = compute_H(lattice, J);
            avg_e += e;
            avg_e2 += pow(e, 2.);
            
            for (int k = 0; k < num_sweeps; k++) {
                sweep_metropolis(lattice, J, beta);
            }
        }
        
        avg_e /= num_snapshots;
        avg_e2 /= num_snapshots;
        
        Es[i] = avg_e / (size_x * size_y);
        Cs[i] = pow(beta, 2) * (avg_e2 - pow(avg_e, 2)) / (size_x * size_y);
        acceptance[i] = num_accepted / (num_snapshots * num_sweeps * size_x * size_y);
    }

    ofstream out_file;
    out_file.open("metropolis_sweep_c.txt");
    for (int i = 0; i < num_T; i++) {
        out_file << Ts[i] << " " << Es[i] << " " << Cs[i] << " " << acceptance[i] << "\n";
    }
    out_file.close();

}

int main() {
    srand(time(NULL));

    int lattice[size_x][size_y];

    initialize_lattice(lattice);

    metropolis_simulation(lattice, 1, 0.01, 5, 20); 
}
