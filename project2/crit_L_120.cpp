#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
using namespace std;

const int size_x = 120;
const int size_y = 120;
const int num_nn = 4;

const int num_therm_sweeps = 10000;
const int num_sweeps = 150;
const int num_snapshots = 10000;

long num_accepted = 0;

double probs[10] = {0};

double compute_H(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J) {
    double H = 0;
    int i, j, left, right, top, bottom, sigma_ij, neighbor_left, neighbor_right, neighbor_top, neighbor_bottom;

    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            sigma_ij = lattice[i][j];

            left = neighbors[i][j][0];
            right = neighbors[i][j][1];
            bottom = neighbors[i][j][2];
            top = neighbors[i][j][3];

            neighbor_left = lattice[left][j];
            neighbor_right = lattice[right][j];
            neighbor_top = lattice[i][top];
            neighbor_bottom = lattice[i][bottom];

            H += -J * sigma_ij * ( neighbor_left + neighbor_right + neighbor_top + neighbor_bottom );
        }
    }

    return H / 2;
}

double compute_M(int lattice[size_x][size_y], double J) {
    int i, j;
    double M = 0;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            M += lattice[i][j];
        }
    }
    return fabs(M);
}

void initialize_lattice(int lattice[size_x][size_y]) {
    int i, j;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            lattice[i][j] = 2 * (rand() % 2) - 1;
        }
    }
}

void initialize_neighbors(int neighbors[size_x][size_y][num_nn]) {
    int curr_x, curr_y, left, right, bottom, top;

    for (curr_x = 0; curr_x < size_x; curr_x++) {
        for (curr_y = 0; curr_y < size_y; curr_y++) {
            left = curr_x - 1;
            right = curr_x + 1;
            top = curr_y - 1;
            bottom = curr_y + 1;

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

            neighbors[curr_x][curr_y][0] = left;
            neighbors[curr_x][curr_y][1] = right;
            neighbors[curr_x][curr_y][2] = bottom;
            neighbors[curr_x][curr_y][3] = top;
        }
    }
}

void metropolis_step(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], int i, int j, double J, double beta) {
    double H_mol = 0;
    int left, right, top, bottom;

    left = neighbors[i][j][0];
    right = neighbors[i][j][1];
    bottom = neighbors[i][j][2];
    top = neighbors[i][j][3];

    int prob_idx;
    int sum_spins = lattice[left][j] + lattice[right][j] + lattice[i][top] + lattice[i][bottom];

    prob_idx = (sum_spins + 4) / 2;
    
    if (lattice[i][j] == -1) {
        prob_idx += 5;
    }

    if (drand48() < probs[prob_idx]) {
        lattice[i][j] = -lattice[i][j];
        num_accepted += 1;
    }
}

void sweep_metropolis(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double beta) {
    int i, j;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            metropolis_step(lattice, neighbors, i, j, J, beta);
        }
    }
}

void thermalize(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double beta) {
    int i;
    for (i = 0; i < num_therm_sweeps; i++) {
        sweep_metropolis(lattice, neighbors, J, beta);
    }
}

void crit_exp_simulation(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double min_T, double max_T, int num_T, string filename) {
    double Ts[num_T];
    double Es[num_T];
    double Cs[num_T];
    double Ms[num_T];
    double Xs[num_T];
    double Bs[num_T];

    for (int i = 0; i < num_T; i++) {
        Ts[i] = min_T + (i / (double) (num_T - 1)) * (max_T - min_T);
    }

    double beta;

    double e, avg_e, avg_e2, m, avg_m, avg_m2, avg_m4;

    long num_calls = 1;
    num_calls *= num_sweeps;
    num_calls *= num_snapshots;
    num_calls *= size_x;
    num_calls *= size_y;

    for (int i = 0; i < num_T; i++) {

        beta = 1. / Ts[i];

        // Pre-compute the acceptance probabilities for all energy levels
        probs[0] = exp( -beta * 2 * (1) * J * (-4) );
        probs[1] = exp( -beta * 2 * (1) * J * (-2) );
        probs[2] = exp( -beta * 2 * (1) * J * (0) );
        probs[3] = exp( -beta * 2 * (1) * J * (2) );
        probs[4] = exp( -beta * 2 * (1) * J * (4) );
        probs[5] = exp( -beta * 2 * (-1) * J * (-4) );
        probs[6] = exp( -beta * 2 * (-1) * J * (-2) );
        probs[7] = exp( -beta * 2 * (-1) * J * (0) );
        probs[8] = exp( -beta * 2 * (-1) * J * (2) );
        probs[9] = exp( -beta * 2 * (-1) * J * (4) );

        thermalize(lattice, neighbors, J, beta);

        num_accepted = 0;

        avg_e = 0;
        avg_e2 = 0;

        avg_m = 0;
        avg_m2 = 0;
        avg_m4 = 0;

        for (int j = 0; j < num_snapshots; j++) {
            e = compute_H(lattice, neighbors, J);
            avg_e += e;
            avg_e2 += e*e;

            m = compute_M(lattice, J);
            avg_m += m;
            avg_m2 += m*m;
            avg_m4 += m*m*m*m;

            if (j % 100 == 0) {
                printf("> %0.2f%% done       \r", 100. * (double) (i * num_snapshots + j) / (num_T * num_snapshots));
                fflush(stdout);
            }

            for (int k = 0; k < num_sweeps; k++) {
                sweep_metropolis(lattice, neighbors, J, beta);
            }
        }

        avg_e /= num_snapshots;
        avg_e2 /= num_snapshots;
        
        avg_m /= num_snapshots;
        avg_m2 /= num_snapshots;
        avg_m4 /= num_snapshots;

        Es[i] = avg_e / (size_x * size_y);
        Cs[i] = (beta*beta) * (avg_e2 - avg_e*avg_e) / (size_x * size_y);
        Ms[i] = avg_m / (size_x * size_y);
        Xs[i] = beta * (avg_m2 - avg_m*avg_m) / (size_x * size_y);
        Bs[i] = 1. - avg_m4 / (3 * avg_m2*avg_m2);
    }

    ofstream out_file;
    out_file.open(filename.c_str());
    for (int i = 0; i < num_T; i++) {
        out_file << Ts[i] << " " << Es[i] << " " << Cs[i] << " " << Ms[i] << " " << Xs[i] << " " << Bs[i] <<  "\n";
    }
    out_file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s low_T high_T num_T outfilename \n", argv[0]);
        return 1;
    }

    double T_low = atof(argv[1]);
    double T_high = atof(argv[2]);
    int n_T = atoi(argv[3]);
    string filename = argv[4];

    srand(time(NULL));

    int lattice[size_x][size_y];
    int neighbors[size_x][size_y][num_nn];

    initialize_lattice(lattice);
    initialize_neighbors(neighbors);

    double J = 1;

    crit_exp_simulation(lattice, neighbors, J, T_low, T_high, n_T, filename); 
}
