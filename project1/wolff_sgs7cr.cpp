#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <algorithm>
using namespace std;

const int size_x = 30;
const int size_y = 30;
const int num_nn = 4;

long cluster_size = 0;

void print_vector(vector<int> vec) {
    if (vec.size() == 0) {
        printf("\n");
        return;
    }
    int i, x, y, val;
    for (i = 0; i < vec.size() - 1; i++) {
        val = vec[i];
        y = val % size_x;
        x = (val - y) / size_x;
        printf("(%d, %d), ", x, y);
    }
    val = vec[vec.size() - 1];
    y = val % size_x;
    x = (val - y) / size_x;
    printf("(%d, %d)\n", x, y);
}

double compute_H(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], int J) {
    double H = 0;
    int i, j, left, right, top, bottom, sigma_ij, neighbor_left, neighbor_right, neighbor_top, neighbor_bottom;

    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {

            sigma_ij = lattice[i][j];

            left = neighbors[i][j][0];
            right = neighbors[i][j][1];
            top = neighbors[i][j][2];
            bottom = neighbors[i][j][3];

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

void wolff_step(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], int x, int y, double J, double beta) {

    vector<int> F_old = {x * size_x + y};

    int i, left, right, top, bottom, sigma_xy, x_flip, y_flip, curr_x, curr_y, curr_point, tag;

    double p_connect = 1 - exp(-2 * beta * J);

    int C[size_x*size_y] = {0};
    
    tag = x * size_x + y;
    C[tag] = 1;

    while (F_old.size() != 0) {
        vector<int> F_new;

        for (i = 0; i < F_old.size(); i++) {
            curr_point = F_old[i];
            curr_y = curr_point % size_x;
            curr_x = (curr_point - curr_y) / size_x;

            sigma_xy = lattice[curr_x][curr_y];
           
            
            left = neighbors[curr_x][curr_y][0];
            right = neighbors[curr_x][curr_y][1];
            bottom = neighbors[curr_x][curr_y][2];
            top = neighbors[curr_x][curr_y][3];
            
            tag = left * size_x + curr_y;
            if (C[tag] == 0 && lattice[left][curr_y] == sigma_xy) {
                if (drand48() < p_connect) {
                    F_new.push_back(tag);
                    C[tag] = 1; 
                }
            }

            tag = right * size_x + curr_y;
            if (C[tag] == 0 && lattice[right][curr_y] == sigma_xy) {
                if (drand48() < p_connect) {
                    F_new.push_back(tag);
                    C[tag] = 1; 
                }
            }

            tag = curr_x * size_x + top;
            if (C[tag] == 0 && lattice[curr_x][top] == sigma_xy) {
                if (drand48() < p_connect) {
                    F_new.push_back(tag);
                    C[tag] = 1; 
                }
            }
            
            tag = curr_x * size_x + bottom;
            if (C[tag] == 0 && lattice[curr_x][bottom] == sigma_xy) {
                if (drand48() < p_connect) {
                    F_new.push_back(tag);
                    C[tag] = 1; 
                }
            }
        }

        F_old = F_new;
    }

    for (i = 0; i < size_x * size_y; i++) {
        if (C[i] == 1) {
            y_flip = i % size_x;
            x_flip = (i - y_flip) / size_x;

            lattice[x_flip][y_flip] = -lattice[x_flip][y_flip];
            
            cluster_size += 1;
        }
    }

}

void sweep_wolff(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double beta) {
    int i, j;
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            wolff_step(lattice, neighbors, i, j, J, beta);
        }
    }
}

void thermalize(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double beta) {
    int num_sweeps = 50;
    int i;
    for (i = 0; i < num_sweeps; i++) {
        sweep_wolff(lattice, neighbors, J, beta);
    }
}

void wolff_simulation(int lattice[size_x][size_y], int neighbors[size_x][size_y][num_nn], double J, double min_T, double max_T, int num_T, string filename) {
    double Ts[num_T];
    double Es[num_T];
    double Cs[num_T];
    double C_sizes[num_T];

    for (int i = 0; i < num_T; i++) {
        Ts[i] = min_T + (i / (double) (num_T - 1)) * (max_T - min_T);
    }

    double beta;
    int num_snapshots = 2000;
    int num_sweeps = 20;

    double avg_e, avg_e2, e;

    for (int i = 0; i < num_T; i++) {
        beta = 1. / Ts[i];

        thermalize(lattice, neighbors, J, beta);

        cluster_size = 0;

        avg_e = 0;
        avg_e2 = 0;

        for (int j = 0; j < num_snapshots; j++) {
            e = compute_H(lattice, neighbors, J);
            avg_e += e;
            avg_e2 += e*e;

            for (int k = 0; k < num_sweeps; k++) {
                printf("> %0.2f%% done      \r", 100. * (double) (i * num_snapshots * num_sweeps + j * num_sweeps + k) / (num_T * num_snapshots * num_sweeps));
                fflush(stdout);
                sweep_wolff(lattice, neighbors, J, beta);
            }
        }

        printf("Avg cluster size: %f\n", (double) cluster_size / (num_snapshots * num_sweeps * size_x * size_y));

        avg_e /= num_snapshots;
        avg_e2 /= num_snapshots;

        Es[i] = avg_e / (size_x * size_y);
        Cs[i] = beta*beta * (avg_e2 - avg_e*avg_e) / (size_x * size_y);
        C_sizes[i] = (double) cluster_size / (num_snapshots * num_sweeps * size_x * size_y);
    }
    printf("\n");


    ofstream out_file;
    out_file.open(filename.c_str());
    for (int i = 0; i < num_T; i++) {
        out_file << Ts[i] << " " << Es[i] << " " << Cs[i] << " " << C_sizes[i] << "\n";
    }
    out_file.close();

}

void print_lattice(int lattice[size_x][size_y]) {
    int i, j;

    printf("\n");
    for (i = 0; i < size_x; i++) {
        for (j = 0; j < size_y; j++) {
            if (lattice[i][j] > 0) {
                printf("  %d ", lattice[i][j]);
            } else {
                printf(" %d ", lattice[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s low_T high_T num_T, out_filename \n", argv[0]);
        return 1;
    }
    float low_T = atof(argv[1]);
    float high_T = atof(argv[2]);
    int n = atoi(argv[3]);
    string filename = argv[4];
    
    srand(time(NULL));

    int lattice[size_x][size_y];
    int neighbors[size_x][size_y][num_nn];

    initialize_lattice(lattice);   
    initialize_neighbors(neighbors);

    double J = 1;

    wolff_simulation(lattice, neighbors, J, low_T, high_T, n, filename);
}
