#!/bin/bash
#SBATCH --output=metro-0.0-0.1.out
#SBATCH --job-name="m 0.0-0.1"
./metropolis 0.01 0.09048381818181817 3 metro-0.0-0.1.txt
