#!/bin/bash
#SBATCH --output=metro-1.7-1.8.out
#SBATCH --job-name="m 1.7-1.8"
./metropolis 1.7001612903225807 1.8208870967741937 4 metro-1.7-1.8.txt
