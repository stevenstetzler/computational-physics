#!/bin/bash
#SBATCH --output=metro-3.7-3.8.out
#SBATCH --job-name="m 3.7-3.8"
./metropolis 3.712258064516129 3.832983870967742 4 metro-3.7-3.8.txt
