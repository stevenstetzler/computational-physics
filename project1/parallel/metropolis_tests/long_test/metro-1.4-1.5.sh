#!/bin/bash
#SBATCH --output=metro-1.4-1.5.out
#SBATCH --job-name="m 1.4-1.5"
./metropolis 1.378225806451613 1.4989516129032259 4 metro-1.4-1.5.txt
