#!/bin/bash
#SBATCH --output=metro-3.6-3.7.out
#SBATCH --job-name="m 3.6-3.7"
./metropolis 3.551290322580645 3.672016129032258 4 metro-3.6-3.7.txt
