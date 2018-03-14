#!/bin/bash
#SBATCH --output=metro-2.5-2.7.out
#SBATCH --job-name="m 2.5-2.7"
./metropolis 2.545241935483871 2.665967741935484 4 metro-2.5-2.7.txt
