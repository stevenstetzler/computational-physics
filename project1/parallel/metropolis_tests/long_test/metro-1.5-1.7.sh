#!/bin/bash
#SBATCH --output=metro-1.5-1.7.out
#SBATCH --job-name="m 1.5-1.7"
./metropolis 1.5391935483870969 1.65991935483871 4 metro-1.5-1.7.txt
