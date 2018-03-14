#!/bin/bash
#SBATCH --output=metro-0.7-0.8.out
#SBATCH --job-name="m 0.7-0.8"
./metropolis 0.6941129032258065 0.8148387096774194 4 metro-0.7-0.8.txt
