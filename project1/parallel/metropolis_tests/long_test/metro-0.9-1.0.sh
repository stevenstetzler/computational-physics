#!/bin/bash
#SBATCH --output=metro-0.9-1.0.out
#SBATCH --job-name="m 0.9-1.0"
./metropolis 0.8550806451612903 1.0160483870967743 5 metro-0.9-1.0.txt
