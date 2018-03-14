#!/bin/bash
#SBATCH --output=metro-3.2-3.4.out
#SBATCH --job-name="m 3.2-3.4"
./metropolis 3.2293548387096775 3.3500806451612903 4 metro-3.2-3.4.txt
