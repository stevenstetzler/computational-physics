#!/bin/bash
#SBATCH --output=metro-3.1-3.2.out
#SBATCH --job-name="m 3.1-3.2"
./metropolis 3.0683870967741935 3.1891129032258063 4 metro-3.1-3.2.txt
