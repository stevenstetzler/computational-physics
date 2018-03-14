#!/bin/bash
#SBATCH --output=metro-4.2-4.4.out
#SBATCH --job-name="m 4.2-4.4"
./metropolis 4.235403225806452 4.356129032258065 4 metro-4.2-4.4.txt
