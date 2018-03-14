#!/bin/bash
#SBATCH --output=metro-2.7-2.8.out
#SBATCH --job-name="m 2.7-2.8"
./metropolis 2.7062096774193547 2.826935483870968 4 metro-2.7-2.8.txt
