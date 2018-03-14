#!/bin/bash
#SBATCH --output=metro-0.5-0.7.out
#SBATCH --job-name="m 0.5-0.7"
./metropolis 0.5331451612903226 0.6538709677419355 4 metro-0.5-0.7.txt
