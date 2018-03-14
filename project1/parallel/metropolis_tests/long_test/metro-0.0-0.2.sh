#!/bin/bash
#SBATCH --output=metro-0.0-0.2.out
#SBATCH --job-name="m 0.0-0.2"
./metropolis 0.01 0.1709677419354839 5 metro-0.0-0.2.txt
