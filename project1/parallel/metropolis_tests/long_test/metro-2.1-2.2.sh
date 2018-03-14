#!/bin/bash
#SBATCH --output=metro-2.1-2.2.out
#SBATCH --job-name="m 2.1-2.2"
./metropolis 2.0623387096774195 2.1830645161290323 4 metro-2.1-2.2.txt
