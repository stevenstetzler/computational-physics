#!/bin/bash
#SBATCH --output=metro-3.4-3.5.out
#SBATCH --job-name="m 3.4-3.5"
./metropolis 3.390322580645161 3.5110483870967744 4 metro-3.4-3.5.txt
