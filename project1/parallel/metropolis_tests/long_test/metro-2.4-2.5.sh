#!/bin/bash
#SBATCH --output=metro-2.4-2.5.out
#SBATCH --job-name="m 2.4-2.5"
./metropolis 2.384274193548387 2.505 4 metro-2.4-2.5.txt
