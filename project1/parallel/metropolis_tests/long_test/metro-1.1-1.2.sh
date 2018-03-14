#!/bin/bash
#SBATCH --output=metro-1.1-1.2.out
#SBATCH --job-name="m 1.1-1.2"
./metropolis 1.0562903225806453 1.177016129032258 4 metro-1.1-1.2.txt
