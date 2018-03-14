#!/bin/bash
#SBATCH --output=metro-0.2-0.3.out
#SBATCH --job-name="m 0.2-0.3"
./metropolis 0.21120967741935487 0.3319354838709678 4 metro-0.2-0.3.txt
