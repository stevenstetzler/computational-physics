#!/bin/bash
#SBATCH --output=metro-0.4-0.5.out
#SBATCH --job-name="m 0.4-0.5"
./metropolis 0.3721774193548387 0.49290322580645163 4 metro-0.4-0.5.txt
