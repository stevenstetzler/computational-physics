#!/bin/bash
#SBATCH --output=metro-1.9-2.0.out
#SBATCH --job-name="m 1.9-2.0"
./metropolis 1.8611290322580647 2.0220967741935483 5 metro-1.9-2.0.txt
