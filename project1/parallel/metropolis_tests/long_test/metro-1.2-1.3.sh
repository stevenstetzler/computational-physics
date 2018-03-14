#!/bin/bash
#SBATCH --output=metro-1.2-1.3.out
#SBATCH --job-name="m 1.2-1.3"
./metropolis 1.217258064516129 1.337983870967742 4 metro-1.2-1.3.txt
