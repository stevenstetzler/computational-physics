#!/bin/bash
#SBATCH --output=wolff-1.7-2.1.out
#SBATCH --job-name="1.7-2.1"
./wolff 1.7001612903225807 2.1025806451612903 11 wolff-1.7-2.1.txt
