#!/bin/bash
#SBATCH --output=wolff-0.0-0.4.out
#SBATCH --job-name="0.0-0.4"
./wolff 0.01 0.41241935483870973 11 wolff-0.0-0.4.txt
