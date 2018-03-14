#!/bin/bash
#SBATCH --output=wolff-2.1-2.5.out
#SBATCH --job-name="2.1-2.5"
./wolff 2.142822580645161 2.505 10 wolff-2.1-2.5.txt
