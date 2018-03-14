#!/bin/bash
#SBATCH --output=wolff-2.5-2.9.out
#SBATCH --job-name="2.5-2.9"
./wolff 2.545241935483871 2.9074193548387095 10 wolff-2.5-2.9.txt
