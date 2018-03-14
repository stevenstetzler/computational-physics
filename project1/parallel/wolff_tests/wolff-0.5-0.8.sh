#!/bin/bash
#SBATCH --output=wolff-0.5-0.8.out
#SBATCH --job-name="0.5-0.8"
./wolff 0.4526612903225807 0.8148387096774194 10 wolff-0.5-0.8.txt
