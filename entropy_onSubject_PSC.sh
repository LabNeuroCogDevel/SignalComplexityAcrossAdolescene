#!/bin/bash
#SBATCH --job-name=MSE_MGS
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64

# Bash script for running MSE on the supercomputer for the Entropy project
module load matlab
matlab -nodesktop -r "run_Entropy_PSC('$INPUTFILE','MGS', 'delay', 6)"

