#!/bin/bash
#SBATCH --job-name=MSE_MGS
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=64

# Bash script for running MSE on the supercomputer for the Entropy project

matlab -nodesktop -r "run_Entropy_PSC('$inputfile','MGS', 'delay', 8); run_Entropy_PSC('$inputfile','MGS', 'delay', 6)"

