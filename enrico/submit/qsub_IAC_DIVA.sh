#!/bin/bash

# -----------------------------------------------
# Environment setup
# -----------------------------------------------


#SBATCH --job-name=fermilat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

echo "~~~~~~~~ SETUP ENVIRONMENT ~~~~~~~~ "

conda activate fermi

# This will allow FTOOLS to run on the cluster, whereas if this is not
# set you get an error about '/dev/tty' e.g. when running ftcopy
export HEADASNOQUERY=1
#setenv HEADASNOQUERY 1
echo HEADASNOQUERY: $HEADASNOQUERY
echo PYTHONPATH: $PYTHONPATH 
echo PATH: $PATH
