#!/bin/bash -l

# -----------------------------------------------
# Environment setup
# -----------------------------------------------


#SBATCH --job-name=fermilat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL

echo "~~~~~~~~ SETUP ENVIRONMENT ~~~~~~~~ "

__conda_setup="$('/net/diva/scratch1/mnievas/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/net/diva/scratch1/mnievas/miniconda/etc/profile.d/conda.sh" ]; then
        . "/net/diva/scratch1/mnievas/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/net/diva/scratch1/mnievas/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate fermi

# This will allow FTOOLS to run on the cluster, whereas if this is not
# set you get an error about '/dev/tty' e.g. when running ftcopy
export HEADASNOQUERY=1
#setenv HEADASNOQUERY 1
echo HEADASNOQUERY: $HEADASNOQUERY
echo PYTHONPATH: $PYTHONPATH 
echo PATH: $PATH
