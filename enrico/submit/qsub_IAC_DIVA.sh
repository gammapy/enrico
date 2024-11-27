#!/bin/bash -l

# -----------------------------------------------
# Environment setup
# -----------------------------------------------


#SBATCH --job-name=fermilat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH --partition=batch
#SBATCH --time=3:00:00


echo "~~~~~~~~ SETUP ENVIRONMENT ~~~~~~~~ "

__conda_setup="$('/net/diva/scratch/mnievas/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/net/diva/scratch/mnievas/miniconda/etc/profile.d/conda.sh" ]; then
        . "/net/diva/scratch/mnievas/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/net/diva/scratch/mnievas/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/net/diva/scratch/mnievas/miniconda/etc/profile.d/mamba.sh" ]; then
    . "/net/diva/scratch/mnievas/miniconda/etc/profile.d/mamba.sh"
fi




export PATH=/net/diva/scratch/mnievas/miniconda/bin/:$PATH
conda activate fermi
export FARM='IAC_DIVA'
export QUEUE='std'
export ENRICO_DIR=/net/diva/scratch/mnievas/enrico/
source $ENRICO_DIR/enrico-init.sh
export USE_FULLMISSION_SPACECRAFT=True
export COMPRESS_WEEKLY_FILES=True
