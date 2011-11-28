#!/usr/bin/bash

# -----------------------------------------------
# Environment setup
# -----------------------------------------------

source ~/.bashrc
source $ASTRO_INIT/init.sh
source $HEADAS/headas-init.sh
source $FERMI_DIR/fermi-init.sh
export PYTHONPATH=
source $FERMI_DIR/fermi-init.sh;
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

# This will echo each command before executing it
set -x

# This will allow FTOOLS to run on the cluster, whereas if this is not
# set you get an error about '/dev/tty' e.g. when running ftcopy
export HEADASNOQUERY=1
echo HEADASNOQUERY: $HEADASNOQUERY

echo $PWD
echo $HOSTNAME

echo FERMI_DIR: $FERMI_DIR
echo HEADAS: $HEADAS

export PATH=$PATH:/home/hfm/deil/cvs/Galactic.Population.Study/src/scripts/fermi/survey
export PYTHONPATH=$PYTHONPATH:/home/hfm/deil/cvs/Galactic.Population.Study/src/scripts/fermi/utils

export FERMI=/home/hfm/deil/storage/fermi/public
echo FERMI: $FERMI
echo PATH: $PATH
echo PYTHONPATH: $PYTHONPATH

# -----------------------------------------------
# Checks for debugging cluster processing
# -----------------------------------------------

#which ftcopy
#which gtselect
#ftcopy asdf=asdf
#gtselect asdf=asdf

# -----------------------------------------------
# Actual commands to be executed
# -----------------------------------------------

export HOME={workdir}
cd {workdir}
echo workdir: $PWD
make #-i -k
