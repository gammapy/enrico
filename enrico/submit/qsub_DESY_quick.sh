#!/usr/bin/bash

# -----------------------------------------------
# Environment setup
# -----------------------------------------------

source ~/.bashrc
source ~/.bashrc_loadfermi

echo "~~~~~~~~ SETUP ENVIRONMENT ~~~~~~~~ "

# This will allow FTOOLS to run on the cluster, whereas if this is not
# set you get an error about '/dev/tty' e.g. when running ftcopy
export HEADASNOQUERY=1
echo HEADASNOQUERY: $HEADASNOQUERY
#PBS -q std.q


