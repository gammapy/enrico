#!/usr/local/bin/bash -f

# -----------------------------------------------
# Environment setup
# -----------------------------------------------

if [ ! -d "${PFILES}" ]; then
    mkdir ${PFILES}
fi


echo "~~~~~~~~ SETUP ENVIRONMENT ~~~~~~~~ "

# This will allow FTOOLS to run on the cluster, whereas if this is not
# set you get an error about '/dev/tty' e.g. when running ftcopy
export HEADASNOQUERY=1
#setenv HEADASNOQUERY 1
echo HEADASNOQUERY: $HEADASNOQUERY
echo PYTHONPATH: $PYTHONPATH 
