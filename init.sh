# Example enrico init script for the MPIK cluster
# TODO: Until we find a better solution (e.g. ask user interactively?)
# you have to manually adapt it to your system.
# Maybe the best is to only work with some environment variables
# FERMI_DIR, FERMI_DATA_DIR and enrico is installed, i.e. no ENRICO_DIR?

# This should be done as an extra step I think
#echo "Initializing Fermi Science Tools"
#export FERMI_DIR=/lfs/l1/hess/users/deil/bin/HESSSurveyBin/fermi/latest/x86_64-unknown-linux-gnu-libc2.5
#source $FERMI_DIR/fermi-init.sh

echo "Adding Enrico to PATH and PYTHONPATH"
export ENRICO_DIR=/home/hfm/deil/git/enrico
export PATH=$PATH:$ENRICO_DIR/scripts
export PYTHONPATH=$ENRICO_DIR:$PYTHONPATH

echo "Setting environment variables for data file locations"
export FERMI_DATA_DIR=/home/hfm/deil/storage/fermi/public
export FERMI_CATALOG_DIR=$FERMI_DATA_DIR/catalog
export FERMI_DIFFUSE_DIR=$FERMI_DATA_DIR/diffuse
export FERMI_DOWNLOAD_DIR=$FERMI_DATA_DIR/download
export FERMI_PREPROCESSED_DIR=$FERMI/preprocessed

echo "Run enrico_setup to check if you are set up correctly."
echo "Check your PATH and PYTHONPATH if enrico_setup is not found or gives an error."
