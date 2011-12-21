echo "Setting up Fermi and Enrico on the MPIK cluster"

export FERMI_DIR=/lfs/l1/hess/users/deil/bin/HESSSurveyBin/fermi/latest/x86_64-unknown-linux-gnu-libc2.5
export FERMI=/home/hfm/deil/storage/fermi/public
export FERMI_CATALOG_DIR=$FERMI/catalog
export FERMI_DIFFUSE_DIR=$FERMI/diffuse
export FERMI_DOWNLOAD_DIR=$FERMI/download
export FERMI_DATA_DIR=$FERMI/data

source $FERMI_DIR/fermi-init.sh

export ENRICO_DIR=/nfs/d22/hfm/sanchez/enrico

export PATH=$PATH:$ENRICO_DIR/scripts
export PYTHONPATH=$PYTHONPATH:$ENRICO_DIR
