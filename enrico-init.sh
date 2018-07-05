# enrico-init.sh sets some environment variables for enrico.
# 
# One way is to put this e.g. in ~/.bashrc or ~/.profile:
#
# export FERMI_DIR=< location of your Fermi software installation >
# export FERMI_DATA_DIR=< location of your Fermi weekly and preprocessed data >
# alias init_fermi="source $FERMI_DIR/fermi-init.sh"
#
# export ENRICO_DIR=< location of your Enrico software checkout >
# alias init_enrico="source $ENRICO_DIR/enrico-init.sh"
#
# This way your python, PATH, PYTHONPATH, ... is not set to the Fermi
# or Enrico software when you log in to your system, yet you can
# conveniently set up your shell for Fermi and Enrico by executing the aliases:
# $ init_fermi
# $ init_enrico

#export ROOTSYS=/afs/ifh.de/group/cta/cta/software/root/root_v5.34.32
#source $ROOTSYS/bin/thisroot.sh

echo "Adding Enrico to PATH and PYTHONPATH"


export PATH=$PATH:$ENRICO_DIR/bin
export PYTHONPATH=$ENRICO_DIR:$PYTHONPATH

echo "Setting environment variables for data file locations"
if [ ! -d "${FERMI_DATA_DIR}" ]; then
    export FERMI_DATA_DIR=$ENRICO_DIR/Data
fi

if [ ! -d "${FERMI_CATALOG_DIR}" ]; then
    export FERMI_CATALOG_DIR=$FERMI_DATA_DIR/catalog
fi

if [ ! -d "${FERMI_DIFFUSE_DIR}" ]; then
    export FERMI_DIFFUSE_DIR=$FERMI_DATA_DIR/diffuse
fi

if [ ! -d "${FERMI_DOWNLOAD_DIR}" ]; then
    export FERMI_DOWNLOAD_DIR=$FERMI_DATA_DIR/download
fi

if [ ! -d "${FERMI_PREPROCESSED_DIR}" ]; then
    export FERMI_PREPROCESSED_DIR=$FERMI_DATA_DIR/preprocessed
fi

mkdir -p $FERMI_DATA_DIR
mkdir -p $FERMI_DIFFUSE_DIR
mkdir -p $FERMI_DOWNLOAD_DIR
mkdir -p $FERMI_PREPROCESSED_DIR
mkdir -p $FERMI_CATALOG_DIR

# Only set this environment variable to the default if it isn't set already,
# in order not to overwrite the user's setting.
if [ -z "${FARM}" ]; then
    #Currently supported : LAPP-Annecy (LAPP), MPIK-Heidelberg (MPIK), CCIN2P3
    export FARM=LAPP
fi

echo "Run enrico_setupcheck to check if you are set up correctly."
echo "Check your PATH and PYTHONPATH if enrico_setup is not found or gives an error."
