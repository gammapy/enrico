# enrico-init.csh sets some environment variables for enrico.
# 
# One way is to put this e.g. in ~/.cshrc or ~/.login:
#
# setenv FERMI_DIR < location of your Fermi software installation >
# setenv FERMI_DATA_DIR < location of your Fermi weekly and preprocessed data >
# alias init_fermi "source $FERMI_DIR/fermi-init.csh"
#
# setenv ENRICO_DIR < location of your Enrico software checkout >
# alias init_enrico "source $ENRICO_DIR/enrico-init.csh"
#
# This way your python, PATH, PYTHONPATH, ... is not set to the Fermi
# or Enrico software when you log in to your system, yet you can
# conveniently set up your shell for Fermi and Enrico by executing the aliases:
# $ init_fermi
# $ init_enrico

echo "Adding Enrico to PATH and PYTHONPATH"
setenv PATH ${PATH}:${ENRICO_DIR}/bin
setenv PYTHONPATH ${ENRICO_DIR}:${PYTHONPATH}

echo "Setting environment variables for data file locations"
if ( ! $?FERMI_DATA_DIR ) then
    setenv FERMI_DATA_DIR /sps/hess/prod/fermi
endif

if ( ! $?FERMI_CATALOG_DIR ) then
    setenv FERMI_CATALOG_DIR $FERMI_DATA_DIR/catalog
endif

if ( ! $?FERMI_DIFFUSE_DIR ) then
    setenv FERMI_DIFFUSE_DIR $FERMI_DATA_DIR/diffuse
endif

if ( ! $?FERMI_DOWNLOAD_DIR ) then
    setenv FERMI_DOWNLOAD_DIR $FERMI_DATA_DIR/download
endif

if ( ! $?FERMI_PREPROCESSED_DIR ) then
    setenv FERMI_PREPROCESSED_DIR $FERMI_DATA_DIR/preprocessed
endif

mkdir -p $FERMI_DATA_DIR
mkdir -p $FERMI_DIFFUSE_DIR
mkdir -p $FERMI_DOWNLOAD_DIR
mkdir -p $FERMI_PREPROCESSED_DIR
mkdir -p $FERMI_CATALOG_DIR

# Only set this environment variable to the default if it isn't set already,
# in order not to overwrite the user's setting.
if ( ! $?FARM ) then
    #Currently supported : LAPP-Annecy (LAPP), MPIK-Heidelberg (MPIK), CCIN2P3
    setenv FARM CCIN2P3
endif

echo "Run enrico_setupcheck to check if you are set up correctly."
echo "Check your PATH and PYTHONPATH if enrico_setup is not found or gives an error."
