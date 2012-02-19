# enrico init script for my Mac

echo "Adding Enrico to PATH and PYTHONPATH"
export ENRICO_DIR=/Users/deil/code/enrico
export PATH=$PATH:$ENRICO_DIR/scripts
export PYTHONPATH=$ENRICO_DIR:$PYTHONPATH

echo "Setting environment variables for data file locations"
export FERMI_DATA_DIR=/Users/deil/data/fermi
export FERMI_CATALOG_DIR=$FERMI_DATA_DIR/catalog
export FERMI_DIFFUSE_DIR=$FERMI_DATA_DIR/diffuse
export FERMI_DOWNLOAD_DIR=$FERMI_DATA_DIR/download
export FERMI_PREPROCESSED_DIR=$FERMI_DATA_DIR/preprocessed

echo "Run enrico_setup to check if you are set up correctly."
echo "Check your PATH and PYTHONPATH if enrico_setup is not found or gives an error."

