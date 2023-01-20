#!/bin/bash

# Set environmental variables:
# ==============================================================================
# Parent directory:
REPO_DIR=$PWD

# export CC=gcc
# export CPP=cpp
# export CXX=g++
# export FC=gfortran

# Directory where tarball is located:
HDF5_TAR_DIR=$REPO_DIR"/tarfiles"

# HDF5 installation directory:
HDF5_INSTALL_DIR=$REPO_DIR"/lib_HDF5"

# HDF5 version to install:
HDF5_VERSION="hdf5-1.10.9"

# Clean up directories:
# ==============================================================================
# Delete any previous installation directories:
rm -r $HDF5_INSTALL_DIR

# Create new folders for local installation of HDF5:
mkdir $HDF5_INSTALL_DIR

# Delete existing HDF5 build directories:
rm -r $HDF5_VERSION

# Start local installation of HDF5 library:
# ==============================================================================
# Extract the source code of HDF5 library:
tar -xvf $HDF5_TAR_DIR"/"$HDF5_VERSION$".tar.gz"

# Enter source code directory of HDF5:
cd $HDF5_VERSION

# Set the prefix for installation folder:
PREFIX=$HDF5_INSTALL_DIR

./configure --prefix=$PREFIX --enable-cxx --enable-build-mode=production

make

if [ $? -eq 0 ] ; then
make install

cd ../

rm -r $HDF5_VERSION
else
echo 'ERROR: Uh-oh! Something went wrong in HDF5 installation!'
return
fi
