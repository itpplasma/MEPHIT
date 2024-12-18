#!/bin/bash

module load cmake
module load ninja
module load gcc
module load openmpi
module load mkl
module load hdf5-serial
module load netcdf-serial
module load fftw-serial
module load boost
module load gsl
module load metis
export METIS_DIR=$METIS_HOME
module load anaconda/3/2023.03

export CC=`which gcc`
export CXX=`which g++`
export FC=`which gfortran`
