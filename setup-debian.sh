#!/bin/bash

sudo apt update
sudo apt install cmake make ninja-build gcc gfortran \
    libopenblas-dev libsuitesparse-dev libhdf5-dev \
    libnetcdf-dev libnetcdff-dev libboost-dev \
    libfftw3-dev
