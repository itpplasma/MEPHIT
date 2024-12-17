#!/bin/bash

sudo apt-get update -y && sudo apt-get install -y -q --no-install-recommends \
    cmake make ninja-build git gfortran gcc patch libgsl-dev \
    libopenblas-dev libsuitesparse-dev libmetis-dev libhdf5-dev \
    libnetcdf-dev libnetcdff-dev libboost-dev libfftw3-dev \
    openmpi-bin libopenmpi-dev python3 python3-dev python3-pip python3-numpy
