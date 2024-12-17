#!/bin/bash

sudo apt-get update -y && sudo apt-get install -y -q --no-install-recommends \
    libopenblas-dev libsuitesparse-dev libmetis-dev libhdf5-dev \
    libnetcdf-dev libnetcdff-dev libboost-dev \
    libfftw3-dev openmpi-bin libopenmpi-dev
