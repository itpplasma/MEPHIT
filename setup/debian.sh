#!/bin/bash

sudo apt-get update -y && sudo apt-get install -y -q --no-install-recommends \
    cmake make ninja-build git gfortran gcc patch \
    libopenblas-dev libsuitesparse-dev libmetis-dev zlib1g-dev \
    libboost-dev \
    openmpi-bin libopenmpi-dev python3 python3-dev python3-pip python3-numpy
