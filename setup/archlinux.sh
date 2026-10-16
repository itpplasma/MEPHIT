#!/bin/bash

sudo -Sy --needed git cmake make ninja patch gcc gcc-fortran \
    openblas zlib suitesparse boost openmpi \
    fftw gsl python python-numpy debugedit fakeroot
