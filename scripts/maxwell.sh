#!/bin/bash
FreeFem++-mpi -nw $(dirname "$0")/maxwell.edp "$@" 1>> freefem.out 2>&1
