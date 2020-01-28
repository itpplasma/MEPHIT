#!/bin/bash
FreeFem++-mpi -nw $(dirname "$0")/../FEM/maxwell.edp "$@" 1>> freefem.out 2>&1
