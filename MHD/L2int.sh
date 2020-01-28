#!/bin/bash
FreeFem++-nw $(dirname "$0")/../FEM/L2int.edp "$@" 1>> freefem.out 2>&1
