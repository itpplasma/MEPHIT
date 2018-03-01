#!/bin/sh

gfortran -o tri_mesh.x SRC/tri_mesh.f90
make -f Axis.mk
make -f RC.mk
make -f Writeout.mk
./axis.x && ./tri_mesh.x && ./readcarre_m.x && ./writeout.x
python extmesh.py inputformaxwell

