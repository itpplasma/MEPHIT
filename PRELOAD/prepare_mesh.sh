#!/bin/sh

make -f Axis.mk && ./axis.x && ./tri_mesh.x && ./readcarre_m.x && ./writeout.x
python extmesh.py inputformaxwell

