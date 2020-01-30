#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 09:14:00 2018

@author: calbert
"""


from readmsh import readmsh

infile = '../PRELOAD/inputformaxwell_ext.msh'
[node, tri, edge, nlab, tlab, elab] = readmsh(infile)

# for meshio
#import meshio
#cells = {
#    'triangle': tri-1
#    }

#mesh = meshio.Mesh(node, cells)

#meshio.write('../PRELOAD/inputformaxwell_ext.exo', mesh)
#meshio.write('../PRELOAD/inputformaxwell_ext_gmsh.msh', mesh, write_binary=False)
#meshio.write('../PRELOAD/inputformaxwell_ext.vtk', mesh)

ntri = len(tri)
nnode = len(node)
nedge = len(edge)

with open('../PRELOAD/inputformaxwell_ext.mesh', 'w') as f:
    f.write('MFEM mesh v1.0\n\n')
    f.write('dimension 2\n')
    
    f.write('\nelements\n{}\n'.format(ntri))
    for kt in range(ntri):
        f.write('{} {} {} {} {}\n'.format(tlab[kt]+1, 2, 
                tri[kt][0]-1, tri[kt][1]-1, tri[kt][2]-1))

    f.write('\nboundary\n{}\n'.format(nedge))
    for ke in range(nedge):
        f.write('{} {} {} {}\n'.format(nlab[ke]+1, 1, 
                edge[ke][0]-1, edge[ke][1]-1))

    f.write('\nvertices\n{}\n2\n'.format(nnode))
    for kn in range(nnode):
        f.write('{} {}\n'.format(node[kn,0], node[kn,1]))

#%%