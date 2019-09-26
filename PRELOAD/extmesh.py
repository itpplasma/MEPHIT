#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 18:45:11 2016

@author: Christopher Albert
"""

import numpy as np
from io import StringIO
import sys
import codecs

from meshpy import triangle 

infile = sys.argv[1]

# read input
with codecs.open(infile + '.msh', mode = 'r', encoding = 'utf-8') as f:
    data = f.readlines()

# generate arrays for nodes, triangles, boundary edges    
[NN,NT,NE] = np.genfromtxt(StringIO(data[0]),dtype=int)

node = np.genfromtxt(StringIO(''.join(data[1:NN+1])),dtype=float)
nlab = np.array(node[:,2],dtype=int)
node = node[:,0:2]

tri = np.genfromtxt(StringIO(''.join(data[NN+1:NN+NT+1])),dtype=int)
tlab = tri[:,3]
tri = tri[:,0:3]

edge = np.genfromtxt(StringIO(''.join(data[NN+NT+1:])),dtype=int)
elab = edge[:,2]
edge = edge[:,0:2]

# swap node order in reversely indexed triangles
for kt in range(NT):
    v = node[tri[kt,:]-1,:]
    det = np.linalg.det(np.c_[v,np.ones(3)])
    if det<0:
        h = tri[kt,0]
        tri[kt,0] = tri[kt,1]
        tri[kt,1] = h
 
# write fixed mesh with det>0 for all triangles       
#with open(infile+'_fixed.msh', 'w') as f:        
#    f.write('{} {} {}\n'.format(NN,NT,NE))
#    for kn in range(NN):
#        nu=0
#        if (kn+1) in edge:
#            nu=1
#        f.write('{:.17f} {:.17f} {}\n'.format(node[kn,0],node[kn,1],nu))
#    for kt in range(NT):
#        f.write('{} {} {} {}\n'.format(
#        tri[kt,0],tri[kt,1],tri[kt,2],0))
#    for ke in range(NE):
#        f.write('{} {} {}\n'.format(edge[ke,0],edge[ke,1],1))


# find boundary nodes (label is 1 there)
cond = nlab > 0
eind = np.where(cond)
epoints = node[cond,:]

# construct center of mesh
r0 = np.mean(node,0)
r0[0] = r0[0]+30. 

# boundary nodes for inner and outer boundary (equal number)
points = np.zeros([2*len(epoints),2])
points[0:len(epoints)] = epoints # inner

# construct outer boundary (ellipse with aspect ratio 1.5)
phi = np.linspace(0, 2*np.pi, len(epoints), endpoint=False)
points[len(epoints):,:] = \
    np.array([r0[0]+(r0[0]-1e-5)*np.cos(phi),\
    r0[1]+1.5*r0[0]*np.sin(phi)]).transpose() # outer
facets = []

# construct inner boundary edge connectivity
for e in edge:
    n1 = np.nonzero(eind==e[0]-1)[1][0]
    n2 = np.nonzero(eind==e[1]-1)[1][0]
    facets.append((n1,n2))

# construct outer boundary edge connectivity
for k in range(len(edge)-1):
    n1 = len(epoints)+k
    n2 = len(epoints)+k+1
    facets.append((n1,n2))
facets.append((len(points)-1,len(points)-len(epoints)))

# initialize MeshPy 2D (triangle) for outer mesh
info = triangle.MeshInfo()
info.set_points([(p[0],p[1]) for p in points])
info.set_facets([(fa[0],fa[1]) for fa in facets])
info.set_holes([(r0[0], r0[1])]) # one point inside inner region

# how to refine:
# the further away from the center, the coarser the mesh
# max_area: how big is a triangle allowed to be
def needs_refinement(vertices, area):
    bary = np.sum(np.array(vertices), axis=0)/3.0
    max_area = .003*((bary[0]-(r0[0]-50))**2+((bary[1]-(r0[1])+25)/1.5)**2)
    return bool(area > max_area)

# build outer mesh
# allow_boundary_steiner=False : do not split boundary edges
mesh = triangle.build(info, refinement_func=needs_refinement,
                      allow_boundary_steiner=False)

mesh_points = np.array(mesh.points)
mesh_tris = np.array(mesh.elements)

# plot outer mesh
import matplotlib.pyplot as plt
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)

# add new mesh parts
NNE = NN + len(mesh_points) - len(epoints)
NTE = NT + len(mesh_tris)
NEE = len(facets[len(epoints):])

for kt in range(len(mesh_tris)):
    for k in range(3):
        if(mesh_tris[kt,k] < len(epoints)):
            mesh_tris[kt,k] = eind[0][mesh_tris[kt,k]] + 1
        else:
            mesh_tris[kt,k] = mesh_tris[kt,k] - len(epoints) + NN + 1
            
node = np.concatenate((node,mesh_points[len(epoints):,:]))
tri = np.concatenate((tri,mesh_tris))
edge = np.array(facets[len(epoints):]) - len(epoints) + NN + 1

# plot full mesh
plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=0.05)
plt.axis('equal')
plt.show()
plt.savefig('mesh.png', dpi=300)


for kt in range(NTE):
    v = node[tri[kt,:]-1,:]
    det = np.linalg.det(np.c_[v,np.ones(3)])
    if det<0:
        h = tri[kt,0]
        tri[kt,0] = tri[kt,1]
        tri[kt,1] = h
        #print("{}: det<0 fixed".format(kt+1))
        
# write extended mesh    
with open(infile+'_ext.msh', 'w') as f:        
    f.write('{} {} {}\n'.format(NNE,NTE,NEE))
    for kn in range(NNE):
        nu=0
        if (kn+1) in edge:
            nu=1
        f.write('{:.17f} {:.17f} {}\n'.format(node[kn,0],node[kn,1],nu))
    for kt in range(NT):
        f.write('{} {} {} {}\n'.format(
        tri[kt,0],tri[kt,1],tri[kt,2],0))
    for kt in range(NT,NTE):
        f.write('{} {} {} {}\n'.format(
        tri[kt,0],tri[kt,1],tri[kt,2],3))
    for ke in range(NEE):
        f.write('{} {} {}\n'.format(edge[ke,0],edge[ke,1],1))
        
# write extended mesh for GMSH
        
