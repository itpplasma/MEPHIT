#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:01:58 2018

@author: calbert
"""

import numpy as np
import matplotlib.pyplot as plt
import dolfin as df
from dolfin import div

from readmsh import readmsh

infile = '../PRELOAD/inputformaxwell_ext.msh'
#currnfile = '../FEM/currtest.dat'
currnfile = '../MHD/Bn.dat'
#currnfile = '../MHD/Bn_fenics.dat'
#currnfile = '../MHD/currn.dat'

[node, tri, edge, nlab, tlab, elab] = readmsh(infile)
currn = np.loadtxt(currnfile)

mesh = df.Mesh()
meshgen = df.MeshEditor()
meshgen.open(mesh, 'triangle', 2, 2)
meshgen.init_vertices(len(node))
meshgen.init_cells(len(tri))
for kn in range(len(node)):
    meshgen.add_vertex(kn, node[kn])
for kt in range(len(tri)):
    meshgen.add_cell(kt, tri[kt,:]-1)
meshgen.close()

# Map for re-ordering of indexes and orientations
map2facing = np.array([1,2,0]) # index of edge facing vertex
map2lex = np.argsort(tri[:,map2facing],1) # map to lexical order in each triangle cell
assert np.max(mesh.cells()[:]+1 - np.sort(tri)) == 0, "triangle sort order"

# Function spaces
P0 = df.FunctionSpace(mesh, 'DG', 0)
P1 = df.FunctionSpace(mesh, 'Lagrange', 1)
Hdiv = df.FunctionSpace(mesh, 'RT', 1)
Hcurl = df.FunctionSpace(mesh, 'N1curl', 1)
Hdivmap = Hdiv.dofmap()
Hcurlmap = Hdiv.dofmap()

# Read data
cells = mesh.cells()
ncells = len(cells)
nodes = mesh.coordinates()
ndofs = len(Hdivmap.dofs())

Jvec = [np.zeros(ndofs), np.zeros(ndofs)]
J = [df.Function(Hdiv), df.Function(Hdiv)]
Jphi = [df.Function(P0), df.Function(P0)]
acell = np.zeros(ncells)

for kcell in range(len(currn)):
    nodeind = cells[kcell]
    nodecoords = nodes[nodeind]
    cellarea = .5*np.linalg.det(np.concatenate(([[1,1,1]],nodecoords.T)))
    cellorient = np.sign(cellarea)
    acell[kcell] = cellarea*cellorient
    orient = np.sign([nodeind[2]+(-1)*nodeind[1], 
                      nodeind[0]+(-1)*nodeind[2], 
                      nodeind[1]+(-1)*nodeind[0]]) # (-1) for signed int
    dofs = Hdivmap.cell_dofs(kcell) # DOFs in cell k
    for k in range(2): # real and imaginary part
        Jvec[k][dofs] = cellorient*orient*currn[kcell,2*map2lex[kcell]+k]
        Jphi[k].vector()[kcell] = currn[kcell, 6+k]/np.abs(cellarea)
        
    #print(np.sum(orient*Jvec[0][dofs])-n*Jphi[1].vector()[kcell]*cellarea)
    #print(np.sum(orient*Jvec[1][dofs])+n*Jphi[0].vector()[kcell]*cellarea)
        
J[0].vector()[:] = Jvec[0]
J[1].vector()[:] = Jvec[1]
    
    
#%% plot mesh
plt.figure()
plt.subplot(1,2,1)
plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=.1)
plt.title('mesh (triplot)')
plt.subplot(1,2,2)
df.plot(mesh, linewidth=.1)
plt.title('mesh (fenics)')

#%%
plt.figure()
p = df.plot(J[0])
#p.scale = 5e12
#p.scaleunits = 'xy'
#p.minshaft = .01
plt.colorbar(p)
plt.title('Re poloidal')

#%%
plt.figure()
p = df.plot(J[1])
plt.colorbar(p)
plt.title('Im poloidal')
plt.tight_layout()

#%%
plt.figure()
p = df.plot(div(J[0]))
plt.colorbar(p)
plt.title('Re div poloidal')
plt.figure()
p = df.plot(div(J[1]))
plt.colorbar(p)
plt.title('Im div poloidal')

plt.figure()
p=df.plot(Jphi[0])
plt.colorbar(p)
plt.title('Re toroidal')
plt.figure()
p=df.plot(Jphi[1])
plt.colorbar(p)
plt.title('Im toroidal')
    