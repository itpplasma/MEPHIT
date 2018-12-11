#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:01:58 2018

@author: calbert
"""

import numpy as np
import matplotlib.pyplot as plt
import codecs
from io import StringIO
import dolfin as df
from dolfin import grad, div, curl, inner, dx

df.parameters["reorder_dofs_serial"] = False

# read msh file
def readmsh(infile):
    with codecs.open(infile, mode = 'r', encoding = 'utf-8') as f:
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
    
    return [node, tri, edge, nlab, tlab, elab]

def eps_permut(v):
    if (v == [0,1,2] or v == [1,2,0] or v == [2,0,1]): return 1
    if (v == [0,2,1] or v == [2,1,0] or v == [1,0,2]): return -1
    return 0
    

#infile = sys.argv[1]
#infile = '../PRELOAD/inputformaxwell.msh'
infile = '../PRELOAD/inputformaxwell_ext.msh'
currnfile = '../MHD/currn.dat'
plotcurrnfile = '../MHD/plot_currn.dat'

[node, tri, edge, nlab, tlab, elab] = readmsh(infile)
currn = np.loadtxt(currnfile)
plotcurrn = np.loadtxt(plotcurrnfile)

# plot mesh - without dolfin
#plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=.1)
#plt.axis('equal')
#plt.show()

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

df.plot(mesh, linewidth=.1)
#%%
map2facing = [1,2,0] # index of edge facing vertex
map2lex = np.argsort(tri,1) # map to lexical order in each triangle cell
assert np.max(mesh.cells()[:]+1 - np.sort(tri)) == 0, "triangle sort order"

#%%

P0 = df.FunctionSpace(mesh, 'DG', 0)
P1 = df.FunctionSpace(mesh, 'Lagrange', 1)
Hdiv = df.FunctionSpace(mesh, 'RT', 1)
Hcurl = df.FunctionSpace(mesh, 'N1curl', 1)
Hdivmap = Hdiv.dofmap()
Hcurlmap = Hdiv.dofmap()

# Test: plot p0 functions
ftest = np.zeros(len(tri))
ftest[0:len(plotcurrn)] = plotcurrn[:,0]

#plot for comparison:
#plt.figure()
#plt.tripcolor(node[:, 0], node[:, 1], tri-1, ftest)

plt.figure()
f = df.Function(P0)
f.vector()[:] = ftest
df.plot(f)

#%%

cells = mesh.cells()
nodes = mesh.coordinates()

Jvec0 = np.zeros(len(Hdivmap.dofs()))
J = df.Function(Hdiv)
Jphi = df.Function(P0)

for kcell in range(len(currn)):
    nodeind = cells[kcell]
    nodecoords = nodes[nodeind]
    cellorient = np.sign(np.linalg.det(np.concatenate(([[1,1,1]],nodecoords.T))))
    inddiff = np.roll(nodeind,1)+(-1)*np.roll(nodeind,-1) # (-1) to convert uint to int
    orient = np.sign(inddiff)
    dofs = Hdivmap.cell_dofs(kcell) # DOFs in cell k
    Jvec0[dofs] = cellorient*orient*currn[kcell,2*map2lex[kcell,map2facing]]
    Jphi.vector()[kcell] = currn[kcell, 6]
    
J.vector()[:] = Jvec0

#%%

plt.figure()
df.plot(Jphi)

plt.figure()
df.plot(J)

#%%

# Ampere
def boundary(_, on_boundary):
    return on_boundary

u = df.TrialFunction(Hcurl)
v = df.TestFunction(Hcurl)

A_D = df.Constant((0, 0))
bc = df.DirichletBC(Hcurl, A_D, boundary)

#Tanalyt = df.Expression(('exp(-(pow(x[0]-150,2)+pow(x[1],2))/pow(100,2))'), degree=2)
#T = df.interpolate(Tanalyt, P1)
#J = curl(T)

#plt.figure()
#p = df.plot(T)
#plt.colorbar(p)

a = (inner(curl(u), curl(v)) + .1*inner(u,v)) * dx
b = inner(J, v) * dx

A = df.Function(Hcurl)
df.solve(a == b, A, bc)

plt.figure()
p = df.plot(div(J))
plt.colorbar(p)

plt.figure()
p = df.plot(A)
plt.colorbar(p)

plt.figure()
p = df.plot(curl(A))
plt.colorbar(p)

plt.show()

#solve Ampere([ax, ay], [wx, wy], solver = UMFPACK) = // defines the PDE
#  int2d(Th)(x * (dx(wy) - dy(wx)) * (dx(ay) - dy(ax)))
#  + int2d(Th)(n^2 / x * [wx, wy] '* [ax, ay])
#  - int2d(Th)(4.0 * pi / c * [wx, wy] '* [Jr,Jz])
#  + on(1, ax = 0, ay = 0);