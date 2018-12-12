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
from dolfin import derivative, div, curl, inner, dx

df.parameters["reorder_dofs_serial"] = False

doplot = True

c = 29979245800.0

n = 2.0

#infile = sys.argv[1]
#infile = '../PRELOAD/inputformaxwell.msh'
infile = '../PRELOAD/inputformaxwell_ext.msh'
currnfile = '../MHD/currn.dat'
plotcurrnfile = '../MHD/plot_currn.dat'
#currnfile = '../MHD/Bn_nonres.dat'
#plotcurrnfile = '../MHD/plot_Bn_nonres.dat'

#currnfile = '../PRELOAD/VACFIELD/Bn_flux.dat'

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
map2facing = np.array([2,0,1]) # index of edge facing vertex
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
nodes = mesh.coordinates()
ndofs = len(Hdivmap.dofs())

Jvec = [np.zeros(ndofs), np.zeros(ndofs)]
J = [df.Function(Hdiv), df.Function(Hdiv)]
Jphi = [df.Function(P0), df.Function(P0)]

for kcell in range(len(currn)):
    nodeind = cells[kcell]
    nodecoords = nodes[nodeind]
    cellarea = .5*np.linalg.det(np.concatenate(([[1,1,1]],nodecoords.T)))
    cellorient = np.sign(cellarea)
    inddiff = np.roll(nodeind,1)+(-1)*np.roll(nodeind,-1) # (-1) to convert uint to int
    orient = np.sign(inddiff)
    dofs = Hdivmap.cell_dofs(kcell) # DOFs in cell k
    for k in range(2): # real and imaginary part
        Jvec[k][dofs] = cellorient*orient*currn[kcell,2*map2lex[kcell]+k]
        Jphi[k].vector()[kcell] = currn[kcell, 6+k]/np.abs(cellarea)
        
    #print(np.sum(orient*Jvec[0][dofs])-n*Jphi[1].vector()[kcell]*cellarea)
    #print(np.sum(orient*Jvec[1][dofs])+n*Jphi[0].vector()[kcell]*cellarea)
        
J[0].vector()[:] = Jvec[0]
J[1].vector()[:] = Jvec[1]
    
#%% Solve Ampere
def boundary(_, on_boundary):
    return on_boundary

u = df.TrialFunction(Hcurl)
v = df.TestFunction(Hcurl)

# Homogeneous Dirichlet boundary condition
A_D = df.Constant((0, 0))
bc = df.DirichletBC(Hcurl, A_D, boundary)

rweight = df.Expression('x[0]', degree=1)

# variational bilinear form for 2D curl-curl harmonics (->stiffness matrix)
a = (rweight*inner(curl(u), curl(v)) + n**2/rweight*inner(u,v)) * dx

# linear form (->right-hand side vector)
b = []
b.append(4.0*np.pi/c*inner(J[0], v) * dx) # real part
b.append(4.0*np.pi/c*inner(J[1], v) * dx) # imaginary part

A = [df.Function(Hcurl), df.Function(Hcurl)] # poloidal vector potential
df.solve(a == b[0], A[0], bc) # real part
df.solve(a == b[1], A[1], bc) # imaginary part

# TODO: conversion to B field and output of fluxes

#%%   ------------------------ 
###   Plotting and comparisons
###   ------------------------
if (doplot):
    plotcurrn = np.loadtxt(plotcurrnfile)
    
    #%% plot mesh
    plt.subplot(1,2,1)
    plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=.1)
    plt.title('mesh (triplot)')
    plt.subplot(1,2,2)
    df.plot(mesh, linewidth=.1)
    plt.title('mesh (fenics)')
    
    #%%
    plt.figure()
    plt.subplot(1,2,1)
    df.plot(J[0])
    plt.title('Re J')
    plt.subplot(1,2,2)
    df.plot(J[1])
    plt.title('Im J')
    plt.tight_layout()
    
    #%%
    plt.figure()
    plt.subplot(1,2,1)
    p = df.plot(div(J[0])-n*Jphi[1])
    plt.colorbar(p)
    plt.title('Re div J')
    plt.subplot(1,2,2)
    p = df.plot(div(J[1])+n*Jphi[0])
    plt.colorbar(p)
    plt.title('Im div J')
    
    
    #%%
    
    Jphiplot = [np.zeros(len(tri)),np.zeros(len(tri))]
    Jphiplot[0][0:len(plotcurrn)] = plotcurrn[:,6]
    Jphiplot[1][0:len(plotcurrn)] = plotcurrn[:,7]
    
    plt.figure()
    plt.tripcolor(node[:, 0], node[:, 1], tri-1, Jphiplot[0])
    plt.colorbar()
    plt.title('Re Jphi (tripcolor)')
    
    plt.figure()
    plt.tripcolor(node[:, 0], node[:, 1], tri-1, Jphiplot[1])
    plt.colorbar()
    plt.title('Im Jphi (tripcolor)')
    
    plt.figure()
    p=df.plot(Jphi[0])
    plt.colorbar(p)
    plt.title('Re Jphi')
    plt.figure()
    p=df.plot(Jphi[1])
    plt.colorbar(p)
    plt.title('Im Jphi')
    
    #%%
    plt.figure()
    plt.subplot(1,2,1)
    p = df.plot(A[0])
    plt.colorbar(p)
    plt.title('Re A')
    plt.subplot(1,2,2)
    p = df.plot(A[1])
    plt.colorbar(p)
    plt.title('Im A')
    plt.tight_layout()
    
    plt.figure()
    plt.subplot(1,2,1)
    p = df.plot(curl(A[0]))
    plt.colorbar(p)
    plt.title('Re Bphi')
    plt.subplot(1,2,2)
    p = df.plot(curl(A[1]))
    plt.colorbar(p)
    plt.title('Im Bphi')
    plt.tight_layout()
    
    plt.show()
