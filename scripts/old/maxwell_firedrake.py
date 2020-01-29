#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:01:58 2018

@author: calbert
"""

import numpy as np
import matplotlib.pyplot as plt
import firedrake as fd
from firedrake import div, curl, inner, dx
from pyop2.mpi import COMM_WORLD

from readmsh import readmsh

fd.parameters['reorder_meshes'] = False

doplot = True

c = 29979245800.0

n = 2.0

#infile = sys.argv[1]
#infile = '../PRELOAD/inputformaxwell.msh'
infile = '../PRELOAD/inputformaxwell_ext.msh'
infile2 = '../PRELOAD/inputformaxwell_ext.exo'
currnfile = '../MHD/currn.dat'
plotcurrnfile = '../MHD/plot_currn.dat'
#currnfile = '../MHD/Bn_nonres.dat'
#plotcurrnfile = '../MHD/plot_Bn_nonres.dat'

#currnfile = '../PRELOAD/VACFIELD/Bn_flux.dat'

outfile = '../MHD/Bn_fenics.dat'

currn = np.loadtxt(currnfile)


[node, tri, edge, nlab, tlab, elab] = readmsh(infile)
mesh0 = fd.mesh._from_cell_list(2, tri-1, node, COMM_WORLD)
mesh = fd.Mesh(mesh0)
#mesh = fd.Mesh(infile2)
p = fd.plot(mesh)

# Function spaces
P0 = fd.FunctionSpace(mesh, 'DG', 0)
P1 = fd.FunctionSpace(mesh, 'Lagrange', 1)
Hdiv = fd.FunctionSpace(mesh, 'RT', 1)
Hcurl = fd.FunctionSpace(mesh, 'N1curl', 1)
Hdivmap = Hdiv.cell_node_map()

# Map for re-ordering of indexes and orientations
map2facing = np.array([2,0,1]) # index of edge facing vertex
map2lex = np.argsort(tri[:,map2facing],1) # map to lexical order in each triangle cell
assert np.max(P1.cell_node_list+1 - np.sort(tri)) == 0, "triangle sort order"
#%%

#%%
# Read data
ncells = len(tri)
ndofs = Hdiv.dof_count

Jvec = [np.zeros(ndofs), np.zeros(ndofs)]
J = [fd.Function(Hdiv), fd.Function(Hdiv)]
Jphi = [fd.Function(P0), fd.Function(P0)]
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

#%% TODO: output bflux
Bflux = np.zeros((len(currn),8))
Avec = [A[0].vector()[:], A[1].vector()[:]]
for kcell in range(len(currn)):
    nodeind = cells[kcell]
    nodecoords = nodes[nodeind]
    cellarea = .5*np.linalg.det(np.concatenate(([[1,1,1]],nodecoords.T)))
    cellorient = np.sign(cellarea)
    orient = np.sign([nodeind[2]+(-1)*nodeind[1], 
                      nodeind[0]+(-1)*nodeind[2], 
                      nodeind[1]+(-1)*nodeind[0]]) # (-1) for signed int
    dofs = Hdivmap.cell_dofs(kcell) # DOFs in cell k
    Bflux[kcell,2*map2lex[kcell]] = -n * Avec[1][dofs]*cellorient*orient
    Bflux[kcell,2*map2lex[kcell]+1] = n * Avec[0][dofs]*cellorient*orient
    Bflux[kcell,6] = -sum(Avec[0][dofs]*cellorient*orient) # Re R B^phi
    Bflux[kcell,7] = -sum(Avec[1][dofs]*cellorient*orient) # Im R B^phi

np.savetxt(outfile, Bflux)

#%%   ------------------------ 
###   Plotting and comparisons
###   ------------------------
if (doplot):
    plotcurrn = np.loadtxt(plotcurrnfile)
    
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
    p = df.plot(div(J[0]))
    plt.colorbar(p)
    plt.title('Re div J')
    plt.subplot(1,2,2)
    p = df.plot(div(J[1]))
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
    
    #%% Plot resulting magnetic field
    B = [df.Function(Hdiv), df.Function(Hdiv)]
    Bphi = [df.Function(P0), df.Function(P0)]
    B[0].vector()[:] = -n * Avec[1]
    B[1].vector()[:] = n * Avec[0]
    Bphivec = [np.zeros(ncells), np.zeros(ncells)]
    Bphivec[0][0:len(Bflux)] = Bflux[:,6]/acell[0:len(Bflux)]
    Bphivec[1][0:len(Bflux)] = Bflux[:,7]/acell[0:len(Bflux)]
    Bphi[0].vector()[:] = Bphivec[0]
    Bphi[1].vector()[:] = Bphivec[1]
    
    plt.figure()
    df.plot(B[0])
    plt.title('Re B')
    plt.figure()
    df.plot(B[1])
    plt.title('Im B')
    plt.figure()
    df.plot(Bphi[0])
    plt.title('Re Bphi')
    plt.figure()
    df.plot(Bphi[1])
    plt.title('Im Bphi')
    
    plt.show()
