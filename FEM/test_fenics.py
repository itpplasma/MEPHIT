#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 08:59:56 2018

@author: calbert
"""

import numpy as np
import dolfin as df
from dolfin import inner, grad, div, curl, dx, ds, solve

# Create mesh and define function space
mesh = df.UnitSquareMesh(32, 32)
V = df.FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < df.DOLFIN_EPS or x[0] > 1.0 - df.DOLFIN_EPS

# Define boundary condition
u0 = df.Constant(0.0)
bc = df.DirichletBC(V, u0, boundary)

# Define variational problem
u = df.TrialFunction(V)
v = df.TestFunction(V)
f = df.Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = df.Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = df.Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = df.File("poisson.pvd")
file << u

# Plot solution
import matplotlib.pyplot as plt
df.plot(u)
plt.show()

#%%
mesh = df.UnitSquareMesh(10, 10)
P1 = df.FunctionSpace(mesh, "Lagrange", 1)
Hcurl = df.FunctionSpace(mesh, "N1curl", 1)
Hdiv = df.FunctionSpace(mesh, "RT", 1)
f  = df.Expression("x[0]", degree=2)
fh = df.interpolate(f, P1)
B  = df.Expression(('1','0'), degree=2)
Bh = df.interpolate(B, Hdiv)

df.plot(Bh)

#%%
#
# Summary: fenics returns Nedelec and Raviart-Thomas edge DOFs
# opposite to vertex that is numbered in triangle
#
# Edge orientation is given uniquely pointing from nodes with
# lower index to the one with higher index. The same is done
# for triangles, that are ordered with increasing node indexes.
# Thus, triangles can also be numbered clockwise, but the
# edge orientation matches always between two adjacent ones.
#
#
#
# see also: https://arxiv.org/pdf/1205.3085.pdf
plt.figure()
df.plot(mesh)
Hdivmap = Hdiv.dofmap() # returns a map that yields 3 DOFs for each triangle

kcell = 1
nodeind = mesh.cells()[kcell]
nodecoords = mesh.coordinates()[nodeind]
cellorient = np.sign(np.linalg.det(np.concatenate(([[1,1,1]],nodecoords.T))))
print('Cell {}'.format(kcell))
print('Indexes: {}'.format(nodeind))
print('Coordinates:\n {}'.format(nodecoords))
print('Cell orientation: {}'.format(cellorient))
inddiff = np.roll(nodeind,1)+(-1)*np.roll(nodeind,-1) # (-1) to convert uint to int
orient = np.sign(inddiff)
print('Edge orientations: {}'.format(orient))
dofs = Hdivmap.cell_dofs(kcell) # DOFs in cell k
Bflux = cellorient*orient*Bh.vector()[dofs] # outward fluxes
print('Bflux = {}'.format(Bflux))

nodes = mesh.cells()[0]

#%%

def boundary(_, on_boundary):
    return on_boundary

J = df.Function(Hdiv)
u = df.TrialFunction(Hcurl)
v = df.TestFunction(Hcurl)

A_D = df.Constant((0, 1))
bc = df.DirichletBC(Hcurl, A_D, boundary)
f = df.Constant((0, 0))

a = (inner(curl(u), curl(v)) + .001*inner(u,v)) * dx
b = inner(f, v) * dx

A = df.Function(Hcurl)
df.solve(a == b, A, bc)

plt.figure()
df.plot(A)

plt.show()
# %%
mesh = df.UnitSquareMesh(10, 10)

Hdiv = df.FunctionSpace(mesh, 'RT', 1)
Hcurl = df.FunctionSpace(mesh, 'N1curl', 1)

u = df.TrialFunction(Hcurl)
v = df.TestFunction(Hcurl)

J = df.Expression(('x[1]', '-x[0]'), degree=1)


a = (inner(curl(u), curl(v)) + inner(u,v)) * dx
b = inner(J, v) * dx

A = df.Function(Hcurl)

df.solve(a == b, A)

plt.figure()
df.plot(A)

plt.show()

#%%

