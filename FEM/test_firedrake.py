#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 08:59:56 2018

@author: calbert
"""

import firedrake as fd
from firedrake import as_vector, inner, curl, dx
import matplotlib.pyplot as plt

mesh = fd.UnitSquareMesh(10, 10)

P0 = fd.FunctionSpace(mesh, 'DG', 0)
P1 = fd.FunctionSpace(mesh, 'Lagrange', 1)
Hdiv = fd.FunctionSpace(mesh, 'RT', 1)
Hcurl = fd.FunctionSpace(mesh, 'N1curl', 1)

u = fd.TrialFunction(Hcurl)
v = fd.TestFunction(Hcurl)

r, z = fd.SpatialCoordinate(mesh)
J = fd.Function(Hdiv).project(as_vector((z, -r)))

fd.plot(J)

a = (r*inner(curl(u), curl(v)) + 1.0/r*inner(u,v)) * dx
b = inner(J, v) * dx

A = fd.Function(Hcurl)

fd.solve(a == b, A)

fd.plot(A)

plt.show()

#%%