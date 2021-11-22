#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 12:22:12 2021

@author: patrick
"""

import h5py
from os import path
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from magdifplot import run_dir
import numpy as np


def circumcircle(x, y):
    D = 2.0 * (x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]))
    Cx = ((x[0] ** 2 + y[0] ** 2) * (y[1] - y[2]) + (x[1] ** 2 + y[1] ** 2) * (y[2] - y[0]) + (
                x[2] ** 2 + y[2] ** 2) * (y[0] - y[1])) / D
    Cy = ((x[0] ** 2 + y[0] ** 2) * (x[2] - x[1]) + (x[1] ** 2 + y[1] ** 2) * (x[0] - x[2]) + (
                x[2] ** 2 + y[2] ** 2) * (x[1] - x[0])) / D
    r = np.hypot(x[0] - Cx, y[0] - Cy)
    return (Cx, Cy), r


def wspace(interval, stretch):
    return [0.5 * ((1.0 - stretch) * interval[1] + (1.0 + stretch) * interval[0]),
            0.5 * ((1.0 + stretch) * interval[1] + (1.0 - stretch) * interval[0])]


def remove_pop(lst):
    for elem in lst:
        elem.remove()
    lst.clear()


work_dir = run_dir + '/33353_2325'
testcase = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')
ktri_lo = testcase['/mesh/kt_low'][5] + 11
ktri_hi = ktri_lo + 8
tris = testcase['/mesh/tri_node'][ktri_lo:ktri_hi, :] - 1
tris = tris[:, [0, 1, 2, 0]]
kpoi = np.sort(np.unique(np.ravel(tris)))
R = testcase['/mesh/node_R'][()] - testcase['/mesh/R_O'][()]
Z = testcase['/mesh/node_Z'][()] - testcase['/mesh/Z_O'][()]
range_R = wspace([np.min(R[kpoi]), np.max(R[kpoi])], 1.2)
range_Z = wspace([np.min(Z[kpoi]), np.max(Z[kpoi])], 1.2)
intermediate = []
fig = Figure(figsize=(7, 3.5))
ax = fig.subplots()
pdf = PdfPages(path.join(work_dir, 'poloidal_resolution.pdf'))
ax.set_aspect('equal')
ax.set_xlabel(r'$(R - R_{O})$ / \si{\centi\meter}')
ax.set_ylabel(r'$(Z - Z_{O})$ / \si{\centi\meter}')
ax.set_xlim(range_R)
ax.set_ylim(range_Z)
ax.plot(R[kpoi], Z[kpoi], 'xk')
ax.plot(R[kpoi[:5]], Z[kpoi[:5]], '-k', lw=0.375)
ax.plot(R[kpoi[5:]], Z[kpoi[5:]], '-k', lw=0.375)
ax.plot(R[np.ravel(tris[0, :])], Z[np.ravel(tris[0, :])], '-k')
ax.plot(R[np.ravel(tris[1, :])], Z[np.ravel(tris[1, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
intermediate.append(ax.plot(R[kpoi[[1, 7, 2]]], Z[kpoi[[1, 7, 2]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[kpoi[[1, 6, 7]]], Z[kpoi[[1, 6, 7]]])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='r', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
intermediate.append(ax.plot(R[kpoi[[6, 2, 7]]], Z[kpoi[[6, 2, 7]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[np.ravel(tris[2, :])], Z[np.ravel(tris[2, :])])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='g', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
ax.plot(R[np.ravel(tris[2, :])], Z[np.ravel(tris[2, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
intermediate.append(ax.plot(R[kpoi[[2, 7, 3]]], Z[kpoi[[2, 7, 3]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[np.ravel(tris[3, :])], Z[np.ravel(tris[3, :])])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='g', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
ax.plot(R[np.ravel(tris[3, :])], Z[np.ravel(tris[3, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
intermediate.append(ax.plot(R[kpoi[[2, 8, 3]]], Z[kpoi[[2, 8, 3]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[np.ravel(tris[4, :])], Z[np.ravel(tris[4, :])])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='g', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
ax.plot(R[np.ravel(tris[4, :])], Z[np.ravel(tris[4, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
intermediate.append(ax.plot(R[kpoi[[2, 9, 3]]], Z[kpoi[[2, 9, 3]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[kpoi[[2, 8, 9]]], Z[kpoi[[2, 8, 9]]])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='r', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
intermediate.append(ax.plot(R[kpoi[[8, 3, 9]]], Z[kpoi[[8, 3, 9]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[np.ravel(tris[5, :])], Z[np.ravel(tris[5, :])])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='g', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
ax.plot(R[np.ravel(tris[5, :])], Z[np.ravel(tris[5, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
intermediate.append(ax.plot(R[kpoi[[3, 9, 4]]], Z[kpoi[[3, 9, 4]]], ':k', lw=0.25)[0])
C, r = circumcircle(R[np.ravel(tris[6, :])], Z[np.ravel(tris[6, :])])
intermediate.append(ax.add_patch(Circle(C, r, fill=False, color='g', lw=0.5)))
fig.savefig(pdf, format='pdf', dpi=600)
remove_pop(intermediate)
ax.plot(R[np.ravel(tris[6, :])], Z[np.ravel(tris[6, :])], '-k')
fig.savefig(pdf, format='pdf', dpi=600)
pdf.close()
