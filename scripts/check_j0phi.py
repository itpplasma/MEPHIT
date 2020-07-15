#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:34:35 2020

@author: patrick
"""

import numpy as np
import matplotlib
from matplotlib import rcParams, use, pyplot as plt
from os import path

rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'
use('Agg')
scifmt = matplotlib.ticker.ScalarFormatter()
scifmt.set_powerlimits((-3, 3))
canvas = (6.6, 3.6)
res = 300
thin = 0.5

dyn_per_cm3_to_N_per_m3 = 10.0
c_cgs = 2.9979246e+10
statA_per_cm2_to_A_per_m2 = 1.0e+05 / c_cgs

# work_dir = '/home/patrick/git/NEO-EQ/run/ED6_30835_3200'
# work_dir = '/home/patrick/itp-temp/NEO-EQ/run/30835_3200'
work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
data_gs = np.loadtxt(path.join(work_dir, 'j0_gs.dat'))
data_amp = np.loadtxt(path.join(work_dir, 'j0_amp.dat'))
data_prof = np.loadtxt(path.join(work_dir, 'cmp_prof.dat'))

nflux = 91
kt_max = np.empty((nflux + 1), dtype=int)
kt_max[0] = 300
kt_max[1:] = 600
kt_low = np.empty((nflux + 1), dtype=int)
kt_low[0] = 0
kt_low[1:] = np.cumsum(kt_max[:-1]) + kt_low[0]

kf = 50

theta = np.linspace(0.0, 360.0, kt_max[kf])
plt.figure(figsize=canvas)
plt.plot(theta, data_prof[kt_low[kf]:kt_low[kf]+kt_max[kf], 0] *
         dyn_per_cm3_to_N_per_m3, '-k', lw=thin, label=r'$\nabla p_{0}$')
plt.plot(theta, data_prof[kt_low[kf]:kt_low[kf]+kt_max[kf], 1] *
         dyn_per_cm3_to_N_per_m3, '--r', lw=thin, label='Ampere')
plt.plot(theta, data_prof[kt_low[kf]:kt_low[kf]+kt_max[kf], 2] *
         dyn_per_cm3_to_N_per_m3, '--b', lw=thin, label='Grad-Shafranov')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend()
plt.xlabel(r'$\theta$ / \textdegree')
plt.ylabel(r'$f$ / N m\textsuperscript{-3}')
plt.title('iMHD force balance')
plt.savefig(path.join(work_dir, "check_iMHD.png"), dpi=res)
plt.close()

plt.figure(figsize=canvas)
plt.plot(theta, data_gs[kt_low[kf]:kt_low[kf]+kt_max[kf], 0] *
         statA_per_cm2_to_A_per_m2, '-k', lw=thin, label='Grad-Shafranov')
plt.plot(theta, data_amp[kt_low[kf]:kt_low[kf]+kt_max[kf], 0] *
         statA_per_cm2_to_A_per_m2, '--r', lw=thin, label='Ampere')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend()
plt.xlabel(r'$\theta$ / \textdegree')
plt.ylabel(r'$J_{0}^{R}$ / A m\textsuperscript{-2}')
plt.title('Comparison of equilibrium current')
plt.savefig(path.join(work_dir, "check_j0R.png"), dpi=res)
plt.close()

plt.figure(figsize=canvas)
plt.plot(theta, data_gs[kt_low[kf]:kt_low[kf]+kt_max[kf], 1] *
         statA_per_cm2_to_A_per_m2, '-k', lw=thin, label='Grad-Shafranov')
plt.plot(theta, data_amp[kt_low[kf]:kt_low[kf]+kt_max[kf], 1] *
         statA_per_cm2_to_A_per_m2, '--r', lw=thin, label='Ampere')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend()
plt.xlabel(r'$\theta$ / \textdegree')
plt.ylabel(r'$R J_{0}^{\varphi}$ / A m\textsuperscript{-2}')
plt.title('Comparison of equilibrium current')
plt.savefig(path.join(work_dir, "check_j0phi.png"), dpi=res)
plt.close()

plt.figure(figsize=canvas)
plt.plot(theta, data_gs[kt_low[kf]:kt_low[kf]+kt_max[kf], 2] *
         statA_per_cm2_to_A_per_m2, '-k', lw=thin, label='Grad-Shafranov')
plt.plot(theta, data_amp[kt_low[kf]:kt_low[kf]+kt_max[kf], 2] *
         statA_per_cm2_to_A_per_m2, '--r', lw=thin, label='Ampere')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend()
plt.xlabel(r'$\theta$ / \textdegree')
plt.ylabel(r'$J_{0}^{Z}$ / A m\textsuperscript{-2}')
plt.title('Comparison of equilibrium current')
plt.savefig(path.join(work_dir, "check_j0Z.png"), dpi=res)
plt.close()
