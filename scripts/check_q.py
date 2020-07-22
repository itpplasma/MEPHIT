#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:48:52 2020

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

work_dir = '/home/patrick/git/NEO-EQ/run/ED6_30835_3200'
# work_dir = '/home/patrick/git/NEO-EQ/run/EQH_30835_3200'
# work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'

data_step = np.loadtxt(path.join(work_dir, 'check_q_step.dat'))
data_cont = np.loadtxt(path.join(work_dir, 'check_q_cont.dat'))

plt.figure(figsize=canvas)
plt.plot(data_cont[:, 0], data_cont[:, 3], '-k', lw=thin, label='gfile')
plt.plot(data_cont[:, 0], data_cont[:, 2], '--r', lw=thin, label='field line')
plt.step(data_step[:, 0], data_step[:, 2], where='mid', lw=thin,
         label='triangle grid')
plt.gca().yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
plt.gca().grid(which='major', axis='y', lw=0.5*thin)
### TCFP
# plt.xlim([55.6, 56.2])
# plt.ylim([0.00448, 0.00452])
### ED6
# plt.xlim([74.0, 76.5])
# plt.ylim([-3.05, -2.95])
plt.gca().legend()
plt.xlabel(r'$\hat{\psi}$')
plt.ylabel(r'$q$')
plt.title('Comparison of safety factor approximations')
plt.savefig(path.join(work_dir, "check_q.png"), dpi=res)
plt.close()
