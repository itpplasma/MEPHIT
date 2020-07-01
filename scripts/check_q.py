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

work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
data_step = np.loadtxt(path.join(work_dir, 'check_q_step.dat'))
data_cont = np.loadtxt(path.join(work_dir, 'check_q_cont.dat'))

plt.figure(figsize=canvas)
for m in range(3, 12):
    plt.axhline(m / 2000, lw=0.25 * thin, color='k')
plt.plot(data_cont[:, 0], data_cont[:, 2], '-k', lw=thin, label='gfile')
plt.plot(data_cont[:, 0], data_cont[:, 1], '--r', lw=thin, label='field line')
plt.step(data_step[:, 0], data_step[:, 1], where='mid', lw=thin, label='triangle grid')
# plt.xlim([55.6, 56.2])
# plt.ylim([0.00448, 0.00452])
plt.gca().legend(loc='lower right')
plt.xlabel(r'$d$ / cm')
plt.ylabel(r'$q$')
plt.title('Comparison of safety factor approximations')
plt.savefig(path.join(work_dir, "check_q.png"), dpi=res)
plt.close()
