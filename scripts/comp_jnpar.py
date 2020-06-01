#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:27:00 2020

@author: patrick
"""

import copy
import numpy as np
from scipy import interpolate as interp
import matplotlib
from matplotlib import rcParams, use, pyplot as plt
import h5py
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

cm2_to_m2 = 1.0e-04
statA_per_cm2_to_A_per_m2 = 1.0 / 2.9979246e+05
c1_statA_per_cm2_to_A_per_m2 = 1.0e+05
c1_statA_to_A = 1.0e+01

test_dir = '/home/patrick/itp-temp/NEO-EQ/run'
work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'

nlabel = 3000
m_min = 3
m_max = 9

# values for high pressure test case from Philipp's handwritten notes
r_res = np.array([16.8, 27.5, 35.1, 41.2, 46.6, 51.4, 55.8])
dr_res = np.array([0.5, 0.6, 0.9, 1.3, 1.8, 2.3, 2.8])
I_res = np.array([4., 17., 38., 69., 116., 173., 269.])

magdif_r = np.empty((nlabel, m_max - m_min + 1))
magdif_jnpar = np.empty((nlabel, m_max - m_min + 1), dtype=complex)
for m in np.arange(m_min, m_max + 1):
    k = m - m_min
    data = np.loadtxt(path.join(test_dir,
                                'TCFP_Ipar_{}/currmn_par.dat'.format(m)))
    magdif_r[:, k] = data[:, 0].copy()
    # h_z = data[:, 1].copy()
    magdif_jnpar[:, k].real = data[:, 2] * statA_per_cm2_to_A_per_m2
    magdif_jnpar[:, k].imag = data[:, 3] * statA_per_cm2_to_A_per_m2

plt.figure(figsize=canvas)
ax = plt.gca()
for m in np.arange(m_min, m_max + 1):
    k = m - m_min
    plt.plot(magdif_r[:, k], np.abs(magdif_jnpar[:, k]), '-', lw=thin,
             label='m = {}'.format(m))
ax.get_yaxis().set_major_formatter(scifmt)
plt.legend(loc='upper right')
full_r_range = plt.xlim()
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax_ins_1 = ax.inset_axes([3.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
ax_ins_1.plot(magdif_r[:, 0], np.abs(magdif_jnpar[:, 0]), '-', lw=thin, c=c[0])
ax_ins_1.get_yaxis().set_major_formatter(copy.copy(scifmt))
ax_ins_1.set_xlim(r_res[0] - 0.5, r_res[0] + 0.5)
ax_ins_1.set_ylim(0.0, 0.8e+04)
ax.indicate_inset_zoom(ax_ins_1)
ax_ins_2 = ax.inset_axes([19.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
ax_ins_2.plot(magdif_r[:, 6], np.abs(magdif_jnpar[:, 6]), '-', lw=thin, c=c[6])
ax_ins_2.get_yaxis().set_major_formatter(copy.copy(scifmt))
ax_ins_2.set_xlim(r_res[6] - 0.5, r_res[6] + 0.5)
ax_ins_2.set_ylim(0.0, 0.8e+04)
ax.indicate_inset_zoom(ax_ins_2)
plt.xlabel(r'$r$ / cm')
plt.ylabel(r'$|J_{m n}^{\parallel}|$ / A m\textsuperscript{-2}')
plt.title('Parallel current density of poloidal modes from magdif')
plt.savefig(path.join(work_dir, 'plot_Jnpar_magdif.png'), dpi=res)
plt.close()

kilca_r = dict()
kilca_jnpar = dict()
kilca_rres = dict()
kilca_d = dict()
kilca_Inpar = dict()
data = h5py.File(path.join(work_dir, 'TCFP_flre_hip.hdf5'), 'r')
for m in range(m_min, m_max + 1):
    k = m - m_min + 1
    loc = '/output/postprocessor{}'.format(k)
    kilca_r[m] = data[loc + '/r'][0, :].copy()
    kilca_jnpar[m] = np.empty(data[loc + '/Jpar'].shape[1:], dtype=complex)
    kilca_jnpar[m].real = data[loc + '/Jpar'][0, ...] * \
        c1_statA_per_cm2_to_A_per_m2
    if data[loc + '/Jpar'].shape[0] == 2:
        kilca_jnpar[m].imag = data[loc + '/Jpar'][1, ...] * \
            c1_statA_per_cm2_to_A_per_m2
    else:
        kilca_jnpar[m].imag = np.zeros(data[loc + '/Jpar'].shape[1:])
    kilca_rres[m] = data[loc + '/rres'][0, 0].copy()
    kilca_d[m] = data[loc + '/d'][0, 0].copy()
    kilca_Inpar[m] = data[loc + '/Ipar'][0, 0] * c1_statA_to_A * 0.5 / np.pi

plt.figure(figsize=canvas)
for m in np.arange(m_min, m_max + 1):
    plt.axvline(kilca_rres[m], lw=0.5*thin, color='k')
    plt.axvline(kilca_rres[m] - kilca_d[m], lw=0.5*thin, color='k', ls='--')
    plt.axvline(kilca_rres[m] + kilca_d[m], lw=0.5*thin, color='k', ls='--')
    plt.plot(kilca_r[m], np.abs(kilca_jnpar[m]), '-', lw=thin,
             label='m = {}'.format(m))
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend(loc='upper left')
plt.xlim(full_r_range)
plt.xlabel(r'$r$ / cm')
plt.ylabel(r'$|J_{m n}^{\parallel}|$ / A m\textsuperscript{-2}')
plt.title('Parallel current density of poloidal modes from KiLCA')
plt.savefig(path.join(work_dir, 'plot_Jnpar_kilca.png'), dpi=res)
plt.close()

I_res_magdif = np.zeros((m_max - m_min + 1), dtype=complex)
I_res_kilca = np.zeros((m_max - m_min + 1), dtype=complex)
for m in range(m_min, m_max + 1):
    k = m - m_min
    interval = np.linspace(kilca_rres[m] - kilca_d[m],
                           kilca_rres[m] + kilca_d[m], 512)
    sampling_magdif = interp.interp1d(magdif_r[:, k], magdif_r[:, k] * 
                                      magdif_jnpar[:, k], kind='cubic',
                                      fill_value='extrapolate',
                                      assume_sorted=True)
    I_res_magdif[k] = np.trapz(sampling_magdif(interval), interval) * cm2_to_m2
    sampling_kilca = interp.interp1d(kilca_r[m], kilca_r[m] * kilca_jnpar[m],
                                     kind='cubic', fill_value='extrapolate',
                                     assume_sorted=True)
    I_res_kilca[k] = np.trapz(sampling_kilca(interval), interval) * cm2_to_m2
    print(f"I_res_magdif({m}) = {np.abs(I_res_magdif[k]):.15f}, "
          f"I_res_kilca({m}) = {np.abs(I_res_kilca[k]):.15f}, "
          f"kilca_Ipar({m}) = {kilca_Inpar[m]:.15f}")
