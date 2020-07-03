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

c_cgs = 2.9979246e+10
cm_to_m = 1.0e-02
c1_statA_to_A = 1.0e+01
c1_statA_per_cm2_to_A_per_m2 = c1_statA_to_A / cm_to_m ** 2
statA_to_A = c1_statA_to_A / c_cgs
statA_per_cm2_to_A_per_m2 = c1_statA_per_cm2_to_A_per_m2 / c_cgs

test_dir = '/home/patrick/itp-temp/NEO-EQ/run'
work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
rad_resolution = 2048
m_min = 3
m_max = 9

magdif_r = dict()
magdif_jnpar = dict()
magdif_part_int = dict()
magdif_bndry = dict()
for m in np.arange(m_min, m_max + 1):
    data = np.loadtxt(path.join(test_dir, f"TCFP_Ipar_{m}/currn_par.dat"))
    magdif_r[m] = data[:, 0].copy()
    magdif_jnpar[m] = np.empty((rad_resolution), dtype=complex)
    magdif_jnpar[m].real = data[:, 1].copy()
    magdif_jnpar[m].imag = data[:, 2].copy()
    magdif_part_int[m] = np.empty((rad_resolution), dtype=complex)
    magdif_part_int[m].real = data[:, 5].copy()
    magdif_part_int[m].imag = data[:, 6].copy()
    magdif_bndry[m] = np.empty((rad_resolution), dtype=complex)
    magdif_bndry[m].real = data[:, 9].copy()
    magdif_bndry[m].imag = data[:, 10].copy()

kilca_r = dict()
kilca_jnpar = dict()
kilca_rres = dict()
kilca_d = dict()
kilca_jnpar_int = dict()
data = h5py.File(path.join(work_dir, 'TCFP_flre_hip.hdf5'), 'r')
for m in range(m_min, m_max + 1):
    loc = f"/output/postprocessor{m - m_min + 1}"
    kilca_r[m] = data[loc + '/r'][0, :].copy()
    kilca_jnpar[m] = np.empty(data[loc + '/Jpar'].shape[1:], dtype=complex)
    kilca_jnpar[m].real = data[loc + '/Jpar'][0, ...].copy()
    if data[loc + '/Jpar'].shape[0] == 2:
        kilca_jnpar[m].imag = data[loc + '/Jpar'][1, ...].copy()
    else:
        kilca_jnpar[m].imag = np.zeros(data[loc + '/Jpar'].shape[1:])
    kilca_jnpar_int[m] = interp.interp1d(kilca_r[m], kilca_jnpar[m],
                   kind='cubic', fill_value='extrapolate', assume_sorted=True)
    kilca_rres[m] = data[loc + '/rres'][0, 0].copy()
    kilca_d[m] = data[loc + '/d'][0, 0].copy()
kilca_hz = (data['/output/background/b0z'][0, :] /
            data['/output/background/b0'][0, :])
kilca_r[0] = data['/output/background/R'][0, :].copy()
kilca_hz_int = interp.interp1d(kilca_r[0], kilca_hz, kind='cubic',
                               fill_value='extrapolate', assume_sorted=True)
full_r_range = (0.0, np.amax(kilca_r[0]))

plt.figure(figsize=canvas)
ax = plt.gca()
for m in np.arange(m_min, m_max + 1):
    plt.plot(magdif_r[m], np.abs(magdif_jnpar[m]) * statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
ax.get_yaxis().set_major_formatter(scifmt)
plt.legend(loc='upper right')
plt.xlim(full_r_range)
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax_ins_1 = ax.inset_axes([3.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
ax_ins_1.plot(magdif_r[3], np.abs(magdif_jnpar[3]) * statA_per_cm2_to_A_per_m2,
              '-', lw=thin, c=c[0])
ax_ins_1.get_yaxis().set_major_formatter(copy.copy(scifmt))
ax_ins_1.set_xlim(kilca_rres[3] - 0.5, kilca_rres[3] + 0.5)
ax_ins_1.set_ylim(0.0, 0.8e+04)
ax.indicate_inset_zoom(ax_ins_1)
ax_ins_2 = ax.inset_axes([19.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
ax_ins_2.plot(magdif_r[9], np.abs(magdif_jnpar[9]) * statA_per_cm2_to_A_per_m2,
              '-', lw=thin, c=c[6])
ax_ins_2.get_yaxis().set_major_formatter(copy.copy(scifmt))
ax_ins_2.set_xlim(kilca_rres[9] - 0.5, kilca_rres[9] + 0.5)
ax_ins_2.set_ylim(0.0, 0.8e+04)
ax.indicate_inset_zoom(ax_ins_2)
plt.xlabel(r'$r$ / cm')
plt.ylabel(r'$|J_{m n}^{\parallel}|$ / A m\textsuperscript{-2}')
plt.title('Parallel current density of poloidal modes from magdif')
plt.savefig(path.join(work_dir, 'plot_Jnpar_magdif.png'), dpi=res)
plt.close()

plt.figure(figsize=canvas)
for m in np.arange(m_min, m_max + 1):
    plt.axvline(kilca_rres[m], lw=0.25 * thin, color='k')
    plt.axvline(kilca_rres[m] - 0.5 * kilca_d[m], lw=0.25 * thin, color='k',
                ls='--')
    plt.axvline(kilca_rres[m] + 0.5 * kilca_d[m], lw=0.25 * thin, color='k',
                ls='--')
    plt.plot(kilca_r[m], np.abs(kilca_jnpar[m]) * c1_statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend(loc='upper left')
plt.xlim(full_r_range)
plt.xlabel(r'$r$ / cm')
plt.ylabel(r'$|J_{m n}^{\parallel}|$ / A m\textsuperscript{-2}')
plt.title('Parallel current density of poloidal modes from KiLCA')
plt.savefig(path.join(work_dir, 'plot_Jnpar_kilca.png'), dpi=res)
plt.close()

magdif_Imnpar = dict()
kilca_Imnpar = dict()
for m in range(m_min, m_max + 1):
    # magdif
    magdif_min = np.searchsorted(magdif_r[m],
                                 kilca_rres[m] - kilca_d[m], 'left')
    magdif_max = np.searchsorted(magdif_r[m],
                                 kilca_rres[m] + kilca_d[m], 'right') - 1
    magdif_mid = np.searchsorted(magdif_r[m], kilca_rres[m])
    magdif_half = min(magdif_max - magdif_mid + 1, magdif_mid - magdif_min)
    magdif_width = np.empty(magdif_half)
    magdif_intJ = np.empty(magdif_half, dtype=complex)
    magdif_intB = np.empty(magdif_half, dtype=complex)
    magdif_bndryB = np.empty(magdif_half, dtype=complex)
    for w in range(0, magdif_half):
        magdif_lo = magdif_mid - 1 - w
        magdif_hi = magdif_mid + w
        magdif_width[w] = magdif_r[m][magdif_hi] - magdif_r[m][magdif_lo]
        magdif_intJ[w] = np.trapz(magdif_jnpar[m][magdif_lo:magdif_hi] *
                   magdif_r[m][magdif_lo:magdif_hi],
                   magdif_r[m][magdif_lo:magdif_hi]) * 2.0 * np.pi * statA_to_A
        magdif_intB[w] = np.trapz(magdif_part_int[m][magdif_lo:magdif_hi],
                   magdif_r[m][magdif_lo:magdif_hi]) * 0.5 * c_cgs * statA_to_A
        magdif_bndryB[w] = (magdif_bndry[m][magdif_hi] -
                     magdif_bndry[m][magdif_lo]) * 0.5 * c_cgs * statA_to_A
    magdif_lo = np.searchsorted(magdif_r[m],
                                kilca_rres[m] - 0.5 * kilca_d[m], 'left')
    magdif_hi = np.searchsorted(magdif_r[m],
                                 kilca_rres[m] + 0.5 * kilca_d[m], 'right') - 1
    magdif_Imnpar[m] = np.trapz(magdif_jnpar[m][magdif_lo:magdif_hi] *
                 magdif_r[m][magdif_lo:magdif_hi],
                 magdif_r[m][magdif_lo:magdif_hi]) * 2.0 * np.pi * statA_to_A
    # KiLCA
    interval = np.linspace(kilca_rres[m] - kilca_d[m],
                           kilca_rres[m] + kilca_d[m], rad_resolution)
    kilca_mid = rad_resolution // 2
    kilca_half = rad_resolution // 2 - 1
    kilca_width = np.empty(kilca_half)
    kilca_intJ = np.empty(kilca_half, dtype=complex)
    ### kilca_intB = np.empty(kilca_half, dtype=complex)
    ### kilca_bndryB = np.empty(kilca_half, dtype=complex)
    for w in range(0, kilca_half):
        kilca_lo = kilca_mid - 1 - w
        kilca_hi = kilca_mid + w
        kilca_width[w] = interval[kilca_hi] - interval[kilca_lo]
        kilca_intJ[w] = np.trapz(kilca_jnpar_int[m](interval[kilca_lo:kilca_hi]) *
                  interval[kilca_lo:kilca_hi] *
                  kilca_hz_int(interval[kilca_lo:kilca_hi]),
                  interval[kilca_lo:kilca_hi]) * 2.0 * np.pi * c1_statA_to_A
    interval = np.linspace(kilca_rres[m] - 0.5 * kilca_d[m],
                           kilca_rres[m] + 0.5 * kilca_d[m], rad_resolution // 2)
    kilca_Imnpar[m] = np.trapz(kilca_jnpar_int[m](interval) * interval *
                  kilca_hz_int(interval), interval) * 2.0 * np.pi * c1_statA_to_A
    plt.figure(figsize=canvas)
    plt.axvline(kilca_d[m], lw=0.25 * thin, color='k')
    plt.axhline(np.abs(kilca_Imnpar[m]), lw=0.25 * thin, color='k')
    plt.plot(kilca_width, np.abs(kilca_intJ), '-k', lw=thin, label='KiLCA, J')
    plt.plot(magdif_width, np.abs(magdif_intJ), '-b', lw=thin,
             label='magdif, J')
    plt.plot(magdif_width, np.abs(magdif_intB + magdif_bndryB), '-r', lw=thin,
             label='magdif, B')
    plt.plot(magdif_width, np.abs(magdif_intB), '--r', lw=thin,
             label='magdif, B, part.int.')
    plt.plot(magdif_width, np.abs(magdif_bndryB), ':r', lw=thin,
             label='magdif, B, bndry.')
    plt.gca().legend(loc='lower right')
    plt.xlabel(r'$d$ / cm')
    plt.ylabel(r'$\vert I_{m n}^{\parallel} \vert$ / A')
    plt.title(f"Parallel current depending on assumed layer width (m = {m})")
    plt.savefig(path.join(work_dir, f"cmp_Imnpar_{m}.png"), dpi=res)
    plt.close()
