#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:19:58 2020

@author: patrick
"""

import numpy as np
import matplotlib
from matplotlib import rcParams, use, pyplot as plt
import h5py
import netCDF4
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

work_dir = '/home/patrick/git/NEO-EQ/run/30835_3200_ed6'
n = 2  ### read from config

gpec = netCDF4.Dataset(path.join(work_dir, 'gpec_profile_output_n2.nc'), 'r')
gpec_q = np.array(gpec.variables['q_rational'][:])
gpec_m_min = int(np.rint(n * np.amin(np.abs(gpec_q))))
gpec_m_max = int(np.rint(n * np.amax(np.abs(gpec_q))))
gpec_I_res = dict()
gpec_w_isl = dict()
for m in range(gpec_m_min, gpec_m_max + 1):
    k = m - gpec_m_min
    gpec_I_res[m] = np.complex(gpec.variables['I_res'][0, k],
                               gpec.variables['I_res'][1, k])
    gpec_w_isl[m] = gpec.variables['w_isl'][k]

magdif = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')
magdif_m_min = magdif['/mesh/m_res_min'][()]
magdif_m_max = magdif['/mesh/m_res_max'][()] - 1  ### mode 11 is too far out
magdif_psi = magdif['/cache/fs/psi'][()]
magdif_rad = magdif['/cache/fs/rad'][()]
### fix mesh%psi_res in Fortran!
magdif_rad_res = magdif['/mesh/rad_norm_res'][()] * magdif_rad[-1]
magdif_nrad = 2048  ### read from config

mephit_psi = dict()
mephit_r = dict()
mephit_jnpar = dict()
mephit_part_int = dict()
mephit_bndry = dict()
mephit_Ichar = dict()
mephit_Delta = dict()
for m in np.arange(magdif_m_min, magdif_m_max + 1):
    data = np.loadtxt(path.join(work_dir, f"currn_par_{m}.dat"))
    mephit_psi[m] = data[:, 0].copy()
    mephit_r[m] = data[:, 1].copy()
    mephit_jnpar[m] = np.empty((magdif_nrad), dtype=complex)
    mephit_jnpar[m].real = data[:, 2].copy()
    mephit_jnpar[m].imag = data[:, 3].copy()
    mephit_part_int[m] = np.empty((magdif_nrad), dtype=complex)
    mephit_part_int[m].real = data[:, 6].copy()
    mephit_part_int[m].imag = data[:, 7].copy()
    mephit_bndry[m] = np.empty((magdif_nrad), dtype=complex)
    mephit_bndry[m].real = data[:, 10].copy()
    mephit_bndry[m].imag = data[:, 11].copy()
    mephit_Ichar[m] = data[:, 14].copy()
    mephit_Delta[m] = np.empty((magdif_nrad), dtype=complex)
    mephit_Delta[m].real = data[:, 15].copy()
    mephit_Delta[m].imag = data[:, 16].copy()

mephit_Imnpar = dict()
for m in range(magdif_m_min, min(gpec_m_max, magdif_m_max) + 1):
    k = m - magdif_m_min
    # MEPHIT
    mephit_mid = np.searchsorted(mephit_r[m], magdif_rad_res[k])
    mephit_half = min(magdif_nrad - mephit_mid, mephit_mid)
    mephit_width = np.empty(mephit_half)
    mephit_intJ = np.empty(mephit_half, dtype=complex)
    mephit_intB = np.empty(mephit_half, dtype=complex)
    mephit_bndryB = np.empty(mephit_half, dtype=complex)
    mephit_GPEC = np.empty(mephit_half, dtype=complex)
    for w in range(0, mephit_half):
        mephit_lo = mephit_mid - 1 - w
        mephit_hi = mephit_mid + w
        mephit_width[w] = mephit_r[m][mephit_hi] - mephit_r[m][mephit_lo]
        mephit_intJ[w] = np.trapz(mephit_jnpar[m][mephit_lo:mephit_hi],
                   mephit_psi[m][mephit_lo:mephit_hi]) * statA_to_A
        mephit_intB[w] = np.trapz(mephit_part_int[m][mephit_lo:mephit_hi],
                   mephit_psi[m][mephit_lo:mephit_hi]) * 2.0 * np.pi * statA_to_A
        mephit_bndryB[w] = (mephit_bndry[m][mephit_hi] -
                     mephit_bndry[m][mephit_lo]) * 2.0 * np.pi * statA_to_A
        mephit_GPEC[w] = -1j / m * mephit_Ichar[m][mephit_mid] * (
            mephit_Delta[m][mephit_hi] - mephit_Delta[m][mephit_lo]) * statA_to_A
    mephit_lo = mephit_mid - 1 - mephit_half // 2
    mephit_hi = mephit_mid + mephit_half // 2
    mephit_Imnpar[m] = np.trapz(mephit_jnpar[m][mephit_lo:mephit_hi],
                   mephit_psi[m][mephit_lo:mephit_hi]) * statA_to_A
    # plots
    plt.figure(figsize=canvas)
    plt.axhline(np.abs(gpec_I_res[m]), lw=0.25 * thin, color='k', label='GPEC')
    plt.plot(mephit_width, 2.0 * np.pi * np.abs(mephit_intJ), '--b', lw=thin,
             label='J, symfluxcoord.')
# =============================================================================
#     plt.plot(mephit_width, np.abs(mephit_intB + mephit_bndryB), '-.c', lw=thin,
#              label='B, symfluxcoord.')
#     plt.plot(mephit_width, np.abs(mephit_intB), '--c', lw=thin,
#              label='B, symfluxcoord., part.int.')
# =============================================================================
    plt.plot(mephit_width, np.abs(mephit_bndryB), ':c', lw=thin,
             label='B, symfluxcoord, bndry.')
# =============================================================================
#     plt.plot(mephit_width, np.abs(mephit_GPEC), '--g', lw=thin,
#              label=r'$\Delta_{mn}$')
# =============================================================================
    plt.gca().legend(loc='lower left')
    plt.xlabel(r'$d$ / cm')
    plt.ylabel(r'$\vert I_{m n}^{\parallel} \vert$ / A')
    plt.title(f"Parallel current depending on assumed layer width (m = {m})")
    plt.savefig(path.join(work_dir, f"cmp_Imnpar_{m}.png"), dpi=res)
    plt.close()
