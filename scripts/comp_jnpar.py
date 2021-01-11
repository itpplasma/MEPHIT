#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:27:00 2020

@author: patrick
"""

from magdifplot import scifmt, parcurr, statA_per_cm2_to_A_per_m2, \
    c1_statA_per_cm2_to_A_per_m2
import numpy as np
import h5py
from matplotlib import pyplot as plt
from os import path

canvas = (6.6, 3.6)
res = 300
thin = 0.5

test_dir = '/home/patrick/itp-temp/git/NEO-EQ/run'
work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
nrad = 2048  ### read from config
m_min = 3
m_max = 9

magdif = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')
rad_max = magdif['/cache/fs/rad'][-1]
rad_res = {}
for k, m in enumerate(range(magdif['/mesh/rad_norm_res'].attrs['lbounds'][0],
                            magdif['/mesh/rad_norm_res'].attrs['ubounds'][0])):
    rad_res[m] = magdif['/mesh/rad_norm_res'][k] * rad_max
cylcoor = parcurr()
cylcoor.process_magdif(lambda m: path.join(test_dir, f"TCFP_Ipar_{m}/currn_par.dat"),
                       range(m_min, m_max + 1), nrad, rad_res, symfluxcoord=False)
fluxcoor = parcurr()
fluxcoor.process_magdif(lambda m: path.join(test_dir, f"TCFP_Ipar_{m}/currn_par_{m}.dat"),
                        range(m_min, m_max + 1), nrad, rad_res, symfluxcoord=True)
kilca_hic = parcurr()
kilca_hic.process_KiLCA(path.join(work_dir, 'TCFP_flre_hic.hdf5'), nrad)
kilca_mec = parcurr()
kilca_mec.process_KiLCA(path.join(work_dir, 'TCFP_flre_mec.hdf5'), nrad)
kilca_loc = parcurr()
kilca_loc.process_KiLCA(path.join(work_dir, 'TCFP_flre_loc.hdf5'), nrad)

full_r_range = (0.0, np.amax(kilca_mec.rad[0]))
plt.figure(figsize=canvas)
ax = plt.gca()
for m in range(m_min, m_max + 1):
    plt.plot(cylcoor.rad[m], np.abs(cylcoor.jnpar[m]) * statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
ax.get_yaxis().set_major_formatter(scifmt())
plt.legend(loc='upper right')
plt.xlim(full_r_range)
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
# =============================================================================
# ax_ins_1 = ax.inset_axes([3.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
# ax_ins_1.plot(cylcoor.rad[3], np.abs(cylcoor.jnpar[3]) * statA_per_cm2_to_A_per_m2,
#               '-', lw=thin, c=c[0])
# ax_ins_1.get_yaxis().set_major_formatter(scifmt())
# ax_ins_1.set_xlim(kilca_mec.rres[3] - 0.5, kilca_mec.rres[3] + 0.5)
# ax_ins_1.set_ylim(0.0, 0.8e+04)
# ax.indicate_inset_zoom(ax_ins_1)
# ax_ins_2 = ax.inset_axes([19.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
# ax_ins_2.plot(cylcoor.rad[9], np.abs(cylcoor.jnpar[9]) * statA_per_cm2_to_A_per_m2,
#               '-', lw=thin, c=c[6])
# ax_ins_2.get_yaxis().set_major_formatter(scifmt())
# ax_ins_2.set_xlim(kilca_mec.rres[9] - 0.5, kilca_mec.rres[9] + 0.5)
# ax_ins_2.set_ylim(0.0, 0.8e+04)
# ax.indicate_inset_zoom(ax_ins_2)
# =============================================================================
plt.xlabel(r'$r$ / \si{\centi\meter}')
plt.ylabel(r'$\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
plt.title('Parallel current density of poloidal modes from MEPHIT (cylcoord)')
plt.savefig(path.join(work_dir, 'plot_Jnpar_cycloor.png'), dpi=res)
plt.close()

plt.figure(figsize=canvas)
ax = plt.gca()
for m in range(m_min, m_max + 1):
    plt.plot(fluxcoor.rad[m], np.abs(fluxcoor.jnpar[m]) * statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
ax.get_yaxis().set_major_formatter(scifmt())
plt.legend(loc='upper right')
plt.xlim(full_r_range)
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
# =============================================================================
# ax_ins_1 = ax.inset_axes([3.0, 2.5e+03, 10.0, 3.5e+03], transform=ax.transData)
# ax_ins_1.plot(fluxcoor.rad[3], np.abs(fluxcoor.jnpar[3]) * statA_per_cm2_to_A_per_m2,
#               '-', lw=thin, c=c[0])
# ax_ins_1.get_yaxis().set_major_formatter(scifmt())
# ax_ins_1.set_xlim(kilca_mec.rres[3] - 0.5, kilca_mec.rres[3] + 0.5)
# ax_ins_1.set_ylim(0.0, 3.5e+02)
# ax.indicate_inset_zoom(ax_ins_1)
# ax_ins_2 = ax.inset_axes([19.0, 2.5e+03, 10.0, 3.5e+03], transform=ax.transData)
# ax_ins_2.plot(fluxcoor.rad[9], np.abs(fluxcoor.jnpar[9]) * statA_per_cm2_to_A_per_m2,
#               '-', lw=thin, c=c[6])
# ax_ins_2.get_yaxis().set_major_formatter(scifmt())
# ax_ins_2.set_xlim(kilca_mec.rres[9] - 0.5, kilca_mec.rres[9] + 0.5)
# ax_ins_2.set_ylim(0.0, 3.5e+02)
# ax.indicate_inset_zoom(ax_ins_2)
# =============================================================================
plt.xlabel(r'$r$ / \si{\centi\meter}')
plt.ylabel(r'$\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
plt.title('Parallel current density of poloidal modes from MEPHIT (symfluxcoord)')
plt.savefig(path.join(work_dir, 'plot_Jnpar_fluxcoor.png'), dpi=res)
plt.close()

plt.figure(figsize=canvas)
for m in range(m_min, m_max + 1):
    plt.axvline(kilca_hic.rres[m], lw=0.25 * thin, color='k')
    plt.axvline(kilca_hic.rres[m] - 0.5 * kilca_hic.Bnvac_R_Im[m],
                lw=0.25 * thin, color='k', ls='--')
    plt.axvline(kilca_hic.rres[m] + 0.5 * kilca_hic.Bnvac_R_Im[m],
                lw=0.25 * thin, color='k', ls='--')
    plt.plot(kilca_hic.rad[m], np.abs(kilca_hic.jnpar[m]) * c1_statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
plt.gca().get_yaxis().set_major_formatter(scifmt())
plt.gca().legend(loc='upper left')
plt.xlim(full_r_range)
plt.xlabel(r'minor radius $r$ / \si{\centi\meter}')
plt.ylabel(r'parallel current density $\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
plt.title('Parallel current density from KiLCA ($c = 100$)')
plt.savefig(path.join(work_dir, 'plot_Jnpar_kilca_hic.png'), dpi=res)
plt.close()
plt.figure(figsize=canvas)
for m in range(m_min, m_max + 1):
    plt.axvline(kilca_loc.rres[m], lw=0.25 * thin, color='k')
    plt.axvline(kilca_loc.rres[m] - 0.5 * kilca_loc.Bnvac_R_Im[m],
                lw=0.25 * thin, color='k', ls='--')
    plt.axvline(kilca_loc.rres[m] + 0.5 * kilca_loc.Bnvac_R_Im[m],
                lw=0.25 * thin, color='k', ls='--')
    plt.plot(kilca_loc.rad[m], np.abs(kilca_loc.jnpar[m]) * c1_statA_per_cm2_to_A_per_m2,
             '-', lw=thin, label=f"m = {m}")
plt.gca().get_yaxis().set_major_formatter(scifmt())
plt.gca().legend(loc='upper left')
plt.xlim(full_r_range)
plt.xlabel(r'minor radius $r$ / \si{\centi\meter}')
plt.ylabel(r'parallel current density $\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
plt.title('Parallel current density from KiLCA ($c = 0.01$)')
plt.savefig(path.join(work_dir, 'plot_Jnpar_kilca_loc.png'), dpi=res)
plt.close()

for m in range(m_min, m_max + 1):
    # plot: compare KiLCA collisionality
    plt.figure(figsize=canvas)
    plt.axvline(kilca_mec.Bnvac_R_Im[m], lw=thin, color='k', ls='-')
    plt.axhline(np.abs(kilca_mec.Imnpar[m]), lw=thin, color='k', ls='-')
    plt.plot(kilca_mec.width[m], np.abs(kilca_mec.intJ[m]), '-k', lw=thin,
             label='KiLCA ($c = 1$)')
    plt.axvline(kilca_hic.Bnvac_R_Im[m], lw=thin, color='b', ls='--')
    plt.axhline(np.abs(kilca_hic.Imnpar[m]), lw=thin, color='b', ls='--')
    plt.plot(kilca_hic.width[m], np.abs(kilca_hic.intJ[m]), '--b', lw=thin,
             label='KiLCA ($c = 100$)')
    plt.axvline(kilca_loc.Bnvac_R_Im[m], lw=thin, color='r', ls='--')
    plt.axhline(np.abs(kilca_loc.Imnpar[m]), lw=thin, color='r', ls='--')
    plt.plot(kilca_loc.width[m], np.abs(kilca_loc.intJ[m]), '--r', lw=thin,
             label='KiLCA ($c = 0.01$)')
    plt.gca().legend(loc='lower right', fontsize='x-small')
    plt.xlabel(r'layer width $Bnvac_R_Im$ / cm')
    plt.ylabel(r'parallel current $\vert I_{m n}^{\parallel} \vert$ / A')
    plt.title(f"Comparison: different collisionalities for KiLCA (m = {m})")
    plt.savefig(path.join(work_dir, f"cmp_KiLCA_{m}.png"), dpi=res)
    plt.close()
    # plot: compare KiLCA and MEPHIT
    plt.figure(figsize=canvas)
    plt.axvline(kilca_mec.Bnvac_R_Im[m], lw=0.25 * thin, color='k')
    plt.axhline(np.abs(kilca_mec.Imnpar[m]), lw=0.25 * thin, color='k', ls='-')
    plt.plot(kilca_mec.width[m], np.abs(kilca_mec.intJ[m]), '-k', lw=thin,
             label='KiLCA')
    plt.plot(cylcoor.width[m], np.abs(cylcoor.intJ[m]), '--r', lw=thin,
             label='J, cylcoord.')
    plt.plot(cylcoor.width[m], np.abs(cylcoor.intB[m] + cylcoor.bndryB[m]), '-.m', lw=thin,
             label='B, cylcoord.')
    plt.plot(cylcoor.width[m], np.abs(cylcoor.intB[m]), '--m', lw=thin,
             label='B, cylcoord., part.int.')
    plt.plot(cylcoor.width[m], np.abs(cylcoor.bndryB[m]), ':m', lw=thin,
             label='B, cylcoord, bndry.')
    plt.plot(cylcoor.width[m], np.abs(fluxcoor.intJ[m]), '--b', lw=thin,
             label='J, symfluxcoord.')
# =============================================================================
#     plt.plot(cylcoor.width[m], np.abs(fluxcoor.intB[m] + fluxcoor.bndryB[m]), '-.c', lw=thin,
#              label='B, symfluxcoord.')
#     plt.plot(cylcoor.width[m], np.abs(fluxcoor.intB[m]), '--c', lw=thin,
#              label='B, symfluxcoord., part.int.')
# =============================================================================
    plt.plot(cylcoor.width[m], np.abs(fluxcoor.bndryB[m]), ':c', lw=thin,
             label='B, symfluxcoord, bndry.')
# =============================================================================
#     plt.plot(cylcoor.width[m], np.abs(fluxcoor.GPEC[m]), '--g', lw=thin,
#              label=r'$\Delta_{mn}$')
# =============================================================================
    plt.gca().legend(loc='lower right', fontsize='x-small')
    plt.xlabel(r'layer width $Bnvac_R_Im$ / \si{\centi\meter}')
    plt.ylabel(r'parallel current $\lvert I_{m n}^{\parallel} \rvert$ / \si{\ampere}')
    plt.title(f"Comparison: KiLCA and MEPHIT (m = {m})")
    plt.savefig(path.join(work_dir, f"cmp_Imnpar_{m}.png"), dpi=res)
    plt.close()
