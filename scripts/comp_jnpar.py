#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:27:00 2020

@author: patrick
"""

from magdifplot import parcurr, statA_per_cm2_to_A_per_m2, c1_statA_per_cm2_to_A_per_m2
import numpy as np
import h5py
from matplotlib import rcParams
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from os import path

figsize = (6.6, 3.6)
res = 300
thin = 0.5

test_dir = '/home/patrick/itp-temp/git/NEO-EQ/run'
work_dir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
m_min = 3
m_max = 9
m_range = range(m_min, m_max + 1)

magdif = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')
nrad = magdif['/config/nrad_Ipar'][()]
rad_max = magdif['/cache/fs/rad'][-1]
rad_res = {}
for k, m in enumerate(range(magdif['/mesh/rad_norm_res'].attrs['lbounds'][0],
                            magdif['/mesh/rad_norm_res'].attrs['ubounds'][0])):
    rad_res[m] = magdif['/mesh/rad_norm_res'][k] * rad_max
h5_files = {}
for m in m_range:
    h5_files[m] = h5py.File(path.join(test_dir, f"TCFP_Ipar_{m}/magdif.h5"), 'r')
cylcoor = parcurr()
cylcoor.process_magdif(lambda m: h5_files[m], m_range, nrad, rad_res, symfluxcoord=False)
fluxcoor = parcurr()
fluxcoor.process_magdif(lambda m: h5_files[m], m_range, nrad, rad_res, symfluxcoord=True)
kilca_hic = parcurr()
kilca_hic.process_KiLCA(path.join(work_dir, 'TCFP_flre_hic.hdf5'), nrad)
kilca_mec = parcurr()
kilca_mec.process_KiLCA(path.join(work_dir, 'TCFP_flre_mec.hdf5'), nrad)
kilca_loc = parcurr()
kilca_loc.process_KiLCA(path.join(work_dir, 'TCFP_flre_loc.hdf5'), nrad)

full_r_range = (0.0, np.amax(kilca_mec.rad[0]))
fig = Figure(figsize=figsize)
ax = fig.subplots()
for m in range(m_min, m_max + 1):
    ax.plot(cylcoor.rad[m], np.abs(cylcoor.jnpar[m]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, label=f"m = {m}")
ax.legend(loc='upper right')
ax.set_xlim(full_r_range)
c = rcParams['axes.prop_cycle'].by_key()['color']
# =============================================================================
# ax_ins_1 = ax.inset_axes([3.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
# ax_ins_1.plot(cylcoor.rad[3], np.abs(cylcoor.jnpar[3]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, c=c[0])
# ax_ins_1.set_xlim(kilca_mec.rres[3] - 0.5, kilca_mec.rres[3] + 0.5)
# ax_ins_1.set_ylim(0.0, 0.8e+04)
# ax.indicate_inset_zoom(ax_ins_1)
# ax_ins_2 = ax.inset_axes([19.0, 0.6e+05, 10.0, 0.8e+05], transform=ax.transData)
# ax_ins_2.plot(cylcoor.rad[9], np.abs(cylcoor.jnpar[9]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, c=c[6])
# ax_ins_2.set_xlim(kilca_mec.rres[9] - 0.5, kilca_mec.rres[9] + 0.5)
# ax_ins_2.set_ylim(0.0, 0.8e+04)
# ax.indicate_inset_zoom(ax_ins_2)
# =============================================================================
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
ax.set_title('Parallel current density of poloidal modes from MEPHIT (cylcoord)')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'plot_Jnpar_cycloor.png'), dpi=res)

fig = Figure(figsize=figsize)
ax = fig.subplots()
for m in range(m_min, m_max + 1):
    ax.plot(fluxcoor.rad[m], np.abs(fluxcoor.jnpar[m]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, label=f"m = {m}")
ax.legend(loc='upper right')
ax.set_xlim(full_r_range)
c = rcParams['axes.prop_cycle'].by_key()['color']
# =============================================================================
# ax_ins_1 = ax.inset_axes([3.0, 2.5e+03, 10.0, 3.5e+03], transform=ax.transData)
# ax_ins_1.plot(fluxcoor.rad[3], np.abs(fluxcoor.jnpar[3]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, c=c[0])
# ax_ins_1.set_xlim(kilca_mec.rres[3] - 0.5, kilca_mec.rres[3] + 0.5)
# ax_ins_1.set_ylim(0.0, 3.5e+02)
# ax.indicate_inset_zoom(ax_ins_1)
# ax_ins_2 = ax.inset_axes([19.0, 2.5e+03, 10.0, 3.5e+03], transform=ax.transData)
# ax_ins_2.plot(fluxcoor.rad[9], np.abs(fluxcoor.jnpar[9]) * statA_per_cm2_to_A_per_m2, '-', lw=thin, c=c[6])
# ax_ins_2.set_xlim(kilca_mec.rres[9] - 0.5, kilca_mec.rres[9] + 0.5)
# ax_ins_2.set_ylim(0.0, 3.5e+02)
# ax.indicate_inset_zoom(ax_ins_2)
# =============================================================================
ax.set_xlabel(r'$r$ / \si{\centi\meter}')
ax.set_ylabel(r'$\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
ax.set_title('Parallel current density of poloidal modes from MEPHIT (symfluxcoord)')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'plot_Jnpar_fluxcoor.png'), dpi=res)

fig = Figure(figsize=figsize)
ax = fig.subplots()
for m in range(m_min, m_max + 1):
    ax.axvline(kilca_hic.rres[m], lw=0.25 * thin, color='k')
    ax.axvline(kilca_hic.rres[m] - 0.5 * kilca_hic.d[m], lw=0.25 * thin, color='k', ls='--')
    ax.axvline(kilca_hic.rres[m] + 0.5 * kilca_hic.d[m], lw=0.25 * thin, color='k', ls='--')
    ax.plot(kilca_hic.rad[m], np.abs(kilca_hic.jnpar[m]) * c1_statA_per_cm2_to_A_per_m2, '-', lw=thin, label=f"m = {m}")
ax.legend(loc='upper left')
ax.set_xlim(full_r_range)
ax.set_xlabel(r'minor radius $r$ / \si{\centi\meter}')
ax.set_ylabel(r'parallel current density $\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
ax.set_title('Parallel current density from KiLCA ($c = 100$)')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'plot_Jnpar_kilca_hic.png'), dpi=res)

fig = Figure(figsize=figsize)
ax = fig.subplots()
for m in range(m_min, m_max + 1):
    ax.axvline(kilca_loc.rres[m], lw=0.25 * thin, color='k')
    ax.axvline(kilca_loc.rres[m] - 0.5 * kilca_loc.d[m], lw=0.25 * thin, color='k', ls='--')
    ax.axvline(kilca_loc.rres[m] + 0.5 * kilca_loc.d[m], lw=0.25 * thin, color='k', ls='--')
    ax.plot(kilca_loc.rad[m], np.abs(kilca_loc.jnpar[m]) * c1_statA_per_cm2_to_A_per_m2, '-', lw=thin, label=f"m = {m}")
ax.legend(loc='upper left')
ax.set_xlim(full_r_range)
ax.set_xlabel(r'minor radius $r$ / \si{\centi\meter}')
ax.set_ylabel(r'parallel current density $\lvert J_{m n}^{\parallel} \rvert$ / \si{\ampere\per\meter\squared')
ax.set_title('Parallel current density from KiLCA ($c = 0.01$)')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'plot_Jnpar_kilca_loc.png'), dpi=res)

for m in range(m_min, m_max + 1):
    # plot: compare KiLCA collisionality
    fig = Figure(figsize=figsize)
    ax = fig.subplots()
    ax.axvline(kilca_mec.d[m], lw=thin, color='k', ls='-')
    ax.axhline(np.abs(kilca_mec.Imnpar[m]), lw=thin, color='k', ls='-')
    ax.plot(kilca_mec.width[m], np.abs(kilca_mec.intJ[m]), '-k', lw=thin, label='KiLCA ($c = 1$)')
    ax.axvline(kilca_hic.d[m], lw=thin, color='b', ls='--')
    ax.axhline(np.abs(kilca_hic.Imnpar[m]), lw=thin, color='b', ls='--')
    ax.plot(kilca_hic.width[m], np.abs(kilca_hic.intJ[m]), '--b', lw=thin, label='KiLCA ($c = 100$)')
    ax.axvline(kilca_loc.d[m], lw=thin, color='r', ls='--')
    ax.axhline(np.abs(kilca_loc.Imnpar[m]), lw=thin, color='r', ls='--')
    ax.plot(kilca_loc.width[m], np.abs(kilca_loc.intJ[m]), '--r', lw=thin, label='KiLCA ($c = 0.01$)')
    ax.legend(loc='lower right', fontsize='x-small')
    ax.set_xlabel(r'layer width $d$ / cm')
    ax.set_ylabel(r'parallel current $\vert I_{m n}^{\parallel} \vert$ / A')
    ax.set_title(f"Comparison: different collisionalities for KiLCA (m = {m})")
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"cmp_KiLCA_{m}.png"), dpi=res)
    # plot: compare KiLCA and MEPHIT
    fig = Figure(figsize=figsize)
    ax = fig.subplots()
    ax.axvline(kilca_mec.d[m], lw=0.25 * thin, color='k')
    ax.axhline(np.abs(kilca_mec.Imnpar[m]), lw=0.25 * thin, color='k', ls='-')
    ax.plot(kilca_mec.width[m], np.abs(kilca_mec.intJ[m]), '-k', lw=thin, label='KiLCA')
    ax.plot(cylcoor.width[m], np.abs(cylcoor.intJ[m]), '--r', lw=thin, label='J, cylcoord.')
    ax.plot(cylcoor.width[m], np.abs(cylcoor.intB[m] + cylcoor.bndryB[m]), '-.m', lw=thin, label='B, cylcoord.')
    ax.plot(cylcoor.width[m], np.abs(cylcoor.intB[m]), '--m', lw=thin, label='B, cylcoord., part.int.')
    ax.plot(cylcoor.width[m], np.abs(cylcoor.bndryB[m]), ':m', lw=thin, label='B, cylcoord, bndry.')
    ax.plot(cylcoor.width[m], np.abs(fluxcoor.intJ[m]), '--b', lw=thin, label='J, symfluxcoord.')
# =============================================================================
#     ax.plot(cylcoor.width[m], np.abs(fluxcoor.intB[m] + fluxcoor.bndryB[m]), '-.c', lw=thin, label='B, symfluxcoord.')
#     ax.plot(cylcoor.width[m], np.abs(fluxcoor.intB[m]), '--c', lw=thin, label='B, symfluxcoord., part.int.')
# =============================================================================
    ax.plot(cylcoor.width[m], np.abs(fluxcoor.bndryB[m]), ':c', lw=thin, label='B, symfluxcoord, bndry.')
# =============================================================================
#     ax.plot(cylcoor.width[m], np.abs(fluxcoor.GPEC[m]), '--g', lw=thin, label=r'$\Delta_{mn}$')
# =============================================================================
    ax.legend(loc='lower right', fontsize='x-small')
    ax.set_xlabel(r'layer width $d$ / \si{\centi\meter}')
    ax.set_ylabel(r'parallel current $\lvert I_{m n}^{\parallel} \rvert$ / \si{\ampere}')
    ax.set_title(f"Comparison: KiLCA and MEPHIT (m = {m})")
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"cmp_Imnpar_{m}.png"), dpi=res)
