#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:19:58 2020

@author: patrick
"""

from magdifplot import run_dir, parcurr
import numpy as np
import h5py
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from os import path

figsize = (6.6, 3.6)
res = 300
thin = 0.5

work_dir = run_dir + '/30835_3200_ed6'
m_min = 3
m_max = 9

magdif = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')
nrad = magdif['/config/nrad_Ipar'][()]
rad_max = magdif['/cache/fs/rad'][-1]
rad_res = {}
for k, m in enumerate(range(magdif['/mesh/rad_norm_res'].attrs['lbounds'][0],
                            magdif['/mesh/rad_norm_res'].attrs['ubounds'][0])):
    rad_res[m] = magdif['/mesh/rad_norm_res'][k] * rad_max
mephit = parcurr()
mephit.process_magdif(lambda m: magdif, range(m_min, m_max + 1), nrad, rad_res, symfluxcoord=True)
gpec = parcurr()
gpec.process_GPEC(path.join(work_dir, 'gpec_profile_output_n2.nc'))

for m in range(m_min, m_max + 1):
    # plots
    fig = Figure(figsize=figsize)
    ax = fig.subplots()
    ax.axhline(np.abs(gpec.I_res[m]), lw=0.25 * thin, color='k', label='GPEC')
    ax.plot(mephit.width[m], np.abs(mephit.intJ[m]), '--b', lw=thin, label='J, symfluxcoord.')
# =============================================================================
#     ax.plot(mephit.width[m], np.abs(mephit.intB[m] + mephit.bndryB[m]), '-.c', lw=thin, label='B, symfluxcoord.')
#     ac.plot(mephit.width[m], np.abs(mephit.intB[m]), '--c', lw=thin, label='B, symfluxcoord., part.int.')
# =============================================================================
    ax.plot(mephit.width[m], np.abs(mephit.bndryB[m]), ':c', lw=thin, label='B, symfluxcoord, bndry.')
# =============================================================================
#     ac.plot(mephit.width[m], np.abs(mephit.GPEC[m]), '--g', lw=thin, label=r'$\Delta_{mn}$')
# =============================================================================
    ax.legend(loc='lower left', fontsize='x-small')
    ax.set_xlabel(r'layer width $d$ / \si{\centi\meter}')
    ax.set_ylabel(r'parallel current $\lvert I_{m n}^{\parallel} \rvert$ / \si{\ampere}')
    ax.set_title(f"Comparison: GPEC and MEPHIT (m = {m})")
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"cmp_Imnpar_{m}.png"), dpi=res)
