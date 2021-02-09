#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib import use
import matplotlib.pyplot as plt
import h5py
from os import path

use('Agg')
work_dir = '/temp/lainer_p/git/NEO-EQ/run/Bvac_ImBm_g33353.2325'
magdif = h5py.File(path.join(work_dir, 'magdif.h5'), 'r')

R = magdif['/debug_coordinates/R'][()]
Z = magdif['/debug_coordinates/Z'][()] + 1.8
R_GPEC = magdif['/debug_coordinates/R_GPEC'][()]
Z_GPEC = magdif['/debug_coordinates/Z_GPEC'][()]
psi_divisor = 109
theta_divisor = 32

plt.figure(figsize=(3.3, 4.4))
plt.plot(R_GPEC[:, ::psi_divisor], Z_GPEC[:, ::psi_divisor], 'k-', lw=0.0125)
plt.plot(R[:, ::psi_divisor], Z[:, ::psi_divisor], 'r--', lw=0.0125)
plt.gca().set_aspect('equal')
plt.savefig(path.join(work_dir, 'debug_coord_psi.pdf'), dpi=600)
plt.close()

plt.figure(figsize=(3.3, 4.4))
for kpol in range(0, R.shape[0] + 1, R.shape[0] // 32):
    plt.plot(R_GPEC[kpol, :], Z_GPEC[kpol, :], 'k-', lw=0.0125)
    plt.plot(R[kpol, :], Z[kpol, :], 'r--', lw=0.0125)
plt.gca().set_aspect('equal')
plt.savefig(path.join(work_dir, 'debug_coord_theta.pdf'), dpi=600)
plt.close()

magdif.close()
