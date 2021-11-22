#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from magdifplot import run_dir
import matplotlib.pyplot as plt
import h5py
from os import path

work_dir = run_dir + '/Bvac_ImBm_g33353.2325'
magdif = h5py.File(path.join(work_dir, 'mephit.h5'), 'r')

ZMID = -1.8
R_O = magdif['/mesh/R_O'][()]
Z_O = magdif['/mesh/Z_O'][()]
R = magdif['/debug_coordinates/R'][()]
Z = magdif['/debug_coordinates/Z'][()]
R_GPEC = magdif['/debug_coordinates/R_GPEC'][()]
Z_GPEC = magdif['/debug_coordinates/Z_GPEC'][()]
xi_n_R = magdif['/debug_coordinates/xi_n_R'][()]
xi_n_Z = magdif['/debug_coordinates/xi_n_Z'][()]
psi_divisor = 109
theta_divisor = 32

plt.figure(figsize=(3.3, 4.4))
plt.plot(R_O, Z_O, 'xk', ms=0.1)
h_GPEC = plt.plot(R_GPEC[:, ::psi_divisor], Z_GPEC[:, ::psi_divisor] + ZMID, '-k', lw=0.0125)
h_MEPHIT = plt.plot(R[:, ::psi_divisor], Z[:, ::psi_divisor], '--r', lw=0.0125)
plt.gca().set_aspect('equal')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')
plt.title(r'$\psi$ coordinate lines')
h_GPEC[0].set_label('GPEC')
h_MEPHIT[0].set_label('MEPHIT')
plt.legend()
plt.savefig(path.join(work_dir, 'debug_coord_psi.pdf'), dpi=600)
plt.close()

plt.figure(figsize=(3.3, 4.4))
plt.plot(R_O, Z_O, 'xk', ms=0.1)
for kpol in range(0, R.shape[0] + 1, R.shape[0] // 32):
    h_GPEC = plt.plot(R_GPEC[kpol, :], Z_GPEC[kpol, :] + ZMID, '-k', lw=0.0125)
    h_MEPHIT = plt.plot(R[kpol, :], Z[kpol, :], '--r', lw=0.0125)
plt.gca().set_aspect('equal')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')
plt.title(r'$\vartheta$ coordinate lines')
h_GPEC[0].set_label('GPEC')
h_MEPHIT[0].set_label('MEPHIT')
plt.gca().legend()
plt.savefig(path.join(work_dir, 'debug_coord_theta.pdf'), dpi=600)
plt.close()

plt.figure(figsize=(3.3, 4.4))
plt.plot(R[:, -1], Z[:, -1], '-k', label='equilibrium')
plt.plot(R[:, -1] + 10 * xi_n_R.real, Z[:, -1] + 10 * xi_n_Z.real, '--r', lw=0.75, label='+ perturbation x10')
plt.gca().set_aspect('equal')
plt.xlabel(r'$R$ / cm')
plt.ylabel(r'$Z$ / cm')
# plt.title('Flux surface displacement')
plt.gca().legend(fontsize='x-small')
plt.savefig(path.join(work_dir, 'xi_n.pdf'), dpi=600)
plt.close()

magdif.close()
