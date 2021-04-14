#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

from os import path
from glob import iglob
import numpy as np
from copy import deepcopy
from functools import partial
from magdifplot import magdif, fslabel, polmodes, magdif_poloidal_plots, magdif_2d_rectplots, cm_to_m, G_to_T
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.path import Path as Polygon
import netCDF4

# ylabel_abs = r'$\lvert \mathcal{B}_{mn}^{\psi} \rvert$ / \si{\weber}'
# ylabel_arg = r'$\arg \mathcal{B}_{mn}^{\psi}$ / \si{\degree}'
# ylabel_grad_abs = r'$\partial_{\hat{\psi}} \lvert \mathcal{B}_{mn}^{\psi} \rvert$ / \si{\weber}'
ylabel_abs = r'$\lvert B_{mn}^{\hat{\psi}} \rvert$ / \si{\tesla}'
ylabel_arg = r'$\arg B_{mn}^{\hat{\psi}}$ / \si{\degree}'
ylabel_grad_abs = r'$\partial_{\hat{\psi}} \lvert B_{mn}^{\hat{\psi}} \rvert$ / \si{\tesla}'
cutoff = 0.2

for workdir in iglob('/temp/lainer_p/git/NEO-EQ/run/Bvac_ImBm_g33353.2325'):
    testcase = magdif(workdir)
    testcase.read_configfile()
    testcase.read_datafile()
    sgn_dpsi = np.sign(testcase.data['/cache/fs/psi'][-1] -
                       testcase.data['/cache/fs/psi'][0])
    helicity = -np.sign(testcase.data['/cache/fs/q'][-1])

    vac = polmodes('vacuum perturbation (GPEC)', 'k--')
    vac.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi_x')
    pert = polmodes('full perturbation (GPEC)', 'b--')
    pert.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi')
    mephit_vac = polmodes('vacuum perturbation (MEPHIT)', 'r--')
    mephit_vac.read_magdif(testcase.data, fslabel.psi_norm, '/Bmnvac/comp_psi_contravar_dens')
    for m in mephit_vac.var.keys():
        mephit_vac.var[m] /= testcase.data['/mesh/gpec_jarea'][()] * cm_to_m ** 2
    mephit_pert = polmodes('full perturbation (MEPHIT)', 'g--')
    mephit_pert.read_magdif(testcase.data, fslabel.psi_norm, '/postprocess/Bmn/coeff_rad')
    for m in mephit_pert.var.keys():
        mephit_pert.var[m] /= testcase.data['/mesh/gpec_jarea'][()] * cm_to_m ** 2
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_abs, np.abs, vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psin_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_arg, partial(np.angle, deg=True), vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_abs, np.abs, vac, mephit_vac
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_arg, partial(np.angle, deg=True), vac, mephit_vac
    ))
    vac_grad = deepcopy(vac)
    for m in vac_grad.var.keys():
        vac_grad.var[m] = np.gradient(np.abs(vac_grad.var[m]), vac_grad.rho[m])
        vac_grad.var[m] = np.delete(vac_grad.var[m], np.nonzero(vac_grad.rho[m] > cutoff))
        vac_grad.rho[m] = np.delete(vac_grad.rho[m], np.nonzero(vac_grad.rho[m] > cutoff))
    mephit_vac_grad = deepcopy(mephit_vac)
    for m in mephit_vac.var.keys():
        mephit_vac_grad.var[m] = np.gradient(np.abs(mephit_vac_grad.var[m]), mephit_vac_grad.rho[m])
        mephit_vac_grad.var[m] = np.delete(mephit_vac_grad.var[m], np.nonzero(mephit_vac_grad.rho[m] > cutoff))
        mephit_vac_grad.rho[m] = np.delete(mephit_vac_grad.rho[m], np.nonzero(mephit_vac_grad.rho[m] > cutoff))
    vac_zoom = deepcopy(vac)
    for m in vac_zoom.var.keys():
        vac_zoom.var[m] = np.delete(vac_zoom.var[m], np.nonzero(vac_zoom.rho[m] > cutoff))
        vac_zoom.rho[m] = np.delete(vac_zoom.rho[m], np.nonzero(vac_zoom.rho[m] > cutoff))
    mephit_vac_zoom = deepcopy(mephit_vac)
    for m in mephit_vac.var.keys():
        mephit_vac_zoom.var[m] = np.delete(mephit_vac_zoom.var[m], np.nonzero(mephit_vac_zoom.rho[m] > cutoff))
        mephit_vac_zoom.rho[m] = np.delete(mephit_vac_zoom.rho[m], np.nonzero(mephit_vac_zoom.rho[m] > cutoff))
    hack = magdif_poloidal_plots(
        workdir, 'zoom_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_abs, np.abs, vac_zoom, mephit_vac_zoom
    )
    hack.omit_res = True
    testcase.plots.append(hack)
    hack = magdif_poloidal_plots(
        workdir, 'grad_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_grad_abs, lambda v: v, vac_grad, mephit_vac_grad
    )
    hack.omit_res = True
    testcase.plots.append(hack)

    rootgrp = netCDF4.Dataset(path.join(workdir, 'gpec_cylindrical_output_n2.nc'), 'r')
    R = [testcase.data['/Bnvac/rect_R'][()] * cm_to_m, np.array(rootgrp['R'])]
    Z = [testcase.data['/Bnvac/rect_Z'][()] * cm_to_m, np.array(rootgrp['z'])]
    gpec_dataset = {'R': 'b_r', 'Z': 'b_z', 'phi': 'b_t'}
    latex_coord = {'R': 'R', 'Z': 'Z', 'phi': r'\varphi'}
    latex_cmplx = {'Re': r'\Real', 'Im': r'\Imag'}
    name_coord = {'R': 'radial', 'Z': 'axial', 'phi': 'azimuthal'}
    name_cmplx = {'Re': 'real', 'Im': 'imaginary'}
    kp_max_lcfs = testcase.data['/mesh/kp_max'][-1]
    separatrix = Polygon(cm_to_m * np.column_stack((testcase.data['/mesh/node_R'][-kp_max_lcfs:],
                                                    testcase.data['/mesh/node_Z'][-kp_max_lcfs:])))
    RR, ZZ = np.meshgrid(R[0], Z[0])
    nan_mask = np.logical_not(separatrix.contains_points(np.column_stack((RR.ravel(), ZZ.ravel())))).reshape(RR.shape)
    RR[nan_mask] = np.nan
    ZZ[nan_mask] = np.nan
    kR = 83
    kZ = 150
    Bnvac = {}
    for coord, dataset in gpec_dataset.items():
        covar = np.ones(RR.shape) if coord != 'phi' else RR
        Bnvac['Re'] = [covar * testcase.data['/Bnvac/rect_comp_' + coord][()].real * G_to_T,
                       covar * np.array(rootgrp[dataset][0, :, :] - rootgrp[dataset + '_plasma'][0, :, :]) * 0.5]
        Bnvac['Im'] = [covar * testcase.data['/Bnvac/rect_comp_' + coord][()].imag * G_to_T,
                       covar * np.array(rootgrp[dataset][1, :, :] - rootgrp[dataset + '_plasma'][1, :, :]) * -0.5]
        B0 = [covar * testcase.data['/equil/B0_' + coord][()] * G_to_T, covar * rootgrp[dataset + '_equil'][:, :]]
        for part in latex_cmplx.keys():
            Bnvac[part][0][nan_mask] = np.nan
            Bnvac[part][1][nan_mask] = np.nan
            label = '$' + latex_cmplx[part] + r' B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}'
            title = [f"Vacuum field from MEPHIT\n({name_cmplx[part]} {name_coord[coord]} component)",
                     f"Vacuum field from GPEC\n({name_cmplx[part]} {name_coord[coord]} component)"]
            testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac[part], label, title=title,
                                                      filename=path.join(workdir, f'debug_Bnvac_{coord}_{part}.pdf')))
        B0[0][nan_mask] = np.nan
        B0[1][nan_mask] = np.nan
        centered = coord != 'phi'
        label = '$' + latex_cmplx[part] + r' B_{0 ' + latex_coord[coord] + r'}$ / \si{\tesla}'
        title = [f"Equilibrium field from MEPHIT\n({name_coord[coord]} component)",
                 f"Equilibrium field from GPEC\n({name_coord[coord]} component)"]
        testcase.plots.append(magdif_2d_rectplots(R, Z, B0, label, title=title, centered=centered,
                                                  filename=path.join(workdir, f'debug_B0_{coord}.pdf')))

        Bnvac_diff = [Bnvac['Re'][1] - Bnvac['Re'][0], Bnvac['Im'][1] - Bnvac['Im'][0]]
        label = (r'$\Delta B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        title = [f"Vacuum field: GPEC - MEPHIT\n(real {name_coord[coord]} component)",
                 f"Vacuum field: GPEC - MEPHIT\n(imaginary {name_coord[coord]} component)"]
        testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac_diff, label, title=title,
                                                  filename=path.join(workdir, f'debug_Bnvac_{coord}_absdiff.pdf')))
        B0_diff = [B0[1] - B0[0], B0[1] - B0[0]]
        label = (r'$\Delta B_{0 ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        title = [f"Equilibrium field: GPEC - MEPHIT\n({name_coord[coord]} component)",
                 f"Equilibrium field: GPEC - MEPHIT\n({name_coord[coord]} component)"]
        testcase.plots.append(magdif_2d_rectplots(R, Z, B0_diff, label, title=title,
                                                  filename=path.join(workdir, f'debug_B0_{coord}_absdiff.pdf')))

        fig = Figure()
        axs = fig.subplots(1, 2, sharex='all', sharey='all')
        axs[0].plot(Z[1], Bnvac['Re'][1][:, kR], '-k', label='GPEC')
        axs[0].plot(Z[0], Bnvac['Re'][0][:, kR], '--r', label='MEPHIT')
        axs[0].legend()
        axs[0].set_xlabel(r"$Z$ / \si{\meter}")
        axs[0].set_ylabel(r'$\Real B_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].plot(Z[1], Bnvac['Im'][1][:, kR], '-k', label='GPEC')
        axs[1].plot(Z[0], Bnvac['Im'][0][:, kR], '--r', label='MEPHIT')
        axs[1].legend()
        axs[1].set_xlabel(r"$Z$ / \si{\meter}")
        axs[1].set_ylabel(r'$\Imag B_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_Z_Bnvac_{coord}.pdf"))
        fig = Figure()
        axs = fig.subplots(1, 2, sharex='all', sharey='all')
        axs[0].plot(R[1], Bnvac['Re'][1][kZ, :], '-k', label='GPEC')
        axs[0].plot(R[0], Bnvac['Re'][0][kZ, :], '--r', label='MEPHIT')
        axs[0].legend()
        axs[0].set_xlabel(r"$R$ / \si{\meter}")
        axs[0].set_ylabel(r'$\Real B_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].plot(R[1], Bnvac['Im'][1][kZ, :], '-k', label='GPEC')
        axs[1].plot(R[0], Bnvac['Im'][0][kZ, :], '--r', label='MEPHIT')
        axs[1].legend()
        axs[1].set_xlabel(r"$R$ / \si{\meter}")
        axs[1].set_ylabel(r'$\Imag B_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_R_Bnvac_{coord}.pdf"))

        B0[0][nan_mask] = np.nan
        B0[1][nan_mask] = np.nan
        fig = Figure()
        axs = fig.subplots(1, 2)
        axs[0].plot(Z[1], B0[1][:, kR], '-k', label='GPEC')
        axs[0].plot(Z[0], B0[0][:, kR], '--r', label='MEPHIT')
        axs[0].legend()
        axs[0].set_xlabel(r"$Z$ / \si{\meter}")
        axs[0].set_ylabel(r'$B_{0}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].plot(R[1], B0[1][kZ, :] * R[1], '-k', label='GPEC')
        axs[1].plot(R[0], B0[0][kZ, :] * R[0], '--r', label='MEPHIT')
        axs[1].legend()
        axs[1].set_xlabel(r"$R$ / \si{\meter}")
        axs[1].set_ylabel(r'$B_{0}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_B0_{coord}.pdf"))

    testcase.dump_plots()
