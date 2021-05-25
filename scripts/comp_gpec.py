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
from magdifplot import magdif, fslabel, polmodes, magdif_poloidal_plots, magdif_2d_rectplots, cm_to_m, G_to_T, Mx_to_Wb
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.path import Path as Polygon
import netCDF4

psi_abs = r'$\abs [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\tesla}'
psi_arg = r'$\arg [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\degree}'
psin_abs = r'$\abs [B_{n}^{\perp}]_{m}$ / \si{\tesla}'
psin_arg = r'$\arg [B_{n}^{\perp}]_{m}$ / \si{\degree}'

for workdir in iglob('/temp/lainer_p/git/NEO-EQ/run/Bvac_ImBm_g33353.2325'):
    testcase = magdif(workdir)
    testcase.read_configfile()
    testcase.read_datafile()
    sgn_dpsi = np.sign(testcase.data['/cache/fs/psi'][-1] -
                       testcase.data['/cache/fs/psi'][0])
    helicity = -np.sign(testcase.data['/cache/fs/q'][-1])
    nemov = magdif(workdir, datafile='magdif_nemov.h5')
    nemov.read_configfile()
    nemov.read_datafile()

    vac = polmodes('vacuum perturbation (GPEC)', 'k--')
    vac.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi_x')
    pert = polmodes('full perturbation (GPEC)', 'b--')
    pert.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi')
    mephit_vac = polmodes('vacuum perturbation (MEPHIT)', 'r--')
    mephit_vac.read_magdif(testcase.data, fslabel.psi_norm, '/Bmnvac/comp_psi_contravar_dens')
    for m in mephit_vac.var.keys():
        mephit_vac.var[m] /= testcase.data['/mesh/gpec_jacfac'][:, 16] * cm_to_m ** 2
    mephit_pert = polmodes('full perturbation (MEPHIT)', 'g--')
    mephit_pert.read_magdif(testcase.data, fslabel.psi_norm, '/postprocess/Bmn/coeff_rad')
    for m in mephit_pert.var.keys():
        mephit_pert.var[m] /= testcase.data['/mesh/gpec_jacfac'][:, 16] * cm_to_m ** 2
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psi_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psi_abs, np.abs, vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psi_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psi_arg, partial(np.angle, deg=True), vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psi_abs, np.abs, vac, mephit_vac
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psi_arg, partial(np.angle, deg=True), vac, mephit_vac
    ))
    vacn = polmodes('vacuum perturbation (GPEC)', 'k--')
    vacn.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'b_n_x')
    pertn = polmodes('full perturbation (GPEC)', 'b--')
    pertn.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'b_n')
    mephit_vacn = polmodes('vacuum perturbation (MEPHIT)', 'r--')
    mephit_vacn.read_magdif(testcase.data, fslabel.psi_norm, '/Bmnvac/comp_n')
    for m in mephit_vacn.var.keys():
        mephit_vacn.var[m] /= cm_to_m ** 2
    mephit_pertn = polmodes('full perturbation (MEPHIT)', 'g--')
    mephit_pertn.read_magdif(testcase.data, fslabel.psi_norm, '/postprocess/Bmn/coeff_n')
    for m in mephit_pertn.var.keys():
        mephit_pertn.var[m] /= cm_to_m ** 2
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psin_abs, np.abs, vacn, pertn, mephit_vacn, mephit_pertn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psin_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psin_arg, partial(np.angle, deg=True), vacn, pertn, mephit_vacn, mephit_pertn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psin_abs, np.abs, vacn, mephit_vacn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, psin_arg, partial(np.angle, deg=True), vacn, mephit_vacn
    ))
    jacfac = deepcopy(vac)
    data = np.loadtxt(path.join(workdir, 'gpec_diagnostics_jacfac_1.out'), skiprows=2)
    m_range = np.unique(data[:, 1].astype(int))
    npsi = data.shape[0] // m_range.size
    for m in m_range:
        jacfac.rho[m] = np.empty((npsi,), dtype='d')
        jacfac.var[m] = np.empty((npsi,), dtype='D')
    k = 0
    for kpsi in np.arange(npsi):
        for m in m_range:
            jacfac.rho[m][kpsi] = data[k, 0]
            jacfac.var[m].real[kpsi] = data[k, 2] * sgn_dpsi
            jacfac.var[m].imag[kpsi] = data[k, 3] * sgn_dpsi * helicity
            k += 1
    mephit_jacfac = polmodes('weighting factor (MEPHIT)', 'r--')
    mephit_jacfac.read_magdif(testcase.data, fslabel.psi_norm, '/mesh/gpec_jacfac')
    for m in mephit_jacfac.var.keys():
        mephit_jacfac.var[m] /= Mx_to_Wb * testcase.data['/mesh/gpec_jacfac'][:, 16]
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_jacfac_abs.pdf', testcase.config, testcase.data, fslabel.psi_norm,
        r'$\abs [\sqrt{g} \lVert \nabla \psi \rVert]_{m}$ / \si{\square\meter}',
        np.abs, jacfac, mephit_jacfac
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_jacfac_arg.pdf', testcase.config, testcase.data, fslabel.psi_norm,
        r'$\arg [\sqrt{g} \lVert \nabla \psi \rVert]_{m} $ / \si{\square\meter}',
        partial(np.angle, deg=True), jacfac, mephit_jacfac
    ))

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

        Nemov_Re = covar * nemov.data['/Bnvac/rect_comp_' + coord][()].real * G_to_T
        Nemov_Re[nan_mask] = np.nan
        Nemov_Im = covar * nemov.data['/Bnvac/rect_comp_' + coord][()].imag * G_to_T
        Nemov_Im[nan_mask] = np.nan
        fig = Figure()
        axs = fig.subplots(1, 2, sharex='all', sharey='all')
        axs[0].plot(Z[1], Bnvac['Re'][1][:, kR], label='GPEC')
        axs[0].plot(Z[0], Bnvac['Re'][0][:, kR], label='MEPHIT')
        axs[0].plot(Z[0], Nemov_Re[:, kR], label='Nemov')
        axs[0].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[0].legend()
        axs[0].set_xlabel(r'$Z$ / \si{\meter}')
        axs[0].set_ylabel(r'$\Real B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[0].set_title(f"Vacuum field at $R = \\SI{{{R[0][kR]:1.3}}}{{\\meter}}$\n" +
                         f"(real {name_coord[coord]} component)")
        axs[1].plot(Z[1], Bnvac['Im'][1][:, kR], label='GPEC')
        axs[1].plot(Z[0], Bnvac['Im'][0][:, kR], label='MEPHIT')
        axs[1].plot(Z[0], Nemov_Im[:, kR], label='Nemov')
        axs[1].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[1].legend()
        axs[1].set_xlabel(r'$Z$ / \si{\meter}')
        axs[1].set_ylabel(r'$\Imag B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].set_title(f"Vacuum field at $R = \\SI{{{R[0][kR]:1.3}}}{{\\meter}}$\n" +
                         f"(imaginary {name_coord[coord]} component)")
        axs[1].yaxis.set_tick_params(labelleft=True)
        axs[1].yaxis.offsetText.set_visible(True)
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_Z_Bnvac_{coord}.pdf"))
        fig = Figure()
        axs = fig.subplots(1, 2, sharex='all', sharey='all')
        axs[0].plot(R[1], Bnvac['Re'][1][kZ, :], label='GPEC')
        axs[0].plot(R[0], Bnvac['Re'][0][kZ, :], label='MEPHIT')
        axs[0].plot(R[0], Nemov_Re[kZ, :], label='Nemov')
        axs[0].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[0].legend()
        axs[0].set_xlabel(r'$R$ / \si{\meter}')
        axs[0].set_ylabel(r'$\Real B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[0].set_title(f"Vacuum field at $Z = \\SI{{{Z[0][kZ]:1.3}}}{{\\meter}}$\n" +
                         f"(real {name_coord[coord]} component)")
        axs[1].plot(R[1], Bnvac['Im'][1][kZ, :], label='GPEC')
        axs[1].plot(R[0], Bnvac['Im'][0][kZ, :], label='MEPHIT')
        axs[1].plot(R[0], Nemov_Im[kZ, :], label='Nemov')
        axs[1].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[1].legend()
        axs[1].set_xlabel(r'$R$ / \si{\meter}')
        axs[1].set_ylabel(r'$\Imag B_{n \mathrm{v} ' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].set_title(f"Vacuum field at $Z = \\SI{{{Z[0][kZ]:1.3}}}{{\\meter}}$\n" +
                         f"(imaginary {name_coord[coord]} component)")
        axs[1].yaxis.set_tick_params(labelleft=True)
        axs[1].yaxis.offsetText.set_visible(True)
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_R_Bnvac_{coord}.pdf"))

        B0[0][nan_mask] = np.nan
        B0[1][nan_mask] = np.nan
        fig = Figure()
        axs = fig.subplots(1, 2)
        axs[0].plot(Z[1], B0[1][:, kR], '-k', label='GPEC')
        axs[0].plot(Z[0], B0[0][:, kR], '--r', label='MEPHIT')
        axs[0].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[0].legend()
        axs[0].set_xlabel(r'$Z$ / \si{\meter}')
        axs[0].set_ylabel(r'$B_{0}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[0].set_title(f"Equilibrium field at $R = \\SI{{{R[0][kR]:1.3}}}{{\\meter}}$\n" +
                         f"({name_coord[coord]} component)")
        axs[1].plot(R[1], B0[1][kZ, :] * R[1], '-k', label='GPEC')
        axs[1].plot(R[0], B0[0][kZ, :] * R[0], '--r', label='MEPHIT')
        axs[1].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
        axs[1].legend()
        axs[1].set_xlabel(r'$R$ / \si{\meter}')
        axs[1].set_ylabel(r'$B_{0}^{' + latex_coord[coord] + r'}$ / \si{\tesla}')
        axs[1].set_title(f"Equilibrium field at $Z = \\SI{{{Z[0][kZ]:1.3}}}{{\\meter}}$\n" +
                         f"({name_coord[coord]} component)")
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"cut_B0_{coord}.pdf"))

    testcase.dump_plots()
