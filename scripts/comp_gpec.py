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
from magdifplot import run_dir, magdif, fslabel, polmodes, magdif_poloidal_plots, magdif_2d_rectplots, \
    cm_to_m, G_to_T, Mx_to_Wb
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.path import Path as Polygon
from matplotlib.colors import Normalize
import colorcet
from scipy.interpolate import UnivariateSpline
import netCDF4

psi_abs = r'$\abs [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\tesla}'
psi_arg = r'$\arg [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\degree}'
psin_abs = r'$\abs [B_{n}^{\perp}]_{m}$ / \si{\tesla}'
psin_arg = r'$\arg [B_{n}^{\perp}]_{m}$ / \si{\degree}'

psi_n = r'$\hat{\psi}$'
theta_n = r'$\hat{\vartheta}$'
sqrt_g = r'$\sqrt{g}$ / \si{\meter\per\tesla}'

for workdir in iglob(run_dir + '/Bvac_ImBm_g33353.2325'):
    testcase = magdif(workdir)
    testcase.read_datafile()
    sgn_dpsi = np.sign(testcase.data['/cache/fs/psi'][-1] -
                       testcase.data['/cache/fs/psi'][0])
    helicity = -np.sign(testcase.data['/cache/fs/q'][-1])
    nemov = magdif(workdir, datafile='magdif_nemov.h5')
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
        workdir, 'GPEC_Bmn_psi_abs.pdf', testcase.data,
        fslabel.psi_norm, psi_abs, np.abs, vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psi_arg.pdf', testcase.data,
        fslabel.psi_norm, psi_arg, partial(np.angle, deg=True), vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_abs.pdf', testcase.data,
        fslabel.psi_norm, psi_abs, np.abs, vac, mephit_vac
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_arg.pdf', testcase.data,
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
        workdir, 'GPEC_Bmn_psin_abs.pdf', testcase.data,
        fslabel.psi_norm, psin_abs, np.abs, vacn, pertn, mephit_vacn, mephit_pertn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psin_arg.pdf', testcase.data,
        fslabel.psi_norm, psin_arg, partial(np.angle, deg=True), vacn, pertn, mephit_vacn, mephit_pertn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_abs.pdf', testcase.data,
        fslabel.psi_norm, psin_abs, np.abs, vacn, mephit_vacn
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psin_arg.pdf', testcase.data,
        fslabel.psi_norm, psin_arg, partial(np.angle, deg=True), vacn, mephit_vacn
    ))

    psi = sgn_dpsi * Mx_to_Wb * testcase.data['/cache/fs/psi'][()]
    psi -= psi[0]
    psipol_max = psi[-1]
    psi /= psipol_max
    jac_psi = Mx_to_Wb / psipol_max * testcase.data['/debug_GPEC/jac/psi'][()]
    jac_nrad = jac_psi.size
    jac_theta = 0.5 / np.pi * testcase.data['/debug_GPEC/jac/theta'][()]
    jac_npol = jac_theta.size
    jac_gpec = cm_to_m / G_to_T * np.transpose(testcase.data['/debug_GPEC/jac/jac'][()])
    jac_mephit = cm_to_m / G_to_T * np.transpose(testcase.data['/debug_GPEC/jac/sqrt_g'][()])
    debug_psi = Mx_to_Wb / psipol_max * testcase.data['/debug_GPEC/psi'][()]
    debug_nrad = debug_psi.size
    debug_theta = 0.5 / np.pi * testcase.data['/debug_GPEC/theta'][()]
    debug_npol = debug_theta.size
    debug_R = cm_to_m * np.transpose(testcase.data['/debug_GPEC/R_GPEC'][()])
    R = cm_to_m * np.transpose(testcase.data['/debug_GPEC/R'][()])
    rootgrp = netCDF4.Dataset(path.join(workdir, 'dcon_output_n2.nc'), 'r')
    spline_F = UnivariateSpline(jac_psi, np.array(rootgrp['f']), s=0)
    spline_q = UnivariateSpline(jac_psi, np.array(rootgrp['q']), s=0)
    rootgrp.close()
    jac_pest = np.empty(debug_R.shape)
    for k in range(jac_pest.shape[0]):
        jac_pest[k, :] = debug_R[k, :] ** 2 * np.abs(spline_q(debug_psi[k]) / spline_F(debug_psi[k]))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(jac_psi, jac_gpec[:, 0], label=r'GPEC, $\hat{\vartheta} = 0$')
    ax.plot(jac_psi, jac_mephit[:, 0], label=r'MEPHIT, $\hat{\vartheta} = 0$')
    ax.plot(jac_psi, jac_gpec[:, jac_npol // 2], label=r'GPEC, $\hat{\vartheta} = 0.5$')
    ax.plot(jac_psi, jac_mephit[:, jac_npol // 2], label=r'MEPHIT, $\hat{\vartheta} = 0.5$')
    ax.set_xlabel(psi_n)
    ax.set_ylabel(sqrt_g)
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'comp_jac_psi.pdf'))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(debug_psi, jac_pest[:, 0], label=r'PEST, $\hat{\vartheta} = 0$')
    ax.plot(jac_psi, jac_mephit[:, 0], label=r'MEPHIT, $\hat{\vartheta} = 0$')
    ax.plot(debug_psi, jac_pest[:, debug_npol // 2], label=r'PEST, $\hat{\vartheta} = 0.5$')
    ax.plot(jac_psi, jac_mephit[:, jac_npol // 2], label=r'MEPHIT, $\hat{\vartheta} = 0.5$')
    ax.set_xlabel(psi_n)
    ax.set_ylabel(sqrt_g)
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'debug_jac_psi.pdf'))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(jac_theta, jac_gpec[jac_nrad // 2, :], label=r'GPEC, $\hat{\psi} = 0.5$')
    ax.plot(jac_theta, jac_mephit[jac_nrad // 2, :], label=r'MEPHIT, $\hat{\psi} = 0.5$')
    ax.plot(jac_theta, jac_gpec[jac_nrad // 4 * 3, :], label=r'GPEC, $\hat{\psi} = 0.75$')
    ax.plot(jac_theta, jac_mephit[jac_nrad // 4 * 3, :], label=r'MEPHIT, $\hat{\psi} = 0.75$')
    ax.set_xlabel(theta_n)
    ax.set_ylabel(sqrt_g)
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'comp_jac_theta.pdf'))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(jac_psi, jac_gpec[:, 0] / jac_mephit[:, 0], label=r'$\hat{\vartheta} = 0$')
    ax.plot(jac_psi, jac_gpec[:, jac_npol // 2] / jac_mephit[:, jac_npol // 2], label=r'$\hat{\vartheta} = 0.5$')
    ax.set_xlabel(psi_n)
    ax.set_ylabel('ratio of Jacobians (GPEC / MEPHIT)')
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'comp_jac_psi_fac.pdf'))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(jac_theta, jac_gpec[jac_nrad // 2, :] / jac_mephit[jac_nrad // 2, :], label=r'$\hat{\psi} = 0.5$')
    ax.plot(jac_theta, jac_gpec[jac_nrad // 4 * 3, :] / jac_mephit[jac_nrad // 4 * 3, :], label=r'$\hat{\psi} = 0.75$')
    ax.set_xlabel(theta_n)
    ax.set_ylabel('ratio of Jacobians (GPEC / MEPHIT)')
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'comp_jac_theta_fac.pdf'))
    fig = Figure()
    ax = fig.subplots()
    ax.plot(debug_psi, debug_R[:, debug_npol // 4], label=r'GPEC, $\hat{\vartheta} = 0.25$')
    ax.plot(debug_psi, R[:, debug_npol // 4], label=r'MEPHIT, $\hat{\vartheta} = 0.25$')
    ax.plot(debug_psi, debug_R[:, debug_npol // 4 * 3], label=r'GPEC, $\hat{\vartheta} = 0.75$')
    ax.plot(debug_psi, R[:, debug_npol // 4 * 3], label=r'MEPHIT, $\hat{\vartheta} = 0.75$')
    ax.set_xlabel(psi_n)
    ax.set_ylabel(r'$R$ / \si{\meter}')
    ax.legend()
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'comp_R.pdf'))

    fig = Figure(figsize=(8.0, 5.0))
    axs = fig.subplots(1, 2, sharex='all', sharey='all')
    images = [axs[0].pcolormesh(jac_theta, jac_psi, jac_mephit, cmap=colorcet.cm.CET_L3,
                                shading='gouraud', rasterized=True),
              axs[1].pcolormesh(debug_theta, debug_psi, jac_pest, cmap=colorcet.cm.CET_L3,
                                shading='gouraud', rasterized=True)]
    axs[0].set_xlabel(theta_n)
    axs[0].set_ylabel(psi_n)
    axs[0].set_title(r'$\sqrt{g}$ from MEPHIT')
    axs[1].set_xlabel(theta_n)
    axs[1].set_ylabel(psi_n)
    axs[1].set_title(r'$\sqrt{g}$ from GPEC (PEST)')
    axs[1].yaxis.set_tick_params(labelleft=True)
    axs[1].yaxis.offsetText.set_visible(True)
    clim = (min(np.amin(image.get_array()) for image in images), max(np.amax(image.get_array()) for image in images))
    norm = Normalize(vmin=clim[0], vmax=clim[1])
    for im in images:
        im.set_norm(norm)
    axs[0].contour(jac_theta, jac_psi, jac_mephit, np.linspace(clim[0], clim[1], 10), colors='w', linewidths=0.1)
    axs[1].contour(debug_theta, debug_psi, jac_pest, np.linspace(clim[0], clim[1], 10), colors='w', linewidths=0.1)
    cbar = fig.colorbar(images[0], ax=axs)
    cbar.set_label(sqrt_g, rotation=90)
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'debug_jac.pdf'))

    fig = Figure(figsize=(8.0, 5.0))
    axs = fig.subplots(1, 2, sharex='all', sharey='all')
    images = [axs[0].pcolormesh(debug_theta, debug_psi, R, cmap=colorcet.cm.CET_L3,
                                shading='gouraud', rasterized=True),
              axs[1].pcolormesh(debug_theta, debug_psi, debug_R, cmap=colorcet.cm.CET_L3,
                                shading='gouraud', rasterized=True)]
    axs[0].set_xlabel(psi_n)
    axs[0].set_ylabel(theta_n)
    axs[0].set_title(r'$R$ from MEPHIT')
    axs[1].set_xlabel(psi_n)
    axs[1].set_ylabel(theta_n)
    axs[1].set_title(r'$R$ from GPEC')
    axs[1].yaxis.set_tick_params(labelleft=True)
    axs[1].yaxis.offsetText.set_visible(True)
    clim = (min(np.amin(image.get_array()) for image in images), max(np.amax(image.get_array()) for image in images))
    norm = Normalize(vmin=clim[0], vmax=clim[1])
    for im in images:
        im.set_norm(norm)
    axs[0].contour(debug_theta, debug_psi, R, np.linspace(clim[0], clim[1], 10), colors='w', linewidths=0.1)
    axs[1].contour(debug_theta, debug_psi, debug_R, np.linspace(clim[0], clim[1], 10), colors='w', linewidths=0.1)
    cbar = fig.colorbar(images[0], ax=axs)
    cbar.set_label(r'$R$ / \si{\meter}', rotation=90)
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(workdir, 'debug_R.pdf'))

    jac_modes_gpec = np.fft.rfft(jac_gpec[:, :-1]) / (jac_npol - 1)
    jac_modes_mephit = np.fft.rfft(jac_mephit[:, :-1]) / (jac_npol - 1)
    for m in range(0, 9):
        fig = Figure()
        ax = fig.subplots()
        ax.plot(jac_psi, jac_modes_gpec[:, m].real, label=fr"$\Real \sqrt{{g}}_{{m = {m}}}$ (GPEC)")
        ax.plot(jac_psi, jac_modes_mephit[:, m].real, label=fr"$\Real \sqrt{{g}}_{{m = {m}}}$ (MEPHIT)")
        ax.plot(jac_psi, jac_modes_gpec[:, m].imag, label=fr"$\Imag \sqrt{{g}}_{{m = {m}}}$ (GPEC)")
        ax.plot(jac_psi, jac_modes_mephit[:, m].imag, label=fr"$\Imag \sqrt{{g}}_{{m = {m}}}$ (MEPHIT)")
        ax.set_xlabel(psi_n)
        ax.set_ylabel(sqrt_g)
        ax.legend()
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"comp_jac_psi_{m}.pdf"))
    delpsi = G_to_T * cm_to_m * np.transpose(testcase.data['/debug_GPEC/delpsi'][()])
    grad_psi = G_to_T * cm_to_m * np.transpose(testcase.data['/debug_GPEC/grad_psi'][()])
    modes_delpsi = np.fft.rfft(delpsi[:, :-1]) / (jac_npol - 1)
    modes_grad_psi = np.fft.rfft(grad_psi[:, :-1]) / (jac_npol - 1)
    for m in range(0, 9):
        fig = Figure()
        ax = fig.subplots()
        ax.plot(debug_psi[::10], modes_delpsi[:, m].real,
                label=fr"$\Real \lVert \nabla \psi \rVert_{{m = {m}}}$ (GPEC)")
        ax.plot(debug_psi[::10], modes_grad_psi[:, m].real,
                label=fr"$\Real \lVert \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)")
        ax.plot(debug_psi[::10], modes_delpsi[:, m].imag,
                label=fr"$\Imag \lVert \nabla \psi \rVert_{{m = {m}}}$ (GPEC)")
        ax.plot(debug_psi[::10], modes_grad_psi[:, m].imag,
                label=fr"$\Imag \lVert \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)")
        ax.set_xlabel(psi_n)
        ax.set_ylabel(r'$\lVert \nabla \psi \rVert$ / \si{\tesla\meter}')
        ax.legend()
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"comp_grad_psi_{m}.pdf"))
    jacfac = cm_to_m ** 2 * np.transpose(testcase.data['/debug_GPEC/jacfac'][()])
    contradenspsi = cm_to_m ** 2 * np.transpose(testcase.data['/debug_GPEC/contradenspsi'][()])
    modes_jacfac = np.fft.fft(jacfac[:, :-1]) / (jac_npol - 1)
    modes_contradenspsi = np.fft.rfft(contradenspsi[:, :-1]) / (jac_npol - 1)
    for m in range(0, 9):
        fig = Figure()
        ax = fig.subplots()
        ax.plot(debug_psi[::10], modes_jacfac[:, m].real,
                label=fr"$\Real \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (GPEC)")
        ax.plot(debug_psi[::10], modes_contradenspsi[:, m].real,
                label=fr"$\Real \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)")
        ax.plot(debug_psi[::10], modes_jacfac[:, m].imag,
                label=fr"$\Imag \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (GPEC)")
        ax.plot(debug_psi[::10], modes_contradenspsi[:, m].imag,
                label=fr"$\Imag \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)")
        ax.set_xlabel(psi_n)
        ax.set_ylabel(r'$\sqrt{g} \lVert \nabla \psi \rVert$ / \si{\meter\squared}')
        ax.legend()
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(workdir, f"comp_jacfac_{m}.pdf"))
    continue
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
