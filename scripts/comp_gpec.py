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
    vac.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), 1, 'b_n_x')  # sgn_dpsi, 'Jbgradgpsi_x'
    pert = polmodes('full perturbation (GPEC)', 'b--')
    pert.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), 1, 'b_n')  # sgn_dpsi, 'Jbgradpsi'
    mephit_vac = polmodes('vacuum perturbation (MEPHIT)', 'r--')
    mephit_vac.read_magdif(testcase.data, fslabel.psi_norm,
                           '/Bmnvac/comp_n')
    for m in mephit_vac.var.keys():
        mephit_vac.var[m] /= cm_to_m ** 2  # Gauss to Tesla instead of Maxwell to Weber
    mephit_pert = polmodes('full perturbation (MEPHIT)', 'g--')
    mephit_pert.read_magdif(testcase.data, fslabel.psi_norm,
                            '/postprocess/Bmn/coeff_n')
    for m in mephit_pert.var.keys():
        mephit_pert.var[m] /= cm_to_m ** 2  # Gauss to Tesla instead of Maxwell to Weber
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_n_abs.pdf', testcase.config, testcase.data,  # _psi_
        fslabel.psi_norm, ylabel_abs, np.abs, vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_n_arg.pdf', testcase.config, testcase.data,  # _psi_
        fslabel.psi_norm, ylabel_arg, partial(np.angle, deg=True), vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_n_abs.pdf', testcase.config, testcase.data,  # _psi_
        fslabel.psi_norm, ylabel_abs, np.abs, vac, mephit_vac
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_n_arg.pdf', testcase.config, testcase.data,  # _psi_
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
        workdir, 'zoom_Bmn_n_abs.pdf', testcase.config, testcase.data,  # _psi_
        fslabel.psi_norm, ylabel_abs, np.abs, vac_zoom, mephit_vac_zoom
    )
    hack.omit_res = True
    testcase.plots.append(hack)
    hack = magdif_poloidal_plots(
        workdir, 'grad_Bmn_n_abs.pdf', testcase.config, testcase.data,  # _psi_
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
    scaling = {'R': (0.1, 0.1), 'Z': (0.1, 0.1), 'phi': (0.001, 0.001)}
    Bnvac = {}
    for coord, dataset in gpec_dataset.items():
        Bnvac['Re'] = [testcase.data['/Bnvac/rect_comp_' + coord][()].real * G_to_T,
                       np.array(rootgrp[dataset][0, :, :] - rootgrp[dataset + '_plasma'][0, :, :])]
        Bnvac['Im'] = [testcase.data['/Bnvac/rect_comp_' + coord][()].imag * G_to_T,
                       np.array(rootgrp[dataset][1, :, :] - rootgrp[dataset + '_plasma'][1, :, :]) * -helicity]
        for part in latex_cmplx.keys():
            label = '$' + latex_cmplx[part] + r' \mathcal{B}_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\tesla}'
            title = [f"Vacuum field from MEPHIT\n({name_cmplx[part]} {name_coord[coord]} component)",
                     f"Vacuum field from GPEC\n({name_cmplx[part]} {name_coord[coord]} component)"]
            testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac[part], label, title=title, clim_scale=scaling[coord],
                                                      filename=path.join(workdir, f'debug_Bnvac_{coord}_{part}.pdf')))

    testcase.dump_plots()
