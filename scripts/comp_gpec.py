#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

from os import path
from glob import iglob
import numpy as np
from functools import partial
from magdifplot import magdif, fslabel, polmodes, magdif_poloidal_plots, magdif_2d_rectplots, cm_to_m, G_to_T
import netCDF4

ylabel_abs = r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\tesla}'
ylabel_arg = r'$\arg \sqrt{g} B_{mn}^{\psi}$ / \si{\degree}'

for workdir in iglob('/temp/lainer_p/git/NEO-EQ/run/Bvac_*'):
    testcase = magdif(workdir)
    testcase.read_configfile()
    testcase.read_datafile()
    sgn_dpsi = np.sign(testcase.data['/cache/fs/psi'][-1] -
                       testcase.data['/cache/fs/psi'][0])
    helicity = -np.sign(testcase.data['/cache/fs/q'][-1])

    vac = polmodes('vacuum perturbation (GPEC)', 'r--')
    vac.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi_x')
    pert = polmodes('full perturbation (GPEC)', 'k--')
    pert.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_dpsi, 'Jbgradpsi')
    mephit_vac = polmodes('vacuum perturbation (MEPHIT)', 'b--')
    mephit_vac.read_magdif(testcase.data, fslabel.psi_norm,
                           '/Bmnvac/comp_psi_contravar_dens')
    mephit_pert = polmodes('full perturbation (MEPHIT)', 'g--')
    mephit_pert.read_magdif(testcase.data, fslabel.psi_norm,
                            '/postprocess/Bmn/coeff_rad')
    fouriermodes = polmodes('vacuum perturbation (amn.dat)', 'c--')
    fouriermodes.read_fouriermodes(testcase.data)
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psi_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_abs, np.abs, vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'GPEC_Bmn_psi_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_arg, partial(np.angle, deg=True), vac, pert, mephit_vac, mephit_pert
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_abs.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_abs, np.abs, vac, mephit_vac, fouriermodes
    ))
    testcase.plots.append(magdif_poloidal_plots(
        workdir, 'debug_Bmn_psi_arg.pdf', testcase.config, testcase.data,
        fslabel.psi_norm, ylabel_arg, partial(np.angle, deg=True), vac, mephit_vac, fouriermodes
    ))

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
            label = '$' + latex_cmplx[part] + r' \mathcal{B}_{n \mathrm{v}}^{' + latex_coord[coord] + r'}$ / \si{\gauss}'
            title = [f"Vacuum field from MEPHIT\n({name_cmplx[part]} {name_coord[coord]} component)",
                     f"Vacuum field from GPEC\n({name_cmplx[part]} {name_coord[coord]} component)"]
            testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac[part], label, title=title, clim_scale=scaling[coord],
                                                      filename=path.join(workdir, f'debug_Bnvac_{coord}_{part}.pdf')))

    testcase.dump_plots()
    del testcase, rootgrp
