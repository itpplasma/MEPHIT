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
from magdifplot import magdif, fslabel, polmodes, magdif_poloidal_plots, magdif_2d_rectplots
import netCDF4

ylabel_abs = r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\maxwell}'
ylabel_arg = r'$\arg \sqrt{g} B_{mn}^{\psi}$ / \si{\degree}'

for workdir in iglob('/temp/lainer_p/git/NEO-EQ/run/Bvac_*'):
    testcase = magdif(workdir)
    testcase.read_configfile()
    testcase.read_datafile()
    sgn_Itor = -np.sign(testcase.data['/cache/fs/psi'][-1] -
                        testcase.data['/cache/fs/psi'][0])

    vac = polmodes('vacuum perturbation (GPEC)', 'r--')
    vac.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_Itor, 'Jbgradpsi_x')
    pert = polmodes('full perturbation (GPEC)', 'k--')
    pert.read_GPEC(path.join(workdir, 'gpec_profile_output_n2.nc'), sgn_Itor, 'Jbgradpsi')
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
    R = [testcase.data['/Bnvac/rect_R'][()], np.array(rootgrp['R']) * 1.0e+02]
    Z = [testcase.data['/Bnvac/rect_Z'][()], np.array(rootgrp['z']) * 1.0e+02]
    Bnvac_R_Re = [testcase.data['/Bnvac/rect_comp_R'][()].real, np.array(rootgrp['b_r'][0, :, :]) * 1.0e+04]
    Bnvac_R_Im = [testcase.data['/Bnvac/rect_comp_R'][()].imag, np.array(rootgrp['b_r'][1, :, :]) * 1.0e+04]
    label_Re = r'$\Real B_{n \mathrm{v}}^{R}$ / \si{\gauss}'
    label_Im = r'$\Imag B_{n \mathrm{v}}^{R}$ / \si{\gauss}'
    title_Re = ["Vacuum field from MEPHIT\n(real radial component)",
                "Vacuum field from GPEC\n(real radial component)"]
    title_Im = ["Vacuum field from MEPHIT\n(imaginary radial component)",
                "Vacuum field from GPEC\n(imaginary radial component)"]
    testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac_R_Re, label_Re, title=title_Re, clim_scale=(0.1, 0.1),
                                              filename=path.join(workdir, 'debug_Bnvac_R_Re.pdf')))
    testcase.plots.append(magdif_2d_rectplots(R, Z, Bnvac_R_Im, label_Im, title=title_Im, clim_scale=(0.1, 0.1),
                                              filename=path.join(workdir, 'debug_Bnvac_R_Im.pdf')))

    testcase.dump_plots()
    del testcase, rootgrp
