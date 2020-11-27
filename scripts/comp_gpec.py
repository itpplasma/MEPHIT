#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

from os import path
import numpy as np
from magdifplot import magdif, fslabel, polmodes, magdif_poloidal_plots

workdir = '/home/patrick/git/NEO-EQ/run/30835_3200_ed6'
ncfile = 'gpec_profile_output_n2.nc'
ylabel = r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\maxwell}'

testcase = magdif(workdir)
testcase.read_configfile()
testcase.read_datafile()

vac = polmodes('vacuum perturbation (GPEC)', 'r--')
vac.read_GPEC(path.join(workdir, ncfile), 'b_n_x')
pert = polmodes('full perturbation (GPEC)', 'k--')
pert.read_GPEC(path.join(workdir, ncfile), 'b_n')
mephit_vac = polmodes('vacuum perturbation (MEPHIT)', 'b--')
mephit_vac.read_magdif(testcase.data, fslabel.psi_norm,
                       '/postprocess/Bmn_vac/coeff_rad')
mephit_pert = polmodes('full perturbation (MEPHIT)', 'g--')
mephit_pert.read_magdif(testcase.data, fslabel.psi_norm,
                        '/postprocess/Bmn/coeff_rad')
testcase.plots.append(magdif_poloidal_plots(
    workdir, 'GPEC_Bmn_psi.pdf', testcase.config, testcase.data,
    fslabel.psi_norm, ylabel, np.abs, vac, pert, mephit_vac, mephit_pert
))

testcase.dump_plots()
