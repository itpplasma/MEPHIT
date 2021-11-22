#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:04:25 2020

@author: patrick
"""

from os import path
import numpy as np
from magdifplot import run_dir, magdif, fslabel, polmodes, magdif_poloidal_plots

workdir = run_dir + '/geomint_TCFP'
kilcafile = 'TCFP_flre_hic.hdf5'
ylabel = r'perpendicular perturbation field $\lvert B_{mn}^{r} \rvert$ / \si{\gauss}'

testcase = magdif(workdir)
testcase.read_datafile()

pert = polmodes('full perturbation (MEPHIT)', 'k-')
pert.read_magdif(testcase.data, fslabel.r, '/postprocess/Bmn/coeff_rad')
vac = polmodes('vacuum perturbation', 'r--')
vac.read_magdif(testcase.data, fslabel.r, '/postprocess/Bmn_vac/coeff_rad')
kilca = polmodes('full perturbation (KiLCA)', 'b-.')
kilca.read_KiLCA(path.join(workdir, kilcafile), 'Br')

testcase.plots.append(magdif_poloidal_plots(
    workdir, 'KiLCA_Bmn_r.pdf', testcase.data,
    fslabel.r, ylabel, np.abs, vac, pert, kilca
))

testcase.dump_plots()
