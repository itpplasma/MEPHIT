#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

import magdifplot
from os import path

workdir = '/temp/lainer_p/NEO-EQ/GPEC'
ncfile = '/temp/ulbl_p/GPEC/30835_3200_EQH_Boozer/gpec_profile_output_n2.nc'

testcase = magdifplot.magdif(workdir, path.join(workdir, 'magdif.in'),
                             'inputformaxwell.msh')
testcase.read_configfile()
testcase.read_fluxvar()
testcase.calc_resonances()

pert = magdifplot.magdif_GPEC_bnormal(workdir, ncfile,
                                      'b_n', 'full perturbation')
pert.process()
vac = magdifplot.magdif_GPEC_bnormal(workdir, ncfile,
                                     'b_n_x', 'vacuum perturbation')
vac.process()
testcase.plots.append(magdifplot.magdif_poloidal_plots(
        testcase.config['n'], testcase.psi_norm, testcase.q, testcase.resonance,
        r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
        pert, vac
))

testcase.dump_plots()
