#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

import magdifplot

testcase = magdifplot.magdif(
    '/afs/itp.tugraz.at/user/lainer_p/git/NEO-EQ/MHD',
    'magdif.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.read_configfile()
testcase.read_fluxvar()
testcase.calc_resonances()

pert = magdifplot.magdif_GPEC_bnormal('/temp/lainer_p/NEO-EQ/GPEC',
                                      '/temp/ulbl_p/GPEC/30835_3200_EQH/gpec_xbnormal_n2.out',
                                      7, 'full perturbation')
pert.process()
vac = magdifplot.magdif_GPEC_bnormal('/temp/lainer_p/NEO-EQ/GPEC',
                                     '/temp/ulbl_p/GPEC/30835_3200_EQH/gpec_vbnormal_n2.out',
                                     6, 'vacuum perturbation')
vac.process()
testcase.plots.append(magdifplot.magdif_poloidal_plots(
        testcase.config['n'], testcase.psi_norm, testcase.q, testcase.resonance,
        r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
        pert, vac
))

testcase.dump_plots()
