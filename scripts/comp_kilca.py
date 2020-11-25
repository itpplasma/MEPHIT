#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 23:04:25 2020

@author: patrick
"""

import magdifplot

workdir = '/home/patrick/git/NEO-EQ/run/geomint_TCFP'
kilcafile = 'TCFP_flre_hic.hdf5'

testcase = magdifplot.magdif(workdir)
testcase.read_configfile()
testcase.read_datafile()

vac = magdifplot.magdif_vec_polmodes(
    testcase.data, '/postprocess/Bmn_vac/coeff_rad',
    magdifplot.fslabel.r, 'vacuum perturbation')
vac.process()
pert = magdifplot.magdif_vec_polmodes(
    testcase.data, '/postprocess/Bmn/coeff_rad',
    magdifplot.fslabel.r, 'full perturbation (MEPHIT)')
pert.process()
kilca = magdifplot.magdif_KiLCA_Bmn_r(workdir, kilcafile,
                                      'full perturbation (KiLCA)')
kilca.process()

testcase.plots.append(magdifplot.magdif_poloidal_plots(
    workdir, 'KiLCA_Bmn_r.pdf', testcase.config, testcase.data, magdifplot.fslabel.r,
    r'perpendicular perturbation field $\lvert B_{mn}^{r} \rvert$ / \si{\gauss}',
    pert, kilca
))

testcase.dump_plots()
