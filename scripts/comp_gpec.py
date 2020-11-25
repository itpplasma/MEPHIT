#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:50:35 2020

@author: lainer_p
"""

import magdifplot

workdir = '/home/patrick/git/NEO-EQ/run/33353_2325'
ncfile = 'gpec_profile_output_n2.nc'

testcase = magdifplot.magdif(workdir)
testcase.read_configfile()
testcase.read_datafile()

pert = magdifplot.magdif_GPEC_bnormal(workdir, ncfile,
                                      'b_n', 'full perturbation')
pert.process()
vac = magdifplot.magdif_GPEC_bnormal(workdir, ncfile,
                                     'b_n_x', 'vacuum perturbation')
vac.process()
testcase.plots.append(magdifplot.magdif_poloidal_plots(
    workdir, 'GPEC_Bmn_psi.pdf', testcase.config, testcase.data,
    magdifplot.fslabel.psi_norm, poldata=pert, refdata=vac,
    ylabel=r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\maxwell}'))

testcase.dump_plots()
