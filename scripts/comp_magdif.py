#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:32:15 2019

@author: lainer_p
"""

import os
from magdifplot import magdif, magdif_1d_cutplot, magdif_mnDat, fslabel, magdif_poloidal_plots

testcase = magdif(
    '/temp/lainer_p/NEO-EQ/magdif',
    'magdif.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.read_configfile()
testcase.read_fluxvar()
testcase.load_mesh()
testcase.calc_resonances()
testcase.plots.append(magdif_1d_cutplot(
        testcase.r, r'$r$ / cm', testcase.psi, r'$\psi$ / Mx',
        'disc poloidal flux', os.path.join(testcase.datadir, 'plot_psi.pdf')
))
testcase.plots.append(magdif_1d_cutplot(
        testcase.r, r'$r$ / cm', testcase.q, r'$q$',
        'safety factor', os.path.join(testcase.datadir, 'plot_q.pdf')
))
testcase.plots.append(magdif_1d_cutplot(
        testcase.r, r'$r$ / cm', testcase.pres0,
        r'$p_{0}$ / dyn cm\textsuperscript{-2}',
        'pressure', os.path.join(testcase.datadir, 'plot_pres0.pdf')
))
pert = magdif_mnDat(testcase.datadir, 'Bmn_psi.dat', 0,
                    fslabel.psi_norm, 'full perturbation')
pert.process()
vac = magdif_mnDat(testcase.datadir, 'Bmn_vac_psi.dat', 0,
                   fslabel.psi_norm, 'vacuum perturbation')
vac.process()
testcase.plots.append(magdif_poloidal_plots(
        testcase.config['n'], testcase.psi_norm, testcase.q,
        testcase.resonance,
        r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
        pert, vac
))
pert = magdif_mnDat(testcase.datadir, 'currmn_000_theta.dat', 0,
                    fslabel.psi_norm, 'initial perturbation')
pert.process()
testcase.plots.append(magdif_poloidal_plots(
        testcase.config['n'], testcase.psi_norm, testcase.q,
        testcase.resonance,
        r'$\left\vert J_{mn \theta}^{(0)} \right\vert$'
        + r' / statA cm\textsuperscript{-1}', pert
))
testcase.dump_plots_parallel()
