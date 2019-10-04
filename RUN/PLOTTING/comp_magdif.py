#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:32:15 2019

@author: lainer_p
"""

import os
import magdifplot

testcase = magdifplot.magdif(
    '/temp/lainer_p/NEO-EQ/magdif',
    'magdif.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.read_configfile()
testcase.read_fluxvar()
testcase.load_mesh()
testcase.plots.append(magdifplot.magdif_1d_cutplot(
        testcase.rho, r'$r$ / cm', testcase.psi, r'$\psi$ / Mx',
        'disc poloidal flux', os.path.join(testcase.datadir, 'plot_psi.pdf')
))
testcase.plots.append(magdifplot.magdif_1d_cutplot(
        testcase.rho, r'$r$ / cm', testcase.q, r'$q$',
        'safety factor', os.path.join(testcase.datadir, 'plot_q.pdf')
))
testcase.plots.append(magdifplot.magdif_1d_cutplot(
        testcase.rho, r'$r$ / cm', testcase.pres0,
        r'$p_{0}$ / dyn cm\textsuperscript{-2}',
        'pressure', os.path.join(testcase.datadir, 'plot_pres0.pdf')
))
if (os.path.isfile(os.path.join(testcase.datadir, 'Bmn_psi.dat')) and
        os.path.isfile(os.path.join(testcase.datadir, 'Bmn_vac_psi.dat'))):
    testcase.plots.append(magdifplot.magdif_poloidal_modes(
            testcase.config['n'], testcase.s, testcase.q,
            testcase.datadir, 'Bmn_psi.dat',
            r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
            'Bmn_vac_psi.dat'
    ))
if os.path.isfile(os.path.join(testcase.datadir, 'currmn_000_theta.dat')):
    testcase.plots.append(magdifplot.magdif_poloidal_modes(
            testcase.config['n'], testcase.s, testcase.q,
            testcase.datadir, 'currmn_000_theta.dat',
            r'$\left\vert J_{mn \theta}^{(0)} \right\vert$'
            + r' / statA cm\textsuperscript{-1}'
    ))
testcase.dump_plots_parallel()
