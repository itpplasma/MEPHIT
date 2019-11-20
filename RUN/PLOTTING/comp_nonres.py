#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:29:38 2019

@author: lainer_p
"""

import numpy as np
from magdifplot import magdif, magdif_2d_triplot, magdif_mnDat, fslabel, magdif_poloidal_plots

testcase = magdif(
    '/temp/lainer_p/NEO-EQ/test_nonres_precon_lowdens',
    'test_nonres_precon_lowdens.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.load_mesh()
presn = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                   + 'test_nonres_precon_lowdens/presn.dat')
presn = presn[:, 0] + 1j * presn[:, 1]
testcase.plots.append(magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(presn),
    title=r'$\mathrm{Re} \: p_{n}$ / dyn cm\textsuperscript{-2}',
    filename='/temp/lainer_p/NEO-EQ/nonres_presn_Re.png'
))
testcase.plots.append(magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.imag(presn),
    title=r'$\mathrm{Im} \: p_{n}$ / dyn cm\textsuperscript{-2}',
    filename='/temp/lainer_p/NEO-EQ/nonres_presn_Im.png'
))
presn_000 = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                       + 'test_nonres_precon_lowdens/presn_000.dat')
presn_000 = presn_000[:, 0] + 1j * presn_000[:, 1]
testcase.plots.append(magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(presn_000),
    title=r'$\mathrm{Re} \: p_{n}^{(0)}$ / dyn cm\textsuperscript{-2}',
    filename='/temp/lainer_p/NEO-EQ/nonres_presn_000_Re.png'
))
testcase.plots.append(magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.imag(presn_000),
    title=r'$\mathrm{Im} \: p_{n}^{(0)}$ / dyn cm\textsuperscript{-2}',
    filename='/temp/lainer_p/NEO-EQ/nonres_presn_000_Im.png'
))
testcase.dump_plots_parallel()

testcase = magdif(
    '/temp/lainer_p/NEO-EQ/test_nonres_precon_highdens',
    'test_nonres_precon_highdens.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.read_configfile()
testcase.read_fluxvar()
testcase.load_mesh()
testcase.calc_resonances()
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
