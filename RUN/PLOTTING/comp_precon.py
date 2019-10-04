#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:22:56 2019

@author: lainer_p
"""

import numpy as np
import magdifplot

testcase = magdifplot.magdif(
    '/temp/lainer_p/NEO-EQ/test_res_precon_lowdens',
    'test_res_precon_lowdens.in',
    '../PRELOAD/inputformaxwell.msh'
)
testcase.read_configfile()
testcase.load_mesh()
testcase.plots.append(magdifplot.magdif_conv_plot(
        testcase.datadir, 'convergence.dat', testcase.config['tol'])
)
magdifplot.magdif_conv_plot(
        '/temp/lainer_p/NEO-EQ/test_res_direct_lowdens',
        'convergence.dat',
        np.hypot(-6.5150263298574818E-002, -1.4910530480484588E-002)
).dump_plot()
vacuum = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                    + 'test_res_precon_highdens/plot_Bn_vac.dat')
precon = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                    + 'test_res_precon_lowdens/plot_Bn.dat')
direct = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                    + 'test_res_direct_lowdens/plot_Bn.dat')
highdens = np.loadtxt('/temp/lainer_p/NEO-EQ/'
                      + 'test_res_precon_highdens/plot_Bn.dat')
vacuum = vacuum[:, 4] + 1j * vacuum[:, 5]
precon = precon[:, 4] + 1j * precon[:, 5]
direct = direct[:, 4] + 1j * direct[:, 5]
highdens = highdens[:, 4] + 1j * highdens[:, 5]
testcase.plots.append(magdifplot.magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(vacuum),
    title=r'$\mathrm{Re} \: B_{\mathrm{v} n}^{Z}$ / G',
    filename='/temp/lainer_p/NEO-EQ/Bn_vac_Z_Re.png'
))
testcase.plots.append(magdifplot.magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(precon - vacuum),
    title=r'$\mathrm{Re} \:  B_{\mathrm{p} n}^{Z}$ / G',
    filename='/temp/lainer_p/NEO-EQ/Bn_Z_Re.png'
))
testcase.plots.append(magdifplot.magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(precon - direct),
    title=r'$\mathrm{Re} \: \Delta B_{n}^{Z}$ / G',
    filename='/temp/lainer_p/NEO-EQ/Bn_diff_Z_Re.png'
))
testcase.plots.append(magdifplot.magdif_2d_triplot(
    node=testcase.node, tri=testcase.tri,
    data=np.real(highdens - vacuum),
    title=r'$\mathrm{Re} \: B_{\mathrm{p} n}^{Z}$ / G',
    filename='/temp/lainer_p/NEO-EQ/Bn_highdens_Z_Re.png'
))
testcase.dump_plots_parallel()
