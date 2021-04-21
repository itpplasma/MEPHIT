#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:34:35 2020

@author: patrick
"""

from magdifplot import magdif, statA_per_cm2_to_A_per_m2
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from os import path
from glob import iglob

dyn_per_cm3_to_N_per_m3 = 10.0
latex_coord = {'R': 'R', 'Z': 'Z', 'phi': r'\varphi'}
for work_dir in iglob("/temp/lainer_p/git/NEO-EQ/run/g*"):

    testcase = magdif(work_dir)
    testcase.read_datafile()
    psi = testcase.data['/debug_equil/psi'][()]
    kf = psi.size // 2
    theta = testcase.data['/debug_equil/theta'][()]

    fig = Figure()
    ax = fig.subplots()
    ax.plot(theta, testcase.data['/debug_equil/grad_p0'][kf, :] * dyn_per_cm3_to_N_per_m3, label=r'$\nabla p_{0}$')
    ax.plot(theta, testcase.data['/debug_equil/Ampere_Lorentz'][kf, :] * dyn_per_cm3_to_N_per_m3, label='Ampere')
    ax.plot(theta, testcase.data['/debug_equil/GS_Lorentz'][kf, :] * dyn_per_cm3_to_N_per_m3, label='Grad--Shafranov')
    ax.legend()
    ax.set_xlabel(r'$\theta$ / \si{rad}')
    ax.set_ylabel(r'$f$ / \si{\newton\per\meter\cubed')
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, "check_iMHD.pdf"))

    for coord, symbol in latex_coord.items():
        fig = Figure()
        ax = fig.subplots()
        ax.plot(theta, testcase.data[f"/debug_equil/GS_j0_{coord}"][kf, :] * statA_per_cm2_to_A_per_m2,
                label='Grad--Shafranov')
        ax.plot(theta, testcase.data[f"/debug_equil/Ampere_j0_{coord}"][kf, :] * statA_per_cm2_to_A_per_m2,
                label='Ampere')
        ax.legend()
        ax.set_xlabel(r'$\theta$ / \si{\degree}')
        ax.set_ylabel(fr'$J_{{0 {symbol}}}$ / \si{{\ampere\per\meter\squared}}')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(work_dir, f"check_j0{coord}.pdf"))
