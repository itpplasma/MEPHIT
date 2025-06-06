#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import clip
from mephit_plot import run_dir, set_matplotlib_defaults, Mephit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
from os import path
from glob import iglob

set_matplotlib_defaults()
q_max = 10.0
for work_dir in iglob(run_dir + '/GSE/g3*'):
    try:
        testcase = Mephit(work_dir)
        testcase.read_datafile()
        fig = Figure()
        ax = fig.subplots()
        ax.plot(testcase.data['/debug_q/GEQDSK_psi_norm'][()],
                clip(testcase.data['/debug_q/GEQDSK_q'][()], -q_max, q_max), label='gfile')
        ax.plot(testcase.data['/debug_q/RK_psi_norm'][()],
                testcase.data['/debug_q/RK_q'][()], label='field line')
        ax.step(testcase.data['/debug_q/step_psi_norm'][()],
                testcase.data['/debug_q/step_q'][()], where='mid', label='triangle grid')
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.grid(which='major', axis='y', lw=0.25)
        ax.legend()
        ax.set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
        ax.set_ylabel(r'safety factor $q$')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(work_dir, "check_q.pdf"))
    except:
        print(f"error in processing {work_dir}")
        continue
