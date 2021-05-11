#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from magdifplot import magdif, statA_per_cm2_to_A_per_m2
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from os import path
from numpy import clip

dyn_per_cm3_to_N_per_m3 = 10.0
latex_coord = {'R': 'R', 'Z': 'Z', 'phi': r'\varphi'}
work_dir = "/temp/lainer_p/git/NEO-EQ/run/g33353.2325"

equil_base = magdif(work_dir)
equil_base.read_datafile()
equil_IDE = magdif(work_dir + '_AUGD_IDE_ed1')
equil_IDE.read_datafile()
equil_EQH = magdif(work_dir + '_AUGD_EQH_ed1')
equil_EQH.read_datafile()
equil_EQB = magdif(work_dir + '_MICDU_EQB_ed3')
equil_EQB.read_datafile()

kf_base = equil_base.data['/debug_equil/psi'].size // 2
kf_IDE = equil_IDE.data['/debug_equil/psi'].size // 2
kf_EQH = equil_EQH.data['/debug_equil/psi'].size // 2
kf_EQB = equil_EQB.data['/debug_equil/psi'].size // 2
print(f"base: psi = {equil_base.data['/debug_equil/psi'][kf_base]}")
print(f"IDE:  psi = {equil_IDE.data['/debug_equil/psi'][kf_IDE]}")
print(f"EQH:  psi = {equil_EQH.data['/debug_equil/psi'][kf_EQH]}")
print(f"EQB:  psi = {equil_EQB.data['/debug_equil/psi'][kf_EQB]}")

theta = equil_base.data['/debug_equil/theta'][()]
fig = Figure(figsize=(6.5, 5.2))
axs = fig.subplots(2, 2, sharex='all', sharey='all')
axs[0, 0].plot(theta, 2.0 * equil_IDE.data['/debug_equil/Ampere_Lorentz'][kf_IDE, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Ampere')
axs[0, 0].plot(theta, 2.0 * equil_IDE.data['/debug_equil/GS_Lorentz'][kf_IDE, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Grad--Shafranov')
axs[0, 1].plot(theta, equil_EQH.data['/debug_equil/Ampere_Lorentz'][kf_EQH, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Ampere')
axs[0, 1].plot(theta, equil_EQH.data['/debug_equil/GS_Lorentz'][kf_EQH, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Grad--Shafranov')
axs[1, 0].plot(theta, equil_EQB.data['/debug_equil/Ampere_Lorentz'][kf_EQB, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Ampere')
axs[1, 0].plot(theta, equil_EQB.data['/debug_equil/GS_Lorentz'][kf_EQB, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Grad--Shafranov')
axs[1, 1].plot(theta, -equil_base.data['/debug_equil/Ampere_Lorentz'][kf_base, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Ampere')
axs[1, 1].plot(theta, -equil_base.data['/debug_equil/GS_Lorentz'][kf_base, :] * dyn_per_cm3_to_N_per_m3,
               label=r'Grad--Shafranov')
axs[0, 0].set_title(r'g33353.2325\_AUGD\_IDE\_ed1', fontfamily='monospace', fontsize='small')
axs[0, 1].set_title(r'g33353.2325\_AUGD\_EQH\_ed1', fontfamily='monospace', fontsize='small')
axs[1, 0].set_title(r'g33353.2325\_MICDU\_EQB\_ed3', fontfamily='monospace', fontsize='small')
axs[1, 1].set_title('g33353.2325', fontfamily='monospace', fontsize='small')
for ax in axs.flatten().tolist():
    ax.legend()
    ax.set_xlabel(r'symmetry flux poloidal angle $\vartheta$ / \si{rad}')
    ax.set_ylabel(r'force density $f$ / \si{\newton\per\meter\cubed')
for ax in axs[:, 1:].flatten().tolist():
    ax.yaxis.set_tick_params(labelleft=True)
    ax.yaxis.offsetText.set_visible(True)
for ax in axs[:-1, :].flatten().tolist():
    ax.xaxis.set_tick_params(labelbottom=True)
for ax in axs.flatten().tolist():
    ax.yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, "comp_iMHD.pdf"))

q_max = 6.0
fig = Figure(figsize=(6.5, 5.2))
axs = fig.subplots(2, 2, sharex='all', sharey='all')
axs[0, 0].plot(equil_IDE.data['/debug_q/GEQDSK_psi_norm'][()],
               clip(equil_IDE.data['/debug_q/GEQDSK_q'][()], 0.0, q_max),
               label=r'GEQDSK')
axs[0, 0].plot(equil_IDE.data['/debug_q/RK_psi_norm'][()], equil_IDE.data['/debug_q/RK_q'][()],
               label=r'Runge--Kutta')
axs[0, 1].plot(equil_EQH.data['/debug_q/GEQDSK_psi_norm'][()],
               clip(equil_EQH.data['/debug_q/GEQDSK_q'][()], 0.0, q_max),
               label=r'GEQDSK')
axs[0, 1].plot(equil_EQH.data['/debug_q/RK_psi_norm'][()], equil_EQH.data['/debug_q/RK_q'][()],
               label=r'Runge--Kutta')
axs[1, 0].plot(equil_EQB.data['/debug_q/GEQDSK_psi_norm'][()],
               clip(equil_EQB.data['/debug_q/GEQDSK_q'][()], 0.0, q_max),
               label=r'GEQDSK')
axs[1, 0].plot(equil_EQB.data['/debug_q/RK_psi_norm'][()], equil_EQB.data['/debug_q/RK_q'][()],
               label=r'Runge--Kutta')
axs[1, 1].plot(equil_base.data['/debug_q/GEQDSK_psi_norm'][()],
               clip(-equil_base.data['/debug_q/GEQDSK_q'][()], 0.0, q_max),
               label=r'GEQDSK')
axs[1, 1].plot(equil_base.data['/debug_q/RK_psi_norm'][()], -equil_base.data['/debug_q/RK_q'][()],
               label=r'Runge--Kutta')
axs[0, 0].set_title(r'g33353.2325\_AUGD\_IDE\_ed1', fontfamily='monospace', fontsize='small')
axs[0, 1].set_title(r'g33353.2325\_AUGD\_EQH\_ed1', fontfamily='monospace', fontsize='small')
axs[1, 0].set_title(r'g33353.2325\_MICDU\_EQB\_ed3', fontfamily='monospace', fontsize='small')
axs[1, 1].set_title('g33353.2325', fontfamily='monospace', fontsize='small')
for ax in axs.flatten().tolist():
    ax.legend()
    ax.set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
    ax.set_ylabel(r'safety factor $q$')
for ax in axs[:, 1:].flatten().tolist():
    ax.yaxis.set_tick_params(labelleft=True)
for ax in axs[:-1, :].flatten().tolist():
    ax.xaxis.set_tick_params(labelbottom=True)
for ax in axs.flatten().tolist():
    ax.yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, "comp_q.pdf"))
