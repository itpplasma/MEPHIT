# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     language: python
#     name: python3
# ---

# %%
from os import environ, path
from numpy import abs, pi, sign
from scipy.constants import c as clight
import matplotlib.pyplot as plt
from mephit_plot import Mephit, Gpec, set_matplotlib_defaults


# %%
set_matplotlib_defaults()
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')
testcase = Mephit(work_dir, 'mephit.h5')
testcase.open_datafile()
testcase.postprocess()
reference = Gpec(testcase.work_dir, 2)
reference.open_datafiles()

# %%
m_res_min = testcase.data['/mesh/m_res_min'][()]
kf_res = testcase.data['/mesh/res_ind'][()]
nflux = testcase.data['/mesh/nflux'][()]
res = testcase.normalize_psi(testcase.data['/mesh/psi_res'][()])
delta_mn = testcase.normalize_psi_diff(testcase.data['/mesh/delta_psi_mn'][()])
q = testcase.data['/cache/fs_half/q'][()]
sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])
conversion = 4.0 * pi / clight
jmnpar_Bmod = [
    testcase.get_polmodes('total', '/debug_KiLCA/jmnpar_Bmod_total/coeff', conversion, L1=True),
    testcase.get_polmodes('gyrokinetic', '/debug_KiLCA/jmnpar_Bmod_KiLCA/coeff', conversion, L1=True),
    testcase.get_polmodes('iMHD', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversion, L1=True),
    testcase.get_polmodes('iMHD (undamped)', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
]
conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][()]
Bmn = [
    testcase.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    testcase.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
    reference.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi'),
    reference.get_polmodes('GPEC vacuum perturbation', sgn_dpsi, 'Jbgradpsi_x'),
]

# %%
for m in testcase.post['m_res']:
    if abs(m) > 7:
        continue
    k = abs(m) - m_res_min
    fig = plt.figure(layout='constrained', dpi=150)
    ax = fig.subplots()
    mask = (jmnpar_Bmod[0]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
        (jmnpar_Bmod[0]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
    ax.axvline(res[k], color='k', lw=0.5)
    ax.axvline(res[k] - 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    ax.axvline(res[k] + 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    for polmode in jmnpar_Bmod:
        ax.semilogy(polmode['rho'][m][mask], abs(polmode['var'][m][mask]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\mu_{0} \lvert (j_{\parallel} / B_{0})_{\V{m}} \rvert$ [\si{\per\meter}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

    fig = plt.figure(layout='constrained', dpi=150)
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5)
    for polmode in Bmn:
        mask = (polmode['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
            (polmode['rho'][m] <= res[k] + 4.5 * delta_mn[k])
        ax.plot(polmode['rho'][m], abs(polmode['var'][m]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\V{m}} \rvert A^{-1}$ [\si{\tesla}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

