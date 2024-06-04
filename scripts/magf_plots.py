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
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %%
# %matplotlib inline
from os import environ, getcwd, path
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import h5py
from mephit_plot import Mephit, Gpec, Mars
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')


# %%
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')  # AUG
# work_dir = path.join(environ['MEPHIT_DIR'], 'run/47046')  # MAST-U
mephit = Mephit(work_dir, 'mephit.h5')
mephit.open_datafile()
mephit.postprocess()
gpec = Gpec(mephit.work_dir, mephit.data['/config/n'][()])
gpec.open_datafiles()
# ref_dir = path.join(environ['HOME'], 'TU/PhD/MARS_MEPHIT/forPatrick')
# mars = Mars(ref_dir)
# mars.open_datafiles()

# %%
m_res_min = mephit.data['/mesh/m_res_min'][()]
res = mephit.normalize_psi(mephit.data['/mesh/psi_res'][()])
delta_mn = mephit.normalize_psi_diff(mephit.data['/mesh/delta_psi_mn'][()])
sgn_dpsi = np.sign(mephit.data['/cache/fs/psi'][-1] - mephit.data['/cache/fs/psi'][0])
conversion = 1.0e-04 / mephit.data['/mesh/gpec_jacfac'][()]  # TODO: jacfac not accounted for in MARS output
Bmn = [
    mephit.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
    gpec.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi'),
    gpec.get_polmodes('GPEC vacuum perturbation', sgn_dpsi, 'Jbgradpsi_x'),
    # mars.get_polmodes('MARS vacuum perturbation', 'VACUUM'),
    # mars.get_polmodes('MARS full perturbation', 'PLASMA'),
]

# %%
for m in mephit.post['m_res']:
    if abs(m) > 12:
        break
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5)
    for polmode in Bmn:
        mask = (polmode['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
            (polmode['rho'][m] <= res[k] + 4.5 * delta_mn[k])
        ax.plot(polmode['rho'][m], np.abs(polmode['var'][m]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ [\si{\tesla}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

# %%
for m in mephit.post['sgn_m_res'] * np.arange(6):
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    for polmode in Bmn:
        grad = np.full(polmode['var'][m].shape, np.nan, dtype='D')
        grad.real = np.gradient(polmode['var'][m].real, polmode['rho'][m])
        grad.imag = np.gradient(polmode['var'][m].imag, polmode['rho'][m])
        mask = (polmode['rho'][m] >= 0.001) & (polmode['rho'][m] <= 0.1)
        ax.plot(polmode['rho'][m][mask], np.abs(polmode['var'][m][mask]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ [\si{\tesla}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

# %%
mephit.close_datafile()
gpec.close_datafiles()
