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
import colorcet as cc
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
n = mephit.data['/config/n'][()]
m_res_min = mephit.data['/mesh/m_res_min'][()]
m_res_max = mephit.data['/mesh/m_res_max'][()]
q = mephit.data['/cache/fs/q'][()]
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
for polmode in Bmn:
    polmode['m_range'] = np.arange(min(polmode['var'].keys()), max(polmode['var'].keys()) + 1)
    polmode['arr'] = np.zeros((polmode['m_range'].size, polmode['rho'][0].size), dtype=complex)
    for i, m in enumerate(polmode['m_range']):
        polmode['arr'][i, :] = polmode['var'][m]


# %%
mephit.close_datafile()
gpec.close_datafiles()

# %%
fig = plt.figure(figsize=(6.6, 3.6), dpi=150)
axs = fig.subplots(1, 2, sharex='all', sharey='all')
ims = [None, None]
for i in range(2):
    # ims[i] = axs[i].pcolormesh(*np.meshgrid(Bmn[i]['m_range'], Bmn[i]['rho'][0], indexing='ij'),
    #                             np.abs(Bmn[i]['arr']), shading='nearest', cmap=cc.m_fire_r)
    ims[i] = axs[i].contourf(*np.meshgrid(Bmn[i]['m_range'], Bmn[i]['rho'][0], indexing='ij'),
                              np.abs(Bmn[i]['arr']), levels=256, cmap=cc.m_fire_r)
    for m in range(m_res_min, m_res_max):  # ignore 'last' resonance
        legend_handle = axs[i].plot(m * mephit.post['sgn_m_res'],
                                    mephit.post['psi_norm'][np.argmin(np.abs((np.abs(q) * n - m)))],
                                    '+', color='tab:cyan', label='rational surfaces ($n = 2$)')
    axs[i].set_xlabel('poloidal mode number $m$')
axs[0].set_ylabel(r'normalized poloidal flux $\hat{\psi}$')
axs[0].set_title('vacuum perturbation')
axs[1].set_title('with plasma response')
cmax = max([im.get_clim()[1] for im in ims])
ims[0].set_clim([0.0, cmax])
ims[1].set_clim([0.0, cmax])
order_of_magnitude = 10.0 ** np.floor(np.log10(cmax))
ticks = np.arange(0.0, cmax, order_of_magnitude)
cbar = fig.colorbar(ims[1], ticks=ticks)
cbar.set_label(r'$\lvert (\sqrt{g} B^{\psi}_{n})_{m} \rvert A^{-1}$ / \si{\tesla}', rotation=90)
axs[0].legend(handles=legend_handle, fontsize='small', loc='lower left')
fig.savefig(path.join(work_dir, 'Bn_spectrum.png'))
plt.show()

# %%
for m in mephit.post['m_res']:
    if abs(m) > 12:
        break
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5, ls=':')
    for polmode in Bmn:
        mask = (polmode['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
            (polmode['rho'][m] <= res[k] + 4.5 * delta_mn[k])
        ax.plot(polmode['rho'][m], np.abs(polmode['var'][m]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi}_{n})_{m} \rvert A^{-1}$ / \si{\tesla}')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    fig.savefig(path.join(work_dir, f'Bmn_psi_{abs(m)}.pdf'), backend='pgf', dpi=150)
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
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi}_{n})_{m} \rvert A^{-1}$ / \si{\tesla}')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    fig.savefig(path.join(work_dir, f'Bmn_psi_zoom_axis_{abs(m)}.pdf'), backend='pgf', dpi=150)
    plt.show()
