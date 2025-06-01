# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# %matplotlib inline
from os import environ, getcwd, path
from matplotlib import rcParams, ticker
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import elementary_charge, c as clight
from scipy.interpolate import CubicSpline
import h5py
from mephit_plot import Mephit, Gpec  # for polmode plots
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')

# %%
work_dir = path.join(environ['MEPHIT_RUN_DIR'], '33353_2900_EQH')

# %%
data = {}
with h5py.File(path.join(work_dir, 'mephit.h5'), 'r') as f:
    data['Delta_E_r'] = f['/resonance_sweep/Delta_E_r'][()]
    data['E_r_zero'] = f['/resonance_sweep/E_r_zero'][()]
    data['Imn_res'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
    E_r_interp = CubicSpline(f['/equil/profiles/E_r/x'][()][::-1],
                             f['/equil/profiles/E_r/y'][()][::-1])
    dens_interp = CubicSpline(f['/equil/profiles/dens_e/x'][()][::-1],
                              f['/equil/profiles/dens_e/y'][()][::-1])
    temp_e_interp = CubicSpline(f['/equil/profiles/temp_e/x'][()][::-1],
                                f['/equil/profiles/temp_e/y'][()][::-1])
    temp_i_interp = CubicSpline(f['/equil/profiles/temp_i/x'][()][::-1],
                                f['/equil/profiles/temp_i/y'][()][::-1])
    nu_e_interp = CubicSpline(f['/equil/profiles/nu_e/x'][()][::-1],
                              f['/equil/profiles/nu_e/y'][()][::-1])
    data['B0'] = np.abs(f['/equil/bcentr'][()])
    data['R0'] = f['/mesh/R_O'][()]
    data['rsmall'] = f['/cache/fs/rsmall'][()]
    data['psi'] = f['/cache/fs/psi'][()]
    data['psi_res'] = f['/mesh/psi_res'][()]
    data['n'] = f['/mesh/n'][()]
    data['m_res_min'] = f['/mesh/m_res_min'][()]
    data['m_res_max'] = f['/mesh/m_res_max'][()]
    sgn_m_res = int(np.sign(-f['/cache/fs/q'][-1]))  # consolidate with mephit_plot.py
    data['Delta_res'] = f['/mesh/Delta_psi_res'][()]  / (data['psi'][-1] - data['psi'][0])
with h5py.File(path.join(work_dir, 'mephit_protium.h5'), 'r') as f:
    data['Imn_res_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight

# %%
data['E_r_res'] = E_r_interp(data['psi_res'])
data['v_ExB'] = -1.0e2 * clight / data['B0'] * (data['E_r_res'] + data['Delta_E_r'][:, np.newaxis])
data['v_ExB_orig'] = -1.0e2 * clight / data['B0'] * data['E_r_res']
data['V_e_diamag'] = 1.0e1 / (-elementary_charge * data['B0'] * dens_interp(data['psi_res'])) * \
    (dens_interp(data['psi_res']) * temp_e_interp(data['psi_res'], nu=1) +
     dens_interp(data['psi_res'], nu=1) * temp_e_interp(data['psi_res']))
data['V_i_diamag'] = 1.0e1 / (elementary_charge * data['B0'] * dens_interp(data['psi_res'])) * \
    (dens_interp(data['psi_res']) * temp_i_interp(data['psi_res'], nu=1) +
     dens_interp(data['psi_res'], nu=1) * temp_i_interp(data['psi_res']))

# %%
m_min = 6
m_max = 7


# %%
def normalize_psi(psi):
    return (psi - data['psi'][0]) / (data['psi'][-1] - data['psi'][0])

psi_norm_to_rsmall = CubicSpline(normalize_psi(data['psi']), data['rsmall'])
rsmall_to_psi_norm = CubicSpline(data['rsmall'], normalize_psi(data['psi']))

# %%
rsmall_res = psi_norm_to_rsmall(normalize_psi(data['psi_res']))
q = np.arange(data['m_res_min'], data['m_res_max'] + 1) / data['n']
k_perp = data['n'] * q / rsmall_res
v_Te = 4.19e7 * np.sqrt(np.abs(temp_e_interp(data['psi_res'])))
nu_e = nu_e_interp(data['psi_res'])
delta_mn = q * data['R0'] * data['v_ExB_orig'] / v_Te * max(1,
    np.sqrt(nu_e / (k_perp * data['v_ExB_orig'])))

# %%
for k in range(data['E_r_res'].size):
    print(f"m = {k + data['m_res_min']:2}, E_r = {data['E_r_res'][k]:+.16e}")
print(f"Delta_E_r(6) at 0: {-data['E_r_res'][3]:+.16e}")
DeltaDelta_E_r = 0.2e5 * -1.0e-2 / clight * data['B0']
print(f"Delta_E_r(6) at +: {-data['E_r_res'][3] + DeltaDelta_E_r:+.16e}")
print(f"Delta_E_r(6) at -: {-data['E_r_res'][3] - DeltaDelta_E_r:+.16e}")

# %%
fig = plt.figure(figsize=(6, 4), dpi=150)
ax = fig.subplots()
# ax.axs[k]hline(0, color='black', linewidth=0.5)
for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    ax.axvline(normalize_psi(data['psi_res'][k]), linewidth=0.5)
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res'][:, k]), label=f'D, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_H'][:, k]), label=f'H, $m = {m}$')
ax.set_xlabel(r'electric zero at $\hat{\psi}$')
ax.set_ylabel(r'$I_{\vec{m} \parallel}$ / \si{\ampere}')
# ax.set_ylim((0.0, 0.03e10))
ax.legend()
plt.show()

# %%
fig = plt.figure(figsize=(6, 4), dpi=150)
ax = fig.subplots()
k_min = 0
k_max = 100
for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res'][k_min:k_max, k]), label=f'D, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_H'][k_min:k_max, k]), label=f'H, $m = {m}$')
ax.set_xlabel(r'electric field offset $\Delta E_r$')
ax.set_ylabel(r'$I_{\vec{m} \parallel}$ / \si{\ampere}')
ax.legend()
plt.show()


# %%
eresoff = np.array([34, 60])
eres = np.array([44, 81])

for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    print(f"Delta_E_r off el.fl.res. {m}: {data['Delta_E_r'][eresoff[m - m_min]]:.16e}")
    print(f"Delta_E_r  on el.fl.res. {m}: {data['Delta_E_r'][eres[m - m_min]]:.16e}")
    print(f"v_ExB     off el.fl.res. {m}: {data['v_ExB'][eresoff[m - m_min], k]:.16e}")
    print(f"v_ExB      on el.fl.res. {m}: {data['v_ExB'][eres[m - m_min], k]:.16e}")


# %%
fig = plt.figure(figsize=(3.6 * (m_max - m_min + 1), 2.9), dpi=150)
axs = fig.subplots(1, m_max - m_min + 1)
for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    # axs[m - m_min].axhline(0.0, color='black', linewidth=0.5)
    axs[m - m_min].axvline(0.0, linewidth=0.5, linestyle=':')
    axs[m - m_min].axvline(data['v_ExB'][eresoff[m - m_min], k], linewidth=0.5, linestyle='-')
    axs[m - m_min].axvline(data['v_ExB'][eres[m - m_min], k], linewidth=0.5, linestyle='--')
    axs[m - m_min].axvline(data['v_ExB_orig'][m - data['m_res_min']], linewidth=0.5, linestyle='-.')
    # axs[m - m_min].axvline(0.5e-6 * data['V_e_diamag'][k], linewidth=0.5, linestyle='--')
    # axs[m - m_min].axvline(0.5e-6 * data['V_i_diamag'][k], linewidth=0.5, linestyle=':')
    mask = data['v_ExB'][:, k] > -data['v_ExB'][0, k]
    axs[m - m_min].semilogy(data['v_ExB'][mask, k],
                            np.abs(data['Imn_res'][mask, k]),
                            color='k', ls='-', label='D')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k],
                            np.abs(data['Imn_res_H'][mask, k]),
                            color='tab:orange', ls='-', label='H')
    xlim = np.array(axs[m - m_min].get_xlim())
    axs[m - m_min].set_xlim(np.array([-1.0, 1.0]) * np.amax(np.abs(xlim)))
    axs[m - m_min].set_xlabel(r'$V_{E \times B}$ / \si{\cm\per\second}')
    axs[m - m_min].set_ylabel(r'$I_{\parallel \vec{m}}$ / a.u.')  # \si{\ampere}
    axs[m - m_min].legend(loc='lower right')
    axs[m - m_min].text(0.975, 0.25, fr"$\vec{{m}} = ({m * sgn_m_res}, 2)$",
                        ha='right', va='bottom', transform=axs[m - m_min].transAxes)
    axs[m - m_min].text(data['v_ExB'][eres[m - m_min], k], 0.95,
                        'on el.fl.res.', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
    axs[m - m_min].text(data['v_ExB'][eresoff[m - m_min], k], 0.95,
                        'off el.fl.res.', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
    axs[m - m_min].text(data['v_ExB_orig'][m - data['m_res_min']], 0.45 if m == 7 else 0.95,
                        r'orig.\ value', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
fig.savefig(path.join(work_dir, 'sweep_resonance.pdf'), backend='pgf')
plt.show()

# %%
mephit = {}
for isotope in ['D', 'H']:
    for shift in ['eresoff', 'eres']:
        for m in [6, 7]:
            sim = f'{isotope}_{shift}_{m}'
            mephit[sim] = Mephit(work_dir, f'mephit_{sim}.h5')
            mephit[sim].open_datafile()
            mephit[sim].postprocess()
mephit['iMHD'] = Mephit(work_dir, 'mephit_iMHD.h5')
mephit['iMHD'].open_datafile()
mephit['iMHD'].postprocess()
gpec = Gpec(work_dir, 2)
gpec.open_datafiles()

# %%
res = mephit['iMHD'].normalize_psi(mephit['iMHD'].data['/mesh/psi_res'][()])
conversion = 1.0e-04 / mephit['iMHD'].data['/mesh/gpec_jacfac'][()]
sgn_dpsi = np.sign(mephit['iMHD'].data['/cache/fs/psi'][-1] - mephit['iMHD'].data['/cache/fs/psi'][0])
Bmn = {'vac': mephit['iMHD'].get_polmodes('vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
       'iMHD': mephit['iMHD'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       'GPEC': gpec.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi')}
Ires = {}
for isotope in ['D', 'H']:
    for shift in ['eresoff', 'eres']:
        for m in [6, 7]:
            sim = f'{isotope}_{shift}_{m}'
            Bmn[sim] = mephit[sim].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion)
            Ires[sim] = mephit[sim].get_Ires()

# %%
for m in [6, 7]:
    for shift in ['eresoff', 'eres']:
        for isotope in ['H', 'D']:
            sim = f'{isotope}_{shift}_{m}'
            print(f"{sim}: {Ires[sim][m * sgn_m_res]}")

# %%
fig = plt.figure(figsize=(7.2, 2.9), dpi=150)
axs = fig.subplots(2, 1, sharey='all')
for k, m in enumerate([6, 7]):
    m_res = sgn_m_res * m
    axs[k].axhline(0.0, color='k', lw=0.5)
    axs[k].axvline(res[m - data['m_res_min']], color='k', lw=0.5)
    axs[k].plot(Bmn['vac']['rho'][m_res][1:], np.abs(Bmn['vac']['var'][m_res][1:]),
                ls=':', color='k', label='vacuum perturbation')  # purple
    axs[k].plot(Bmn['iMHD']['rho'][m_res][1:], np.abs(Bmn['iMHD']['var'][m_res][1:]),
                ls='-.', color='tab:green', label='MEPHIT, iMHD')
    axs[k].plot(Bmn['GPEC']['rho'][m_res], np.abs(Bmn['GPEC']['var'][m_res]),
                ls='-.', color='tab:purple', label='GPEC')
    sim = f'D_eresoff_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='k', label='MEPHIT, D off el.fl.res.')
    sim = f'H_eresoff_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='tab:orange', label='MEPHIT, H off el.fl.res.')
    sim = f'D_eres_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='--', color='k', label='MEPHIT, D on el.fl.res.')
    sim = f'H_eres_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='--', color='tab:orange', label='MEPHIT, H on el.fl.res.')
    axs[k].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    axs[k].yaxis.set_major_locator(ticker.MultipleLocator(0.5e-4))
    axs[k].set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
    axs[k].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
    axs[k].text(0.1, 0.1, fr"$\vec{{m}} = ({m_res}, 2)$",
        ha='center', va='bottom', transform=axs[k].transAxes)
# axs[0].set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ / \si{\tesla}')
axs[0].set_ylabel(r'normal mag. field perturbation / \si{\tesla}')
axs[0].legend(loc='upper left')
fig.savefig(path.join(work_dir, f'Bmn_6_kinetic.pdf'), backend='pgf', dpi=150)
plt.show()

# %%
fig = plt.figure(figsize=(6, 3.6), dpi=150)
axs = fig.subplots(1, 2, sharey='all').ravel()
for k, m in enumerate([6, 7]):
    axs[k].axhline(0.0, color='k', lw=0.5)
    axs[k].axvline(res[m - data['m_res_min']], color='k', lw=0.5, ls=':')
    axs[k].axvline(res[m - data['m_res_min']] - 0.5 *
                   data['Delta_res'][m - data['m_res_min']],
                   color='k', lw=0.5, ls='--')
    axs[k].axvline(res[m - data['m_res_min']] + 0.5 *
                   data['Delta_res'][m - data['m_res_min']],
                   color='k', lw=0.5, ls='--')
    for j, shift in enumerate(['eresoff', 'eres']):
        for i, isotope in enumerate(['D', 'H']):
            sim = f'{isotope}_{shift}_{m}'
            mask = (Bmn[sim]['rho'][m * sgn_m_res] >= res[m - data['m_res_min']] -
                    15 * data['Delta_res'][m - data['m_res_min']]) & \
                   (Bmn[sim]['rho'][m * sgn_m_res] <= res[m - data['m_res_min']] +
                    15 * data['Delta_res'][m - data['m_res_min']])
            axs[k].plot(Bmn[sim]['rho'][m * sgn_m_res][mask],
                        np.abs(Bmn[sim]['var'][m * sgn_m_res][mask]),
                        color=('tab:orange' if i == 1 else 'k'),
                        ls=('--' if j == 1 else '-'),
                        label=isotope + (' on ' if j == 1 else ' off ') + 'el.fl.res.')
    axs[k].set_xlabel(r'normalized ploidal flux $\hat{\psi}$')
    axs[k].set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ / \si{\tesla}')
    axs[k].xaxis.set_major_locator(ticker.MultipleLocator(0.025))
    axs[k].yaxis.set_tick_params(labelleft=True)
    axs[k].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
    axs[k].text(0.2, 0.05, fr"$\vec{{m}} = ({m * sgn_m_res}, 2)$",
                ha='center', va='bottom', transform=axs[k].transAxes)
axs[0].legend(loc='upper left', fontsize='small')
ylims = axs[0].get_ylim()
axs[0].set_ylim((ylims[0], 1.1 * ylims[1]))
fig.savefig(path.join(work_dir, f'zoom_resonance.pdf'), backend='pgf', dpi=150)
plt.show()

# %%
for isotope in ['D', 'H']:
    txt = np.loadtxt(f'/temp/kasilov/MEPHIT-KILCA/VARENNA24/STANDALONE_PERTEFIELD/{isotope}/E_parallel.dat')
    if isotope == 'D':
        data['rad'] = txt[:, 0]
    data[f'Epar_m6_{isotope}'] = np.empty(txt.shape[0], dtype=complex)
    data[f'Epar_m6_{isotope}'].real = txt[:, 1] * clight * 1.0e-4
    data[f'Epar_m6_{isotope}'].imag = txt[:, 2] * clight * 1.0e-4
    data[f'Eperp_m6_{isotope}'] = np.empty(txt.shape[0], dtype=complex)
    data[f'Eperp_m6_{isotope}'].real = txt[:, 4] * clight * 1.0e-4
    data[f'Eperp_m6_{isotope}'].imag = txt[:, 5] * clight * 1.0e-4
    txt = np.loadtxt(f'/temp/kasilov/MEPHIT-KILCA/VARENNA24/STANDALONE_PERTEFIELD/{isotope}/pert_current.dat')
    data[f'Jpar_m6_{isotope}'] = np.empty(txt.shape[0], dtype=complex)
    data[f'Jpar_m6_{isotope}'].real = txt[:, 1] * 1.75e-8  # a.u.
    data[f'Jpar_m6_{isotope}'].imag = txt[:, 2] * 1.75e-8  # a.u.

# %%
fig = plt.figure(figsize=(6, 3.6), dpi=150)
ax = fig.subplots()
ax.plot(data['rad'], np.abs(data[f'Jpar_m6_D']), color='k', ls='-',
        label=r'$\left\lvert J_{\parallel \vec{m}} \right\rvert$, D')
ax.plot(data['rad'], np.abs(data[f'Jpar_m6_H']), color='tab:orange', ls='-',
        label=r'$\left\lvert J_{\parallel \vec{m}} \right\rvert$, H')
ax.plot(data['rad'], np.abs(data[f'Epar_m6_D']), color='k', ls='--',
        label=r'$\bigl\lvert E_{\parallel \vec{m}}^{\text{MA}} \bigr\rvert$, D')
ax.plot(data['rad'], np.abs(data[f'Epar_m6_H']), color='tab:orange', ls='--',
        label=r'$\bigl\lvert E_{\parallel \vec{m}}^{\text{MA}} \bigr\rvert$, H')
ax.set_xlim((48.75, 55.25))
ax.set_xlabel(r'effective radius $r$ / \si{\centi\meter}')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
secax = ax.secondary_xaxis(-0.2, functions=(rsmall_to_psi_norm, psi_norm_to_rsmall))
secax.set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
secax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
ax.set_ylabel(r'$\left\lvert J_{\parallel \vec{m}} \right\rvert$ / a.u., ' +
              r'$\bigl\lvert E_{\parallel \vec{m}}^{\text{MA}} \bigr\rvert$ / \si{\volt\per\meter}')
fig.legend(bbox_to_anchor=(1, 1), bbox_transform=ax.transAxes)
ax.text(0.1, 0.9, fr"$\vec{{m}} = (-6, 2)$", ha='center', va='bottom', transform=ax.transAxes)
fig.savefig(path.join(work_dir, f'Epar_Jpar.pdf'), backend='pgf', dpi=150)
plt.show()
