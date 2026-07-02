# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %%
# %matplotlib inline
from os import environ, getcwd, path, mkdir, remove
from shutil import copy2
from subprocess import Popen, PIPE, STDOUT
from datetime import datetime
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
work_dir = '/temp/zenzphil/MEPHIT/run/33353_2900_EQH_sweep'

# %%
data = {}
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_4_0_D_par_smooth.h5'), 'r') as f:
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
    data['Delta_res'] = f['/mesh/Delta_psi_res_curr'][()]  / (data['psi'][-1] - data['psi'][0])
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_1_H.h5'), 'r') as f:
    data['Imn_res_01_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
"""with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_2.h5'), 'r') as f:
    data['Imn_res_02'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight"""
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_4_D.h5'), 'r') as f:
    data['Imn_res_04_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_4_H.h5'), 'r') as f:
    data['Imn_res_04_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_5_D.h5'), 'r') as f:
    data['Imn_res_05_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_5_H.h5'), 'r') as f:
    data['Imn_res_05_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_6_D.h5'), 'r') as f:
    data['Imn_res_06_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_6_H.h5'), 'r') as f:
    data['Imn_res_06_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
"""with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_0_7.h5'), 'r') as f:
    data['Imn_res_07'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight"""
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_1_0_D.h5'), 'r') as f:
    data['Imn_res_10_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_1_0_H.h5'), 'r') as f:
    data['Imn_res_10_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
"""with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_1_5.h5'), 'r') as f:
    data['Imn_res_15'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight"""
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_2_0_D.h5'), 'r') as f:
    data['Imn_res_20_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_2_0_H.h5'), 'r') as f:
    data['Imn_res_20_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_4_0_D.h5'), 'r') as f:
    data['Imn_res_40_D'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight
with h5py.File(path.join(work_dir, 'mephit_cross_fade_new_int_4_0_H.h5'), 'r') as f:
    data['Imn_res_40_H'] = f['/resonance_sweep/Imn_res'][()] * 10.0 / clight


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
delta_mn = q * data['R0'] * data['v_ExB_orig'] / v_Te * np.maximum(1,
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
for m in range(m_min, m_max):
    k = m - data['m_res_min']
    ax.axvline(normalize_psi(data['psi_res'][k]), linewidth=0.5)
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_0001'][:, k]), label=f'D 0.001, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res'][:, k]), label=f'D 0.1, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_01_H'][:, k]), label=f'H 0.1, $m = {m}$')
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_02'][:, k]), label=f'D 0.2, $m = {m}$')
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_03'][:, k]), label=f'D 0.3, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_04_D'][:, k]), label=f'D 0.4, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_04_H'][:, k]), label=f'H 0.4, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_05_D'][:, k]), label=f'D 0.5, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_05_H'][:, k]), label=f'H 0.5, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_06_D'][:, k]), label=f'D 0.6, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_06_H'][:, k]), label=f'H 0.6, $m = {m}$')
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_07'][:, k]), label=f'D 0.7, $m = {m}$')
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_10'][:, k]), label=f'D 1.0, $m = {m}$')
    #ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_15'][:, k]), label=f'D 1.5, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_20_D'][:, k]), label=f'D 2.0, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_20_H'][:, k]), label=f'H 2.0, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_40_D'][:, k]), label=f'D 4.0, $m = {m}$')
    ax.semilogy(normalize_psi(data['E_r_zero']), np.abs(data['Imn_res_40_H'][:, k]), label=f'H 4.0, $m = {m}$')
ax.set_xlabel(r'electric zero at $\hat{\psi}$')
ax.set_ylabel(r'$I_{\vec{m} \parallel}$ / \si{\ampere}')
# ax.set_ylim((0.0, 0.03e10))
ax.legend()
plt.savefig(path.join(work_dir,'Imn_res_vs_E_r.pdf'))
plt.show()

# %%
fig = plt.figure(figsize=(6, 4), dpi=150)
ax = fig.subplots()
k_min = 0
k_max = 100
for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res'][k_min:k_max, k]), label=f'D 0.1cm, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_01_H'][k_min:k_max, k]), label=f'H 0.1cm, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_04_D'][k_min:k_max, k]), label=f'D 04, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_04_H'][k_min:k_max, k]), label=f'H 04, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_05_D'][k_min:k_max, k]), label=f'D 0.5cm, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_05_H'][k_min:k_max, k]), label=f'H 0.5cm, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_06_D'][k_min:k_max, k]), label=f'D 06, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_06_H'][k_min:k_max, k]), label=f'H 06, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_20_D'][k_min:k_max, k]), label=f'D 20, $m = {m}$')
    #ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_20_H'][k_min:k_max, k]), label=f'H 20, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_40_D'][k_min:k_max, k]), label=f'D 4.0cm, $m = {m}$')
    ax.semilogy(np.arange(k_min, k_max), np.abs(data['Imn_res_40_H'][k_min:k_max, k]), label=f'H 4.0cm, $m = {m}$')
    ax.axvline(35, color='black', linewidth=0.5)
    ax.axvline(43, color='black', linewidth=0.5)
    ax.axvline(60, color='black', linewidth=0.5)
    ax.axvline(80, color='black', linewidth=0.5)
ax.set_xlabel(r'electric field offset $\Delta E_r$')
ax.set_ylabel(r'$I_{\vec{m} \parallel}$ / \si{\ampere}')
ax.legend()
plt.show()


# %%
eresoff = np.array([35, 60])
eres = np.array([43, 80])

for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    print(f"Delta_E_r off el.fl.res. {m}: {data['Delta_E_r'][eresoff[m - m_min]]:.16e}")
    print(f"Delta_E_r  on el.fl.res. {m}: {data['Delta_E_r'][eres[m - m_min]]:.16e}")
    print(f"v_ExB     off el.fl.res. {m}: {data['v_ExB'][eresoff[m - m_min], k]:.16e}")
    print(f"v_ExB      on el.fl.res. {m}: {data['v_ExB'][eres[m - m_min], k]:.16e}")


# %%
fig = plt.figure(figsize=(6, 2.25 * (m_max - m_min + 1)), dpi=150)
axs = fig.subplots(m_max - m_min + 1, 1)
for m in range(m_min, m_max + 1):
    k = m - data['m_res_min']
    # axs[m - m_min].axhline(0.0, color='black', linewidth=0.5)
    axs[m - m_min].axvline(0.0, linewidth=0.5, linestyle=':')
    axs[m - m_min].axvline(data['v_ExB'][eresoff[m - m_min], k]/1e5, linewidth=0.5, linestyle='-')
    axs[m - m_min].axvline(data['v_ExB'][eres[m - m_min], k]/1e5, linewidth=0.5, linestyle='--')
    axs[m - m_min].axvline(data['v_ExB_orig'][m - data['m_res_min']]/1e5, linewidth=0.5, linestyle='-.')
    # axs[m - m_min].axvline(0.5e-6 * data['V_e_diamag'][k], linewidth=0.5, linestyle='--')
    # axs[m - m_min].axvline(0.5e-6 * data['V_i_diamag'][k], linewidth=0.5, linestyle=':')
    mask = data['v_ExB'][:, k] > -data['v_ExB'][0, k]
    """axs[m - m_min].semilogy(data['v_ExB'][mask, k],
                            np.abs(data['Imn_res'][mask, k]),
                            color='k', ls='-', label='D 0.1cm')"""
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_05_D'][mask, k]),
                            color='red', ls='-', label='D 0.5cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_10_D'][mask, k]),
                            color='blue', ls='-', label='D 1.0cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_20_D'][mask, k]),
                            color='green', ls='-', label='D 2.5cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_40_D'][mask, k]),
                            color='k', ls='-', label='D 4.0cm')
    """axs[m - m_min].semilogy(data['v_ExB'][mask, k],
                            np.abs(data['Imn_res_01_H'][mask, k]),
                            color='tab:orange', ls='-', label='H 0.1cm')"""
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_05_H'][mask, k]),
                            color='red', ls='--', label='H 0.5cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_10_H'][mask, k]),
                            color='blue', ls='--', label='H 1.0cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_20_H'][mask, k]),
                            color='green', ls='--', label='H 2.5cm')
    axs[m - m_min].semilogy(data['v_ExB'][mask, k]/1e5,
                            np.abs(data['Imn_res_40_H'][mask, k]),
                            color='k', ls='--', label='H 4.0cm')

    xlim = np.array(axs[m - m_min].get_xlim())
    print(xlim)
    axs[m - m_min].set_xlim(np.array([-1.0, 1.0]) * np.amax(np.abs(xlim)))
    axs[m - m_min].set_xlabel(r'$V_{E \times B}$ / \si{\km\per\second}')
    axs[m - m_min].set_ylabel(r'$I_{\parallel \vec{m}}$ / a.u.')  # \si{\ampere}
    axs[m - m_min].legend(loc='lower right', fontsize='small', ncol=2)
    axs[m - m_min].text(0.1, 0.1, fr"$\vec{{m}} = ({m * sgn_m_res}, 2)$",
                        ha='center', va='bottom', transform=axs[m - m_min].transAxes)
    axs[m - m_min].text(data['v_ExB'][eres[m - m_min], k]/1e5, 0.95,
                        'on el.fl.res.', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
    axs[m - m_min].text(data['v_ExB'][eresoff[m - m_min], k]/1e5, 0.95,
                        'off el.fl.res.', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
    axs[m - m_min].text(data['v_ExB_orig'][m - data['m_res_min']]/1e5, 0.45 if m == 7 else 0.95,
                        r'orig.\ value', rotation=90, ha='right', va='top',
                        transform=axs[m - m_min].get_xaxis_transform())
fig.savefig(path.join(work_dir, 'sweep_resonance.pdf'), backend='pgf')
fig.savefig(path.join(work_dir, 'sweep_resonance.png'), dpi=300)
plt.show()

# %%
# post-process in Fortran before loading HDF5 files below
# see last cell for data returned by Sergei
command = path.join('/afs/itp.tugraz.at/user/zenzphil/code/MEPHIT/build', 'bin', 'mephit_post.x')
out_dir = path.join(work_dir, 'KiLCA-MEPHIT_' + datetime.today().strftime('%Y-%m-%d'))
if not path.exists(out_dir):
    mkdir(out_dir)
for isotope in ['D', 'H']:
    for shift in ['eresoff', 'eres']:
        for m in ['6', '7']:
            sim_in = f'{isotope}_{shift}_{m}'
            with Popen([command, path.join(work_dir, f'mephit_{sim_in}.h5'), m],
                       cwd=work_dir, stdout=PIPE, stderr=STDOUT, text=True) as proc:
                for line in proc.stdout:
                    print(line, end='')
            sim_out = isotope + ('_on' if shift == 'eres' else '_off') + '_el_fl_res_m' + m
            if not path.exists(path.join(out_dir, sim_out)):
                mkdir(path.join(out_dir, sim_out))
            for fname in ['response_current.dat', 'bpsi_over_bphi.dat',
                          'Phi_m.dat', 'Phi_m_aligned.dat']:
                copy2(path.join(work_dir, fname), path.join(out_dir, sim_out, fname))
                remove(path.join(work_dir, fname))

# %%
mephit = {}
for width in ['0_5','2_5','4_0']:
    for isotope in ['D', 'H']:
        for shift in ['eresoff', 'eres']:
            for m in [6, 7]:
                sim = f'{width}_{isotope}_{shift}_{m}'
                mephit[sim] = Mephit(work_dir,f'mephit_cross_fade_new_int_{sim}.h5')
                mephit[sim].open_datafile()
                mephit[sim].postprocess()
mephit['iMHD'] = Mephit(work_dir, 'mephit_imhd.h5')
mephit['05'] = Mephit(work_dir, 'mephit_cross_fade_new_int_0_5_D.h5')
mephit['20'] = Mephit(work_dir, 'mephit_cross_fade_new_int_2_0_D.h5')
mephit['40_D'] = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D.h5')
mephit['40_H'] = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_H.h5')
mephit['vac_resp'] = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_no_response.h5')
mephit['smooth'] = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_par_smooth.h5')
mephit['damp'] = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_par_damp.h5')

#mephit['05_new'] = Mephit('/temp/lainer_p/MEPHIT-run/33353_2900_EQH', 'mephit_cross_fade_new_int_0_5_D.h5')
mephit['iMHD'].open_datafile()
mephit['iMHD'].postprocess()
mephit['05'].open_datafile()
mephit['05'].postprocess()
mephit['20'].open_datafile()
mephit['20'].postprocess()
mephit['40_D'].open_datafile()
mephit['40_D'].postprocess()
mephit['40_H'].open_datafile()
mephit['40_H'].postprocess()
mephit['vac_resp'].open_datafile()
mephit['vac_resp'].postprocess()
mephit['smooth'].open_datafile()
mephit['smooth'].postprocess()
mephit['damp'].open_datafile()
mephit['damp'].postprocess()

#gpec = Gpec(work_dir, 2)
#gpec.open_datafiles()

# %%
res = mephit['iMHD'].normalize_psi(mephit['iMHD'].data['/mesh/psi_res'][()])
conversion = 1.0e-04 / mephit['iMHD'].data['/mesh/gpec_jacfac'][()]
print(len(conversion))
sgn_dpsi = np.sign(mephit['iMHD'].data['/cache/fs/psi'][-1] - mephit['iMHD'].data['/cache/fs/psi'][0])
Bmn = {'vac': mephit['iMHD'].get_polmodes('vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
       'iMHD': mephit['iMHD'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       '20' : mephit['20'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       '05' : mephit['05'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       '40_D': mephit['40_D'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       '40_H': mephit['40_H'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       'vac_resp': mephit['05'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       'smooth': mephit['smooth'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion),
       'damp': mephit['damp'].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion)
       }#'GPEC': gpec.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi')}
Ires = {}
for width in ['0_5','2_5','4_0']:
    for isotope in ['D', 'H']:
        for shift in ['eresoff', 'eres']:
            for m in [6, 7]:
                sim = f'{width}_{isotope}_{shift}_{m}'
                print(f"Loading {sim} ...")
                Bmn[sim] = mephit[sim].get_polmodes('full perturbation', '/iter/Bmn/coeff_rad', conversion)
                Ires[sim] = mephit[sim].get_Ires()

# %%
for m in [6, 7]:
    for width in ['0_5', '2_5', '4_0']:
        for shift in ['eresoff', 'eres']:
            for isotope in ['H', 'D']:
                sim = f'{width}_{isotope}_{shift}_{m}'
                print(f"{sim}: {Ires[sim][m * sgn_m_res]}")

# %%
fig = plt.figure(figsize=(7.2, 7.2), dpi=150)
axs = fig.subplots(2, 1, sharey='all')
for k, m in enumerate([6, 7]):
    m_res = sgn_m_res * m
    axs[k].axhline(0.0, color='k', lw=0.5)
    axs[k].axvline(res[m - data['m_res_min']], color='k', lw=0.5)
    axs[k].plot(Bmn['vac']['rho'][m_res][1:], np.abs(Bmn['vac']['var'][m_res][1:]),
                ls=':', color='k', label='vacuum perturbation')  # purple
    axs[k].plot(Bmn['iMHD']['rho'][m_res][1:], np.abs(Bmn['iMHD']['var'][m_res][1:]),
                ls='-.', color='tab:green', label='MEPHIT, iMHD')
    """axs[k].plot(Bmn['GPEC']['rho'][m_res], np.abs(Bmn['GPEC']['var'][m_res]),
                ls='-.', color='tab:purple', label='GPEC')"""
    sim = f'4_0_D_eres_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='tab:orange', label='MEPHIT, D on el.fl.res. 0.5cm')
    sim = f'4_0_H_eres_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='--', color='tab:orange', label='MEPHIT, H on el.fl.res. 0.5cm')
    sim = f'4_0_D_eresoff_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='k', label='MEPHIT, D off el.fl.res. 0.5cm')
    sim = f'4_0_H_eresoff_{m}'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='--', color='k', label='MEPHIT, H off el.fl.res. 0.5cm')
    sim = 'smooth'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='red', label='MEPHIT, D 4.0cm smooth')
    sim = 'damp'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='green', label='MEPHIT, D 4.0cm damped')
    sim = '40_D'
    axs[k].plot(Bmn[sim]['rho'][m_res][1:], np.abs(Bmn[sim]['var'][m_res][1:]),
                ls='-', color='tab:blue', label='MEPHIT, D 4.0cm')
    #axs[k].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    #axs[k].yaxis.set_major_locator(ticker.MultipleLocator(0.5e-4))
    axs[k].set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
    axs[k].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
    axs[k].text(0.1, 0.1, fr"$\vec{{m}} = ({m_res}, 2)$",
        ha='center', va='bottom', transform=axs[k].transAxes)
    axs[k].set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ / \si{\tesla}')
    #axs[0].set_ylabel(r'normal mag. field perturbation / \si{\tesla}')
    axs[k].legend(loc='upper left')
fig.savefig(path.join(work_dir, f'Bmn_67_pert.pdf'), backend='pgf', dpi=150)
fig.savefig(path.join(work_dir, f'Bmn_67_pert.png'), dpi=300)
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
    colors = ['red', 'blue', 'green', 'k']
    for width in ['4_0']:
        for j, shift in enumerate(['eresoff', 'eres']):
            for i, isotope in enumerate(['D', 'H']):
                sim = f'{width}_{isotope}_{shift}_{m}'
                mask = (Bmn[sim]['rho'][m * sgn_m_res] >= res[m - data['m_res_min']] +
                         data['Delta_res'][m - data['m_res_min']]) & \
                    (Bmn[sim]['rho'][m * sgn_m_res] <= res[m - data['m_res_min']] -
                         data['Delta_res'][m - data['m_res_min']])
                axs[k].plot(Bmn[sim]['rho'][m * sgn_m_res][mask],
                            np.abs(Bmn[sim]['var'][m * sgn_m_res][mask]),
                            color=colors[['0_5', '2_5', '4_0'].index(width)],
                            ls=('--' if i == 1 else '-'),
                            marker='x' if j == 0 else 'none',
                            markevery=0.1,
                            label=isotope + (' on ' if j == 1 else ' off ') + width.replace('_', '.') + 'cm')
    
    for sim in ['smooth', 'damp', '40_D']:
        axs[k].plot(Bmn[sim]['rho'][m * sgn_m_res][mask],
                                np.abs(Bmn[sim]['var'][m * sgn_m_res][mask]),
                                color='k' if sim == 'damp' else ('red' if sim == 'smooth' else 'blue'),
                                ls='-',
                                label= sim)
    axs[k].set_xlabel(r'normalized ploidal flux $\hat{\psi}$')
    axs[k].set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\vec{m}} \rvert A^{-1}$ / \si{\tesla}')
    axs[k].xaxis.set_major_locator(ticker.MultipleLocator(0.025))
    axs[k].yaxis.set_tick_params(labelleft=True)
    axs[k].yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
    axs[k].text(0.2, 0.05, fr"$\vec{{m}} = ({m * sgn_m_res}, 2)$",
                ha='center', va='bottom', transform=axs[k].transAxes)
    axs[k].legend(loc='upper left', fontsize='small', ncol=2)
ylims = axs[0].get_ylim()
axs[0].set_ylim((ylims[0], 1.1 * ylims[1]))
fig.savefig(path.join(work_dir, f'zoom_resonance.pdf'), backend='pgf', dpi=150)
fig.savefig(path.join(work_dir, f'zoom_resonance.png'), dpi=300)
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
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
#secax = ax.secondary_xaxis(-0.2, functions=(rsmall_to_psi_norm, psi_norm_to_rsmall))
#secax.set_xlabel(r'normalized poloidal flux $\hat{\psi}$')
psi_min = normalize_psi(data['psi']).min()
psi_max = normalize_psi(data['psi']).max()
#secax.set_xlim(psi_min, psi_max)
#secax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
ax.set_ylabel(r'$\left\lvert J_{\parallel \vec{m}} \right\rvert$ / a.u., ' +
              r'$\bigl\lvert E_{\parallel \vec{m}}^{\text{MA}} \bigr\rvert$ / \si{\volt\per\meter}')
fig.legend(bbox_to_anchor=(1, 1), bbox_transform=ax.transAxes)
ax.text(0.1, 0.9, fr"$\vec{{m}} = (-6, 2)$", ha='center', va='bottom', transform=ax.transAxes)
fig.savefig(path.join(work_dir, f'Epar_Jpar.pdf'), backend='pgf', dpi=150)
plt.show()
