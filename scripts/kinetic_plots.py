# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
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
from scipy.interpolate import interp1d
import h5py
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')

# %%
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')
data = {}
with h5py.File(path.join(work_dir, 'mephit.h5'), 'r') as f:
    data['E_r_zero'] = f['/resonance_sweep/E_r_zero'][()]
    data['Imn_res'] = f['/resonance_sweep/Imn_res'][()]
    data['rsmall'] = f['/cache/fs/rsmall'][()]
    data['psi'] = f['/cache/fs/psi'][()]
    data['psi_res'] = f['/mesh/psi_res'][()]
    data['m_res_min'] = f['/mesh/m_res_min'][()]
psi_to_rsmall = interp1d(data['psi'], data['rsmall'], 'cubic')
data['rsmall_res'] = psi_to_rsmall(data['psi_res'])

# %%
fig = plt.figure(figsize=(6, 4))
ax = fig.subplots()
# ax.axhline(0, color='black', linewidth=0.5)
for k in range(5):
    m = data['m_res_min'] + k
    ax.axvline(data['rsmall_res'][k], color='black', linewidth=0.5)
    ax.semilogy(data['E_r_zero'], np.abs(data['Imn_res'][:, k]), label=f'$m = {m}$')
ax.set_xlabel(r'$E_{r}$ zero crossing at $r$ / cm')
ax.set_ylabel(r'$I_{mn}^{\parallel}$ / statA')
ax.legend()
plt.show()
