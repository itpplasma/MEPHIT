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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# %matplotlib inline
from os import environ, getcwd, path
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import LogNorm
import numpy as np
import h5py
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')

# %%
data = {}
with h5py.File(path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH/mephit.h5'), 'r') as f:
    data['eigvals'] = f['/iter/eigvals'][()]
    data['rel_err'] = f['/iter/L2int_Bn_diff'][()] / f['/iter/L2int_Bnvac'][()]
    data['triangulation'] = Triangulation(f['/mesh/node_R'][()], f['/mesh/node_Z'][()], f['/mesh/tri_node'][()] - 1)
    if '/iter/debug_kernel' in f:
        data['kernel'] = np.sqrt(
            np.abs(f['/iter/debug_kernel/comp_R'][()]) ** 2 +
            np.abs(f['/iter/debug_kernel/comp_Z'][()]) ** 2 +
            np.abs(f['/iter/debug_kernel/RT0_comp_phi'][()]) ** 2)
data['supremum'] = np.sort(np.abs(data['eigvals'][::-1]))[0]
data['nritz'] = data['eigvals'].size
data['nbad'] = np.sum(np.abs(data['eigvals']) >= 1.0)
data['niter'] = np.sum(np.logical_not(np.isnan(data['rel_err'])))
print(f"{data['nritz']} total eigenvalues and {data['nbad']} bad eigenvalues")
print(f"relative error {data['rel_err'][data['niter'] - 1]} after {data['niter']} iterations")

# %%
fig = plt.figure()
ax = fig.subplots()
iter_range = np.arange(data['niter'])
ax.semilogy(iter_range, data['rel_err'][:data['niter']], '-o')
ax.semilogy(iter_range, data['supremum'] ** (iter_range + 1))
ax.set_xlabel('iteration count')
ax.set_ylabel('relative error after iteration step')
plt.show()

# %%
fig = plt.figure()
ax = fig.subplots(subplot_kw={'projection': 'polar'})
ax.scatter(np.angle(data['eigvals']), np.abs(data['eigvals']))
ax.set_yscale('symlog', linthresh=1.0)
plt.show()

# %%
fig = plt.figure()
ax = fig.subplots()
ax.plot(np.arange(data['nritz']), np.abs(data['eigvals']), '-o')
ax.set_yscale('log')
ax.set_xlabel('eigenvalue index')
ax.set_ylabel('absolute value of eigenvalue')
ax.legend()
plt.show()

# %%
if 'kernel' in data:
    fig = plt.figure()
    ax = fig.subplots()
    im = ax.tripcolor(data['triangulation'], data['kernel'], cmap='magma',
                      norm=LogNorm(vmin=data['kernel'].min(), vmax=data['kernel'].max()))
    ax.set_aspect('equal')
    cbar = fig.colorbar(im)
    plt.show()
