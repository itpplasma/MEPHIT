# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
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
from matplotlib.tri import Triangulation
from matplotlib.colors import LogNorm
import colorcet as cc
import numpy as np
import h5py
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')

# %%
work_dir = path.join(environ['MEPHIT_RUN_DIR'], '33353_2900_EQH')  # 47046
data = {}
with h5py.File(path.join(work_dir, 'mephit_imhd.h5'), 'r') as f:
    data['eigvals'] = f['/iter/eigvals'][()]
    data['supremum'] = f['/config/ritz_threshold'][()]
    data['rel_err'] = f['/iter/L2int_Bn_diff'][()] / f['/iter/L2int_Bnvac'][()]
    data['rel_err_requested'] = f['/config/iter_rel_err'][()]
    data['triangulation'] = Triangulation(f['/mesh/node_R'][()] * 1.0e-2,
                                          f['/mesh/node_Z'][()] * 1.0e-2,
                                          f['/mesh/tri_node'][()] - 1)
    if '/iter/eigvec_001' in f:
        data['eigvec'] = np.sqrt(
            np.abs(f['/iter/eigvec_001/comp_R'][()]) ** 2 +
            np.abs(f['/iter/eigvec_001/comp_Z'][()]) ** 2 +
            np.abs(f['/iter/eigvec_001/RT0_comp_phi'][()]) ** 2)
    if '/iter/debug_kernel' in f:
        data['kernel'] = np.sqrt(
            np.abs(f['/iter/debug_kernel/comp_R'][()]) ** 2 +
            np.abs(f['/iter/debug_kernel/comp_Z'][()]) ** 2 +
            np.abs(f['/iter/debug_kernel/RT0_comp_phi'][()]) ** 2)
data['nritz'] = data['eigvals'].size
data['nbad'] = np.sum(np.abs(data['eigvals']) >= 1.0)
data['niter'] = np.sum(np.logical_not(np.isnan(data['rel_err'])))
print(f"{data['nritz']} total eigenvalues and {data['nbad']} bad eigenvalues")
print(f"relative error {data['rel_err'][data['niter'] - 1]} after {data['niter']} iterations")

# %%
fig = plt.figure(figsize=(6.6, 3.6), dpi=150)
ax = fig.subplots()
iter_range = np.arange(data['niter']) + 1
ax.semilogy(iter_range, data['rel_err'][:data['niter']], '-o', markersize=3)
ax.semilogy(iter_range, data['supremum'] ** iter_range)
ax.axhline(data['rel_err_requested'], linestyle=':', linewidth=0.5, color='k')
ax.set_xlabel('iteration count')
ax.set_ylabel('relative error after iteration step')
fig.savefig(path.join(work_dir, 'convergence.pdf'), backend='pgf')
plt.show()

# %%
fig = plt.figure(figsize=(3.3, 3.3), dpi=150)
ax = fig.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.linspace(0, 2 * np.pi, 100), np.ones(100) * data['supremum'], '--', color='tab:gray')
ax.plot(np.angle(data['eigvals']), np.abs(data['eigvals']), 'o', markersize=3)
ax.set_yscale('symlog', linthresh=1.0)
ax.set_rlabel_position(90)
fig.savefig(path.join(work_dir, 'eigenvalues_polar.pdf'), backend='pgf')
plt.show()

# %%
fig = plt.figure(figsize=(6.6, 3.3), dpi=150)
ax = fig.subplots()
ax.plot(np.arange(data['nritz']), np.abs(data['eigvals']), '-o')
ax.set_yscale('log')
ax.set_xlabel('eigenvalue index')
ax.set_ylabel('absolute value of eigenvalue')
fig.savefig(path.join(work_dir, 'eigenvalues.pdf'), backend='pgf')
plt.show()

# %%
if 'kernel' in data:
    fig = plt.figure(figsize=(3.3, 4.0), dpi=300)
    ax = fig.subplots()
    im = ax.tripcolor(data['triangulation'], data['kernel'], cmap=cc.m_fire, rasterized=True,
                      norm=LogNorm(vmin=data['kernel'].min(), vmax=data['kernel'].max()))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.set_aspect('equal')
    ax.set_xlabel(r'$R$ / \si{\meter}')
    ax.set_ylabel(r'$Z$ / \si{\meter}')
    cbar = fig.colorbar(im)
    fig.savefig(path.join(work_dir, 'preconditioner_kernel.pdf'), backend='pgf')
    plt.show()

# %%
if 'eigvec' in data:
    fig = plt.figure(figsize=(3.3, 4.0), dpi=300)
    ax = fig.subplots()
    im = ax.tripcolor(data['triangulation'], data['eigvec'], cmap=cc.m_fire, rasterized=True)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.set_aspect('equal')
    ax.set_xlabel(r'$R$ / \si{\meter}')
    ax.set_ylabel(r'$Z$ / \si{\meter}')
    cbar = fig.colorbar(im)
    fig.savefig(path.join(work_dir, 'dominant_eigenvector.pdf'), backend='pgf')
    plt.show()
