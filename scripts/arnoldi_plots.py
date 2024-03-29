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
from os import path, getenv
from h5py import File, get_config
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mephit_plot import Mephit, run_dir, set_matplotlib_defaults

# %%
get_config().complex_names = ('real', 'imag')
set_matplotlib_defaults()
run_dir = path.join(getenv('MEPHIT_DIR'), 'run/33353_2900_EQH')

# %%
cases = [
    {'fname': 'mephit.h5', 'title': 'KiLCA, fine'},
    {'fname': 'mephit_coarse.h5', 'title': 'KiLCA, coarse'},
    {'fname': 'mephit_iMHD.h5', 'title': 'iMHD, fine'},
    {'fname': 'mephit_iMHD_coarse.h5', 'title': 'iMHD, coarse'},
    {'fname': 'mephit_iMHD_damped.h5', 'title': 'iMHD, fine, damped'},
    {'fname': 'mephit_iMHD_coarse_damped.h5', 'title': 'iMHD, coarse, damped'},
]

# %%
for c in cases:
    with File(path.join(run_dir, c['fname']), 'r') as f:
        c['eigvals'] = f['/iter/eigvals'][()]
        c['rel_err'] = f['/iter/L2int_Bn_diff'][()] / f['/iter/L2int_Bnvac'][()]
        c['triangulation'] = Triangulation(f['/mesh/node_R'][()],
                                           f['/mesh/node_Z'][()],
                                           f['/mesh/tri_node'][()] - 1)
        nritz = c['eigvals'].size
        nbad = np.sum(np.abs(c['eigvals']) >= 1.0)
        print(f"{c['title']}: {nritz} total eigenvalues, {nbad} bad eigenvalues")
        if '/iter/debug_kernel' in f:
            c['kernel'] = np.sqrt(
                np.abs(f['/iter/debug_kernel/comp_R'][()]) ** 2 +
                np.abs(f['/iter/debug_kernel/comp_Z'][()]) ** 2 +
                np.abs(f['/iter/debug_kernel/RT0_comp_phi'][()]) ** 2)

# %%
fig = plt.figure(layout='constrained')
ax = fig.subplots()
for k, c in enumerate(cases):
    niter = c['rel_err'].size
    ax.semilogy(np.arange(niter), c['rel_err'][:niter], '-', label=c['title'])
ax.legend()
ax.set_xlabel('iteration count')
ax.set_ylabel('relative error after iteration')
fig.show()

# %%
for k, c in enumerate(cases):
    fig = plt.figure(layout='constrained')
    ax = fig.subplots(subplot_kw={'projection': 'polar'})
    ax.scatter(np.angle(c['eigvals']), np.abs(c['eigvals']))
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_title('Eigenvalues: ' + c['title'])
    fig.show()

# %%
fig = plt.figure(layout='constrained')
ax = fig.subplots()
for c in cases:
    ax.plot(np.arange(c['eigvals'].size), np.abs(c['eigvals']),
            '-x', label=c['title'])
ax.set_yscale('log')
ax.set_xlabel('eigenvalue index')
ax.set_ylabel('absolute value of eigenvalue')
ax.legend()
fig.show()

# %%
for c in cases:
    if 'kernel' in c:
        fig = plt.figure(layout='constrained')
        ax = fig.subplots()
        im = ax.tripcolor(c['triangulation'], c['kernel'], cmap='magma')
        ax.set_aspect('equal')
        ax.set_title('Preconditioner kernel: ' + c['title'])
        cbar = fig.colorbar(im)
        # im.set_clim([0.0, max(abs(c['kernel']))])
        fig.show()
