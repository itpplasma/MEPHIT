# %%
# %matplotlib inline
from os import environ, getcwd, path
from copy import copy
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import colorcet as cc
import numpy as np
import h5py
from mephit_plot import Mephit
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')

# %%
cmap_divg_clip = copy(cc.m_CET_D9)
cmap_divg_clip.set_under(color='w', alpha=0.0)
cmap_divg_clip.set_over(color='w', alpha=0.0)

# %%
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')
data = {}
with h5py.File(path.join(work_dir, 'mephit.h5'), 'r') as f:
    data['R'] = f['/mesh/node_R'][()] * 1.0e-2
    data['Z'] = f['/mesh/node_Z'][()] * 1.0e-2
    data['tri'] = f['/mesh/tri_node'][()] - 1
    data['Bn_psi'] = f['/iter/Bn/comp_psi_contravar_dens'][()]  # * 1.0e-5  # Mx to mWb
    data['Bn_vac_psi'] = f['/vac/Bn/comp_psi_contravar_dens'][()]  # * 1.0e-5  # Mx to mWb
    data['R_O'] = f['/mesh/R_O'][()] * 1.0e-2
    data['Z_O'] = f['/mesh/Z_O'][()] * 1.0e-2


# %%
# for configuration see data/mephit_illustration.in
# diverging_q needs to be set to false in the call to mephit_mesh::refine_eqd_partition
# for post-processing, run e.g.
# pdfcrop --luatex mesh.pdf mesh_cropped.pdf
fig = plt.figure(figsize=(3.3, 4.4), dpi=200)  # (5.7, 7.6) for full-page plot
ax = fig.subplots()
ax.triplot(data['R'], data['Z'], data['tri'], lw=0.25)
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / m')
ax.set_ylabel(r'$Z$ / m')
ax.yaxis.set_major_locator(MultipleLocator(0.25))
fig.savefig(path.join(work_dir, 'mesh.pdf'), dpi=600)
plt.show()

# %%
clim = np.amax(data['Bn_psi'].real) * 0.1
for dataset in ['Bn_vac_psi', 'Bn_psi']:
    fig = plt.figure(figsize=(4.0, 4.4), dpi=300)
    ax = fig.subplots()
    im = ax.tripcolor(data['R'], data['Z'], data['tri'], data[dataset].real, cmap=cc.m_CET_D9, rasterized=True)
    # ax.plot(data['R_O'], data['Z_O'], 'k.', markersize=1)
    cbar = fig.colorbar(im)
    # cbar.ax.yaxis.set_offset_position('right')
    # cbar.ax.yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
    cbar.set_label(r'normal magnetic field perturbation $\Real \sqrt{g} B_{n}^{\psi}$ / \si{\gauss\centi\meter\squared}', rotation=90)
    im.set_clim([-clim, clim])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$R$ / m')
    ax.set_ylabel(r'$Z$ / m')
    fig.savefig(path.join(work_dir, f'{dataset}.pdf'), backend='pgf', dpi=300)
    plt.show()
# for postprocessing, run e.g.
# optipng -preserve -clobber -out Bn_psi_optim.png Bn_psi.png
# convert Bn_psi_optim.png -trim +repage Bn_psi.png
