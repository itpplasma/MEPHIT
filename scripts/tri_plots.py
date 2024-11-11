# %%
# %matplotlib inline
from os import environ, getcwd, path
from copy import copy
from matplotlib import rcParams
import matplotlib.pyplot as plt
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
    data['Bn_psi'] = f['/iter/Bn/comp_psi_contravar_dens'][()] * 1.0e-5  # Mx to mWb
    data['R_O'] = f['/mesh/R_O'][()] * 1.0e-2
    data['Z_O'] = f['/mesh/Z_O'][()] * 1.0e-2


# %%
# notes on old configuration: max_Delta_rad = 1.25 cm, refinement = sqrt(2) & deletions = 1 for m = 3, 4, 7
# diverging_q needs to be set to false in the call to mephit_mesh::refine_eqd_partition
# for post-processing, run e.g.
# pdfcrop --luatex mesh.pdf mesh_cropped.pdf
fig = plt.figure(figsize=(3.3, 4.4), dpi=200)
ax = fig.subplots()
ax.triplot(data['R'], data['Z'], data['tri'], lw=0.25)
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / m')
ax.set_ylabel(r'$Z$ / m')
fig.savefig(path.join(work_dir, 'mesh.pdf'), dpi=600)
plt.show()

# %%
clim = np.amax(data['Bn_psi'].real) * 0.1
fig = plt.figure(figsize=(4.0, 4.4), dpi=200)
ax = fig.subplots()
im = ax.tripcolor(data['R'], data['Z'], data['tri'], data['Bn_psi'].real, cmap=cmap_divg_clip)
# ax.plot(data['R_O'], data['Z_O'], 'k.', markersize=1)
cbar = fig.colorbar(im)
# cbar.ax.yaxis.set_offset_position('right')
# cbar.ax.yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
cbar.set_label(r'$\Real \sqrt{g} B_{n}^{\psi}$ / \si{\milli\weber}', rotation=90)
im.set_clim([-clim, clim])
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / m')
ax.set_ylabel(r'$Z$ / m')
fig.savefig(path.join(work_dir, 'Bn_psi.png'), dpi=600)
plt.show()
# for postprocessing, run e.g.
# optipng -preserve -clobber -out Bn_psi_optim.png Bn_psi.png
# convert Bn_psi_optim.png -trim +repage Bn_psi.png
