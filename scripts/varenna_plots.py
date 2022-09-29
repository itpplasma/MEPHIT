from os import path, remove, replace
from copy import copy
import subprocess
from numpy import abs, amax, sign
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mephit_plot import Gpec, Kilca, Mephit, run_dir, set_matplotlib_defaults
import colorcet

cmap_divg_clip = copy(colorcet.cm.CET_D9)
cmap_divg_clip.set_under(color='w', alpha=0.0)
cmap_divg_clip.set_over(color='w', alpha=0.0)
set_matplotlib_defaults()

work_dir = run_dir + '/illustration'
# notes on configuration: max_Delta_rad = 1.25 cm, refinement = sqrt(2) & deletions = 1 for m = 3, 4, 7
# diverging_q needs to be set to false in the call to mephit_mesh::refine_eqd_partition
testcase = Mephit(work_dir)
testcase.open_datafile()
R = testcase.data['/mesh/node_R'][()] * 1.0e-2
Z = testcase.data['/mesh/node_Z'][()] * 1.0e-2
tri = testcase.data['/mesh/tri_node'][()] - 1
fig = Figure(figsize=(3.6, 4.4))
ax = fig.subplots()
ax.triplot(R, Z, tri, lw=0.25)
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / m')
ax.set_ylabel(r'$Z$ / m')
canvas = FigureCanvas(fig)
filename = path.join(work_dir, 'mesh.pdf')
fig.savefig(filename, dpi=600)
try:
    cropped = path.join(work_dir, 'mesh_cropped.pdf')
    subprocess.run(['pdfcrop', '--luatex', filename, cropped], check=True)
    replace(cropped, filename)
except subprocess.CalledProcessError as err:
    print(f"Error calling pdfcrop/mv: {err=}")
except FileNotFoundError as err:
    print(f"Error calling pdfcrop/mv: {err=}")
finally:
    testcase.close_datafile()

testcase = Mephit(run_dir + '/33353_2325')
testcase.open_datafile()
R = testcase.data['/mesh/node_R'][()] * 1.0e-2
Z = testcase.data['/mesh/node_Z'][()] * 1.0e-2
tri = testcase.data['/mesh/tri_node'][()] - 1
data = testcase.data['/iter/Bn/comp_psi_contravar_dens'][()].real * 1.0e-5  # Mx to mWb
clim = amax(data) * 1.75e-3
fig = Figure(figsize=(4.0, 4.4))
ax = fig.subplots()
im = ax.tripcolor(R, Z, tri, data, cmap=cmap_divg_clip)
cbar = fig.colorbar(im)
# cbar.ax.yaxis.set_offset_position('right')
# cbar.ax.yaxis.offsetText.set(x=-0.02, verticalalignment='top', horizontalalignment='right')
cbar.set_label(r'$\Real \sqrt{g} B_{n}^{\psi}$ / \si{\milli\weber}', rotation=90)
im.set_clim([-clim, clim])
ax.set_aspect('equal')
ax.set_xlabel(r'$R$ / m')
ax.set_ylabel(r'$Z$ / m')
canvas = FigureCanvas(fig)
filename = path.join(work_dir, 'Bn_psi.png')
fig.savefig(filename, dpi=600)
try:
    optim = path.join(work_dir, 'Bn_psi_optim.png')
    subprocess.run(['optipng', '-preserve', '-clobber', '-out', optim, filename], check=True)
    subprocess.run(['convert', optim, '-trim', '+repage', filename], check=True)
    remove(optim)
except subprocess.CalledProcessError as err:
    print(f"Error calling opting/ImageMagick/rm: {err=}")
except FileNotFoundError as err:
    print(f"Error calling opting/ImageMagick/rm: {err=}")

testcase.postprocess()
sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])
conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][()]
mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/postprocess/Bmn/coeff_rad', conversion)
mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/postprocess/Bmn_vac/coeff_rad', conversion)
reference = Gpec(testcase.work_dir, 2)
reference.open_datafiles()
gpec_Bmn = reference.get_polmodes('full perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi')
gpec_Bmn_vac = reference.get_polmodes('vacuum perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi_x')
config = {
    'xlabel': r'$\hat{\psi}$', 'rho': testcase.post['psi_norm'],
    'ylabel': r'$\abs\, [\sqrt{g} B_{n}^{\psi}]_{m} A^{-1}$ / \si{\tesla}',
    'poldata': [mephit_Bmn_vac, gpec_Bmn_vac, mephit_Bmn, gpec_Bmn], 'comp': abs,
    'resonances': testcase.post['psi_norm_res'], 'sgn_m_res': testcase.post['sgn_m_res']
}
m_min = 3
horz_plot = 2
vert_plot = 3
squares = (3.3 * horz_plot, 3.3 * vert_plot)
fig = Figure(figsize=squares)
axs = fig.subplots(vert_plot, horz_plot, sharex='all', sharey='all')
for m_abs in range(m_min, m_min + vert_plot):
    for k in range(horz_plot):
        m = (2 * k - 1) * m_abs * config['sgn_m_res']
        axs[m_abs - m_min, k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
        if m in config['resonances']:
            axs[m_abs - m_min, k].axvline(config['resonances'][m],
                                          color='b', alpha=0.5, lw=0.5, label='resonance position')
        for data in config['poldata']:
            axs[m_abs - m_min, k].plot(data['rho'][m], config['comp'](data['var'][m]), label=data['label'])
        axs[m_abs - m_min, k].set_title(('resonant ' if m in config['resonances'] else 'non-resonant ') + fr"$m = {m}$")
        axs[m_abs - m_min, k].set_xlabel(config['xlabel'])
        axs[m_abs - m_min, k].set_ylabel(config['ylabel'])
        if m_abs > m_min:
            axs[m_abs - m_min, k].xaxis.set_tick_params(labelbottom=True)
    axs[m_abs - m_min, 1].yaxis.set_tick_params(labelleft=True)
handles, labels = axs[0, 1].get_legend_handles_labels()
axs[0, 0].legend(handles, labels, loc='upper left', fontsize='small')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'comp_GPEC.pdf'), dpi=300)
testcase.close_datafile()
reference.close_datafiles()

testcase = Mephit(run_dir + '/TCFP')
testcase.open_datafile()
testcase.postprocess()
mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/postprocess/Bmn/coeff_rad', rad=True)
mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/postprocess/Bmn_vac/coeff_rad', rad=True)
reference = Kilca(testcase.work_dir)
reference.open_datafile('TCFP_flre_hip.hdf5')
kilca_Bmn = reference.get_polmodes('full perturbation (KiLCA)')
reference.open_datafile('TCFP_vac_hip.hdf5')
kilca_Bmn_vac = reference.get_polmodes('vacuum perturbation (KiLCA)')
reference.close_datafile()
config = {
    'xlabel': r'$r$ / \si{\centi\meter}', 'rho': testcase.data['/cache/fs/rad'][()],
    'ylabel': r'$\abs\, B_{m, n}^{r}$ / \si{\tesla}', 'comp': abs,
    'resonances': testcase.post['rad_res'], 'sgn_m_res': testcase.post['sgn_m_res'],
    'poldata': [mephit_Bmn_vac, kilca_Bmn_vac, mephit_Bmn, kilca_Bmn]
}
fig = Figure(figsize=squares)
axs = fig.subplots(vert_plot, horz_plot, sharex='all', sharey='all')
for m_abs in range(m_min, m_min + vert_plot):
    for k in range(horz_plot):
        m = (2 * k - 1) * m_abs * config['sgn_m_res']
        axs[m_abs - m_min, k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
        if m in config['resonances']:
            axs[m_abs - m_min, k].axvline(config['resonances'][m],
                                          color='b', alpha=0.5, lw=0.5, label='resonance position')
        for data in config['poldata']:
            if m in data['var']:
                axs[m_abs - m_min, k].plot(data['rho'][m], config['comp'](data['var'][m]), label=data['label'])
            else:
                next(axs[m_abs - m_min, k]._get_lines.prop_cycler)
        axs[m_abs - m_min, k].set_title(('resonant ' if m in config['resonances'] else 'non-resonant ') + fr"$m = {m}$")
        axs[m_abs - m_min, k].set_xlabel(config['xlabel'])
        axs[m_abs - m_min, k].set_ylabel(config['ylabel'])
        if m_abs > m_min:
            axs[m_abs - m_min, k].xaxis.set_tick_params(labelbottom=True)
    axs[m_abs - m_min, 1].yaxis.set_tick_params(labelleft=True)
handles, labels = axs[0, 1].get_legend_handles_labels()
axs[0, 0].legend(handles, labels, loc='upper left', fontsize='small')
canvas = FigureCanvas(fig)
fig.savefig(path.join(work_dir, 'comp_KiLCA.pdf'), dpi=300)
testcase.close_datafile()
