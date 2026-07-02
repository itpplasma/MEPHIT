# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
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
from scipy.interpolate import RectBivariateSpline
from mephit_plot import Mephit, Gpec, Mars
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')


# %%
work_dir = path.join('/temp/zenzphil/MEPHIT/run', '33353_2900_EQH_sweep')  # AUG
# work_dir = path.join(environ['MEPHIT_RUN_DIR'], '47046')  # MAST-U
mephit = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_par.h5')
mephit_D_40 = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_eresoff_7.h5')
mephit_H_40 = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_H_eresoff_7.h5')
mephit_ref = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_D_eres_7.h5')
mephit_thin = Mephit(work_dir, 'mephit_cross_fade_new_int_4_0_H_eres_7.h5')
mephit_D_05_off = Mephit(work_dir, 'mephit_cross_fade_new_int_0_5_D_eresoff_7.h5')
mephit_H_05_off = Mephit(work_dir, 'mephit_cross_fade_new_int_0_5_H_eresoff_7.h5')
mephit_D_05_on  = Mephit(work_dir, 'mephit_cross_fade_new_int_0_5_D_eres_7.h5')
mephit_H_05_on  = Mephit(work_dir, 'mephit_cross_fade_new_int_0_5_H_eres_7.h5')

mephit.open_datafile()
mephit.postprocess()
mephit_D_40.open_datafile()
mephit_D_40.postprocess()
mephit_H_40.open_datafile()
mephit_H_40.postprocess()
mephit_ref.open_datafile()
mephit_ref.postprocess()
mephit_thin.open_datafile()
mephit_thin.postprocess()
mephit_D_05_off.open_datafile()
mephit_D_05_off.postprocess()
mephit_H_05_off.open_datafile()
mephit_H_05_off.postprocess()
mephit_D_05_on.open_datafile()
mephit_D_05_on.postprocess()
mephit_H_05_on.open_datafile()
mephit_H_05_on.postprocess()

# gpec = Gpec(mephit.work_dir, mephit.data['/config/n'][()])
# gpec.open_datafiles()
# ref_dir = path.join(environ['HOME'], 'TU/PhD/MARS_MEPHIT/forPatrick')
# mars = Mars(ref_dir)
# mars.open_datafiles()

# %%
n = mephit.data['/config/n'][()]
m_res_min = mephit.data['/mesh/m_res_min'][()]
m_res_max = mephit.data['/mesh/m_res_max'][()]
q = mephit.data['/cache/fs/q'][()]
res = mephit.normalize_psi(mephit.data['/mesh/psi_res'][()])
delta_mn = np.abs(mephit.normalize_psi_diff(mephit.data['/mesh/Delta_psi_res_curr'][()]))
sgn_dpsi = np.sign(mephit.data['/cache/fs/psi'][-1] - mephit.data['/cache/fs/psi'][0])
conversion = 1.0e-04 / mephit.data['/mesh/gpec_jacfac'][()]  # TODO: jacfac not accounted for in MARS output
Bmn = [
    mephit.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
    mephit.get_polmodes('MEPHIT full perturbation covariant pol', '/iter/Bmn/coeff_pol', 1.0e-04),
    mephit.get_polmodes('MEPHIT full perturbation tor', '/iter/Bmn/coeff_tor', 1.0e-04),
    mephit.get_polmodes('MEPHIT full perturbation pol', '/iter/Bmn/coeff_pol', 1.0e-04),
    mephit.get_polmodes('MEPHIT full perturbation par', '/iter/Bmn/coeff_par', 1.0e-04),
    # gpec.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi'),
    # gpec.get_polmodes('GPEC vacuum perturbation', sgn_dpsi, 'Jbgradpsi_x'),
    # mars.get_polmodes('MARS vacuum perturbation', 'VACUUM'),
    # mars.get_polmodes('MARS full perturbation', 'PLASMA'),
]
B_equil = {
    "B0_R": mephit.data['/equil/B0_R'][:],
    "B0_Z": mephit.data['/equil/B0_Z'][:],
    "B0_phi": mephit.data['/equil/B0_phi'][:],
    "R_rect": mephit.data['/equil/rect_R'][:],
    "Z_rect": mephit.data['/equil/rect_Z'][:],
    "psirz":  mephit.data['/equil/psirz'][:,:],
    "Z_eqdsk": mephit.data['/equil/Z_eqd'][:],
    "R_eqdsk": mephit.data['/equil/R_eqd'][:],
    "R_mesh": mephit.data['/mesh/node_R'][:],
    "Z_mesh": mephit.data['/mesh/node_Z'][:],
    "tri": mephit.data['/mesh/tri_node'][:]-1,
    "Bn_theta": mephit.data['/iter/Bn/comp_theta_covar'][:],
    "Bn_phi": mephit.data['/iter/Bn/RT0_comp_phi'][:],
    "psi_mesh": mephit.data['/cache/fs/psi'][:],
    "rad": mephit.data['/cache/fs_half/rad'][:],
    "q": mephit.data['/cache/fs/q'][:],
    "R0": mephit.data['/equil/rcentr'][()],
}

Bmn_D_40 = [
    mephit_D_40.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit_D_40.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]

Bmn_H_40 = [
    #mephit_H_40.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit_H_40.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_ref = [
    #mephit_ref.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit_ref.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_thin = [
    #mephit_thin.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit_thin.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_D_05_off = [
    mephit_D_05_off.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_H_05_off = [
    mephit_H_05_off.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_D_05_on = [
    mephit_D_05_on.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]
Bmn_H_05_on = [
    mephit_H_05_on.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
]

for dataset in [Bmn, Bmn_D_40, Bmn_H_40, Bmn_ref, Bmn_thin, Bmn_D_05_off, Bmn_H_05_off, Bmn_D_05_on, Bmn_H_05_on]:
    for polmode in dataset:
        polmode['m_range'] = np.arange(min(polmode['var'].keys()), max(polmode['var'].keys()) + 1)
        polmode['arr'] = np.zeros((polmode['m_range'].size, polmode['rho'][0].size), dtype=complex)
        for i, m in enumerate(polmode['m_range']):
            polmode['arr'][i, :] = polmode['var'][m]


# %%
mephit.close_datafile()
mephit_H_40.close_datafile()
mephit_ref.close_datafile()
mephit_thin.close_datafile()
mephit_D_05_off.close_datafile()
mephit_H_05_off.close_datafile()
mephit_D_05_on.close_datafile()
mephit_H_05_on.close_datafile()
# gpec.close_datafiles()

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
#plt.show(block=False)

# %%
"""fig = plt.figure(figsize=(6.6, 3.6), dpi=150)
axs = fig.subplots(1, 2, sharex='all', sharey='all')
ims = [None, None]
for i in range(2):
    # ims[i] = axs[i].pcolormesh(*np.meshgrid(Bmn[i]['m_range'], Bmn[i]['rho'][0], indexing='ij'),
    #                             np.abs(Bmn[i]['arr']), shading='nearest', cmap=cc.m_fire_r)
    ims[i] = axs[i].contourf(*np.meshgrid(Bmn_ref[i]['m_range'], Bmn_ref[i]['rho'][0], indexing='ij'),
                              np.abs(Bmn_ref[i]['arr']), levels=256, cmap=cc.m_fire_r)
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
fig.savefig(path.join(work_dir, 'Bn_spectrum_ref.png'))
plt.show(block=False)"""

# %%
for m in mephit.post['m_res']:
    if abs(m) > 12:
        break
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5, ls=':')

    # Define all simulations you want to plot with colors, styles, and labels
    polmode_list = [
        (Bmn_D_40[0], 'yellow', '-', '4.0cm D off'),
        (Bmn_D_40[1], 'k', '-', '4.0cm D off'),
        *[(pm, 'tab:orange', '-', '4.0cm H off') for pm in Bmn_H_40],
        *[(pm, 'k', '--', '4.0cm D on') for pm in Bmn_ref],
        *[(pm, 'tab:orange', '--', '4.0cm D on') for pm in Bmn_thin],
        *[(pm, 'lightgreen', '-', '0.5cm D off') for pm in Bmn_D_05_off],
        *[(pm, 'tab:red', '-', '0.5cm H off') for pm in Bmn_H_05_off],
        *[(pm, 'lightgreen', '--', '0.5cm D on') for pm in Bmn_D_05_on],
        *[(pm, 'tab:red', '--', '0.5cm H on') for pm in Bmn_H_05_on],
    ]

    for polmode, color, style, label_suffix in polmode_list:
        ax.plot(polmode['rho'][m], np.abs(polmode['var'][m]), label=polmode['label'] + ' ' + label_suffix, ls=style, color=color)
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi}_{n})_{m} \rvert A^{-1}$ / \si{\tesla}')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    fig.savefig(path.join(work_dir, f'Bmn_psi_{abs(m)}.pdf'), dpi=150)

    psi = Bmn[0]['rho'][m]
    n_points = len(psi)
    ydata = Bmn[1]['var'][m]
    np.savetxt(f'Bmn_psi_{np.abs(m)}_4_0_data.txt', np.vstack(psi, ydata.real, ydata.imag),
               fmt='%.8e', header='psi B_psi_real B_psi_imag')

    """psi = polmode_list[0][0]['rho'][m]
    n_points = len(psi)
    ydata = [np.abs(pm['var'][m]) for pm, _, _, _ in polmode_list]
    with open(path.join(work_dir, f'Bmn_psi_{abs(m)}_data.txt'), 'w') as f:
        # Header
        header = "'psi' " + " ".join([f"'{pm}'" for _, _, _, pm in polmode_list]) + "\n"
        f.write(header)
        # Data rows
        for i in range(n_points):
            row = [f"{psi[i]:.8e}"] + [f"{y[i]:.8e}" for y in ydata]
            f.write(" ".join(row) + "\n")

    ax.set_xlim(res[k]-0.3 * delta_mn[k], res[k] + 0.3 * delta_mn[k])
    xlim = ax.get_xlim()
    ymin, ymax = np.inf, -np.inf
    for polmode in [Bmn[1]] + Bmn_H_40 + Bmn_ref + Bmn_thin + Bmn_D_05_off + Bmn_H_05_off + Bmn_D_05_on + Bmn_H_05_on:
        x = polmode['rho'][m]
        y = np.abs(polmode['var'][m])
        mask = (x >= xlim[0]) & (x <= xlim[1])
        if np.any(mask):
            ymin = min(ymin, np.min(y[mask]))
            ymax = max(ymax, np.max(y[mask]))
    ax.set_ylim(ymin, 1.1*ymax)

    ylims = ax.get_ylim()
    ax.set_ylim(ymin,ymax)
    ax.legend(loc='upper left', fontsize='small')

    fig.savefig(path.join(work_dir, f'Bmn_psi_{abs(m)}_zoom.pdf'), dpi=150)"""

    plt.show(block=False)

# %%
for m in mephit.post['sgn_m_res'] * np.arange(6):
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    for polmode in [Bmn[0], Bmn[1]]:
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

# %%
B_equil["mod_B"] = np.sqrt(B_equil["B0_R"]**2 + B_equil["B0_Z"]**2 + B_equil["B0_phi"]**2)
B_equil["h_R"] = B_equil["B0_R"] / B_equil["mod_B"]
B_equil["h_Z"] = B_equil["B0_Z"] / B_equil["mod_B"]
B_equil["h_phi"] = B_equil["B0_phi"] / B_equil["mod_B"]
B_equil["B0_phi_contrav"] = B_equil["B0_phi"] / (B_equil["R_rect"])

interp_psi = RectBivariateSpline(B_equil["Z_eqdsk"], B_equil["R_eqdsk"], B_equil["psirz"])
B_equil["psi_rect"] = interp_psi(B_equil["Z_rect"], B_equil["R_rect"])

B_equil["q_RZ_eqdsk"] = np.interp(B_equil["psirz"].ravel(),
                                  np.flip(B_equil["psi_mesh"]),
                                  np.flip(B_equil["q"])).reshape(B_equil["psirz"].shape)
interp_q = RectBivariateSpline(B_equil["Z_eqdsk"], B_equil["R_eqdsk"], B_equil["q_RZ_eqdsk"])
B_equil["q_RZ_rect"] = interp_q(B_equil["Z_rect"], B_equil["R_rect"])

B_equil["B0_theta_contrav"] = B_equil["B0_phi_contrav"] / B_equil["q_RZ_rect"]
B_equil["h_theta_contrav"] = B_equil["B0_theta_contrav"] / B_equil["mod_B"]

R, Z = np.meshgrid(B_equil["R_rect"], B_equil["Z_rect"])
fig, ax = plt.subplots(figsize=(7, 6))
c = ax.pcolormesh(R, Z, B_equil["mod_B"], shading='nearest', cmap=cc.m_fire_r)
ax.set_xlabel(r'R / \si{\cm}')
ax.set_ylabel(r'Z / \si{\cm}')
ax.set_aspect('equal', 'box')
cbar = fig.colorbar(c)
cbar.set_label(r'$\lvert \vec{B}_0 \rvert$ / \si{\tesla}', rotation=90)
fig.savefig(path.join(work_dir, 'B_equil.png'), dpi=150)
plt.show(block=False)

q_levels = np.array([1.0, 2.0, 3.0, 4.0, 7.0])
q_fmt = {level: f'$q = {int(level)}$' for level in q_levels}
R_eq, Z_eq = np.meshgrid(B_equil["R_eqdsk"], B_equil["Z_eqdsk"])
fig, ax = plt.subplots(figsize=(7, 6))
#ax.contour(R_eq, Z_eq, B_equil["psirz"], levels=20, colors='k')
ax.contour(R, Z, B_equil["psi_rect"], levels=50, colors='r')
#ax.contour(R_eq, Z_eq, B_equil["psirz"], levels=np.flip(B_equil["psi_eqdsk"]), colors='b', linewidths=2.0)
q_contour = ax.contour(R_eq, Z_eq, B_equil["q_RZ_eqdsk"], levels=q_levels, colors='b', linestyles='dashed')
ax.clabel(q_contour, fmt=q_fmt, inline=True, fontsize=8)
# q_contour_rect = ax.contour(R, Z, B_equil["q_RZ_rect"], levels=q_levels, colors='cyan', linestyles='dashed')
# ax.clabel(q_contour_rect, fmt=q_fmt, inline=True, fontsize=8)
ax.set_xlabel(r'R / \si{\cm}')
ax.set_ylabel(r'Z / \si{\cm}')
ax.set_aspect('equal', 'box')

# %%
B0_R_spline     = RectBivariateSpline(B_equil["R_rect"], B_equil["Z_rect"], B_equil["h_R"].T)
B0_theta_spline = RectBivariateSpline(B_equil["R_rect"], B_equil["Z_rect"], B_equil["h_theta_contrav"].T)
B0_phi_spline   = RectBivariateSpline(B_equil["R_rect"], B_equil["Z_rect"], B_equil["h_phi"].T)

def RZ_from_psi_theta(psi, theta, R0=2.0):
    r = np.sqrt(psi)
    R = R0 + r*np.cos(theta)
    Z = r*np.sin(theta)
    return R, Z

# Poloidal mode numbers
m_modes = np.arange(-24, 25)
n_modes = 2

# 1D flux surface grid
print(Bmn[1]['rho'][0].shape)
psi_grid = Bmn[1]['rho'][0]
psi_grid_norm = (psi_grid - psi_grid[0]) / (psi_grid[-1] - psi_grid[0])

# Poloidal angle grid for reconstruction
Ntheta = 2048
theta_grid = np.linspace(0, 2*np.pi, Ntheta, endpoint=False)

B0_R_theta = np.zeros((len(psi_grid), Ntheta))
B0_theta_theta = np.zeros_like(B0_R_theta)
B0_phi_theta = np.zeros_like(B0_R_theta)

for i_psi, psi in enumerate(psi_grid_norm):
    for i_theta, theta in enumerate(theta_grid):
        R, Z = RZ_from_psi_theta(psi, theta)
        B0_R_theta[i_psi, i_theta] = B0_R_spline.ev(R, Z)
        B0_theta_theta[i_psi, i_theta] = B0_theta_spline.ev(R, Z)
        B0_phi_theta[i_psi, i_theta] = B0_phi_spline.ev(R, Z)

# Fourier transform B0 along theta
B0_R_modes = np.fft.fft(B0_R_theta, axis=1) / Ntheta
B0_theta_modes = np.fft.fft(B0_theta_theta, axis=1) / Ntheta
B0_phi_modes = np.fft.fft(B0_phi_theta, axis=1) / Ntheta

# Map FFT frequencies to poloidal mode numbers
fft_freqs = np.fft.fftfreq(Ntheta, d=1/Ntheta)

for key in Bmn[2]['var'].keys():
    Bmn[4]['var'][key] = Bmn[2]['var'][key] / (B_equil["R0"] + B_equil["rad"])

def dict_to_mode_array(mode_dict, m_modes, psi_grid):
    """
    Convert dict of {m: array(n_psi)} into a 2D NumPy array (n_psi, n_m).

    Parameters
    ----------
    mode_dict : dict
        Keys are integer poloidal mode numberswith open(path.join(work_dir, f'Bmn_par_{np.abs(m_modes[i])}_0_1_data.txt'), 'w') as f:
        # Header
        header = "psi     B_par_real     B_par_imag" + "\n"
        f.write(header)
        # Data rows
        for i in range(n_points):
            row = [f"{psi[i]:.8e}"] + [f"{np.real(ydata[i]):.8e}"] + [f"{np.imag(ydata[i]):.8e}"]
            f.write(" ".join(row) + "\n")
plt.show() (m), values are arrays of length n_psi.
    m_modes : list or array
        List of mode numbers (ordering you want for columns).
    psi_grid : array
        The flux surface grid, only used to set n_psi.

    Returns
    -------
    arr : ndarray, shape (n_psi, len(m_modes))
        Array of modes ordered by m_modes.
    """
    n_psi = len(psi_grid)
    arr = np.zeros((n_psi, len(m_modes)), dtype=complex)  # perturbations are usually complex
    for j, m in enumerate(m_modes):
        if m not in mode_dict:
            raise KeyError(f"Mode m={m} not found in dictionary keys {list(mode_dict.keys())}")
        vals = np.asarray(mode_dict[m])
        if vals.shape[0] != n_psi:
            raise ValueError(f"Mode m={m} has length {vals.shape[0]}, expected {n_psi}")
        arr[:, j] = vals
    return arr

B_theta_modes = dict_to_mode_array(Bmn[4]['var'], m_modes, psi_grid)
B_phi_modes   = dict_to_mode_array(Bmn[3]['var'], m_modes, psi_grid)


B_parallel_modes = np.zeros((len(psi_grid), len(m_modes)), dtype=complex)
print(B_parallel_modes.shape)
for i_psi in range(len(psi_grid)):
    for i_m, m in enumerate(m_modes):
        # convolution: sum over m' (B0^{m-m'} * B1^{m'})
        B_par = 0.0+0.0j
        for i_mp, mp in enumerate(m_modes):
            # find closest index in FFT corresponding to m-mp
            idx = np.argmin(np.abs(fft_freqs - (m-mp)))
            B_par += (B0_theta_modes[i_psi,idx]*B_theta_modes[i_psi,i_mp] +
                      B0_phi_modes[i_psi,idx]*B_phi_modes[i_psi,i_mp])
        B_parallel_modes[i_psi,i_m] = B_par

for i in range (15,22):
    fig=plt.figure()
    plt.plot(Bmn[3]['rho'][i], np.abs(B_parallel_modes[:,i]), label=f'm={m_modes[i]}')
    plt.legend()
    plt.xlabel(r'$\hat{\psi}$')
    plt.ylabel(r'$\lvert (B^{\parallel}_{n})_{m} \rvert $ / \si{\tesla}')
    plt.legend(loc='upper left')
    #fig.savefig(path.join(work_dir, f'Bmn_par_{np.abs(m_modes[i])}_4_0.pdf'), dpi=150)

    """psi = Bmn[0]['rho'][i]
    n_points = len(psi)
    ydata = B_parallel_modes[:,i]
    with open(path.join(work_dir, f'Bmn_par_{np.abs(m_modes[i])}_4_0_data.txt'), 'w') as f:
        # Header
        header = "psi     B_par_real     B_par_imag" + "\n"
        f.write(header)
        # Data rows
        for i in range(n_points):
            row = [f"{psi[i]:.8e}"] + [f"{np.real(ydata[i]):.8e}"] + [f"{np.imag(ydata[i]):.8e}"]
            f.write(" ".join(row) + "\n")"""
plt.show()

# %%
print(len(Bmn[2]['var']))
for key in Bmn[2]['var'].keys():
    print(Bmn[2]['var'][key].shape)
    print(B_equil["R0"], B_equil["rad"][1:])
    Bmn[4]['var'][key] = Bmn[2]['var'][key] / (B_equil["R0"] + B_equil["rad"])

for m in mephit.post['m_res']:
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    polmode_list = [
        (Bmn[4], 'orange', '-', '0.1cm D off'),
        (Bmn[3], 'k', '-', '0.1cm D off'),
    ]

    for polmode, color, style, label_suffix in polmode_list:
        ax.plot(polmode['rho'][m], np.abs(polmode['var'][m]),
                label=polmode['label'] + ' ' + label_suffix, ls=style, color=color)

    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (B^{\theta, \phi}_{n})_{m} \rvert $ / \si{\tesla}')
    ax.set_title(f"$m = {m}$")
    ax.set_ylim(0,0.0005)
    ax.legend(loc='upper left')
    fig.savefig(path.join(work_dir, f'Bmn_theta_phi_{np.abs(m)}_0_1.pdf'), dpi=150)
    plt.show()


# %%
for m in mephit.post['m_res']:
    fig = plt.figure()
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    polmode_list = [
        (Bmn[5], 'orange', '-', '0.1cm D off'),
    ]

    for polmode, color, style, label_suffix in polmode_list:
        ax.plot(polmode['rho'][m], np.abs(polmode['var'][m]),
                label=polmode['label'] + ' ' + label_suffix, ls=style, color=color)

    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (B^{\theta, \phi}_{n})_{m} \rvert $ / \si{\tesla}')
    ax.set_title(f"$m = {m}$")
    ax.set_ylim(0, 0.0005)
    ax.legend(loc='upper left')
    plt.savefig(path.join(work_dir, f'Bmn_par_{np.abs(m)}_0_1.pdf'), dpi=150)
    plt.show()

