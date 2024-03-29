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
#     language: python
#     name: python3
# ---

# %%
from os import environ, path
from numpy import abs, pi, sign
from numpy.fft import fft
from scipy.constants import c as clight
from scipy.io import loadmat
import matplotlib.pyplot as plt
from mephit_plot import Mephit, Gpec, set_matplotlib_defaults


# %%
set_matplotlib_defaults()
work_dir = path.join(environ['MEPHIT_DIR'], 'run/47046')
mephit = Mephit(work_dir, 'mephit.h5')
mephit.open_datafile()
mephit.postprocess()
ref_dir = path.join(environ['HOME'], 'TU/PhD/MARS_MEPHIT/forPatrick')
mars = loadmat(ref_dir + '/OUTPUTS.mat') | loadmat(ref_dir + '/INPUTS.mat', squeeze_me=True) | loadmat(ref_dir + '/Normalizations.mat')
gpec = Gpec(work_dir, 2)
gpec.open_datafiles()

# %%
m_res_min = mephit.data['/mesh/m_res_min'][()]
kf_res = mephit.data['/mesh/res_ind'][()]
nflux = mephit.data['/mesh/nflux'][()]
res = mephit.normalize_psi(mephit.data['/mesh/psi_res'][()])
delta_mn = mephit.normalize_psi_diff(mephit.data['/mesh/delta_psi_mn'][()])
q = mephit.data['/cache/fs_half/q'][()]
sgn_dpsi = sign(mephit.data['/cache/fs/psi'][-1] - mephit.data['/cache/fs/psi'][0])
conversion = 4.0 * pi / clight
jmnpar_Bmod = [
    mephit.get_polmodes('total', '/debug_KiLCA/jmnpar_Bmod_total/coeff', conversion, L1=True),
    mephit.get_polmodes('gyrokinetic', '/debug_KiLCA/jmnpar_Bmod_KiLCA/coeff', conversion, L1=True),
    mephit.get_polmodes('iMHD', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversion, L1=True),
    mephit.get_polmodes('iMHD (undamped)', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
]
print(mars.keys())
print(mars['VACUUM']['B'][0, 0].shape)
print(mars['Mag_field'][0, 0].shape)
mars_vacuum = {'m_max': 0, 'label': 'MARS vacuum perturbation', 'rho': dict(), 'var': dict()}
mars_vacuum_Bmn = fft(mars['VACUUM']['B'][0, 0] * mars['Jacobian'], norm='forward')
for m in mars['Mm']:
    mars_vacuum['rho'][-m] = mars['sp']
    mars_vacuum['var'][-m] = mars_vacuum_Bmn[0, :mars['sp'].size, m] * mars['Mag_field'][0, 0] * 1.0e4
mars_plasma = {'m_max': 0, 'label': 'MARS full perturbation', 'rho': dict(), 'var': dict()}
mars_plasma_Bmn = fft(mars['PLASMA']['B'][0, 0] * mars['Jacobian'], norm='forward')
for m in mars['Mm']:
    mars_plasma['rho'][-m] = mars['sp']
    mars_plasma['var'][-m] = mars_plasma_Bmn[0, :mars['sp'].size, m] * mars['Mag_field'][0, 0] * 1.0e4
conversion = 1.0e-04 / mephit.data['/mesh/gpec_jacfac'][()]  # jacfac not accounted for in MARS output
Bmn = [
    mephit.get_polmodes('MEPHIT vacuum perturbation', '/iter/Bmn_vac/coeff_rad', conversion),
    mephit.get_polmodes('MEPHIT full perturbation', '/iter/Bmn/coeff_rad', conversion),
    # mars_vacuum,
    # mars_plasma,
    gpec.get_polmodes('GPEC full perturbation', sgn_dpsi, 'Jbgradpsi'),
    gpec.get_polmodes('GPEC vacuum perturbation', sgn_dpsi, 'Jbgradpsi_x'),
]

# %%
for m in mephit.post['m_res']:
    if abs(m) > 10:
        continue
    k = abs(m) - m_res_min
    fig = plt.figure(layout='constrained')
    ax = fig.subplots()
    mask = (jmnpar_Bmod[0]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
        (jmnpar_Bmod[0]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
    ax.axvline(res[k], color='k', lw=0.5)
    ax.axvline(res[k] - 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    ax.axvline(res[k] + 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    for polmode in jmnpar_Bmod:
        ax.semilogy(polmode['rho'][m][mask], abs(polmode['var'][m][mask]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\mu_{0} \lvert (j_{\parallel} / B_{0})_{\V{m}} \rvert$ [\si{\per\meter}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

    fig = plt.figure(layout='constrained')
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5)
    for polmode in Bmn:
        mask = (polmode['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
            (polmode['rho'][m] <= res[k] + 4.5 * delta_mn[k])
        ax.plot(polmode['rho'][m], abs(polmode['var'][m]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (\sqrt{g} B^{\psi})_{\V{m}} \rvert A^{-1}$ [\si{\tesla}]')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    plt.show()

