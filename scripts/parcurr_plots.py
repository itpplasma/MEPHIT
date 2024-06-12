# %%
# %matplotlib inline
from os import environ, getcwd, path
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c as clight
import h5py
from mephit_plot import Mephit, Gpec, Mars
plt.style.use('./mephit.mplstyle')
rcParams['text.latex.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pdflatex.tex}'
rcParams['pgf.preamble'] = r'\usepackage{import}\import{' + getcwd() + r'}{mephit-pgf.tex}'
h5py.get_config().complex_names = ('real', 'imag')


# %%
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')  # AUG
# work_dir = path.join(environ['MEPHIT_DIR'], 'run/47046')  # MAST-U
mephit = Mephit(work_dir, 'mephit.h5')
mephit.open_datafile()
mephit.postprocess()
gpec = Gpec(mephit.work_dir, mephit.data['/config/n'][()])
gpec.open_datafiles()
# ref_dir = path.join(environ['HOME'], 'TU/PhD/MARS_MEPHIT/forPatrick')
# mars = Mars(ref_dir)
# mars.open_datafiles()

# %%
m_res_min = mephit.data['/mesh/m_res_min'][()]
res = mephit.normalize_psi(mephit.data['/mesh/psi_res'][()])
delta_mn = mephit.normalize_psi_diff(mephit.data['/mesh/delta_psi_mn'][()])
conversion = 4.0 * np.pi / clight
jmnpar_Bmod = [
    mephit.get_polmodes('total', '/debug_KiLCA/jmnpar_Bmod_total/coeff', conversion, L1=True),
    mephit.get_polmodes('gyrokinetic', '/debug_KiLCA/jmnpar_Bmod_KiLCA/coeff', conversion, L1=True),
    mephit.get_polmodes('iMHD', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversion, L1=True),
    mephit.get_polmodes('iMHD (undamped)', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
]
mephit_Ires = mephit.get_Ires()
gpec_Ires = gpec.get_Ires()

# %%
mephit.close_datafile()
gpec.close_datafiles()

# %%
for m in mephit.post['m_res']:
    if abs(m) > 12:
        break
    k = abs(m) - m_res_min
    fig = plt.figure()
    ax = fig.subplots()
    mask = (jmnpar_Bmod[0]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
        (jmnpar_Bmod[0]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
    ax.axvline(res[k], color='k', lw=0.5)
    ax.axvline(res[k] - 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    ax.axvline(res[k] + 0.5 * delta_mn[k], color='k', lw=0.25, ls='--')
    for polmode in jmnpar_Bmod:
        ax.semilogy(polmode['rho'][m][mask], np.abs(polmode['var'][m][mask]), label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\mu_{0} \lvert (J_{n}^{\parallel} B_{0}^{-1})_{m} \rvert$ / \si{\per\meter}')
    ax.set_title(f"$m = {m}$")
    ax.legend(loc='upper left')
    fig.savefig(path.join(work_dir, f'parcurrmn_{abs(m)}.pdf'), backend='pgf', dpi=150)
    plt.show()

# %%
fig = plt.figure(figsize=(6.6, 3.6), dpi=150)
ax = fig.subplots()
ax.semilogy(mephit_Ires.keys(), mephit_Ires.values(), 'o', label='MEPHIT')
ax.semilogy(gpec_Ires.keys(), gpec_Ires.values(), 'x', label='GPEC')
ax.set_xlabel('resonant poloidal mode number $m$')
ax.set_ylabel(r'$\lvert I_{mn}^{\parallel} \rvert$ / \si{\ampere}')
ax.legend(fontsize='small')
fig.savefig(path.join(work_dir, 'Ipar.pdf'), backend='pgf', dpi=150)
plt.show()
