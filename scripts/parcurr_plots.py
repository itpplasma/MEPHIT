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
work_dir = path.join(environ['MEPHIT_RUN_DIR'], '33353_2900_EQH')  # AUG
# work_dir = path.join(environ['MEPHIT_DIR'], 'run/47046')  # MAST-U
mephit = Mephit(work_dir, 'mephit_imhd.h5')
mephit.open_datafile()
mephit.postprocess()
#gpec = Gpec(mephit.work_dir, mephit.data['/config/n'][()])
#gpec.open_datafiles()
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
    # mephit.get_polmodes('iMHD (undamped)', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
]
mephit_Ires = mephit.get_Ires()
#gpec_Ires = gpec.get_Ires()

# %%
mephit.close_datafile()
#gpec.close_datafiles()

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
    ax.axvline(res[k] - delta_mn[k], color='k', lw=0.25, ls='--')
    ax.axvline(res[k] + delta_mn[k], color='k', lw=0.25, ls='--')
    for polmode in jmnpar_Bmod:
        ax.plot(polmode['rho'][m][mask], np.abs(polmode['var'][m][mask]), label=polmode['label'])
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
#ax.semilogy(gpec_Ires.keys(), gpec_Ires.values(), 'x', label='GPEC')
ax.set_xlabel('resonant poloidal mode number $m$')
ax.set_ylabel(r'$\lvert I_{mn}^{\parallel} \rvert$ / \si{\ampere}')
ax.legend(fontsize='small')
fig.savefig(path.join(work_dir, 'Ipar.pdf'), backend='pgf', dpi=150)
plt.show()

#%%
plt.plot(psi, rad, 'o', label='rad')
plt.plot(psi, rsmall, 'o', label='rsmall')
plt.legend()
plt.show()

# %%
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def psi_to_rad_interpolator(psi_vals, rad_vals, kind='linear'):
    return interp1d(psi_vals, rad_vals, kind=kind)

psi_to_rad = psi_to_rad_interpolator(psi/psi[-1], rad)

m=-3
k = abs(m) - m_res_min
mask = (jmnpar_Bmod[1]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
    (jmnpar_Bmod[1]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
x = psi_to_rad(jmnpar_Bmod[1]['rho'][m][mask])
y = np.abs(jmnpar_Bmod[1]['var'][m][mask])
x_cut = np.concatenate([x[:6], x[-6:]])
y_cut = np.concatenate([y[:6], y[-6:]])

x0 = psi_to_rad(res[k])
def gauss(x, a, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
p0 = [y.max(), delta_mn[k]]
popt, _ = curve_fit(gauss, x, y, p0=p0)
print(f"m = {m}, a = {popt[0]}, sigma = {popt[1]}")
popt_cut, _ = curve_fit(gauss, x_cut, y_cut, p0=p0)
print(f"m = {m}, a_cut = {popt_cut[0]}, sigma_cut = {popt_cut[1]}")
x_fit = np.linspace(x.min(), x.max(), 300)
plt.plot(x_fit, gauss(x_fit, *popt), 'r--', label='Gauss fit')
plt.plot(x_fit, gauss(x_fit, *popt_cut), 'b--', label='Gauss fit cut')
plt.plot(x, y, label=f"{jmnpar_Bmod[1]['label']}")
plt.scatter(x_cut, y_cut, label=f"{jmnpar_Bmod[1]['label']} data")
plt.legend()
plt.show()

for m in mephit.post['m_res'][1:]:
    if abs(m) > 12:
        break
    k = abs(m) - m_res_min
    mask = (jmnpar_Bmod[1]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
        (jmnpar_Bmod[1]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
    x = psi_to_rad(jmnpar_Bmod[1]['rho'][m][mask])
    y = np.abs(jmnpar_Bmod[1]['var'][m][mask])

    x0 = psi_to_rad(res[k])
    def gauss(x, a, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
    p0 = [y.max(), delta_mn[k]]
    popt, _ = curve_fit(gauss, x, y, p0=p0)
    print(f"m = {m}, a = {popt[0]}, sigma = {popt[1]}")
    x_fit = np.linspace(x.min(), x.max(), 300)
    plt.plot(x_fit, gauss(x_fit, *popt), 'r--', label='Gauss fit')
    plt.plot(x, y, label=f"{jmnpar_Bmod[1]['label']}")
    plt.axvline(x0+3*popt[1], color='k', lw=0.5, ls='--')
    plt.axvline(x0-3*popt[1], color='k', lw=0.5, ls='--')
    plt.scatter(x, y, label=f"{jmnpar_Bmod[1]['label']} data")
    plt.legend()
    plt.show()

m=-11
k = abs(m) - m_res_min
mask = (jmnpar_Bmod[1]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
    (jmnpar_Bmod[1]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
x = psi_to_rad(jmnpar_Bmod[1]['rho'][m][mask])
y = np.abs(jmnpar_Bmod[1]['var'][m][mask])
x_cut = x[22:]
y_cut = y[22:]

x0 = psi_to_rad(res[k])
def gauss(x, a, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
p0 = [y.max(), delta_mn[k]]
popt, _ = curve_fit(gauss, x_cut, y_cut, p0=p0)
print(f"m = {m}, a = {popt[0]}, sigma = {popt[1]}")
x_fit = np.linspace(x.min(), x.max(), 300)
plt.plot(x_fit, gauss(x_fit, *popt), 'r--', label='Gauss fit')
plt.plot(x, y, label=f"{jmnpar_Bmod[1]['label']}")

plt.scatter(x_cut, y_cut, label=f"{jmnpar_Bmod[1]['label']} data")
plt.legend()
plt.show()