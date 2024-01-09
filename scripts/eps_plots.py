from os import environ, path
from numpy import abs, pi
from scipy.constants import c as clight
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mephit_plot import Mephit, set_matplotlib_defaults


set_matplotlib_defaults()
work_dir = path.join(environ['MEPHIT_DIR'], 'run/33353_2900_EQH')
testcase = Mephit(work_dir, 'mephit.h5')
testcase.open_datafile()
testcase.postprocess()
m_res_min = testcase.data['/mesh/m_res_min'][()]
kf_res = testcase.data['/mesh/res_ind'][()]
nflux = testcase.data['/mesh/nflux'][()]
res = testcase.normalize_psi(testcase.data['/mesh/psi_res'][()])
delta_mn = testcase.normalize_psi_diff(testcase.data['/mesh/delta_psi_mn'][()])
q = testcase.data['/cache/fs_half/q'][()]
conversion = 4.0 * pi / clight
jmnpar_Bmod = [
    testcase.get_polmodes('total', '/debug_KiLCA/jmnpar_Bmod_total/coeff', conversion, L1=True),
    testcase.get_polmodes('gyrokinetic', '/debug_KiLCA/jmnpar_Bmod_KiLCA/coeff', conversion, L1=True),
    testcase.get_polmodes('iMHD', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversion, L1=True),
    testcase.get_polmodes('iMHD (undamped)', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
]
Bmn = [
    testcase.get_polmodes('vacuum perturbation', '/iter/Bmn_vac/coeff_rad'),
    testcase.get_polmodes('full perturbation', '/iter/Bmn/coeff_rad')
]

for m in testcase.post['m_res']:
    k = abs(m) - m_res_min
    fig = Figure(figsize=(3, 2.5))
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
    ax.legend(loc='lower left')
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"plot_jmnpar_Bmod_abs_m{abs(m)}.pdf"))

    fig = Figure(figsize=(3, 2.5))
    ax = fig.subplots()
    ax.axhline(0.0, color='k', lw=0.5)
    ax.axvline(res[k], color='k', lw=0.5)
    for polmode in Bmn:
        ax.plot(polmode['rho'][m], abs(polmode['var'][m]) / q, label=polmode['label'])
    ax.set_xlabel(r'$\hat{\psi}$')
    ax.set_ylabel(r'$\lvert (B^{\psi} / B_{0}^{\varphi})_{\V{m}} \rvert$ [\si{\tesla\meter\squared}]')
    ax.legend(loc='upper left')
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"plot_Bmn_psi_abs_m{abs(m)}.pdf"))

