from os import environ, path
from numpy import abs, amin, amax, arange, real, imag, pi
from scipy.constants import c as clight
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mephit_plot import Mephit, PolmodePlots, set_matplotlib_defaults


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
label = r'\mu_{0} \, \delta J^{\parallel} B_{0}^{-1}'
conversion = 4.0 * pi / clight
units = r'\per\meter'
jmnpar_Bmod = [
    testcase.get_polmodes('P--S excl. damping', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversion, L1=True),
    testcase.get_polmodes('P--S incl. damping', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversion, L1=True),
    testcase.get_polmodes('KiLCA', '/debug_KiLCA/jmnpar_Bmod_KiLCA/coeff', conversion, L1=True),
    testcase.get_polmodes('total', '/debug_KiLCA/jmnpar_Bmod_total/coeff', conversion, L1=True)
]

k = 0
for m in testcase.post['m_res']:
    fig = Figure(figsize=(6.6, 3.3))
    ax = fig.subplots()
    mask = (jmnpar_Bmod[0]['rho'][m] >= res[k] - 4.5 * delta_mn[k]) & \
        (jmnpar_Bmod[0]['rho'][m] <= res[k] + 4.5 * delta_mn[k])
    for polmode in jmnpar_Bmod:
        ax.semilogy(polmode['rho'][m][mask], abs(polmode['var'][m][mask]), label=polmode['label'])
        ax.set_xlabel(r'$\hat{\psi}$')
        ax.set_ylabel(r'$\abs ' + label + r'$ / \si{' + units + '}')
        ax.axhline(0.0, color='k', alpha=0.5, lw=0.5)
        ax.axvline(res[k] - 0.5 * delta_mn[k], color='r', alpha=0.5, lw=0.5)
        ax.axvline(res[k] + 0.5 * delta_mn[k], color='r', alpha=0.5, lw=0.5)
        for pos in polmode['rho'][m][mask]:
            ax.axvline(pos, color='k', alpha=0.5, lw=0.25)
    ax.legend(loc='upper left')
    fig.suptitle(f"poloidal mode $m = {m}$")
    canvas = FigureCanvas(fig)
    fig.savefig(path.join(work_dir, f"plot_debug_jmnpar_Bmod_abs_{abs(m)}.pdf"))
    k += 1

mephit_Bmn = testcase.get_polmodes('full perturbation', '/postprocess/Bmn/coeff_rad')
mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation', '/postprocess/Bmn_vac/coeff_rad')
config = {
    'xlabel': r'$\hat{\psi}$',
    'ylabel': r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m}$ / \si{\tesla\meter\squared}',
    'rho': testcase.post['psi_norm'],
    'q': testcase.data['/cache/fs/q'][()],
    'resonances': testcase.post['psi_norm_res'],
    'sgn_m_res': testcase.post['sgn_m_res'],
    'omit_res': False,
    'poldata': [mephit_Bmn_vac, mephit_Bmn],
    'comp': abs,
    'res_neighbourhood': testcase.post['psi_norm_res_neighbourhood']
}
PolmodePlots(work_dir, 'debug_Bmn_psi_abs.pdf', config).do_plot()
