from mephit_plot import HLine, Id, IterationPlots, Mephit, ParallelPlotter, \
    run_dir, set_matplotlib_defaults, XTicks
from numpy import abs, arange, argmax, full, ix_, nan, pi
from scipy.constants import c as clight


set_matplotlib_defaults()
plotter = ParallelPlotter()
plotter.start()

work_dir = run_dir + '/33353_2900'
testcase = Mephit(work_dir)
testcase.open_datafile()
testcase.postprocess()

conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][()]
mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/iter/Bmn/coeff_rad', conversion)
mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/iter/Bmn_vac/coeff_rad', conversion)
n = testcase.data['/config/n'][()]
m_res_min = testcase.data['/mesh/m_res_min'][()]
m_res_max = testcase.data['/mesh/m_res_max'][()]
niter = testcase.data['/iter/niter'][()]
nflux = testcase.data['/mesh/nflux'][()]
abs_pmn_iter = full((2 * m_res_max + 1, niter, nflux + 1), nan)
abs_jmnpar_Bmod_iter = full((2 * m_res_max + 1, niter, nflux + 1), nan)
abs_Bmn_rad_iter = full((2 * m_res_max + 1, niter, nflux), nan)
kf_res = testcase.data['/mesh/res_ind'][()]
shielding = full((niter, m_res_max - m_res_min + 1), nan)
penetration = full((niter, m_res_max - m_res_min + 1), nan)
Ires = full((niter, m_res_max - m_res_min + 1), nan)
for kiter in range(0, niter):
    pmn = testcase.get_polmodes(None, f"/iter/pmn_{kiter:03}/coeff", 0.1, L1=True)
    jmnpar_Bmod = testcase.get_polmodes(None, f"/iter/jmnpar_Bmod_{kiter:03}/coeff", 4.0 * pi / clight, L1=True)
    Bmn_rad = testcase.get_polmodes(None, f"/iter/Bmn_{kiter:03}/coeff_rad", conversion)
    Ires[kiter, :] = abs(testcase.data[f"/iter/Ires_{kiter:03}"][()]) * 0.1 / clight
    for m in range(-m_res_max, m_res_max + 1):
        abs_pmn_iter[m + m_res_max, kiter, :] = abs(pmn['var'][m])
        abs_jmnpar_Bmod_iter[m + m_res_max, kiter, :] = abs(jmnpar_Bmod['var'][m])
        abs_Bmn_rad_iter[m + m_res_max, kiter, :] = abs(Bmn_rad['var'][m])
    for m in range(m_res_min, m_res_max + 1):
        m_res = m * testcase.post['sgn_m_res']
        kf_max = argmax(abs_Bmn_rad_iter[m_res + m_res_max, kiter, :kf_res[m - m_res_min] - 1]) + 1
        penetration[kiter, m - m_res_min] = abs_Bmn_rad_iter[m_res + m_res_max, kiter, kf_max - 1] / \
            abs(mephit_Bmn_vac['var'][m_res][kf_max - 1])
        shielding[kiter, m - m_res_min] = abs_Bmn_rad_iter[m_res + m_res_max, kiter, kf_res[m - m_res_min] - 1] / \
            abs(mephit_Bmn_vac['var'][m_res][kf_res[m - m_res_min] - 1])
config = {
    'xlabel': r'$\hat{\psi}$',
    'ylabel': [
        r'$\abs\, p_{mn}$ / \si{\pascal}',
        r'$\abs\, \bigl[ \mu_{0} J_{n}^{\parallel} B_{0}^{-1} \bigr]_{m}$ / \si{\per\meter}',
        r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}'
    ],
    'rho': [testcase.post['psi_norm'], testcase.post['psi_norm'], testcase.post['psi_half_norm']],
    'yscale': ['log', 'log', 'linear'], 'global_ylims': True, 'plotargs': {'lw': 0.5}, 'niter': niter,
    'postprocess': [[Id(), Id(), HLine(0.0, color='k', alpha=0.5, lw=0.5)]]
}
m_res_range = arange(m_res_min, m_res_max + 1)
for m in m_res_range:
    m_res = m * testcase.post['sgn_m_res']
    config['title'] = f"Poloidal Modes in Preconditioned Iterations $k$ for $n = {n}$, $m = \pm {m}$"
    config['res_pos'] = testcase.post['psi_norm_res'][m_res]
    config['res_neighbourhood'] = testcase.post['psi_norm_res_neighbourhood'][m_res]
    config['zoom_x'] = config['res_neighbourhood'][ix_([0, -1])]
    config['plotdata'] = [
        abs_pmn_iter[m_res + m_res_max, :, :],
        abs_jmnpar_Bmod_iter[m_res + m_res_max, :, :],
        abs_Bmn_rad_iter[m_res + m_res_max, :, :],
        full(abs_pmn_iter.shape[1:], nan),  # abs_pmn_iter[-m_res + m_res_max, :, :],
        abs_jmnpar_Bmod_iter[-m_res + m_res_max, :, :],
        # abs_Bmn_rad_iter[-m_res + m_res_max, :, :]
    ]
    plotter.plot_objects.put(IterationPlots(work_dir, f"plot_iter_{m}.pdf", config))
config = {
    'title': f"Shielding \& max. penetration vs. resonant current for $n = {n}$",
    'xlabel': r'poloidal mode $m$',
    'ylabel': ['shielding', 'max. penetration', r'$\abs\, I_{mn}^{\parallel}$ / \si{\ampere}'],
    'rho': [m_res_range, m_res_range, m_res_range], 'plotdata': [shielding, penetration, Ires],
    'yscale': ['log', 'log', 'log'], 'global_ylims': True, 'plotargs': {'ls': '', 'marker': 'x'}, 'niter': niter,
    'postprocess': [[XTicks(m_res_range), XTicks(m_res_range), XTicks(m_res_range)],
                    [Id(), Id(), HLine(0.0, color='k', alpha=0.5, lw=0.5)]]
}
plotter.plot_objects.put(IterationPlots(work_dir, f"plot_shielding.pdf", config))

plotter.finish()
testcase.close_datafile()
