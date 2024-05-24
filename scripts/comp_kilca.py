from mephit_plot import ComplexPlot, HLine, Id, IterationPlots, Kilca, LogY, Mephit, ParallelPlotter, Plot1D, \
    PolmodePlots, run_dir, set_matplotlib_defaults, XTicks, YTicks
from numpy import abs, angle, arange, arctan2, argmax, full, ix_, nan, pi, sign, sum
from scipy.constants import c as clight
from functools import partial


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/TCFP'
    testcase = Mephit(work_dir)
    testcase.open_datafile()
    testcase.postprocess()
    sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])

    # MEPHIT validation
    kfs = [testcase.data['/mesh/res_ind'][0], testcase.data['/mesh/res_ind'][0] // 2]
    res_nonres = ['res', 'nonres']
    iters = ['initial', 'final']
    q = []
    q_half = []
    theta = []
    theta_half = []
    grad_pn = dict(zip(iters, [[], []]))
    grad_pn_half = dict(zip(iters, [[], []]))
    lorentz = dict(zip(iters, [[], []]))
    lhs = dict(zip(iters, [[], []]))
    rhs = dict(zip(iters, [[], []]))
    for kf in kfs:
        q.append(f"{testcase.data['/cache/fs/q'][kf]:.3f}")
        q_half.append(f"{testcase.data['/cache/fs_half/q'][kf - 1]:.3f}")
        # poloidal edge index = point index - 1
        ke_min = testcase.data['/mesh/kp_low'][kf - 1]
        ke_max = testcase.data['/mesh/kp_low'][kf - 1] + testcase.data['/mesh/kp_max'][kf - 1] - 1
        kt_min = testcase.data['/mesh/kt_low'][kf - 1] + 1
        kt_max = testcase.data['/mesh/kt_low'][kf - 1] + testcase.data['/mesh/kt_max'][kf - 1]
        # rad to deg
        theta.append(testcase.data['/mesh/node_theta_geom'][ke_min:ke_max] * 180.0 / pi)  # point index
        theta_half.append(arctan2(testcase.data['/mesh/cntr_Z'][kt_min - 1:kt_max - 1] -
                                  testcase.data['/mesh/Z_O'][()],
                                  testcase.data['/mesh/cntr_R'][kt_min - 1:kt_max - 1] -
                                  testcase.data['/mesh/R_O'][()]) * 180 / pi)
        theta_half[-1][theta_half[-1] < 0.0] += 360.0
        # dyn cm^-3 Mx^-1 to Pa Wb^-1
        dp0_dpsi = testcase.data['/cache/fs/dp_dpsi'][kf - 1] * 1.0e+7
        # G to T
        Bmod = testcase.data['/cache/mid_fields/B0'][ke_min - 1:ke_max - 1] * 1.0e-4
        h = full((3, ke_max - ke_min), nan)
        h[0, :] = testcase.data['/cache/mid_fields/B0_R'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        h[1, :] = testcase.data['/cache/mid_fields/B0_phi'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        h[2, :] = testcase.data['/cache/mid_fields/B0_Z'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        Bn_psi = {}
        for it in iters:
            # dyn cm^-3 to N m^-3
            grad_pn[it].append(testcase.data[f"/debug_MDE_{it}/grad_pn"][ke_min - 1:ke_max - 1, :].T * 10)
            grad_pn_half[it].append(testcase.data[f"/debug_currn_{it}/grad_pn"][kt_min - 1:kt_max - 1, :].T * 10)
            lorentz[it].append(testcase.data[f"/debug_currn_{it}/lorentz"][kt_min - 1:kt_max - 1, :].T * 10)
            # G^2 cm to T^2 m
            Bn_psi[it] = testcase.data[f"/debug_MDE_{it}/Bn_psi_contravar"][ke_min - 1:ke_max - 1] * 1.0e-10
            # LHS vs. RHS
            lhs[it].append(sum(h * grad_pn[it][-1], axis=0))
            rhs[it].append(-Bn_psi[it] / Bmod * dp0_dpsi)
    theta_ticks = XTicks(arange(0, 360 + 1, 45))
    config = {
        'xlabel': r'$\theta$ / \si{\degree}',
        'legend': {'fontsize': 'x-small', 'loc': 'lower right'},
    }
    for it in iters:
        config['postprocess'] = [[theta_ticks, theta_ticks], [HLine(0.0, color='k', alpha=0.5, lw=0.5), Id()]]
        # MDE p_n
        config['ylabel'] = [r'abs MDE $p_{n}$ / \si{\newton\per\cubic\meter}', r'arg MDE $p_{n}$ / \si{\degree}']
        config['title'] = f"MDE for pressure perturbation, {it} iteration"
        config['plotdata'] = []
        for k in range(len(kfs)):
            config['plotdata'].append({'x': theta[k], 'y': lhs[it][k], 'args': {'lw': 0.5},
                                       'label': '$q = ' + q[k] + r', B_{0}^{-1} \V{B}_{0} \cdot \opgrad p_{n}$'})
            config['plotdata'].append({'x': theta[k], 'y': lhs[it][k], 'args': {'lw': 0.5},
                                       'label': '$q = ' + q[k] + r', -B_{0}^{-1} B_{n}^{\psi} \partial_{\psi} p_{0}$'})
        plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_MDE_pn_{it}.pdf", config))
        # grad p_n
        config['postprocess'] = [[theta_ticks, theta_ticks], [LogY(), Id()]]
        config['ylabel'] = [r'$\abs \opgrad p_{n}$ / \si{\newton\per\cubic\meter}', r'$\arg \opgrad p_{n}$ / \si{\degree}']
        config['title'] = f"Solution for pressure perturbation, {it} perturbation"
        config['plotdata'] = []
        comps = [r'\partial_{R} p_{n}', r'\tfrac{\im n}{R} p_{n}', r'\partial_{Z} p_{n}']
        for k in range(len(kfs)):
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': theta[k], 'y': grad_pn[it][k][k_comp, :], 'args': {'lw': 0.5},
                                           'label': f"$q = {q[k]}, {comps[k_comp]}$"})
        plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_grad_pn_{it}.pdf", config))
        # iMHD
        config['ylabel'] = [r'$\abs f$ / \si{\newton\per\cubic\meter}', r'$\arg f$ / \si{\degree}']
        comps = [r'$R$ comp.', r'$\phi$ comp.', r'$Z$ comp.']
        for k in range(len(kfs)):
            config['title'] = f"Linearized iMHD force balance ($q = {q_half[k]}$), {it} iteration"
            config['plotdata'] = []
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': theta_half[k], 'y': lorentz[it][k][k_comp, :], 'args': {'lw': 0.5},
                                           'label': f"Lorentz force, {comps[k_comp]}"})
                config['plotdata'].append({'x': theta_half[k], 'y': grad_pn_half[it][k][k_comp, :],
                                           'label': f"pressure gradient, {comps[k_comp]}", 'args': {'lw': 0.5}})
            plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_MHD_{res_nonres[k]}_{it}.pdf", config))

    # KiLCA comparison
    mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/postprocess/Bmn/coeff_rad', rad=True)
    mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/postprocess/Bmn_vac/coeff_rad', rad=True)
    mephit_Ires = testcase.get_Ires()
    reference = Kilca(work_dir)
    reference.open_datafile('TCFP_flre_hip.hdf5')
    kilca_Bmn = reference.get_polmodes('full perturbation (KiLCA)')
    kilca_Ires = reference.get_Ires()
    reference.open_datafile('TCFP_vac_hip.hdf5')
    kilca_Bmn_vac = reference.get_polmodes('vacuum perturbation (KiLCA)')
    reference.close_datafile()
    config = {
        'xlabel': '$m$', 'ylabel': r'$\abs\, I_{m, n}^{\parallel}$ / \si{\ampere}', 'legend': {'fontsize': 'small'},
        'plotdata': [
            {'x': mephit_Ires.keys(), 'y': mephit_Ires.values(), 'args': {'label': 'MEPHIT', 'marker': 'o', 'ls': ''}},
            {'x': kilca_Ires.keys(), 'y': kilca_Ires.values(), 'args': {'label': 'KiLCA', 'marker': 'x', 'ls': ''}}
        ],
        'postprocess': [XTicks(mephit_Ires.keys()), LogY()]
    }
    plotter.plot_objects.put(Plot1D(work_dir, 'plot_Ires.pdf', config))
    config = {
        'xlabel': r'$r$ / \si{\centi\meter}',
        'ylabel': r'$\abs\, B_{m, n}^{r}$ / \si{\tesla}',
        'rho': testcase.data['/cache/fs/rad'][()], 'q': testcase.data['/cache/fs/q'][()],
        'resonances': testcase.post['rad_res'], 'sgn_m_res': testcase.post['sgn_m_res'], 'omit_res': False,
        'poldata': [mephit_Bmn_vac, kilca_Bmn_vac, mephit_Bmn, kilca_Bmn], 'comp': abs,
        'res_neighbourhood': testcase.post['rad_res_neighbourhood']
    }
    plotter.plot_objects.put(PolmodePlots(work_dir, 'KiLCA_Bmn_r_abs.pdf', config))
    phase_ticks = YTicks(arange(-180, 180 + 1, 45))
    config['postprocess'] = [phase_ticks]
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg\, B_{m, n}^{r}$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'KiLCA_Bmn_r_arg.pdf', config))

    # iteration comparisons
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
        jmnpar_Bmod = testcase.get_polmodes(None, f"/iter/jmnpar_Bmod_{kiter:03}/coeff", 4.0 * pi/ clight, L1=True)
        Bmn_rad = testcase.get_polmodes(None, f"/iter/Bmn_{kiter:03}/coeff_rad")
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
        'xlabel': r'$r$ / \si{\centi\meter}',
        'ylabel': [
            r'$\abs\, p_{mn}$ / \si{\pascal}',
            r'$\abs\, \bigl[ \mu_{0} J_{n}^{\parallel} B_{0}^{-1} \bigr]_{m}$ / \si{\per\meter}',
            r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}'
        ],
        'rho': [testcase.post['rad'], testcase.post['rad'], testcase.post['rad_half']],
        'yscale': ['log', 'log', 'linear'], 'global_ylims': False, 'plotargs': {'lw': 0.5}, 'niter': niter,
        'postprocess': [[Id(), Id(), HLine(0.0, color='k', alpha=0.5, lw=0.5)]]
    }
    m_res_range = arange(m_res_min, m_res_max + 1)
    for m in m_res_range:
        m_res = m * testcase.post['sgn_m_res']
        config['title'] = f"Poloidal Modes in Preconditioned Iterations $k$ for $n = {n}$, $m = \pm {m}$"
        config['res_pos'] = testcase.post['rad_res'][m_res]
        config['res_neighbourhood'] = testcase.post['rad_res_neighbourhood'][m_res]
        config['zoom_x'] = config['res_neighbourhood'][ix_([0, -1])]
        config['plotdata'] = [
            abs_pmn_iter[m_res + m_res_max, :, :],
            abs_jmnpar_Bmod_iter[m_res + m_res_max, :, :],
            abs_Bmn_rad_iter[m_res + m_res_max, :, :],
            abs_pmn_iter[-m_res + m_res_max, :, :],
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
