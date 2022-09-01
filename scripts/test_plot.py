from mephit_plot import ComplexPlot, Gpec, HLine, Id, IterationPlots, LogY, Mephit, ParallelPlotter, Plot1D, \
    PolmodePlots, run_dir, set_matplotlib_defaults, XTicks, YTicks
from numpy import abs, angle, arange, arctan2, argmax, array, full, gradient, ix_, nan, nonzero, pi, sign, sum
from numpy.fft import fft, rfft
from scipy.interpolate import UnivariateSpline
from functools import partial


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/33353_2325'
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
                                       'label': '$q = ' + q[k] + r', B_{0}^{-1} \V{B}_{0} \cdot \grad p_{n}$'})
            config['plotdata'].append({'x': theta[k], 'y': lhs[it][k], 'args': {'lw': 0.5},
                                       'label': '$q = ' + q[k] + r', -B_{0}^{-1} B_{n}^{\psi} \partial_{\psi} p_{0}$'})
        plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_MDE_pn_{it}.pdf", config))
        # grad p_n
        config['postprocess'] = [[theta_ticks, theta_ticks], [LogY(), Id()]]
        config['ylabel'] = [r'$\abs \grad p_{n}$ / \si{\newton\per\cubic\meter}', r'$\arg \grad p_{n}$ / \si{\degree}']
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

    # GPEC comparison
    conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][()]
    mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/postprocess/Bmn/coeff_rad', conversion)
    mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/postprocess/Bmn_vac/coeff_rad', conversion)
    reference = Gpec(work_dir, 2)
    reference.open_datafiles()
    gpec_Bmn = reference.get_polmodes('full perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi')
    gpec_Bmn_vac = reference.get_polmodes('vacuum perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi_x')
    config = {
        'xlabel': r'$\hat{\psi}$',
        'ylabel': r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}',
        'rho': testcase.post['psi_norm'], 'q': testcase.data['/cache/fs/q'][()],
        'resonances': testcase.post['psi_norm_res'], 'sgn_m_res': testcase.post['sgn_m_res'], 'omit_res': False,
        'poldata': [mephit_Bmn_vac, gpec_Bmn_vac, mephit_Bmn, gpec_Bmn], 'comp': abs,
        'res_neighbourhood': testcase.post['psi_norm_res_neighbourhood']
    }
    plotter.plot_objects.put(PolmodePlots(work_dir, 'GPEC_Bmn_psi_abs.pdf', config))
    phase_ticks = YTicks(arange(-180, 180 + 1, 45))
    config['postprocess'] = [phase_ticks]
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'GPEC_Bmn_psi_arg.pdf', config))
    ## zoomed plots
    zoom_x = [0.001, 0.1]
    zoom_m_max = 5
    mephit_dpsi_Bmn = {'m_max': zoom_m_max, 'label': mephit_Bmn['label'], 'rho': dict(), 'var': dict()}
    mephit_dpsi_Bmn_vac = {'m_max': zoom_m_max, 'label': mephit_Bmn_vac['label'], 'rho': dict(), 'var': dict()}
    gpec_dpsi_Bmn = {'m_max': zoom_m_max, 'label': gpec_Bmn['label'], 'rho': dict(), 'var': dict()}
    gpec_dpsi_Bmn_vac = {'m_max': zoom_m_max, 'label': gpec_Bmn_vac['label'], 'rho': dict(), 'var': dict()}
    for m in mephit_Bmn['rho']:
        if abs(m) <= zoom_m_max:
            mephit_dpsi_Bmn['rho'][m] = mephit_Bmn['rho'][m]
            mephit_dpsi_Bmn['var'][m] = full(mephit_Bmn['var'][m].shape, nan, dtype='D')
            mephit_dpsi_Bmn['var'][m].real = gradient(mephit_Bmn['var'][m].real, mephit_Bmn['rho'][m])
            mephit_dpsi_Bmn['var'][m].imag = gradient(mephit_Bmn['var'][m].imag, mephit_Bmn['rho'][m])
    for m in mephit_Bmn_vac['rho']:
        if abs(m) <= zoom_m_max:
            mephit_dpsi_Bmn_vac['rho'][m] = mephit_Bmn_vac['rho'][m]
            mephit_dpsi_Bmn_vac['var'][m] = full(mephit_Bmn_vac['var'][m].shape, nan, dtype='D')
            mephit_dpsi_Bmn_vac['var'][m].real = gradient(mephit_Bmn_vac['var'][m].real, mephit_Bmn_vac['rho'][m])
            mephit_dpsi_Bmn_vac['var'][m].imag = gradient(mephit_Bmn_vac['var'][m].imag, mephit_Bmn_vac['rho'][m])
    for m in gpec_Bmn['rho']:
        if abs(m) <= zoom_m_max:
            gpec_dpsi_Bmn['rho'][m] = gpec_Bmn['rho'][m]
            gpec_dpsi_Bmn['var'][m] = full(gpec_Bmn['var'][m].shape, nan, dtype='D')
            gpec_dpsi_Bmn['var'][m].real = gradient(gpec_Bmn['var'][m].real, gpec_Bmn['rho'][m])
            gpec_dpsi_Bmn['var'][m].imag = gradient(gpec_Bmn['var'][m].imag, gpec_Bmn['rho'][m])
    for m in gpec_Bmn_vac['rho']:
        if abs(m) <= zoom_m_max:
            mask = nonzero((zoom_x[0] <= gpec_Bmn_vac['rho'][m]) & (gpec_Bmn_vac['rho'][m] <= zoom_x[1]))
            gpec_dpsi_Bmn_vac['rho'][m] = gpec_Bmn_vac['rho'][m]
            gpec_dpsi_Bmn_vac['var'][m] = full(gpec_Bmn_vac['var'][m].shape, nan, dtype='D')
            gpec_dpsi_Bmn_vac['var'][m].real = gradient(gpec_Bmn_vac['var'][m].real, gpec_Bmn_vac['rho'][m])
            gpec_dpsi_Bmn_vac['var'][m].imag = gradient(gpec_Bmn_vac['var'][m].imag, gpec_Bmn_vac['rho'][m])
    mephit_Bmn['m_max'] = zoom_m_max
    mephit_Bmn_vac['m_max'] = zoom_m_max
    gpec_Bmn['m_max'] = zoom_m_max
    gpec_Bmn_vac['m_max'] = zoom_m_max
    config = {
        'xlabel': r'$\psi$',
        'ylabel': r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}',
        'rho': testcase.post['psi_norm'], 'q': testcase.data['/cache/fs/q'][()],
        'resonances': testcase.post['psi_norm_res'], 'sgn_m_res': testcase.post['sgn_m_res'], 'omit_res': True,
        'poldata': [mephit_Bmn_vac, gpec_Bmn_vac, mephit_Bmn, gpec_Bmn], 'comp': abs, 'zoom_x': zoom_x
    }
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_Jbgradpsi.pdf', config))
    config['ylabel'] = r'$\abs\, \partial_{\psi} [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}'
    config['poldata'] = [mephit_dpsi_Bmn_vac, gpec_dpsi_Bmn_vac, mephit_dpsi_Bmn, gpec_dpsi_Bmn]
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_dpsi_Jbgradpsi.pdf', config))
    ## normalization debugging
    cm_to_m = 1.0e-02
    G_to_T = 1.0e-04
    Mx_to_Wb = 1.0e-08
    psi = sgn_dpsi * Mx_to_Wb * testcase.data['/cache/fs/psi'][()]
    psi -= psi[0]
    psipol_max = psi[-1]
    psi /= psipol_max
    jac_psi = Mx_to_Wb / psipol_max * testcase.data['/debug_GPEC/jac/psi'][()]
    jac_nrad = jac_psi.size
    jac_theta = 0.5 / pi * testcase.data['/debug_GPEC/jac/theta'][()]
    jac_npol = jac_theta.size
    jac_gpec = cm_to_m / G_to_T * testcase.data['/debug_GPEC/jac/jac'][()]
    jac_mephit = cm_to_m / G_to_T * testcase.data['/debug_GPEC/jac/sqrt_g'][()]
    debug_psi = Mx_to_Wb / psipol_max * testcase.data['/debug_GPEC/psi'][()]
    debug_nrad = debug_psi.size
    debug_theta = 0.5 / pi * testcase.data['/debug_GPEC/theta'][()]
    debug_npol = debug_theta.size
    debug_R = cm_to_m * testcase.data['/debug_GPEC/R_GPEC'][()]
    R = cm_to_m * testcase.data['/debug_GPEC/R'][()]
    spline_F = UnivariateSpline(jac_psi, array(reference.data['dcon']['f']), s=0)
    spline_q = UnivariateSpline(jac_psi, array(reference.data['dcon']['q']), s=0)
    jac_pest = full(debug_R.shape, nan)
    for k in range(jac_pest.shape[0]):
        jac_pest[k, :] = debug_R[k, :] ** 2 * abs(spline_q(debug_psi[k]) / spline_F(debug_psi[k]))
    delpsi = G_to_T * cm_to_m * testcase.data['/debug_GPEC/delpsi'][()]
    grad_psi = G_to_T * cm_to_m * testcase.data['/debug_GPEC/grad_psi'][()]
    jacfac = testcase.data['/debug_GPEC/jacfac'][()]
    contradenspsi = testcase.data['/debug_GPEC/contradenspsi'][()]
    jac_modes_gpec = rfft(jac_gpec[:, :-1]) / (jac_npol - 1)
    jac_modes_mephit = rfft(jac_mephit[:, :-1]) / (jac_npol - 1)
    modes_delpsi = rfft(delpsi[:, :-1]) / (jac_npol - 1)
    modes_grad_psi = rfft(grad_psi[:, :-1]) / (jac_npol - 1)
    modes_jacfac = fft(jacfac[:, :-1]) / (jac_npol - 1)
    modes_contradenspsi = rfft(contradenspsi[:, :-1]) / (jac_npol - 1)
    psi_n = r'$\hat{\psi}$'
    theta_n = r'$\hat{\vartheta}$'
    sqrt_g = r'$\sqrt{g}$ / \si{\meter\per\tesla}'
    config = {'xlabel': psi_n, 'ylabel': sqrt_g, 'legend': {'fontsize': 'small'}}
    config['plotdata'] = [
        {'x': jac_psi, 'y': jac_gpec[:, 0], 'args': {'label': r'GPEC, $\hat{\vartheta} = 0$'}},
        {'x': jac_psi, 'y': jac_mephit[:, 0], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0$'}},
        {'x': jac_psi, 'y': jac_gpec[:, jac_npol // 2], 'args': {'label': r'GPEC, $\hat{\vartheta} = 0.5$'}},
        {'x': jac_psi, 'y': jac_mephit[:, jac_npol // 2], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0.5$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'comp_jac_psi.pdf', config))
    config['plotdata'] = [
        {'x': debug_psi, 'y': jac_pest[:, 0], 'args': {'label': r'PEST, $\hat{\vartheta} = 0$'}},
        {'x': jac_psi, 'y': jac_mephit[:, 0], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0$'}},
        {'x': debug_psi, 'y': jac_pest[:, debug_npol // 2], 'args': {'label': r'PEST, $\hat{\vartheta} = 0.5$'}},
        {'x': jac_psi, 'y': jac_mephit[:, jac_npol // 2], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0.5$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'debug_jac_psi.pdf', config))
    config['xlabel'] = theta_n
    config['plotdata'] = [
        {'x': jac_theta, 'y': jac_gpec[jac_nrad // 2, :], 'args': {'label': r'GPEC, $\hat{\psi} = 0.5$'}},
        {'x': jac_theta, 'y': jac_mephit[jac_nrad // 2, :], 'args': {'label': r'MEPHIT, $\hat{\psi} = 0.5$'}},
        {'x': jac_theta, 'y': jac_gpec[jac_nrad // 4 * 3, :], 'args': {'label': r'GPEC, $\hat{\psi} = 0.75$'}},
        {'x': jac_theta, 'y': jac_mephit[jac_nrad // 4 * 3, :], 'args': {'label': r'MEPHIT, $\hat{\psi} = 0.75$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'comp_jac_theta.pdf', config))
    config['xlabel'] = psi_n
    config['ylabel'] = 'ratio of Jacobians (GPEC / MEPHIT)'
    config['plotdata'] = [
        {'x': jac_psi, 'y': jac_gpec[:, 0] / jac_mephit[:, 0],
         'args': {'label': r'$\hat{\vartheta} = 0$'}},
        {'x': jac_psi, 'y': jac_gpec[:, jac_npol // 2] / jac_mephit[:, jac_npol // 2],
         'args': {'label': r'$\hat{\vartheta} = 0.5$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'comp_jac_psi_fac.pdf', config))
    config['xlabel'] = theta_n
    config['ylabel'] = 'ratio of Jacobians (GPEC / MEPHIT)'
    config['plotdata'] = [
        {'x': jac_theta, 'y': jac_gpec[jac_nrad // 2, :] / jac_mephit[jac_nrad // 2, :],
         'args': {'label': r'$\hat{\psi} = 0.5$'}},
        {'x': jac_theta, 'y': jac_gpec[jac_nrad // 4 * 3, :] / jac_mephit[jac_nrad // 4 * 3, :],
         'args': {'label': r'$\hat{\psi} = 0.75$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'comp_jac_theta_fac.pdf', config))
    config['xlabel'] = psi_n
    config['ylabel'] = r'$R$ / \si{\meter}'
    config['plotdata'] = [
        {'x': debug_psi, 'y': debug_R[:, debug_npol // 4], 'args': {'label': r'GPEC, $\hat{\vartheta} = 0.25$'}},
        {'x': debug_psi, 'y': R[:, debug_npol // 4], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0.25$'}},
        {'x': debug_psi, 'y': debug_R[:, debug_npol // 4 * 3], 'args': {'label': r'GPEC, $\hat{\vartheta} = 0.75$'}},
        {'x': debug_psi, 'y': R[:, debug_npol // 4 * 3], 'args': {'label': r'MEPHIT, $\hat{\vartheta} = 0.75$'}}
    ]
    plotter.plot_objects.put(Plot1D(work_dir, 'comp_R.pdf', config))
    config['xlabel'] = psi_n
    config['ylabel'] = sqrt_g
    for m in range(0, 9):
        config['plotdata'] = [
            {'x': jac_psi, 'y': jac_modes_gpec[:, m].real,
             'args': {'label': fr"$\Real \sqrt{{g}}_{{m = {m}}}$ (GPEC)"}},
            {'x': jac_psi, 'y': jac_modes_mephit[:, m].real,
             'args': {'label': fr"$\Real \sqrt{{g}}_{{m = {m}}}$ (MEPHIT)"}},
            {'x': jac_psi, 'y': jac_modes_gpec[:, m].imag,
             'args': {'label': fr"$\Imag \sqrt{{g}}_{{m = {m}}}$ (GPEC)"}},
            {'x': jac_psi, 'y': jac_modes_mephit[:, m].imag,
             'args': {'label': fr"$\Imag \sqrt{{g}}_{{m = {m}}}$ (MEPHIT)"}}
        ]
        plotter.plot_objects.put(Plot1D(work_dir, f"comp_jac_{m}.pdf", config))
    config['xlabel'] = psi_n
    config['ylabel'] = r'$\lVert \nabla \psi \rVert$ / \si{\tesla\meter}'
    for m in range(0, 9):
        config['plotdata'] = [
            {'x': debug_psi[::10], 'y': modes_delpsi[:, m].real,
             'args': {'label': fr"$\Real \lVert \nabla \psi \rVert_{{m = {m}}}$ (GPEC)"}},
            {'x': debug_psi[::10], 'y': modes_grad_psi[:, m].real,
             'args': {'label': fr"$\Real \lVert \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)"}},
            {'x': debug_psi[::10], 'y': modes_delpsi[:, m].imag,
             'args': {'label': fr"$\Imag \lVert \nabla \psi \rVert_{{m = {m}}}$ (GPEC)"}},
            {'x': debug_psi[::10], 'y': modes_grad_psi[:, m].imag,
             'args': {'label': fr"$\Imag \lVert \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)"}}
        ]
        plotter.plot_objects.put(Plot1D(work_dir, f"comp_grad_psi_{m}.pdf", config))
    config['xlabel'] = psi_n
    config['ylabel'] = r'$\sqrt{g} \lVert \nabla \psi \rVert / \langle \sqrt{g} \lVert \nabla \psi \rVert \rangle$'
    for m in range(0, 9):
        config['plotdata'] = [
            {'x': debug_psi[::10], 'y': abs(modes_jacfac[:, m].real),
             'args': {'label': fr"$\Real \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (GPEC)"}},
            {'x': debug_psi[::10], 'y': abs(modes_contradenspsi[:, m].real),
             'args': {'label': fr"$\Real \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)"}},
            {'x': debug_psi[::10], 'y': abs(modes_jacfac[:, m].imag),
             'args': {'label': fr"$\Imag \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (GPEC)"}},
            {'x': debug_psi[::10], 'y': abs(modes_contradenspsi[:, m].imag),
             'args': {'label': fr"$\Imag \lVert \sqrt{{g}} \nabla \psi \rVert_{{m = {m}}}$ (MEPHIT)"}}
        ]
        plotter.plot_objects.put(Plot1D(work_dir, f"comp_jacfac_{m}.pdf", config))

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
        jmnpar_Bmod = testcase.get_polmodes(None, f"/iter/jmnpar_Bmod_{kiter:03}/coeff", 1.0 / 29.9792458, L1=True)
        Bmn_rad = testcase.get_polmodes(None, f"/iter/Bmn_{kiter:03}/coeff_rad", conversion)
        Ires[kiter, :] = abs(testcase.data[f"/iter/Ires_{kiter:03}"][()]) * 0.1 / 2.99792458e+08
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
            r'$\abs\, \bigl[ J_{n}^{\parallel} B_{0}^{-1} \bigr]_{m}$ / \si{\per\henry}',
            r'$\abs\, [\sqrt{g} \V{B}_{n} \cdot \nabla \psi]_{m} A^{-1}$ / \si{\tesla}'
        ],
        'rho': [testcase.post['psi_norm'], testcase.post['psi_norm'], testcase.post['psi_half_norm']],
        'yscale': ['log', 'log', 'linear'], 'global_ylims': False, 'plotargs': {'lw': 0.5}, 'niter': niter,
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
            abs_pmn_iter[-m_res + m_res_max, :, :],
            abs_jmnpar_Bmod_iter[-m_res + m_res_max, :, :],
            abs_Bmn_rad_iter[-m_res + m_res_max, :, :]
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

    # convergence estimation
    niter = testcase.data['/iter/niter'][()]
    sup_eigval = testcase.data['/config/ritz_threshold'][()]
    L2int_Bnvac = testcase.data['/iter/L2int_Bnvac'][()]
    L2int_Bn_diff = testcase.data['/iter/L2int_Bn_diff'][:niter]
    rel_err = testcase.data['/config/iter_rel_err'][()]
    plotter.plot_objects.put(Plot1D.conv_plot(work_dir, sup_eigval, L2int_Bnvac, L2int_Bn_diff, rel_err))

    plotter.finish()
    testcase.close_datafile()
    reference.close_datafiles()
