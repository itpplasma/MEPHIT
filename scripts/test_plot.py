from mephit_plot import ComplexPlot, Gpec, HLine, Id, IterationPlots, LogY, Mephit, ParallelPlotter, Plot1D, \
    PolmodePlots, run_dir, set_matplotlib_defaults, XTicks, YTicks
from numpy import abs, angle, arange, arctan2, argmax, array, full, gradient, ix_, nan, ndarray, nonzero, pi, r_, sign, sum, where
from numpy.fft import fft, rfft
from numpy.linalg import norm
from scipy.constants import c as clight
from scipy.interpolate import UnivariateSpline
from functools import partial


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/33353_2900'
    testcase = Mephit(work_dir)
    testcase.open_datafile()
    testcase.postprocess()
    sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])
    nflux = testcase.data['/mesh/nflux'][()]
    m_max = 24
    m_res_min = testcase.data['/mesh/m_res_min'][()]

    # MEPHIT validation
    iters = ['initial', 'final']
    samples = dict((key, []) for key in ['kf', 'res', 'q', 'theta'])
    for k in range(3):
        samples['kf'].append((testcase.data['/mesh/res_ind'][k] +
                              (testcase.data['/mesh/res_ind'][k - 1] if k > 0 else 0)) // 2 - 1)
        samples['kf'].append(testcase.data['/mesh/res_ind'][k] - 1)
        samples['res'] += ['nonres', 'res']
    valid = {}
    for qty in ['pn', 'pn_lhs', 'pn_rhs', 'jnpar_Bmod', 'jnpar_Bmod_lhs', 'jnpar_Bmod_rhs']:
        samples[qty] = dict((it, []) for it in iters)
        valid[qty] = dict((it, full((2 * m_max + 1, nflux), nan, dtype='D')) for it in iters)
    for qty in ['grad_pn', 'grad_jnpar_Bmod', 'lorentz']:
        samples[qty] = dict((it, []) for it in iters)
        valid[qty] = dict((it, full((2 * m_max + 1, nflux, 3), nan, dtype='D')) for it in iters)
    labels = {
        'pn': r'\delta p',
        'pn_lhs': r'\V{h} \cdot \nabla \delta p',
        'pn_rhs': r"-p_{0}'(\psi) \delta B^{\psi} B_{0}^{-1}",
        'jnpar_Bmod': r'\mu_{0} \delta J^{\parallel} B_{0}^{-1}',
        'jnpar_Bmod_lhs': r'\V{h} \cdot \nabla (\mu_{0} \delta J^{\parallel} B_{0}^{-1})',
        'jnpar_Bmod_rhs': r'-\mu_{0} B_{0}^{-1} \nabla \cdot \delta \V{J}^{\perp}',
        'lorentz': r'\delta (\V{J} \times \V{B})'
    }
    conversions = {
        'pn': 0.1, 'grad_pn': 10.0, 'Bn_psi_contravar': 1.0e-10, 'lorentz': 10.0,
        'jnpar_Bmod': 4.0 * pi / clight, 'grad_jnpar_Bmod': 400.0 * pi / clight, 'div_jnperp': 0.04 * pi / clight
    }
    units = {
        'pn': r'\pascal', 'pn_lhs': r'\newton\per\cubic\meter',
        'jnpar_Bmod': r'\per\meter', 'jnpar_Bmod_lhs': r'\per\square\meter',
        'lorentz': r'\newton\per\cubic\meter'
    }
    for kf in arange(nflux):
        k_low = testcase.data['/cache/kp_low'][kf] - 1
        k_max = 1 << testcase.data['/cache/log2_kp_max'][kf]
        trunc = r_[0:m_max+1, k_max-m_max:k_max]
        s = testcase.data['/cache/sample_polmodes']
        B0 = full((3, k_max), nan)
        B0[0, :] = s['B0_R'][k_low+1:k_low+k_max+1] * 1.0e-4
        B0[1, :] = s['B0_phi'][k_low+1:k_low+k_max+1] * 1.0e-4
        B0[2, :] = s['B0_Z'][k_low+1:k_low+k_max+1]  * 1.0e-4
        Bmod = norm(B0, axis=0)
        dp0_dpsi = testcase.data['/cache/fs/dp_dpsi'][kf+1] * 1.0e+7
        if kf in samples['kf']:
            samples['q'].append(f"{testcase.data['/cache/fs/q'][kf+1]:.3f}")
            samples['theta'].append(s['theta'][k_low+1:k_low+k_max+1] * 180.0 / pi)
        for it in iters:
            slices = {}
            for qty, factor in conversions.items():
                slices[qty] = testcase.data['/debug_MDE_' + it + '/' + qty][k_low+1:k_low+k_max+1, ...].T * factor
            for qty in ['pn', 'jnpar_Bmod']:
                valid[qty][it][:, kf] = fft(slices[qty])[trunc] / k_max
                valid['grad_' + qty][it][:, kf, :] = fft(slices['grad_' + qty])[:, trunc].T / k_max
                valid[qty + '_lhs'][it][:, kf] = fft(sum(B0 * slices['grad_' + qty], axis=0) / Bmod)[trunc] / k_max
            valid['lorentz'][it][:, kf, :] = fft(slices['lorentz'])[:, trunc].T / k_max
            valid['pn_rhs'][it][:, kf] = -dp0_dpsi * fft(slices['Bn_psi_contravar'] / Bmod)[trunc] / k_max
            valid['jnpar_Bmod_rhs'][it][:, kf] = -fft(slices['div_jnperp'] / Bmod)[trunc] / k_max
            if kf in samples['kf']:
                for qty in ['pn', 'grad_pn', 'jnpar_Bmod', 'grad_jnpar_Bmod', 'lorentz']:
                    samples[qty][it].append(slices[qty])
                samples['pn_lhs'][it].append(sum(B0 * slices['grad_pn'], axis=0) / Bmod)
                samples['pn_rhs'][it].append(-dp0_dpsi * slices['Bn_psi_contravar'] / Bmod)
                samples['jnpar_Bmod_lhs'][it].append(sum(B0 * slices['grad_jnpar_Bmod'], axis=0) / Bmod)
                samples['jnpar_Bmod_rhs'][it].append(-slices['div_jnperp'] / Bmod)
    for qty in ['pn', 'pn_lhs', 'pn_rhs', 'jnpar_Bmod', 'jnpar_Bmod_lhs', 'jnpar_Bmod_rhs']:
        label = '$' + labels[qty] + '$'
        for it in iters:
            valid[qty][it] = {'m_max': m_max, 'rho': {m: testcase.post['psi_norm'][1:] for m in range(-m_max, m_max + 1)},
                              'label': label, 'var': {m: valid[qty][it][m, :] for m in range(-m_max, m_max + 1)}}
    for vqty in ['grad_pn', 'grad_jnpar_Bmod', 'lorentz']:
        for ix, coord in enumerate(['R', 'phi', 'Z']):
            if 'grad_' in vqty:
                qty = vqty.replace('grad_', '', 1)
                label = r'$\partial_{' + coord.replace('phi', r'\varphi') + '} (' + labels[qty] + ')$'
                qty = 'd' + qty + '_d' + coord
            else:
                qty = vqty + "_" + coord
                label = '$' + labels[vqty] + r' \cdot \hat{\V{e}}_{' + coord.replace('phi', r'\varphi') + '}$'
            valid[qty] = {}
            for it in iters:
                valid[qty][it] = {'m_max': m_max, 'rho': {m: testcase.post['psi_norm'][1:] for m in range(-m_max, m_max + 1)},
                                  'label': label, 'var': {m: valid[vqty][it][m, :, ix] for m in range(-m_max, m_max + 1)}}
        del valid[vqty]
    arg = partial(angle, deg=True)
    phase_ticks = YTicks(arange(-180, 180 + 1, 45))
    theta_ticks = XTicks(arange(0, 360 + 1, 45))
    config = {
        'xlabel': r'$\hat{\psi}$', 'rho': testcase.post['psi_norm'][1:],
        'q': testcase.data['/cache/fs/q'][1:], 'sgn_m_res': testcase.post['sgn_m_res'],
        'omit_res': False, 'resonances': testcase.post['psi_norm_res'],
        'res_neighbourhood': testcase.post['psi_norm_res_neighbourhood']
    }
    for it in iters:
        for qty in ['pn', 'jnpar_Bmod']:
            config['poldata'] = [valid[qty + '_lhs'][it], valid[qty + '_rhs'][it]]
            config['ylabel'] = 'abs MDE $' + labels[qty] + r'$ / \si{' + units[qty + '_lhs'] + '}'
            config['comp'] = abs
            config['postprocess'] = [LogY()]
            plotter.plot_objects.put(PolmodePlots(work_dir, f"plot_{it}_polmodes_MDE_{qty}_abs.pdf", config))
            config['postprocess'] = [phase_ticks]
            config['comp'] = arg
            config['ylabel'] = 'arg MDE $' + labels[qty] + r'$ / \si{\degree}'
            plotter.plot_objects.put(PolmodePlots(work_dir, f"plot_{it}_polmodes_MDE_{qty}_arg.pdf", config))
        config['poldata'] = []
        for coord in ['R', 'phi', 'Z']:
                config['poldata'] += [valid['dpn_d' + coord][it], valid['lorentz_' + coord][it]]
        config['ylabel'] = r'abs $f$ / \si{' + units['lorentz'] + '}'
        config['comp'] = abs
        config['postprocess'] = [LogY()]
        plotter.plot_objects.put(PolmodePlots(work_dir, f"plot_{it}_polmodes_MHD_abs.pdf", config))
        config['postprocess'] = [phase_ticks]
        config['comp'] = arg
        config['ylabel'] = r'arg $f$ / \si{\degree}'
        plotter.plot_objects.put(PolmodePlots(work_dir, f"plot_{it}_polmodes_MHD_arg.pdf", config))

    config = {
        'xlabel': r'$\vartheta$ / \si{\degree}',
        'legend': {'fontsize': 'x-small', 'loc': 'lower right'}
    }
    for it in iters:
        config['postprocess'] = [[theta_ticks, theta_ticks], [HLine(0.0, color='k', alpha=0.5, lw=0.5), Id()]]
        # MDE p_n
        config['ylabel'] = [r'abs MDE $p_{n}$ / \si{\newton\per\cubic\meter}', r'arg MDE $p_{n}$ / \si{\degree}']
        config['title'] = f"MDE for pressure perturbation, {it} iteration"
        for k in range(3):
            config['plotdata'] = [{'x': samples['theta'][2*k], 'y': samples['pn_lhs'][it][2*k], 'args': {'lw': 0.5},
                                   'label': '$q = ' + samples['q'][2*k] + r', B_{0}^{-1} \V{B}_{0} \cdot \grad p_{n}$'},
                                  {'x': samples['theta'][2*k], 'y': samples['pn_rhs'][it][2*k], 'args': {'lw': 0.5},
                                   'label': '$q = ' + samples['q'][2*k] + r', -B_{0}^{-1} B_{n}^{\psi} \partial_{\psi} p_{0}$'},
                                  {'x': samples['theta'][2*k+1], 'y': samples['pn_lhs'][it][2*k+1], 'args': {'lw': 0.5},
                                   'label': '$q = ' + samples['q'][2*k+1] + r', B_{0}^{-1} \V{B}_{0} \cdot \grad p_{n}$'},
                                  {'x': samples['theta'][2*k+1], 'y': samples['pn_rhs'][it][2*k+1], 'args': {'lw': 0.5},
                                   'label': '$q = ' + samples['q'][2*k+1] + r', -B_{0}^{-1} B_{n}^{\psi} \partial_{\psi} p_{0}$'}]
            plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_{it}_MDE_pn_{m_res_min + k}.pdf", config))
        # grad p_n
        config['postprocess'] = [[theta_ticks, theta_ticks], [LogY(), Id()]]
        config['ylabel'] = [r'$\abs \grad p_{n}$ / \si{\newton\per\cubic\meter}', r'$\arg \grad p_{n}$ / \si{\degree}']
        config['title'] = f"Solution for pressure perturbation, {it} perturbation"
        comps = [r'\partial_{R} p_{n}', r'\tfrac{\im n}{R} p_{n}', r'\partial_{Z} p_{n}']
        for k in range(3):
            config['plotdata'] = []
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': samples['theta'][2*k], 'y': samples['grad_pn'][it][2*k][k_comp, :],
                                           'label': f"$q = {samples['q'][2*k]}, {comps[k_comp]}$", 'args': {'lw': 0.5}})
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': samples['theta'][2*k+1], 'y': samples['grad_pn'][it][2*k+1][k_comp, :],
                                           'label': f"$q = {samples['q'][2*k+1]}, {comps[k_comp]}$", 'args': {'lw': 0.5}})
            plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_{it}_grad_pn_{m_res_min + k}.pdf", config))
        # iMHD
        config['ylabel'] = [r'$\abs f$ / \si{\newton\per\cubic\meter}', r'$\arg f$ / \si{\degree}']
        comps = [r'$R$ comp.', r'$\phi$ comp.', r'$Z$ comp.']
        for k in range(2):
            config['title'] = f"Linearized iMHD force balance ($q = {samples['q'][k]}$), {it} iteration"
            config['plotdata'] = []
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': samples['theta'][k], 'y': samples['lorentz'][it][k][k_comp, :],
                                           'label': f"Lorentz force, {comps[k_comp]}", 'args': {'lw': 0.5}})
                config['plotdata'].append({'x': samples['theta'][k], 'y': samples['grad_pn'][it][k][k_comp, :],
                                           'label': f"pressure gradient, {comps[k_comp]}", 'args': {'lw': 0.5}})
            plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_{it}_MHD_{samples['res'][k]}.pdf", config))
    jmnpar_Bmod_excl = testcase.get_polmodes('excluding KiLCA current', '/debug_KiLCA/jmnpar_Bmod_excl/coeff', conversions['jnpar_Bmod'], L1=True)
    jmnpar_Bmod_incl = testcase.get_polmodes('including KiLCA current', '/debug_KiLCA/jmnpar_Bmod_incl/coeff', conversions['jnpar_Bmod'], L1=True)
    config = {
        'xlabel': r'$\hat{\psi}$',
        'ylabel': r'$\abs\, ' + labels['jnpar_Bmod'] + r'$ / \si{' + units['jnpar_Bmod'] + '}',
        'rho': testcase.post['psi_norm'], 'q': testcase.data['/cache/fs/q'][()],
        'resonances': testcase.post['psi_norm_res'], 'sgn_m_res': testcase.post['sgn_m_res'], 'omit_res': False,
        'poldata': [jmnpar_Bmod_excl, jmnpar_Bmod_incl], 'comp': abs,
        'res_neighbourhood': testcase.post['psi_norm_res_neighbourhood']
    }
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_debug_polmodes_jnpar_Bmod_abs.pdf', config))
    phase_ticks = YTicks(arange(-180, 180 + 1, 45))
    config['postprocess'] = [phase_ticks]
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg\, ' + labels['jnpar_Bmod'] + r'$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_debug_polmodes_jnpar_Bmod_arg.pdf', config))

    # GPEC comparison
    mephit_Ires = testcase.get_Ires()
    reference = Gpec(work_dir, 2)
    reference.open_datafiles()
    gpec_Ires = reference.get_Ires()
    config = {
        'xlabel': '$m$', 'ylabel': r'$\abs\, I_{m, n}^{\parallel}$ / \si{\ampere}', 'legend': {'fontsize': 'small'},
        'plotdata': [
            {'x': mephit_Ires.keys(), 'y': mephit_Ires.values(), 'args': {'label': 'MEPHIT', 'marker': 'o', 'ls': ''}},
            {'x': gpec_Ires.keys(), 'y': gpec_Ires.values(), 'args': {'label': 'GPEC', 'marker': 'x', 'ls': ''}}
        ],
        'postprocess': [XTicks(mephit_Ires.keys())]
    }
    plotter.plot_objects.put(Plot1D(work_dir, 'plot_Ires.pdf', config))

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
    reference.close_datafiles()
