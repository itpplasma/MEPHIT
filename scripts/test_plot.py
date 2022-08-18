from mephit_plot import Gpec, HLine, Id, LogY, Mephit, ParallelPlotter, Plot1D, PlotObject, PolmodePlots, run_dir, \
    set_matplotlib_defaults, XTicks, YTicks
from numpy import abs, angle, arange, arctan2, argmax, full, nan, pi, sign, sum
from functools import partial


class IterationPlots(PlotObject):
    def do_plot(self):
        from os import path
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.backends.backend_pdf import PdfPages
        from numpy import amin, amax
        super().do_plot()
        horz_plot = 3
        vert_plot = 1
        three_squares = (9.9, 3.3)
        pdf = PdfPages(path.join(self.work_dir, self.filename))
        if self.config['global_ylims']:
            ylims = [(amin(ydata), amax(ydata)) for ydata in self.config['plotdata']]
        else:
            ylims = [None, None, None]
        for kiter in range(0, self.config['niter']):
            print(f"Plotting {self.filename}, k = {kiter + 1} ...")
            fig = Figure(figsize=three_squares)
            axs = fig.subplots(vert_plot, horz_plot)
            for k in range(horz_plot):
                axs[k].set_yscale(self.config['yscale'][k])
                axs[k].plot(self.config['rho'][k], self.config['plotdata'][k][kiter, :], **self.config['plotargs'])
                if not self.config['global_ylims'] and kiter == 0:
                    ylims[k] = axs[k].get_ylim()
                else:
                    axs[k].set_ylim(ylims[k])
                if 'res_pos' in self.config.keys():
                    axs[k].axvline(self.config['res_pos'], color='b', alpha=0.5, lw=0.5)
                if 'postprocess' in self.config.keys():
                    for f in self.config['postprocess']:
                        f[k](fig, axs[k])
                axs[k].set_xlabel(self.config['xlabel'])
                axs[k].set_ylabel(self.config['ylabel'][k])
            fig.suptitle(self.config['title'] + f", $k = {kiter + 1}$")
            canvas = FigureCanvas(fig)
            fig.savefig(pdf, format='pdf', dpi=300)
        pdf.close()


class ComplexPlot(PlotObject):
    def do_plot(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from os import path
        from numpy import abs, angle
        super().do_plot()
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        fig = Figure(figsize=two_squares)
        axs = fig.subplots(vert_plot, horz_plot, sharex='all')
        labels = []
        for plotdata in self.config['plotdata']:
            axs[0].plot(plotdata['x'], abs(plotdata['y']), **plotdata['args'])
            axs[1].plot(plotdata['x'], angle(plotdata['y'], deg=True), **plotdata['args'])
            axs[1].set_yticks(arange(-180, 180+1, 45))
            axs[1].axhline(0.0, color='k', alpha=0.5, lw=0.5)
            labels.append(plotdata['label'])
        for k in range(horz_plot):
            if 'postprocess' in self.config.keys():
                for f in self.config['postprocess']:
                    f[k](fig, axs[k])
            if 'xlabel' in self.config.keys():
                axs[k].set_xlabel(self.config['xlabel'])
            if 'ylabel' in self.config.keys():
                axs[k].set_ylabel(self.config['ylabel'][k])
        if 'legend' in self.config.keys():
            fig.legend(labels=labels, **self.config['legend'])
        if 'title' in self.config.keys():
            fig.suptitle(self.config['title'])
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(self.work_dir, self.filename))


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
    conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][:, 16]
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

    # iteration comparisons
    n = testcase.data['/config/n'][()]
    m_res_min = testcase.data['/mesh/m_res_min'][()]
    m_res_max = testcase.data['/mesh/m_res_max'][()]
    niter = testcase.data['/iter/niter'][()]
    nflux = testcase.data['/mesh/nflux'][()]
    abs_pmn_iter = full((m_res_max - m_res_min + 1, niter, nflux + 1), nan)
    abs_jmnpar_Bmod_iter = full((m_res_max - m_res_min + 1, niter, nflux + 1), nan)
    abs_Bmn_rad_iter = full((m_res_max - m_res_min + 1, niter, nflux), nan)
    kf_res = testcase.data['/mesh/res_ind'][()]
    shielding = full((niter, m_res_max - m_res_min + 1), nan)
    penetration = full((niter, m_res_max - m_res_min + 1), nan)
    Ires = full((niter, m_res_max - m_res_min + 1), nan)
    for kiter in range(0, niter):
        pmn = testcase.get_polmodes(None, f"/iter/pmn_{kiter:03}/coeff", 0.1, L1=True)
        jmnpar_Bmod = testcase.get_polmodes(None, f"/iter/jmnpar_Bmod_{kiter:03}/coeff", 1.0 / 29.9792458, L1=True)
        Bmn_rad = testcase.get_polmodes(None, f"/iter/Bmn_{kiter:03}/coeff_rad", conversion)
        Ires[kiter, :] = abs(testcase.data[f"/iter/Ires_{kiter:03}"][()]) * 0.1 / 2.99792458e+08
        for m in range(m_res_min, m_res_max + 1):
            m_res = m * testcase.post['sgn_m_res']
            abs_pmn_iter[m - m_res_min, kiter, :] = abs(pmn['var'][m_res])
            abs_jmnpar_Bmod_iter[m - m_res_min, kiter, :] = abs(jmnpar_Bmod['var'][m_res])
            abs_Bmn_rad_iter[m - m_res_min, kiter, :] = abs(Bmn_rad['var'][m_res])
            kf_max = argmax(abs_Bmn_rad_iter[m - m_res_min, kiter, :kf_res[m - m_res_min] - 1]) + 1
            penetration[kiter, m - m_res_min] = abs_Bmn_rad_iter[m - m_res_min, kiter, kf_max - 1] / \
                abs(mephit_Bmn_vac['var'][m][kf_max - 1])
            shielding[kiter, m - m_res_min] = abs_Bmn_rad_iter[m - m_res_min, kiter, kf_res[m - m_res_min] - 1] / \
                abs(mephit_Bmn_vac['var'][m][kf_res[m - m_res_min] - 1])
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
    for m in range(m_res_min, m_res_max + 1):
        m_res = m * testcase.post['sgn_m_res']
        config['title'] = f"Poloidal Modes in Preconditioned Iterations $k$ for $n = {n}$, $m = {m_res}$"
        config['res_pos'] = testcase.post['psi_norm_res'][m]
        config['plotdata'] = [
            abs_pmn_iter[m - m_res_min, :, :],
            abs_jmnpar_Bmod_iter[m - m_res_min, :, :],
            abs_Bmn_rad_iter[m - m_res_min, :, :]
        ]
        plotter.plot_objects.put(IterationPlots(work_dir, f"plot_iter_{m}.pdf", config))
    m_res_range = arange(m_res_min, m_res_max + 1)
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
