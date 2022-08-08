from mephit_plot import Gpec, Id, LogY, Mephit, ParallelPlotter, Plot1D, PlotObject, PolmodePlots, run_dir, \
    set_matplotlib_defaults, XTicks
from numpy import abs, angle, arange, arctan2, empty, full, nan, pi, sign, sum
from functools import partial


class IterationPlots(PlotObject):
    def do_plot(self):
        from os import path
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.backends.backend_pdf import PdfPages
        super().do_plot()
        horz_plot = 3
        vert_plot = 1
        three_squares = (9.9, 3.3)
        pdf = PdfPages(path.join(self.work_dir, self.filename))
        ylims = [None, None, None]
        for kiter in range(0, self.config['niter']):
            print(f"Plotting {self.filename}, k = {kiter} ...")
            fig = Figure(figsize=three_squares)
            axs = fig.subplots(vert_plot, horz_plot)
            for k in range(horz_plot):
                axs[k].set_yscale(self.config['yscale'][k])
                axs[k].plot(self.config['rho'][k], self.config['plotdata'][k][kiter, :])
                if kiter == 0:
                    ylims[k] = axs[k].get_ylim()
                else:
                    axs[k].set_ylim(ylims[k])
                axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
                axs[k].axvline(self.config['res_pos'], color='b', alpha=0.5, lw=0.5)
                axs[k].set_xlabel(self.config['xlabel'])
                axs[k].set_ylabel(self.config['ylabel'][k])
            fig.suptitle(self.config['title'] + f", $k = {kiter}$")
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
            labels.append(plotdata['label'])
        for k in range(horz_plot):
            axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
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
    testcase.read_datafile()
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
    grad_pn = []
    grad_pn_half = []
    lorentz = []
    lhs = []
    rhs = []
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
        h = empty((3, ke_max - ke_min))
        h[0, :] = testcase.data['/cache/mid_fields/B0_R'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        h[1, :] = testcase.data['/cache/mid_fields/B0_phi'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        h[2, :] = testcase.data['/cache/mid_fields/B0_Z'][ke_min - 1:ke_max - 1] / Bmod * 1.0e-4
        # dyn cm^-3 to N m^-3
        grad_pn.append([testcase.data['/debug_MDE_000/grad_pn'][ke_min - 1:ke_max - 1, :].T * 10,
                        testcase.data['/debug_MDE/grad_pn'][ke_min - 1:ke_max - 1, :].T * 10])
        grad_pn_half.append([testcase.data['/debug_currn_000/grad_pn'][kt_min - 1:kt_max - 1, :].T * 10,
                             testcase.data['/debug_currn/grad_pn'][kt_min - 1:kt_max - 1, :].T * 10])
        lorentz.append([testcase.data['/debug_currn_000/lorentz'][kt_min - 1:kt_max - 1, :].T * 10,
                        testcase.data['/debug_currn/lorentz'][kt_min - 1:kt_max - 1, :].T * 10])
        # G^2 cm to T^2 m
        Bn_psi = [testcase.data['/debug_MDE_000/Bn_psi_contravar'][ke_min - 1:ke_max - 1] * 1.0e-10,
                  testcase.data['/debug_MDE/Bn_psi_contravar'][ke_min - 1:ke_max - 1] * 1.0e-10]
        # LHS vs. RHS
        lhs.append([])
        rhs.append([])
        for k_iter in range(len(iters)):
            lhs[-1].append(sum(h * grad_pn[-1][k_iter], axis=0))
            rhs[-1].append(-Bn_psi[k_iter] / Bmod * dp0_dpsi)
    theta_ticks = XTicks(arange(0, 360 + 1, 45))
    config = {
        'xlabel': r'$\theta$ / \si{\degree}',
        'legend': {'fontsize': 'x-small', 'loc': 'lower right'},
    }
    for k_iter in range(len(iters)):
        config['postprocess'] = [[theta_ticks, theta_ticks]]
        # MDE p_n
        config['ylabel'] = [r'abs MDE $p_{n}$ / \si{\newton\per\cubic\meter}', r'arg MDE $p_{n}$ / \si{\degree}']
        config['title'] = f"MDE for pressure perturbation, {iters[k_iter]} iteration"
        config['plotdata'] = []
        for k in range(len(kfs)):
            # noinspection PyTypeChecker
            config['plotdata'].append({'x': theta[k], 'y': lhs[k][k_iter], 'args': {'lw': 0.5},
                                       'label': '$q = ' + q[k] + r', B_{0}^{-1} \V{B}_{0} \cdot \grad p_{n}$'})
            config['plotdata'].append({'x': theta[k], 'y': lhs[k][k_iter], 'args': {'lw': 0.5},
                                       'label': '$q = ' + q[k] + r', -B_{0}^{-1} B_{n}^{\psi} \partial_{\psi} p_{0}$'})
        plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_MDE_pn_{iters[k_iter]}.pdf", config))
        # grad p_n
        config['postprocess'].append([LogY(), Id()])
        config['ylabel'] = [r'$\abs \grad p_{n}$ / \si{\newton\per\cubic\meter}', r'$\arg \grad p_{n}$ / \si{\degree}']
        config['title'] = f"Solution for pressure perturbation, {iters[k_iter]} perturbation"
        config['plotdata'] = []
        comps = [r'\partial_{R} p_{n}', r'\tfrac{\im n}{R} p_{n}', r'\partial_{Z} p_{n}']
        for k in range(len(kfs)):
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': theta[k], 'y': grad_pn[k][k_iter][k_comp, :], 'args': {'lw': 0.5},
                                           'label': f"$q = {q[k]}, {comps[k_comp]}$"})
        plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_grad_pn_{iters[k_iter]}.pdf", config))
        # iMHD
        config['ylabel'] = [r'$\abs f$ / \si{\newton\per\cubic\meter}', r'$\arg f$ / \si{\degree}']
        comps = [r'$R$ comp.', r'$\phi$ comp.', r'$Z$ comp.']
        for k in range(len(kfs)):
            config['title'] = f"Linearized iMHD force balance ($q = {q_half[k]}$), {iters[k_iter]} iteration"
            config['plotdata'] = []
            for k_comp in range(len(comps)):
                config['plotdata'].append({'x': theta_half[k], 'y': lorentz[k][k_iter][k_comp, :], 'args': {'lw': 0.5},
                                           'label': f"Lorentz force, {comps[k_comp]}"})
                config['plotdata'].append({'x': theta_half[k], 'y': grad_pn_half[k][k_iter][k_comp, :],
                                           'label': f"pressure gradient, {comps[k_comp]}", 'args': {'lw': 0.5}})
            plotter.plot_objects.put(ComplexPlot(work_dir, f"plot_MHD_{res_nonres[k]}_{iters[k_iter]}.pdf", config))

    # GPEC comparison
    conversion = 1.0e-04 / testcase.data['/mesh/gpec_jacfac'][:, 16]
    mephit_Bmn = testcase.get_polmodes('full perturbation (MEPHIT)', '/postprocess/Bmn/coeff_rad', conversion)
    mephit_Bmn_vac = testcase.get_polmodes('vacuum perturbation (MEPHIT)', '/postprocess/Bmn_vac/coeff_rad', conversion)
    reference = Gpec(work_dir, 2)
    reference.read_datafile()
    gpec_Bmn = reference.get_polmodes('full perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi')
    gpec_Bmn_vac = reference.get_polmodes('vacuum perturbation (GPEC)', sgn_dpsi, 'Jbgradpsi_x')
    config = {
        'xlabel': r'$\hat{\psi}$',
        'ylabel': r'$\abs [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\tesla}',
        'rho': testcase.post['psi_norm'], 'q': testcase.data['/cache/fs/q'][()],
        'resonances': testcase.post['psi_norm_res'], 'sgn_m_res': testcase.post['sgn_m_res'], 'omit_res': False,
        'poldata': [mephit_Bmn_vac, gpec_Bmn_vac, mephit_Bmn, gpec_Bmn], 'comp': abs,
        'res_neighbourhood': testcase.post['psi_norm_res_neighbourhood']
    }
    plotter.plot_objects.put(PolmodePlots(work_dir, 'GPEC_Bmn_psi_abs.pdf', config))
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'GPEC_Bmn_psi_arg.pdf', config))

    # iteration comparisons
    n = testcase.data['/config/n'][()]
    m = 3 * testcase.post['sgn_m_res']
    niter = testcase.data['/iter/niter'][()]
    nflux = testcase.data['/mesh/nflux'][()]
    abs_pmn_iter = full((niter, nflux + 1), nan)
    abs_jmnpar_Bmod_iter = full((niter, nflux + 1), nan)
    abs_Bmn_rad_iter = full((niter, nflux), nan)
    for kiter in range(0, niter):
        pmn = testcase.get_polmodes(None, f"/iter/pmn_{kiter:03}/coeff", 0.1, L1=True)
        jmnpar_Bmod = testcase.get_polmodes(None, f"/iter/jmnpar_Bmod_{kiter:03}/coeff", 1.0 / 29.9792458, L1=True)
        Bmn_rad = testcase.get_polmodes(None, f"/iter/Bmn_{kiter:03}/coeff_rad", conversion)
        abs_pmn_iter[kiter, :] = abs(pmn['var'][m])
        abs_jmnpar_Bmod_iter[kiter, :] = abs(jmnpar_Bmod['var'][m])
        abs_Bmn_rad_iter[kiter, :] = abs(Bmn_rad['var'][m])
    config = {
        'title': f"Poloidal Modes in Preconditioned Iterations $k$ for $n = {n}$, $m = {m}$",
        'xlabel': r'$\hat{\psi}$',
        'ylabel': [
            r'$\abs\, p_{mn}$ / \si{\pascal}',
            r'$\abs\, \bigl[ J_{n}^{\parallel} B_{0}^{-1} \bigr]_{m}$ / \si{\per\henry}',
            r'$\abs\, [B_{n}^{\perp} \sqrt{g} \lVert \nabla \psi \rVert]_{m} A^{-1}$ / \si{\tesla}'
        ],
        'rho': [pmn['rho'][m], jmnpar_Bmod['rho'][m], Bmn_rad['rho'][m]],
        'yscale': ['log', 'log', 'linear'],
        'res_pos': testcase.post['psi_norm_res'][abs(m)],
        'niter': niter,
        'plotdata': [abs_pmn_iter, abs_jmnpar_Bmod_iter, abs_Bmn_rad_iter]
    }
    plotter.plot_objects.put(IterationPlots(work_dir, 'plot_iter.pdf', config))

    # convergence estimation
    niter = testcase.data['/iter/niter'][()]
    sup_eigval = testcase.data['/config/ritz_threshold'][()]
    L2int_Bnvac = testcase.data['/iter/L2int_Bnvac'][()]
    L2int_Bn_diff = testcase.data['/iter/L2int_Bn_diff'][:niter]
    rel_err = testcase.data['/config/iter_rel_err'][()]
    plotter.plot_objects.put(Plot1D.conv_plot(work_dir, sup_eigval, L2int_Bnvac, L2int_Bn_diff, rel_err))

    plotter.finish()
