from mephit_plot import Gpec, Mephit, ParallelPlotter, Plot1D, PlotObject, PolmodePlots, run_dir, \
    set_matplotlib_defaults
from numpy import abs, angle, full, nan, sign
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


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/33353_2325'
    testcase = Mephit(work_dir)
    testcase.read_datafile()
    testcase.postprocess()
    sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])

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
    mephit_jmn_000 = testcase.get_polmodes('initial perturbation', '/postprocess/jmn_000/coeff_pol')
    mephit_jmn = testcase.get_polmodes('full perturbation', '/postprocess/jmn/coeff_pol')
    config['comp'] = abs
    config['ylabel'] = r'$\abs J_{mn\vartheta}$ / \si{\statampere\per\centi\meter\cubed}'
    config['poldata'] = [mephit_jmn_000, mephit_jmn]
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_jmn_theta_abs.pdf', config))
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg J_{mn\vartheta}$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_jmn_theta_arg.pdf', config))
    mephit_pmn_000 = testcase.get_polmodes('initial perturbation', '/postprocess/pmn_000/coeff', L1=True)
    mephit_pmn = testcase.get_polmodes('full perturbation', '/postprocess/pmn/coeff', L1=True)
    config['comp'] = abs
    config['ylabel'] = r'$\abs p_{mn}$ / \si{\dyne\per\centi\meter\squared}'
    config['poldata'] = [mephit_pmn_000, mephit_pmn]
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_pmn_abs.pdf', config))
    config['comp'] = partial(angle, deg=True)
    config['ylabel'] = r'$\arg p_{mn}$ / \si{\degree}'
    plotter.plot_objects.put(PolmodePlots(work_dir, 'plot_pmn_arg.pdf', config))
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
