from mephit_plot import Gpec, Mephit, ParallelPlotter, PolmodePlots, Plot1D, run_dir, set_matplotlib_defaults
from numpy import abs, angle, sign
from functools import partial


if __name__ == "__main__":
    set_matplotlib_defaults()
    plotter = ParallelPlotter()
    plotter.start()

    work_dir = run_dir + '/33353_2325'
    testcase = Mephit(work_dir)
    testcase.read_datafile()
    testcase.postprocess()
    sgn_dpsi = sign(testcase.data['/cache/fs/psi'][-1] - testcase.data['/cache/fs/psi'][0])
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
    niter = testcase.data['/iter/niter'][()]
    sup_eigval = testcase.data['/config/ritz_threshold'][()]
    L2int_Bnvac = testcase.data['/iter/L2int_Bnvac'][()]
    L2int_Bn_diff = testcase.data['/iter/L2int_Bn_diff'][:niter]
    rel_err = testcase.data['/config/iter_rel_err'][()]
    plotter.plot_objects.put(Plot1D.conv_plot(work_dir, sup_eigval, L2int_Bnvac, L2int_Bn_diff, rel_err))

    plotter.finish()