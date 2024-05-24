from mephit_plot import Gpec, Mephit, ParallelPlotter, Plot1D, run_dir, set_matplotlib_defaults
from numpy import abs, array, full, nan, pi, sign
from numpy.fft import fft, rfft
from scipy.interpolate import UnivariateSpline


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
reference = Gpec(work_dir, 2)
reference.open_datafiles()

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

plotter.finish()
testcase.close_datafile()
reference.close_datafiles()
