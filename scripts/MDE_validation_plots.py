from mephit_plot import ComplexPlot, HLine, Id, LogY, Mephit, ParallelPlotter, \
    PolmodePlots, run_dir, set_matplotlib_defaults, XTicks, YTicks
from numpy import abs, angle, arange, full, nan, pi, r_, sum
from numpy.fft import fft
from numpy.linalg import norm
from scipy.constants import c as clight
from functools import partial


set_matplotlib_defaults()
plotter = ParallelPlotter()
plotter.start()

work_dir = run_dir + '/33353_2900'
testcase = Mephit(work_dir)
testcase.open_datafile()
testcase.postprocess()
nflux = testcase.data['/mesh/nflux'][()]
m_max = 24
m_res_min = testcase.data['/mesh/m_res_min'][()]

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

plotter.finish()
testcase.close_datafile()
