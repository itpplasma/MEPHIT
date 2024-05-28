import multiprocessing as mp
from os import path

# supply path to scripts directory
scripts_dir = path.dirname(path.realpath(__file__))
run_dir = path.realpath(scripts_dir + '/../build/run')


def set_matplotlib_defaults():
    from cycler import cycler
    from matplotlib import rcParams
    rcParams['text.usetex'] = True
    # rcParams['font.family'] = 'serif'
    rcParams['mathtext.fontset'] = 'cm'
    rcParams['figure.constrained_layout.use'] = True
    rcParams['xtick.direction'] = 'in'
    rcParams['xtick.top'] = True
    rcParams['ytick.direction'] = 'in'
    rcParams['ytick.right'] = True
    rcParams['lines.linewidth'] = 0.75
    rcParams['axes.formatter.limits'] = (-3, 4)
    dash_dot_dot = (0, (5, 1.2, 1, 1.2, 1, 1.2))
    dash_dash_dot = (0, (5, 1.2, 5, 1.2, 1, 1.2))
    # https://www.w3schools.com/colors/colors_nbs.asp
    # black, vivid orange, strong purple, vivid yellow, vivid light blue, vivid red,
    # grayish yellow, medium gray, vivid green, strong purplish pink, strong blue, strong yellowish pink
    rcParams['axes.prop_cycle'] = (cycler('color', ['#222222', '#f38400', '#875692', '#f3c300', '#a1caf1', '#be0032',
                                                    '#c2b280', '#848482', '#008856', '#e68fac', '#0067a5', '#f99379']) +
                                   cycler('ls', ['-', '--', '-.', ':', dash_dot_dot, dash_dash_dot,
                                                 '-', '--', '-.', ':', dash_dot_dot, dash_dash_dot]))
    rcParams['text.latex.preamble'] = fr"\usepackage{{import}}\import{{{scripts_dir}}}{{mephit-pdflatex.tex}}"


class Mephit:
    def __init__(self, work_dir, h5file='mephit.h5'):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.h5file = h5file
        self.data = None
        self.post = dict()

    def close_datafile(self):
        from h5pickle import File
        if isinstance(self.data, File):
            self.data.close()

    def open_datafile(self):
        from h5py import get_config as h5py_hack
        from h5pickle import File
        from os import path
        self.close_datafile()
        h5py_hack().complex_names = ('real', 'imag')  # complex values are stored as compound types in libneo/hdf5tools
        self.data = File(path.join(self.work_dir, self.h5file), 'r')

    def normalize_psi(self, arr):
        return (arr - self.data['/cache/fs/psi'][0]) / (self.data['/cache/fs/psi'][-1] - self.data['/cache/fs/psi'][0])

    def normalize_psi_diff(self, arr):
        return arr / (self.data['/cache/fs/psi'][-1] - self.data['/cache/fs/psi'][0])

    def postprocess(self):
        from matplotlib.tri import Triangulation
        from numpy import arange, sign
        self.post['triangulation'] = Triangulation(self.data['/mesh/node_R'][()], self.data['/mesh/node_Z'][()],
                                                   self.data['/mesh/tri_node'][()] - 1)
        self.post['psi_norm'] = self.normalize_psi(self.data['/cache/fs/psi'][()])
        self.post['psi_half_norm'] = self.normalize_psi(self.data['/cache/fs_half/psi'][()])
        self.post['rad'] = self.data['/cache/fs/rad'][()]
        self.post['rad_half'] = self.data['/cache/fs_half/rad'][()]
        self.post['sgn_m_res'] = int(sign(-self.data['/cache/fs/q'][-1]))
        self.post['m_res'] = self.post['sgn_m_res'] * arange(self.data['/mesh/m_res_min'][()],
                                                             self.data['/mesh/m_res_max'][()] + 1)
        self.post['psi_norm_res'] = dict(zip(self.post['m_res'], self.normalize_psi(self.data['/mesh/psi_res'][()])))
        self.post['rad_res'] = dict(zip(self.post['m_res'], self.data['/mesh/rad_norm_res'][()] * self.post['rad'][-1]))
        m_res_min = self.data['/mesh/m_res_min'][()]
        res_ind = self.data['/mesh/res_ind'][()] - 1
        nflux = self.data['/mesh/nflux'][()]
        self.post['psi_norm_res_neighbourhood'] = {}
        self.post['rad_res_neighbourhood'] = {}
        for m in self.post['m_res']:
            kf_min = res_ind[abs(m) - m_res_min] - 2
            kf_max = min(res_ind[abs(m) - m_res_min] + 3, nflux - 1)
            self.post['psi_norm_res_neighbourhood'][m] = self.post['psi_norm'][kf_min:kf_max+1]
            self.post['rad_res_neighbourhood'][m] = self.post['rad'][kf_min:kf_max+1]

    def get_polmodes(self, label, var_name='/iter/Bmn/coeff_rad', conversion=1.0e-4, L1=False, rad=False):
        from numpy import array
        polmodes = {'m_max': 0, 'label': label, 'rho': dict(), 'var': dict()}
        polmodes['m_max'] = (self.data[var_name].shape[1] - 1) // 2
        rho = self.post['rad' if rad else 'psi_norm'] if L1 else self.post['rad_half' if rad else 'psi_half_norm']
        for m in range(-polmodes['m_max'], polmodes['m_max'] + 1):
            polmodes['rho'][m] = rho
            polmodes['var'][m] = array(self.data[var_name][:, m + polmodes['m_max']], dtype='D') * conversion
        return polmodes

    def get_Ires(self):
        from numpy import abs
        from scipy.constants import c as clight
        scaling = self.data['/config/kilca_scale_factor'][()]
        scaling = scaling if scaling != 0 else 1
        return dict(zip(self.post['m_res'], abs(self.data['/iter/Ires'][()]) * scaling * 0.1 / clight))


class Kilca:
    def __init__(self, work_dir):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.data = None
        self.post = dict()

    def close_datafile(self):
        from h5pickle import File
        if isinstance(self.data, File):
            self.data.close()

    def open_datafile(self, datafile):
        from h5pickle import File
        from os import path
        self.close_datafile()
        self.data = File(path.join(self.work_dir, datafile), 'r')

    def get_polmodes(self, label, var_name='Br', conversion=1.0e-04):
        from numpy import array, sign, zeros
        polmodes = {'m_max': 0, 'label': label, 'rho': dict(), 'var': dict()}
        sgn_q = int(sign(self.data['/output/background/profiles/q_i'][0, -1]))
        for name, grp in self.data['/output'].items():
            if 'postprocessor' not in name:
                continue
            m = int(grp['mode'][0, 0]) * sgn_q
            polmodes['m_max'] = max(polmodes['m_max'], abs(m))
            polmodes['rho'][m] = array(grp['r'][0, :], dtype='d')
            polmodes['var'][m] = zeros(polmodes['rho'][m].shape, dtype='D')
            polmodes['var'][m].real = grp[var_name][0, :] * conversion
            if grp[var_name].shape[0] == 2:
                polmodes['var'][m].imag = grp[var_name][1, :] * conversion
        return polmodes

    def get_Ires(self):
        from numpy import sign
        sgn_q = int(sign(self.data['/output/background/profiles/q_i'][0, -1]))
        Ires = {}
        for name, grp in self.data['/output'].items():
            if 'postprocessor' not in name:
                continue
            m = int(grp['mode'][0, 0]) * sgn_q
            Ires[m] = grp['Ipar'][0, 0] * 10.0
        return Ires


class Gpec:
    def __init__(self, work_dir, n):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.n = n or 2
        self.data = dict()
        self.post = dict()

    def close_datafiles(self):
        from netCDF4 import Dataset
        for file in self.data.keys():
            if isinstance(self.data[file], Dataset):
                self.data[file].close()

    def open_datafiles(self):
        from os import path
        from netCDF4 import Dataset
        self.close_datafiles()
        for file in ['profile', 'cylindrical', 'control']:
            self.data[file] = Dataset(path.join(self.work_dir, f"gpec_{file}_output_n{self.n}.nc"), 'r')
        self.data['dcon'] = Dataset(path.join(self.work_dir, f"dcon_output_n{self.n}.nc"), 'r')

    def get_polmodes(self, label, sgn_dpsi=1, var_name='Jbgradpsi'):
        from numpy import array, empty
        polmodes = {'m_max': 0, 'label': label, 'rho': dict(), 'var': dict()}
        rho = array(self.data['profile'].variables['psi_n'])
        helicity = int(self.data['profile'].getncattr('helicity'))
        for k, m_out in enumerate(self.data['profile'].variables['m_out'][:]):
            m = m_out * helicity
            polmodes['m_max'] = max(polmodes['m_max'], abs(m))
            polmodes['rho'][m] = rho
            polmodes['var'][m] = empty(rho.shape, dtype='D')
            # GPEC always uses normal vectors pointing outwards
            # and includes the factor 2 for Fourier series of a real function in the coefficient
            polmodes['var'][m].real = self.data['profile'].variables[var_name][0, k, :] * 0.5 * sgn_dpsi
            # GPEC uses clockwise toroidal angle for positive helicity
            # and expands Fourier series in negative toroidal angle
            polmodes['var'][m].imag = self.data['profile'].variables[var_name][1, k, :] * 0.5 * sgn_dpsi * helicity
        return polmodes

    def get_Ires(self):
        from numpy import hypot
        helicity = int(self.data['profile'].getncattr('helicity'))
        Ires = {}
        for k, abs_m in enumerate(self.data['control'].variables['m_rational'][:]):
            m = abs_m * helicity
            Ires[m] = hypot(self.data['profile'].variables['I_res'][0, k],
                            self.data['profile'].variables['I_res'][1, k])
        return Ires


class Mars:
    def __init__(self, work_dir):
        self.work_dir = work_dir
        self.data = dict()
        self.post = dict()

    def open_datafiles(self):
        from os import path
        from scipy.io import loadmat
        self.data = loadmat(path.join(self.work_dir, 'OUTPUTS.mat')) | \
            loadmat(path.join(self.work_dir, 'INPUTS.mat'), squeeze_me=True) | \
            loadmat(path.join(self.work_dir, 'Normalizations.mat'))

    def get_polmodes(self, label, var_name='PLASMA'):
        from numpy.fft import fft
        polmodes = {'m_max': 0, 'label': label, 'rho': dict(), 'var': dict()}
        all_modes = self.data['Mag_field'][0, 0] * \
            fft(self.data[var_name]['B'][0, 0] * self.data['Jacobian'], norm='forward')
        for m in self.data['Mm']:
            polmodes['rho'][-m] = self.data['sp']
            polmodes['var'][-m] = all_modes[0, :self.data['sp'].size, m]
        return polmodes

class PlotObject:
    def __init__(self, work_dir, filename, config):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.filename = filename or 'plot.pdf'
        self.config = config or dict()

    def do_plot(self):
        print(f"Plotting {self.filename} ...")


class ParallelPlotter:
    def __init__(self):
        self.ctx = mp.get_context('fork')
        self.manager = self.ctx.Manager()
        self.plot_objects = self.manager.Queue()
        self.results = self.manager.Queue()
        self.num_processes = min(4, max(1, mp.cpu_count() - 1))  # use 1 to 4 processes
        self.pool = self.ctx.Pool(processes=self.num_processes)
        self.processes = []

    def start(self):
        mp.set_start_method('fork')
        for i in range(self.num_processes):
            new_process = self.ctx.Process(target=ParallelPlotter.plot_worker, args=(self.plot_objects, self.results))
            self.processes.append(new_process)
            new_process.start()

    def finish(self):
        for i in range(self.num_processes):
            self.plot_objects.put(0)
        num_finished_processes = 0
        while True:
            result = self.results.get()
            if not isinstance(result, PlotObject):
                num_finished_processes += 1
                if num_finished_processes == self.num_processes:
                    break

    @classmethod
    def plot_worker(cls, plot_objects_queue, results_queue):
        import traceback
        while True:
            plot_object = plot_objects_queue.get()
            if isinstance(plot_object, PlotObject):
                try:
                    plot_object.do_plot()
                except BaseException:
                    print(f"Error plotting '{plot_object.filename}'")
                    print(traceback.format_exc())
            else:
                results_queue.put(plot_object)
                break
        return


class Id:
    def __call__(self, fig, ax):
        pass


class XTicks:
    def __init__(self, ticks):
        self.ticks = ticks

    def __call__(self, fig, ax):
        ax.set_xticks(self.ticks)


class YTicks:
    def __init__(self, ticks):
        self.ticks = ticks

    def __call__(self, fig, ax):
        ax.set_yticks(self.ticks)


class LogY:
    def __call__(self, fig, ax):
        ax.set_yscale('log')


class HLine:
    def __init__(self, pos, **kwargs):
        self.pos = pos
        self.kwargs = kwargs

    def __call__(self, fig, ax):
        ax.axhline(self.pos, **self.kwargs)


class Plot1D(PlotObject):
    def do_plot(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from os import path
        super().do_plot()
        fig = Figure()
        ax = fig.subplots()
        for plotdata in self.config['plotdata']:
            ax.plot(plotdata['x'], plotdata['y'], **plotdata['args'])
        if 'postprocess' in self.config.keys():
            for f in self.config['postprocess']:
                f(fig, ax)
        if 'xlabel' in self.config.keys():
            ax.set_xlabel(self.config['xlabel'])
        if 'ylabel' in self.config.keys():
            ax.set_ylabel(self.config['ylabel'])
        if 'legend' in self.config.keys():
            ax.legend(**self.config['legend'])
        if 'title' in self.config.keys():
            ax.set_title(self.config['title'])
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(self.work_dir, self.filename))


class PolmodePlots(PlotObject):
    def do_plot(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.backends.backend_pdf import PdfPages
        from numpy import any
        super().do_plot()
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        pdf = PdfPages(path.join(self.work_dir, self.filename))
        # plot symmetric mode and safety factor
        m = 0
        fig = Figure(figsize=two_squares)
        axs = fig.subplots(vert_plot, horz_plot)
        axs[0].axhline(0.0, color='k', alpha=0.5, lw=0.5)
        for data in self.config['poldata']:
            if m in data['var'].keys():
                if 'zoom_x' in self.config.keys():
                    zoom_x = (self.config['zoom_x'][0] <= data['rho'][m]) & (data['rho'][m] <= self.config['zoom_x'][1])
                    axs[0].plot(data['rho'][m][zoom_x], self.config['comp'](data['var'][m][zoom_x]),
                                label=data['label'])
                else:
                    axs[0].plot(data['rho'][m], self.config['comp'](data['var'][m]), label=data['label'])
            else:
                next(axs[0]._get_lines.prop_cycler)
        if 'postprocess' in self.config.keys():
            for f in self.config['postprocess']:
                f(fig, axs[0])
        axs[0].legend(loc='upper left', fontsize='x-small')
        axs[0].set_title(f"$m = {m}$")
        axs[0].set_xlabel(self.config['xlabel'])
        axs[0].set_ylabel(self.config['ylabel'])
        for res in self.config['resonances'].values():
            axs[1].axvline(res, color='b', alpha=0.5, lw=0.5)
        axs[1].plot(self.config['rho'], self.config['q'], '-k')
        if any(self.config['q'] > 0.0):
            axs[1].set_ylim(bottom=0.0)
        else:
            axs[1].set_ylim(top=0.0)
        axs[1].set_xlabel(self.config['xlabel'])
        axs[1].set_ylabel(r'$q$')
        canvas = FigureCanvas(fig)
        fig.savefig(pdf, format='pdf', dpi=300)
        # plot non-symmetric modes
        m_max = max(map(lambda d: d['m_max'], self.config['poldata']))
        for m_abs in range(1, m_max + 1):
            print(f"Plotting {self.filename}, m = ±{m_abs} ...")
            fig = Figure(figsize=two_squares)
            axs = fig.subplots(vert_plot, horz_plot, sharex='all', sharey='all')
            for k in range(horz_plot):
                m = (2 * k - 1) * m_abs * self.config['sgn_m_res']
                axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
                if not self.config['omit_res']:
                    if m in self.config['resonances']:
                        axs[k].axvline(self.config['resonances'][m],
                                       color='b', alpha=0.5, lw=0.5, label='resonance position')
                    if m in self.config['resonances'] and 'res_neighbourhood' in self.config.keys():
                        for pos in self.config['res_neighbourhood'][m]:
                            axs[k].axvline(pos, color='k', alpha=0.5, lw=0.25)
                for data in self.config['poldata']:
                    if m in data['var'].keys():
                        if 'zoom_x' in self.config.keys():
                            zoom_x = (self.config['zoom_x'][0] <= data['rho'][m]) & \
                                     (data['rho'][m] <= self.config['zoom_x'][1])
                            axs[k].plot(data['rho'][m][zoom_x], self.config['comp'](data['var'][m][zoom_x]),
                                        label=data['label'])
                        else:
                            axs[k].plot(data['rho'][m], self.config['comp'](data['var'][m]), label=data['label'])
                    else:
                        next(axs[k]._get_lines.prop_cycler)
                if 'postprocess' in self.config.keys():
                    for f in self.config['postprocess']:
                        f(fig, axs[k])
                axs[k].legend(loc='upper left', fontsize='x-small')
                axs[k].set_title(('resonant ' if m in self.config['resonances'] else 'non-resonant ') + fr"$m = {m}$")
                axs[k].set_xlabel(self.config['xlabel'])
                axs[k].set_ylabel(self.config['ylabel'])
            axs[1].yaxis.set_tick_params(labelleft=True)
            canvas = FigureCanvas(fig)
            fig.savefig(pdf, format='pdf', dpi=300)
        pdf.close()


class IterationPlots(PlotObject):
    def do_plot(self):
        from os import path
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.backends.backend_pdf import PdfPages
        from numpy import fmax, fmin
        super().do_plot()
        horz_plot = 3
        vert_plot = 1
        three_squares = (9.9, 3.3)
        pdf = PdfPages(path.join(self.work_dir, self.filename))
        figs = [Figure(figsize=three_squares) for kiter in range(0, self.config['niter'])]
        for kiter in range(0, self.config['niter']):
            print(f"Plotting {self.filename}, k = {kiter + 1} ...")
            axs = figs[kiter].subplots(vert_plot, horz_plot)
            for k in range(horz_plot):
                axs[k].set_yscale(self.config['yscale'][k])
            for k in range(len(self.config['plotdata'])):
                if 'zoom_x' in self.config.keys():
                    zoom_x = (self.config['zoom_x'][0] <= self.config['rho'][k % horz_plot]) & \
                             (self.config['rho'][k % horz_plot] <= self.config['zoom_x'][1])
                    axs[k % horz_plot].plot(self.config['rho'][k % horz_plot][zoom_x],
                                            self.config['plotdata'][k][kiter, :][zoom_x], **self.config['plotargs'])
                else:
                    axs[k % horz_plot].plot(self.config['rho'][k % horz_plot], self.config['plotdata'][k][kiter, :],
                                            **self.config['plotargs'])
            for k in range(horz_plot):
                if 'res_pos' in self.config.keys():
                    axs[k].axvline(self.config['res_pos'], color='b', alpha=0.5, lw=0.5)
                if 'res_neighbourhood' in self.config.keys():
                    for pos in self.config['res_neighbourhood']:
                        axs[k].axvline(pos, color='k', alpha=0.5, lw=0.25)
                if 'postprocess' in self.config.keys():
                    for f in self.config['postprocess']:
                        f[k](figs[kiter], axs[k])
                axs[k].set_xlabel(self.config['xlabel'])
                axs[k].set_ylabel(self.config['ylabel'][k])
            figs[kiter].suptitle(self.config['title'] + f", $k = {kiter + 1}$")
        ylims = [figs[0].axes[k].get_ylim() for k in range(horz_plot)]
        if self.config['global_ylims']:
            for kiter in range(self.config['niter']):
                for k in range(horz_plot):
                    ylim = figs[kiter].axes[k].get_ylim()
                    ylims[k] = [fmin(ylims[k][0], ylim[0]), fmax(ylims[k][1], ylim[1])]
        for kiter in range(0, self.config['niter']):
            for k in range(horz_plot):
                figs[kiter].axes[k].set_ylim(ylims[k])
            canvas = FigureCanvas(figs[kiter])
            figs[kiter].savefig(pdf, format='pdf', dpi=300)
        pdf.close()


class ComplexPlot(PlotObject):
    def do_plot(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from os import path
        from numpy import abs, angle, arange
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
