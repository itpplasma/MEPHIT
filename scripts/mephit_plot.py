import multiprocessing as mp
from os import path

# supply path to scripts directory
scripts_dir = path.dirname(path.realpath(__file__))
run_dir = path.realpath(scripts_dir + '/../run')


def set_matplotlib_defaults():
    from cycler import cycler
    from matplotlib import rcParams
    rcParams['text.usetex'] = True
    rcParams['font.family'] = 'serif'
    rcParams['mathtext.fontset'] = 'cm'
    rcParams['figure.constrained_layout.use'] = True
    rcParams['xtick.direction'] = 'in'
    rcParams['xtick.top'] = True
    rcParams['ytick.direction'] = 'in'
    rcParams['ytick.right'] = True
    rcParams['lines.linewidth'] = 0.75
    rcParams['axes.formatter.limits'] = (-3, 4)
    rcParams['axes.prop_cycle'] = (cycler('color', ['k', 'tab:orange', 'tab:purple', 'tab:red', 'tab:cyan',
                                                    'tab:brown', 'tab:gray', 'tab:green', 'tab:pink', 'tab:blue']) +
                                   cycler('ls', ['-', '--', '-.', ':', (0, (5, 1.2, 1, 1.2, 1, 1.2)),
                                                 '-', '--', '-.', ':', (0, (5, 1.2, 1, 1.2, 1, 1.2))]))
    latex_preamble = path.join(scripts_dir, 'magdifplot.tex')
    rcParams['text.latex.preamble'] = fr"\input{{{latex_preamble}}}"


class Mephit:
    def __init__(self, work_dir):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.data = None
        self.post = dict()

    def __del__(self):
        from h5pickle import File
        if isinstance(self.data, File):
            self.data.close()

    def read_datafile(self):
        from h5py import get_config as h5py_hack
        from h5pickle import File
        from os import path
        h5py_hack().complex_names = ('real', 'imag')  # complex values are stored as compound types in libneo/hdf5tools
        self.data = File(path.join(self.work_dir, 'mephit.h5'), 'r')

    def normalize_psi(self, arr):
        return (arr - self.data['/cache/fs/psi'][0]) / (self.data['/cache/fs/psi'][-1] - self.data['/cache/fs/psi'][0])

    def postprocess(self):
        from matplotlib.tri import Triangulation
        from scipy import interpolate
        from numpy import arange, sign
        self.post['triangulation'] = Triangulation(self.data['/mesh/node_R'][()], self.data['/mesh/node_Z'][()],
                                                   self.data['/mesh/tri_node'][()] - 1)
        self.post['psi_norm'] = self.normalize_psi(self.data['/cache/fs/psi'][()])
        self.post['sgn_m_res'] = int(sign(-self.data['/cache/fs/q'][-1]))
        m_res = arange(self.data['/mesh/m_res_min'][()], self.data['/mesh/m_res_max'][()] + 1) * self.post['sgn_m_res']
        self.post['psi_norm_res'] = dict(zip(m_res, self.normalize_psi(self.data['/mesh/psi_res'][()])))

    def get_polmodes(self, label, var_name='/postprocess/Bmn/coeff_rad', conversion=1.0):
        from numpy import array
        polmodes = {'m_max': 0, 'label': label, 'rho': dict(), 'var': dict()}
        polmodes['m_max'] = (self.data[var_name].shape[1] - 1) // 2
        rho = self.normalize_psi(self.data['/cache/fs_half/psi'][()])
        for m in range(-polmodes['m_max'], polmodes['m_max'] + 1):
            polmodes['rho'][m] = rho
            polmodes['var'][m] = array(self.data[var_name][:, m + polmodes['m_max']], dtype='D') * conversion
        return polmodes


class Gpec:
    def __init__(self, work_dir, n):
        from os import getcwd
        self.work_dir = work_dir or getcwd()
        self.n = n or 2
        self.data = dict()
        self.post = dict()

    def __del__(self):
        from netCDF4 import Dataset
        for file in self.data.keys():
            if isinstance(self.data[file], Dataset):
                self.data[file].close()

    def read_datafile(self):
        from os import path
        from netCDF4 import Dataset
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
        self.num_processes = max(1, mp.cpu_count() - 1)
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
        while True:
            plot_object = plot_objects_queue.get()
            if isinstance(plot_object, PlotObject):
                plot_object.do_plot()
            else:
                results_queue.put(plot_object)
                break
        return


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
        if 'xlabel' in self.config.keys():
            ax.set_xlabel(self.config['xlabel'])
        if 'ylabel' in self.config.keys():
            ax.set_ylabel(self.config['ylabel'])
        if 'legend' in self.config.keys():
            ax.legend(**self.config['legend'])
        if 'title' in self.config.keys():
            ax.set_title(self.config['title'])
        if 'postprocess' in self.config.keys():
            for f in self.config['postprocess']:
                f(fig, ax)
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
        for data in self.config['poldata']:
            if m in data['var'].keys():
                axs[0].plot(data['rho'][m], self.config['comp'](data['var'][m]), label=data['label'])
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
        m_max = min(map(lambda d: d['m_max'], self.config['poldata']))
        for m_abs in range(1, m_max + 1):
            print(f"Plotting {self.filename}, m = {m_abs} ...")
            fig = Figure(figsize=two_squares)
            axs = fig.subplots(vert_plot, horz_plot, sharex='all', sharey='all')
            for k in range(horz_plot):
                m = (2 * k - 1) * m_abs * self.config['sgn_m_res']
                axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
                if m in self.config['resonances'] and not self.config['omit_res']:
                    axs[k].axvline(self.config['resonances'][m],
                                   color='b', alpha=0.5, lw=0.5, label='resonance position')
                for data in self.config['poldata']:
                    if m in data['var'].keys():
                        axs[k].plot(data['rho'][m], self.config['comp'](data['var'][m]), label=data['label'])
                axs[k].legend(loc='upper left', fontsize='x-small')
                axs[k].set_title(('resonant ' if m in self.config['resonances'] else 'non-resonant ') + fr"$m = {m}$")
                axs[k].set_xlabel(self.config['xlabel'])
                axs[k].set_ylabel(self.config['ylabel'])
            axs[1].yaxis.set_tick_params(labelleft=True)
            canvas = FigureCanvas(fig)
            fig.savefig(pdf, format='pdf', dpi=300)
        pdf.close()
