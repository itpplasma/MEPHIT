#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:39:22 2019

@author: lainer_p
"""

from sys import argv
from os import path, cpu_count
import re
from enum import Enum
import f90nml.parser
import matplotlib
import matplotlib.pyplot as plt
import colorcet
import numpy as np
from scipy import interpolate
from scipy import optimize
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from multiprocessing import Pool

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.use('Agg')
# =============================================================================
# pgf_config = {
#     'pgf.texsystem': 'lualatex',
#     'pgf.rcfonts': False,
#     'pgf.preamble': [
#         r'\usepackage{amsmath}',
#         r'\usepackage[math-style=ISO, bold-style=ISO]{unicode-math}'
#         r'\usepackage{siunitx}',
#     ]
# }
# matplotlib.use('pgf')
# matplotlib.rcParams.update(pgf_config)
# =============================================================================


class magdif_2d_triplot:
    scifmt = matplotlib.ticker.ScalarFormatter()
    scifmt.set_powerlimits((-3, 4))

    def __init__(self, node, tri, data, filename, title=None, clim_scale=None):
        self.node = node
        self.tri = tri
        self.data = data
        self.title = title
        self.filename = filename
        if clim_scale is not None:
            self.clim_scale = clim_scale
        else:
            self.clim_scale = 1.0

    def dump_plot(self):
        print('plotting ', self.filename)
        plt.figure(figsize=(3.3, 4.4))
        plt.tripcolor(
                self.node[:, 0], self.node[:, 1],
                self.tri - 1, self.data, cmap=colorcet.cm.coolwarm
        )
        plt.gca().set_aspect('equal')
        plt.colorbar(format=self.__class__.scifmt)
        plt.clim([-max(abs(self.data)) * self.clim_scale,
                  max(abs(self.data)) * self.clim_scale])
        plt.xlabel(r'$R$ / cm')
        plt.ylabel(r'$Z$ / cm')
        if self.title is not None:
            plt.title(self.title)
        plt.savefig(self.filename, dpi=300)
        plt.close()


class magdif_1d_cutplot:
    def __init__(self, x, xlabel, y, ylabel, title, filename):
        self.x = x
        self.y = y
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.filename = filename

    def dump_plot(self):
        print('plotting ', self.filename)
        plt.figure(figsize=(6.6, 3.6))
        plt.plot(self.x, self.y, '-k')
        plt.ticklabel_format(style='sci', scilimits=(-3, 4))
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.savefig(self.filename)
        plt.close()


class magdif_conv_plot:
    def __init__(self, datadir, conv_file, max_eigval, xlim=None, ylim=None,
                 title=None):
        self.datadir = datadir
        self.conv_file = conv_file
        self.max_eigval = max_eigval
        self.xlim = xlim
        self.ylim = ylim
        self.title = title

    def dump_plot(self):
        try:
            conv = np.loadtxt(path.join(self.datadir, self.conv_file))
        except Exception as err:
            print('Error: {}'.format(err))
            return
        niter = len(conv) - 1
        kiter = np.arange(0, niter + 1)
        plt.figure(figsize=(3.2, 3.2))
        plt.semilogy(
                kiter, conv[0] * self.max_eigval ** kiter, 'r-',
                label=r'$\vert \lambda_{\mathrm{max}} \vert^{k} \Vert'
                + r' \mathbf{B}_{n}^{(0)} \Vert_{2}$'
        )
        plt.semilogy(
                kiter, conv, 'xk',
                label=r'$\Vert \mathbf{B}_{n}^{(k)} \Vert_{2}$'
        )
        if self.xlim is not None:
            plt.xlim(self.xlim)
        if self.ylim is not None:
            plt.ylim(self.ylim)
        plt.gca().legend(loc='upper right')
        plt.xticks(kiter)
        plt.xlabel('iteration step $k$')
        plt.ylabel(r'$\Vert \delta \mathbf{B} \Vert_{2}$ / G cm')
        if self.title is not None:
            plt.title(self.title)
        else:
            plt.title('estimation of convergence')

        plt.tight_layout()
        plt.savefig(path.join(self.datadir, 'convergence.pdf'))
        plt.close()


class fslabel(Enum):
    psi_norm = r'$\hat{\psi}$'
    r = r'$r$ / cm'


class magdif_mnDat:
    def __init__(self, datadir, datafile, q_scaling, rad_coord, label):
        self.datadir = datadir
        self.datafile = datafile
        self.q_scaling = q_scaling
        self.rad_coord = rad_coord
        self.label = label

    def process(self):
        print('reading poloidal modes from', self.datafile)
        try:
            data = np.loadtxt(path.join(self.datadir, self.datafile))
        except Exception as err:
            print('Error: {}'.format(err))
            return
        if self.rad_coord is fslabel.psi_norm:
            # normalize psi
            self.rho = (data[:, 0] - data[0, 0]) / (data[-1, 0] - data[0, 0])
        else:
            self.rho = data[:, 0]
        # sort data
        sorting = np.argsort(self.rho)
        self.rho = self.rho[sorting]
        data = data[sorting, :]
        # process remaining data
        if self.q_scaling == 0:
            self.q = data[:, 1]
        else:
            self.q = data[:, 1] * self.q_scaling
        self.range = (np.shape(data)[1] - 2) // 2
        self.m_max = (self.range - 1) // 2
        self.offset = self.m_max
        self.abs = np.hypot(
                data[:, 2:(2 + self.range)], data[:, (2 + self.range):]
        )
        self.ymax = np.amax(self.abs, axis=0)
        self.ymax = np.fmax(self.ymax[self.offset:-1],
                            self.ymax[self.offset:0:-1])


class magdif_poloidal_plots:
    def __init__(self, n, psi, q, resonance, ylabel, data, ref=None,
                 r_interp=None):
        self.n = n
        self.psi = psi
        self.q = q
        self.resonance = resonance
        self.ylabel = ylabel
        self.data = data
        self.ref = ref
        self.r_interp = r_interp

    def dump_plot(self):
        if self.r_interp is None:
            rho = self.psi
            resonance = self.resonance
            xlabel = fslabel.psi_norm.value
        else:
            rho = self.r_interp(self.psi)
            resonance = {}
            for k, v in self.resonance.items():
                resonance[k] = self.r_interp(v)
            xlabel = fslabel.r.value

        # plot non-symmetric modes
        fmt = path.join(self.data.datadir,
                        path.splitext(self.data.datafile)[0] + '_{}.pdf')
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        for m_abs in range(1, self.data.m_max):
            if self.ref is not None and m_abs <= self.ref.m_max:
                ymax = max(self.data.ymax[m_abs], self.ref.ymax[m_abs])
            else:
                ymax = self.data.ymax[m_abs]
            plt.figure(figsize=two_squares)
            for k in range(horz_plot):
                ax = plt.subplot(vert_plot, horz_plot, k+1)
                m = (2 * k - 1) * m_abs
                if m in resonance:
                    ax.axvline(resonance[m], color='b', alpha=0.5)
                # q is negative in result_spectrum.f90, so index is reversed
                if self.ref is not None and m_abs <= self.ref.m_max:
                    ax.plot(self.ref.rho, self.ref.abs[:, self.ref.offset - m],
                            'r--', label=self.ref.label)
                ax.plot(self.data.rho, self.data.abs[:, self.data.offset - m],
                        label=self.data.label)
                ax.legend()
                ax.ticklabel_format(style='sci', scilimits=(-3, 4))
                ax.set_ylim(0.0, ymax)
                ax.set_title('$m = {}$'.format(m))
                ax.set_ylabel(self.ylabel)
                ax.set_xlabel(xlabel)
            plt.tight_layout()
            plt.savefig(fmt.format(m))
            plt.close()
        # plot symmetric mode and safety factor
        plt.figure(figsize=two_squares)
        ax = plt.subplot(vert_plot, horz_plot, 1)
        if self.ref is not None:
            ax.plot(self.ref.rho, self.ref.abs[:, self.ref.offset],
                    'r--', label=self.ref.label)
        ax.plot(self.data.rho, self.data.abs[:, self.data.offset],
                label=self.data.label)
        ax.legend()
        ax.ticklabel_format(style='sci', scilimits=(-3, 4))
        ax.set_title('$m = 0$')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(self.ylabel)
        ax = plt.subplot(vert_plot, horz_plot, 2)
        # q is negative in result_spectrum.f90
        ax.plot(self.data.rho, np.abs(self.data.q), label='kinetic')
        ax.plot(rho, self.q, 'r--', label='MHD')
        ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$q$')
        plt.tight_layout()
        plt.savefig(fmt.format(0))
        plt.close()


class magdif:

    vector_re = re.compile(r'(?<=plot_)([^/]+)\.dat$')
    scalar_re = re.compile(r'presn([^/]*)\.dat$')
    default_re = re.compile(r'(plot_)?([^/]+)\.dat$')

    def __init__(self, datadir, configfile, meshfile):
        self.plots = []
        self.datadir = datadir
        self.configfile = configfile
        self.meshfile = meshfile

    def read_configfile(self):
        print('reading configuration from ', self.configfile)
        p = f90nml.parser.Parser()
        nml = p.read(self.configfile)
        self.config = nml['settings']

    def read_fluxvar(self):
        print('reading contents of ', self.config['fluxvar_file'])
        fluxvar = np.loadtxt(path.join(self.datadir,
                                       self.config['fluxvar_file']))
        self.r = fluxvar[:, 0]
        self.psi = fluxvar[:, 1]
        self.q = fluxvar[:, 2]
        self.dens = fluxvar[:, 3]
        self.temp = fluxvar[:, 4]
        self.pres0 = fluxvar[:, 5]
        self.psi_norm = (self.psi - self.psi[0]) / (self.psi[-1] - self.psi[0])
        self.r_interp = interpolate.interp1d(self.psi_norm, self.r,
                                             kind='cubic')

    def load_mesh(self):
        with open(self.meshfile, 'r') as f:
            data = f.readlines()
        print('meshfile header: ', data[0])
        [NN, NT, NE] = np.loadtxt(StringIO(data[0]), dtype=int)
        node = np.loadtxt(StringIO(''.join(data[1:NN+1])), dtype=float)
        self.node = node[:, 0:2]
        tri = np.loadtxt(StringIO(''.join(data[NN+1:NN+NT+1])), dtype=int)
        self.tri = tri[:, 0:3]
        edge = np.loadtxt(StringIO(''.join(data[NN+NT+1:])), dtype=int)
        self.edge = edge[:, 0:2]

    def calc_resonances(self):
        n = self.config['n']
        if self.config['kilca_scale_factor'] == 0:
            q = self.q
        else:
            q = self.q * self.config['kilca_scale_factor']
        q_min = np.amin(q)
        q_max = np.amax(q)
        psi_min = np.amin(self.psi_norm)
        psi_max = np.amax(self.psi_norm)
        q_interp = interpolate.interp1d(self.psi_norm, q, kind='cubic')
        m_resonant = -np.arange(
                np.ceil(q_min * n),
                np.floor(q_max * n) + 1,  # +1 to include end point
                dtype=int
        )
        q_resonant = -m_resonant / n
        self.resonance = {}
        for k, m in enumerate(m_resonant):
            def q_resonant_interp(x):
                return q_interp(x) - q_resonant[k]
            self.resonance[m] = optimize.brentq(
                    q_resonant_interp, psi_min, psi_max)
        print(self.resonance.values())
        print(np.array([0.608, 0.760, 0.823, 0.861, 0.891, 0.918])**2)

    @classmethod
    def decorate_filename_vectorplot(cls, datafile, infix):
        return cls.vector_re.sub(r'\1' + infix + '.pdf', datafile)

    @classmethod
    def decorate_filename_presnplot(cls, datafile, infix):
        return cls.scalar_re.sub(r'plot_presn\1' + infix + '.pdf', datafile)

    @classmethod
    def decorate_default(cls, datafile, infix):
        return cls.default_re.sub(r'plot_\2' + infix + '.pdf', datafile)

    def generate_2d_triplots(
            self, datafile, start_column, infix_list, filename_decorator):
        print('reading contents of ', datafile)
        try:
            contents = np.loadtxt(path.join(self.datadir, datafile))
        except Exception as err:
            print('Error: {}'.format(err))
            return
        for column, infix in enumerate(infix_list, start_column):
            self.plots.append(magdif_2d_triplot(
                node=self.node, tri=self.tri,
                data=contents[:, column],
                filename=filename_decorator(datafile, infix)
                ))

    def generate_default_plots(self):
        self.plots.append(magdif_1d_cutplot(
                self.r, r'$r$ / cm', self.psi, r'$\psi$ / Mx',
                'disc poloidal flux', 'plot_psi.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.r, r'$r$ / cm', self.q, r'$q$',
                'safety factor', 'plot_q.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.r, r'$r$ / cm', self.dens,
                r'$n$ / cm\textsuperscript{-3}',
                'particle density', 'plot_dens.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.r, r'$r$ / cm', self.temp, r'$T$ / eV',
                'temperature', 'plot_temp.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.r, r'$r$ / cm', self.pres0,
                r'$p_{0}$ / dyn cm\textsuperscript{-2}',
                'pressure', 'plot_pres0.pdf'
        ))
        self.generate_2d_triplots(
                'j0phi.dat', 3, [''], self.__class__.decorate_default
        )

        Bn_datafiles = ['plot_Bn.dat', 'plot_Bn_000.dat', 'plot_Bn_vac.dat']
        currn_datafiles = ['plot_currn.dat', 'plot_currn_000.dat']
        presn_datafiles = ['presn.dat', 'presn_000.dat']
        vector_infix = [
                '_R_Re', '_R_Im', '_Z_Re', '_Z_Im', '_phi_Re', '_phi_Im',
                '_contradenspsi_Re', '_contradenspsi_Im',
                '_cotheta_Re', '_cotheta_Im'
        ]
        scalar_infix = ['_Re', '_Im']

        for datafile in Bn_datafiles + currn_datafiles:
            self.generate_2d_triplots(
                    datafile, 2, vector_infix,
                    self.__class__.decorate_filename_vectorplot
            )
        for datafile in presn_datafiles:
            self.generate_2d_triplots(
                    datafile, 0, scalar_infix,
                    self.__class__.decorate_filename_presnplot
            )

        self.plots.append(magdif_conv_plot(
                self.datadir, 'convergence.dat', self.config['tol'])
        )

        if self.config['kilca_scale_factor'] == 0:
            pert = magdif_mnDat(self.datadir, 'Bmn_psi.dat', 0,
                                fslabel.psi_norm, 'full perturbation')
            pert.process()
            vac = magdif_mnDat(self.datadir, 'Bmn_vac_psi.dat', 0,
                               fslabel.psi_norm, 'vacuum perturbation')
            vac.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config['n'], self.psi_norm, self.q, self.resonance,
                    r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
                    pert, vac
            ))
            pert = magdif_mnDat(self.datadir, 'currmn_000_theta.dat', 0,
                                fslabel.psi_norm, 'initial perturbation')
            pert.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config['n'], self.psi_norm, self.q, self.resonance,
                    r'$\left\vert J_{mn \theta}^{(0)} \right\vert$'
                    + r' / statA cm\textsuperscript{-1}', pert
            ))

        else:
            pert = magdif_mnDat(self.datadir, 'Bmn_r.dat',
                                self.config['kilca_scale_factor'],
                                fslabel.r, 'full perturbation')
            pert.process()
            vac = magdif_mnDat(self.datadir, 'Bmn_vac_r.dat',
                               self.config['kilca_scale_factor'],
                               fslabel.r, 'vacuum perturbation')
            vac.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config['n'], self.psi_norm, self.q
                    * self.config['kilca_scale_factor'], self.resonance,
                    r'$\left\vert B_{mn}^{r} \right\vert$ / G',
                    pert, vac, self.r_interp
            ))

    def dump_plots(self):
        for p in self.plots:
            p.dump_plot()

    def wrapper(self, index):
        self.plots[index].dump_plot()

    def dump_plots_parallel(self):
        with Pool(max(1, cpu_count() - 1)) as p:
            p.map(self.wrapper, range(len(self.plots)))


if __name__ == '__main__':
    testcase = magdif(argv[1], argv[2], argv[3])
    testcase.read_configfile()
    testcase.read_fluxvar()
    testcase.calc_resonances()
    testcase.load_mesh()
    testcase.generate_default_plots()
    testcase.dump_plots_parallel()
