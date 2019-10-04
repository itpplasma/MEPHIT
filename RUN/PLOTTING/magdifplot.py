#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:39:22 2019

@author: lainer_p
"""

import sys
import os
import re
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

    def __init__(self, node, tri, data, filename, title=None):
        self.node = node
        self.tri = tri
        self.data = data
        self.title = title
        self.filename = filename

    def dump_plot(self):
        print('plotting ', self.filename)
        plt.figure(figsize=(3.3, 4.4))
        plt.tripcolor(
                self.node[:, 0], self.node[:, 1],
                self.tri - 1, self.data, cmap=colorcet.cm.coolwarm
        )
        plt.gca().set_aspect('equal')
        plt.colorbar(format=self.__class__.scifmt)
        plt.clim([-max(abs(self.data)), max(abs(self.data))])
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
    def __init__(self, datadir, conv_file, max_eigval):
        self.datadir = datadir
        self.conv_file = conv_file
        self.max_eigval = max_eigval

    def dump_plot(self):
        try:
            conv = np.loadtxt(os.path.join(self.datadir, self.conv_file))
        except Exception as err:
            print('Error: {}'.format(err))
            return
        niter = len(conv) - 1
        kiter = np.arange(1, niter + 1)
        plt.figure(figsize=(3.2, 3.2))
        plt.semilogy(
                kiter, conv[0] * self.max_eigval ** kiter, 'r-',
                label=r'$\vert \lambda_{\mathrm{max}} \vert^{k} \Vert'
                + r' \delta \mathbf{B}_{\mathrm{v}} \Vert_{2}$'
        )
        plt.semilogy(
                kiter, conv[1:], 'xk',
                label=r'$\Vert \delta \mathbf{B}^{(k)} \Vert_{2}$'
        )
        plt.gca().legend(loc='lower left')
        plt.xticks(kiter)
        plt.xlabel('iteration step $k$')
        plt.ylabel(r'$\Vert \delta \mathbf{B} \Vert_{2}$ / G cm')
        plt.title('estimation of convergence')

        plt.tight_layout()
        plt.savefig(os.path.join(self.datadir, 'convergence.pdf'))
        plt.close()


class magdif_poloidal_modes:
    def __init__(self, n, s, q, datadir, datafile, label, reffile=None,
                 r_s_interp=None):
        self.n = n
        self.s_mhd = s
        self.q_mhd = q
        self.datadir = datadir
        self.datafile = datafile
        self.label = label
        self.reffile = reffile
        self.r_s_interp = r_s_interp

    def dump_plot(self):
        print('plotting poloidal modes from', self.datafile)
        try:
            data = np.loadtxt(os.path.join(self.datadir, self.datafile))
        except Exception as err:
            print('Error: {}'.format(err))
            return
        if self.reffile is not None:
            try:
                ref = np.loadtxt(os.path.join(self.datadir, self.reffile))
            except Exception as err:
                print('Error: {}'.format(err))
                return

        # normalize psi
        s = (data[:, 0] - data[0, 0]) / (data[-1, 0] - data[0, 0])
        q = data[:, 1]
        xlabel = r'$\hat{\psi}$'

        q_min = np.amin(self.q_mhd)
        q_max = np.amax(self.q_mhd)
        q_interp = interpolate.interp1d(self.s_mhd, self.q_mhd, kind='cubic')

        m_resonant = -np.arange(
                np.ceil(q_min * self.n),
                np.floor(q_max * self.n) + 1,  # +1 to include end point
                dtype=int
        )
        q_resonant = -m_resonant / self.n
        s_resonant = np.zeros_like(m_resonant, dtype=float)
        for k, m in enumerate(m_resonant):
            def q_resonant_interp(x):
                return q_interp(x) - q_resonant[k]
            s_resonant[k] = optimize.brentq(
                    q_resonant_interp, np.amin(self.s_mhd), np.amax(self.s_mhd)
            )
        print(s_resonant)
        print(np.array([0.608, 0.760, 0.823, 0.861, 0.891, 0.918])**2)

        data_range = (np.shape(data)[1] - 2) // 2
        abs_data = np.hypot(
                data[:, 2:(2 + data_range)], data[:, (2 + data_range):]
        )
        if self.reffile is not None:
            ref_range = (np.shape(ref)[1] - 2) // 2
            abs_ref = np.hypot(
                    ref[:, 2:(2 + ref_range)], ref[:, (2 + ref_range):]
            )

            if data_range != ref_range:
                print('Different m_max for vacuum and full perturbation field')
                return
        m_max = (data_range - 1) // 2
        offset = m_max

        if self.r_s_interp is not None:
            s = self.r_s_interp(s)
            xlabel = r'$r$'
            s_resonant = self.r_s_interp(s_resonant)

        # plot non-symmetric modes
        fmt = os.path.join(self.datadir,
                           os.path.splitext(self.datafile)[0] + '_{}.pdf')
        horz_plot = 2
        vert_plot = 1
        for m in range(1, m_max + 1):
            if self.reffile is not None:
                yrang = [0, max(
                        np.amax(abs_data[:, offset - m]),
                        np.amax(abs_data[:, offset + m]),
                        np.amax(abs_ref[:, offset - m]),
                        np.amax(abs_ref[:, offset + m])
                )]
            else:
                yrang = [0, max(
                        np.amax(abs_data[:, offset - m]),
                        np.amax(abs_data[:, offset + m])
                )]
            plt.figure(figsize=(6.6, 3.3))
            ax = plt.subplot(vert_plot, horz_plot, 1)
            if self.reffile is not None:
                plt.plot(s, abs_ref[:, offset - m], 'r--',
                         label='vacuum perturbation')
            plt.plot(s, abs_data[:, offset - m], label='full perturbation')
            ax.legend()
            ax.ticklabel_format(style='sci', scilimits=(-3, 4))
            plt.ylim(yrang)
            # q is negative in result_spectrum.f90
            plt.title('$m = {}$'.format(-m))
            plt.ylabel(self.label)
            plt.xlabel(xlabel)
            index = [i for (i, val) in enumerate(m_resonant) if abs(val) == m]
            if len(index) == 1:
                ax.axvline(s_resonant[index], color='b', alpha=0.5)
            ax = plt.subplot(vert_plot, horz_plot, 2)
            if self.reffile is not None:
                plt.plot(s, abs_ref[:, offset + m], 'r--',
                         label='vacuum perturbation')
            plt.plot(s, abs_data[:, offset + m], label='full perturbation')
            ax.legend()
            ax.ticklabel_format(style='sci', scilimits=(-3, 4))
            plt.ylim(yrang)
            # q is negative in result_spectrum.f90
            plt.title('$m = {}$'.format(m))
            plt.xlabel(xlabel)
            plt.ylabel(self.label)
            index = [i for (i, val) in enumerate(m_resonant) if abs(val) == m]
            if len(index) == 1:
                ax.axvline(s_resonant[index], color='b', alpha=0.5)
            plt.tight_layout()
            plt.savefig(fmt.format(m))
            plt.close()
        # plot symmetric mode and safety factor
        plt.figure(figsize=(6.6, 3.3))
        ax = plt.subplot(vert_plot, horz_plot, 1)
        if self.reffile is not None:
            plt.plot(s, abs_ref[:, offset], 'r--',
                     label='vacuum perturbation')
        plt.plot(s, abs_data[:, offset], label='full perturbation')
        ax.legend()
        ax.ticklabel_format(style='sci', scilimits=(-3, 4))
        plt.title('$m = 0$')
        plt.xlabel(r'$s$')
        plt.ylabel(self.label)
        ax = plt.subplot(vert_plot, horz_plot, 2)
        # q is negative in result_spectrum.f90
        plt.plot(s, np.abs(q), label='kinetic')
        if self.r_s_interp is not None:
            plt.plot(self.r_s_interp(self.s_mhd), self.q_mhd, 'r--',
                     label='MHD')
        else:
            plt.plot(self.s_mhd, self.q_mhd, 'r--', label='MHD')
        ax.legend()
        plt.xlabel(xlabel)
        plt.ylabel(r'$q$')
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
        fluxvar = np.loadtxt(os.path.join(self.datadir,
                                          self.config['fluxvar_file']))
        self.rho = fluxvar[:, 0]
        self.psi = fluxvar[:, 1]
        self.q = fluxvar[:, 2]
        self.dens = fluxvar[:, 3]
        self.temp = fluxvar[:, 4]
        self.pres0 = fluxvar[:, 5]
        self.s = (self.psi - self.psi[0]) / (self.psi[-1] - self.psi[0])

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
            contents = np.loadtxt(os.path.join(self.datadir, datafile))
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
                self.rho, r'$\rho$ / cm', self.psi, r'$\psi$ / Mx',
                'disc poloidal flux', 'plot_psi.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.rho, r'$\rho$ / cm', self.q, r'$q$',
                'safety factor', 'plot_q.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.rho, r'$\rho$ / cm', self.dens,
                r'$n$ / cm\textsuperscript{-3}',
                'particle density', 'plot_dens.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.rho, r'$\rho$ / cm', self.temp, r'$T$ / eV',
                'temperature', 'plot_temp.pdf'
        ))
        self.plots.append(magdif_1d_cutplot(
                self.rho, r'$\rho$ / cm', self.pres0,
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

        if (os.path.isfile(os.path.join(self.datadir, 'Bmn_psi.dat')) and
                os.path.isfile(os.path.join(self.datadir, 'Bmn_vac_psi.dat'))):
            self.plots.append(magdif_poloidal_modes(
                    self.config['n'], self.s, self.q,
                    self.datadir, 'Bmn_psi.dat',
                    r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
                    'Bmn_vac_psi.dat'
            ))
        if os.path.isfile(os.path.join(self.datadir, 'currmn_000_theta.dat')):
            self.plots.append(magdif_poloidal_modes(
                    self.config['n'], self.s, self.q,
                    self.datadir, 'currmn_000_theta.dat',
                    r'$\left\vert J_{mn \theta}^{(0)} \right\vert$'
                    + r' / statA cm\textsuperscript{-1}'
            ))
        r_s_interp = interpolate.interp1d(self.s, self.rho, kind='cubic')
        if os.path.isfile(os.path.join(self.datadir, 'Bpmn_r.dat')):
            self.plots.append(magdif_poloidal_modes(
                    self.config['n'] * self.config['kilca_scale_factor'],
                    self.s, self.q, self.datadir, 'Bpmn_r.dat',
                    r'$\left\vert B_{\mathrm{p}mnr} \right\vert$ / G',
                    'Bpmn_vac_r.dat', r_s_interp
            ))

    def dump_plots(self):
        for p in self.plots:
            p.dump_plot()

    def wrapper(self, index):
        self.plots[index].dump_plot()

    def dump_plots_parallel(self):
        with Pool(max(1, os.cpu_count() - 1)) as p:
            p.map(self.wrapper, range(len(self.plots)))


if __name__ == '__main__':
    testcase = magdif(sys.argv[1], sys.argv[2], sys.argv[3])
    testcase.read_configfile()
    testcase.read_fluxvar()
    testcase.load_mesh()
    testcase.generate_default_plots()
    testcase.dump_plots_parallel()
