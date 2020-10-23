#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:39:22 2019

@author: lainer_p
"""

from sys import argv
from os import path, cpu_count
from enum import Enum
from netCDF4 import Dataset
from h5py import get_config as h5py_hack
import h5pickle as h5py
import f90nml.parser
import matplotlib
import matplotlib.pyplot as plt
import colorcet
import numpy as np
from scipy import interpolate
from multiprocessing import Pool

# complex values are stored as compound types in libneo/hdf5tools
h5py_hack().complex_names = ('real', 'imag')
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
latex_preamble = path.join(path.dirname(path.realpath(__file__)),
                           'magdifplot.tex')
matplotlib.rcParams['text.latex.preamble'] = fr"\input{{{latex_preamble}}}"
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

    def __init__(self, mesh, data, label, filename, title=None, clim_scale=None):
        self.mesh = mesh
        self.data = data
        self.label = label
        self.title = title
        self.filename = filename
        if clim_scale is not None:
            self.clim_scale = clim_scale
        else:
            self.clim_scale = 1.0

    def dump_plot(self):
        print(f"plotting {self.filename}")
        plt.figure(figsize=(3.3, 4.4))
        plt.tripcolor(
                self.mesh['node_R'], self.mesh['node_Z'],
                self.mesh['tri_node'][()] - 1, self.data,
                cmap=colorcet.cm.coolwarm
        )
        plt.gca().set_aspect('equal')
        cbar = plt.colorbar(format=self.__class__.scifmt)
        cbar.set_label(self.label, rotation=90)
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
        print(f"plotting {self.filename}")
        plt.figure(figsize=(6.6, 3.6))
        plt.plot(self.x, self.y, '-k')
        plt.ticklabel_format(style='sci', scilimits=(-3, 4))
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.savefig(self.filename)
        plt.close()


class magdif_conv_plot:
    def __init__(self, config, data, filename, xlim=None, ylim=None, title=None):
        self.config = config
        self.data = data
        self.filename = filename
        self.xlim = xlim
        self.ylim = ylim
        self.title = title

    def dump_plot(self):
        sup_eigval = self.config['ritz_threshold']
        L2int_Bnvac = self.data['/iter/L2int_Bnvac']
        L2int_Bn_diff = self.data['/iter/L2int_Bn_diff']
        kiter = np.arange(0, len(L2int_Bn_diff))
        plt.figure(figsize=(3.2, 3.2))
        plt.semilogy(
                kiter, L2int_Bnvac * sup_eigval ** kiter, 'r-',
                label=r'$\vert \lambda_{\mathrm{max}} \vert^{k} \Vert'
                + r' \mathbf{B}_{n}^{(0)} \Vert_{2}$'
        )
        plt.semilogy(
                kiter, L2int_Bn_diff, 'xk',
                label=r'$\Vert \delta \mathbf{B}_{n}^{(k)} \Vert_{2}$'
        )
        if self.xlim is not None:
            plt.xlim(self.xlim)
        if self.ylim is not None:
            plt.ylim(self.ylim)
        plt.gca().legend(loc='upper right')
        plt.xticks(kiter)
        plt.xlabel('iteration step $k$')
        plt.ylabel(r'$\Vert \mathbf{B} \Vert_{2}$ / G cm')
        if self.title is not None:
            plt.title(self.title)
        else:
            plt.title('estimation of convergence')

        plt.tight_layout()
        plt.savefig(self.filename)
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
        print(f"reading poloidal modes from {self.datafile}")
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
        self.arg = np.arctan2(
                data[:, (2 + self.range):], data[:, 2:(2 + self.range)]
        )
        self.ymax = np.amax(self.abs, axis=0)
        self.ymax = np.fmax(self.ymax[self.offset:-1],
                            self.ymax[self.offset:0:-1])


class magdif_GPEC_bnormal:
    def __init__(self, datadir, datafile, variable, label):
        self.datadir = datadir
        self.datafile = datafile
        self.variable = variable
        self.label = label

    def process(self):
        print(f"reading poloidal modes from {self.datafile}")
        rootgrp = Dataset(self.datafile, 'r')
        self.rho = np.array(rootgrp.variables['psi_n'])
        self.q = np.array(rootgrp.variables['q'])
        m = np.array(rootgrp.variables['m_out'])
        self.range = np.size(m)
        self.offset = -m[0]
        self.m_max = np.min([-m[0], m[-1]])
        data = np.array(rootgrp.variables[self.variable])
        self.abs = 1e8 * np.transpose(np.hypot(data[0, :, :], data[1, :, :]))
        self.ymax = np.amax(self.abs, axis=0)
        self.ymax = np.fmax(self.ymax[self.offset:self.offset + self.m_max],
                            self.ymax[self.offset:self.offset - self.m_max:-1])
        rootgrp.close()


class magdif_poloidal_plots:
    def __init__(self, config, data, xlabel, ylabel, poldata, refdata=None):
        self.config = config
        self.data = data
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.poldata = poldata
        self.refdata = refdata

        psi_norm = self.data['/cache/fs/psi']
        psi_norm = (psi_norm - psi_norm[0]) / (psi_norm[-1] - psi_norm[0])
        rad = self.data['/cache/fs/rad']
        self.psi2rad = interpolate.interp1d(psi_norm, rad, kind='cubic')
        self.rad2psi = interpolate.interp1d(rad, psi_norm, kind='cubic')

    def interp_rho(self, data):
        if self.xlabel is fslabel.psi_norm and data.rad_coord is fslabel.r:
            return self.rad2psi(data.rho)
        elif self.xlabel is fslabel.r and data.rad_coord is fslabel.psi_norm:
            return self.psi2rad(data.rho)
        else:
            return data.rho

    def dump_plot(self):
        sgn_m_res = -np.sign(self.data['/cache/fs/q'][-1])
        m_res = range(self.data['/mesh/m_res_min'][()],
                      self.data['/mesh/m_res_max'][()] + 1) * sgn_m_res
        if self.xlabel is fslabel.psi_norm:
            rho = self.data['/cache/fs/psi']
            resonance = dict(zip(m_res, self.data['/mesh/psi_res'] /
                                 (rho[-1] - rho[0])))
            # normalize psi
            rho = (rho - rho[0]) / (rho[-1] - rho[0])
        else:
            rho = self.data['/cache/fs/rad']
            resonance = dict(zip(m_res, self.data['/mesh/rad_norm_res'] *
                                 rho[-1]))
        xlabel = self.xlabel.value

        # plot non-symmetric modes
        fmt = path.join(self.poldata.datadir, path.basename(
                        path.splitext(self.poldata.datafile)[0] + '_{}.pdf'))
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        for m_abs in range(1, self.poldata.m_max):
            if self.refdata is not None and m_abs <= self.refdata.m_max:
                ymax = max(self.poldata.ymax[m_abs], self.refdata.ymax[m_abs])
            else:
                ymax = self.poldata.ymax[m_abs]
            plt.figure(figsize=two_squares)
            for k in range(horz_plot):
                ax = plt.subplot(vert_plot, horz_plot, k+1)
                m = (2 * k - 1) * m_abs
                if m in resonance:
                    ax.axvline(resonance[m], color='b', alpha=0.5)
                if self.refdata is not None and m_abs <= self.refdata.m_max:
                    ax.plot(self.interp_rho(self.refdata),
                            self.refdata.abs[:, self.refdata.offset + m],
                            'r--', label=self.refdata.label)
                ax.plot(self.interp_rho(self.poldata),
                        self.poldata.abs[:, self.poldata.offset + m],
                        label=self.poldata.label)
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
        if self.refdata is not None:
            ax.plot(self.interp_rho(self.refdata),
                    self.refdata.abs[:, self.refdata.offset],
                    'r--', label=self.refdata.label)
        ax.plot(self.interp_rho(self.poldata),
                self.poldata.abs[:, self.poldata.offset],
                label=self.poldata.label)
        ax.legend()
        ax.ticklabel_format(style='sci', scilimits=(-3, 4))
        ax.set_title('$m = 0$')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(self.ylabel)
        ax = plt.subplot(vert_plot, horz_plot, 2)
        for res in resonance.values():
            ax.axvline(res, color='b', alpha=0.5)
        ax.plot(rho, self.data['/cache/fs/q'], 'k-')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$q$')
        plt.tight_layout()
        plt.savefig(fmt.format(0))
        plt.close()


RT0_comp = {'/comp_R': {'file': '_R', 'math': lambda vec: fr"{vec}^{{R}}"},
            '/comp_Z': {'file': '_Z', 'math': lambda vec: fr"{vec}^{{Z}}"},
            '/RT0_comp_phi': {'file': '_phi', 'math': lambda vec: fr"R {vec}^{{\varphi}}"},
            '/comp_psi_contravar_dens': {'file': '_contradenspsi',
                                         'math': lambda vec: fr"\sqrt{{g}} {vec}^{{\psi}}"},
            '/comp_theta_covar': {'file': '_cotheta',
                                  'math': lambda vec: fr"\subscript{{{vec}}}{{\vartheta}}"}}


class magdif:

    def __init__(self, datadir, configfile='magdif.inp', datafile='magdif.h5'):
        self.plots = []
        self.datadir = datadir
        self.configfile = configfile
        self.datafile = datafile

    def read_configfile(self):
        print(f"reading configuration from {self.configfile}")
        p = f90nml.parser.Parser()
        nml = p.read(path.join(self.datadir, self.configfile))
        self.config = {**nml['scalars']['config'], **nml['arrays']}

    def read_datafile(self):
        print(f"reading contents of {self.datafile}")
        self.data = h5py.File(path.join(self.datadir, self.datafile), 'r')

    def generate_RT0_triplots(self, grp, label, filename):
        for dataset, decorator in RT0_comp.items():
            nameparts = path.splitext(filename)
            self.plots.append(magdif_2d_triplot(
                mesh=self.data['/mesh'],
                data=self.data[grp + dataset][()].imag,
                label=r'$\mathrm{Re} \, ' + decorator['math'](label) + r'$',
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_Re' + nameparts[1])))
            self.plots.append(magdif_2d_triplot(
                mesh=self.data['/mesh'],
                data=self.data[grp + dataset][()].imag,
                label=r'$\mathrm{Im} \, ' + decorator['math'](label) + r'$',
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_Im' + nameparts[1])))

    def generate_L1_triplots(self, grp, label, filename):
        nameparts = path.splitext(filename)
        self.plots.append(magdif_2d_triplot(
            mesh=self.data['/mesh'],
            data=self.data[grp + '/L1_DOF'][()].real,
            label=r'$\mathrm{Re} \, ' + label + r'$',
            filename=path.join(self.datadir, nameparts[0] +
                               '_Re' + nameparts[1])))
        self.plots.append(magdif_2d_triplot(
            mesh=self.data['/mesh'],
            data=self.data[grp + '/L1_DOF'][()].imag,
            label=r'$\mathrm{Im} \, ' + label + r'$',
            filename=path.join(self.datadir, nameparts[0] +
                               '_Im' + nameparts[1])))

    def generate_default_plots(self):
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'], r'$r$ / cm',
                self.data['/cache/fs/psi'], r'$\psi$ / Mx',
                'disc poloidal flux', path.join(self.datadir, 'plot_psi.pdf')
        ))
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'], r'$r$ / cm',
                self.data['/cache/fs/q'], r'$q$',
                'safety factor', path.join(self.datadir, 'plot_q.pdf')
        ))
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'], r'$r$ / cm',
                self.data['/cache/fs/p'], r'$p_{0}$ / dyn cm\textsuperscript{-2}',
                'pressure', path.join(self.datadir, 'plot_p0.pdf')
        ))
        # TODO: j0phi edge plot
        self.generate_RT0_triplots('/iter/Bn', r'B_{n}', 'plot_Bn.pdf')
        self.generate_RT0_triplots('/iter/Bn_000', r'B_{n}', 'plot_Bn_000.pdf')
        self.generate_RT0_triplots('/iter/Bn_001', r'B_{n}', 'plot_Bn_001.pdf')
        self.generate_RT0_triplots('/iter/jn', r'J_{n}', 'plot_Jn.pdf')
        self.generate_RT0_triplots('/iter/jn_000', r'J_{n}', 'plot_Jn_000.pdf')
        self.generate_RT0_triplots('/iter/jn_001', r'J_{n}', 'plot_Jn_001.pdf')
        self.generate_L1_triplots('/iter/pn', r'p_{n}', 'plot_pn.pdf')
        self.generate_L1_triplots('/iter/pn_000', r'p_{n}', 'plot_pn_000.pdf')
        self.generate_L1_triplots('/iter/pn_001', r'p_{n}', 'plot_pn_001.pdf')

        self.plots.append(magdif_conv_plot(
            self.config, self.data, path.join(self.datadir, 'convergence.pdf'))
        )

        if 'kilca_scale_factor' in self.config:
            kilca_scale_factor = self.config['kilca_scale_factor']
        else:
            kilca_scale_factor = 0
        if kilca_scale_factor == 0:
            pert = magdif_mnDat(self.datadir, 'Bmn_psi.dat', 0,
                                fslabel.psi_norm, 'full perturbation')
            pert.process()
            vac = magdif_mnDat(self.datadir, 'Bmn_vac_psi.dat', 0,
                               fslabel.psi_norm, 'vacuum perturbation')
            vac.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config, self.data, fslabel.psi_norm,
                    r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx',
                    pert, vac
            ))
            pert = magdif_mnDat(self.datadir, 'currmn_000_theta.dat', 0,
                                fslabel.psi_norm, 'initial perturbation')
            pert.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config, self.data, fslabel.psi_norm,
                    r'$\left\vert J_{mn \theta}^{(0)} \right\vert$'
                    + r' / statA cm\textsuperscript{-1}', pert
            ))

        else:
            pert = magdif_mnDat(self.datadir, 'Bmn_r.dat',
                                kilca_scale_factor,
                                fslabel.r, 'full perturbation')
            pert.process()
            vac = magdif_mnDat(self.datadir, 'Bmn_vac_r.dat',
                               kilca_scale_factor,
                               fslabel.r, 'vacuum perturbation')
            vac.process()
            self.plots.append(magdif_poloidal_plots(
                    self.config, self.data, fslabel.r,
                    r'$\left\vert B_{mn}^{r} \right\vert$ / G',
                    pert, vac
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
    testcase.read_datafile()
    testcase.generate_default_plots()
    testcase.dump_plots_parallel()
