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
#==============================================================================
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
#==============================================================================

class magdif_2d_triplot:
    scifmt = matplotlib.ticker.ScalarFormatter()
    scifmt.set_powerlimits((-3, 4))

    def __init__(self, node, tri, data, filename):
        self.node = node
        self.tri = tri
        self.data = data
        #self.title = title
        self.filename = filename

    def dump_plot(self):
        print('plotting ', self.filename)
        plt.figure(figsize = (3.6, 4.8))
        plt.tripcolor(self.node[:, 0] * 0.01, self.node[:, 1] * 0.01, self.tri - 1,
            self.data, cmap = 'RdBu_r')
        plt.gca().set_aspect('equal')
        plt.colorbar(format = self.__class__.scifmt)
        plt.clim([-max(abs(self.data)), max(abs(self.data))])
        plt.xlabel(r'$R$ / m')
        plt.ylabel(r'$Z$ / m')
        #plt.title(self.title)
        plt.savefig(self.filename)
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
        plt.figure()
        plt.plot(self.x, self.y, '-k')
        plt.ticklabel_format(style = 'sci', scilimits = (-3, 4))
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.savefig(self.filename)
        plt.close()

class magdif_conv_plot:
    def __init__(self, datadir, conv_sum_file, conv_diff_file):
        self.datadir = datadir
        self.conv_sum_file = conv_sum_file
        self.conv_diff_file = conv_diff_file

    def dump_plot(self):
        conv_sum = np.loadtxt(os.path.join(self.datadir, self.conv_sum_file))
        conv_diff = np.loadtxt(os.path.join(self.datadir, self.conv_diff_file))
        niter = len(conv_diff)
        kiter = np.arange(0, niter + 1)
        plt.figure(figsize = (8, 4.5))
        plt.semilogy(kiter[1:], np.sqrt(conv_sum[0]) * 1e-6 / kiter[1:]**2, 'r-',
            label = r'$\frac{1}{k^{2}} \Vert \delta \mathbf{B}_{\mathrm{v}} \Vert_{2}$')
        plt.semilogy(kiter[1:], np.sqrt(conv_diff) * 1e-6, 'xk',
            label = r'$\Vert \delta \mathbf{B}^{(k)} \Vert_{2}$')
        plt.gca().legend(loc = 'lower left', fontsize = 'large')
        plt.xticks(kiter[1:])
        plt.xlabel('iteration step $k$')
        plt.ylabel(r'$\Vert \delta \mathbf{B} \Vert_{2}$ / T m')
        plt.title('magnitude of terms in series expansion of perturbation field')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.datadir, 'convergence.pdf'))
        plt.close()

class magdif_poloidal_modes:
    def __init__(self, n, s, q, datadir, datafile, reffile):
        self.n = n
        self.s_mhd = s
        self.q_mhd = q
        self.datadir = datadir
        self.datafile = datafile
        self.reffile = reffile

    def dump_plot(self):
        print('plotting poloidal modes from', self.datafile)
        data = np.loadtxt(os.path.join(self.datadir, self.datafile))
        ref = np.loadtxt(os.path.join(self.datadir, self.reffile))
        
        # normalize psi
        s = (data[:,0] - data[0,0]) / (data[-1,0] - data[0,0])
        q = data[:,1]
        
        q_min = np.amin(self.q_mhd)
        q_max = np.amax(self.q_mhd)
        q_interp = interpolate.interp1d(self.s_mhd, self.q_mhd, kind = 'cubic')
        
        m_resonant = -np.arange(np.amax([np.ceil(q_min * self.n), self.n + 1]),
            np.floor(q_max * self.n) + 1, dtype = int)  # +1 to include end point
        q_resonant = -m_resonant / self.n
        s_resonant = np.zeros_like(m_resonant, dtype = float)
        for k, m in enumerate(m_resonant):
            def q_resonant_interp(x):
                return q_interp(x) - q_resonant[k]
            s_resonant[k] = optimize.brentq(q_resonant_interp,
                np.amin(self.s_mhd), np.amax(self.s_mhd))
        print(s_resonant)
        print(np.array([0.608, 0.760, 0.823, 0.861, 0.891, 0.918])**2)
        
        data_range = (np.shape(data)[1] - 2) // 2
        abs_data = np.hypot(data[:, 2:(2 + data_range)], data[:, (2 + data_range):])
        ref_range = (np.shape(ref)[1] - 2) // 2
        abs_ref = np.hypot(ref[:, 2:(2 + ref_range)], ref[:, (2 + ref_range):])
        
        if data_range != ref_range:
            raise RuntimeError('Different m_max for vacuum and full perturbation field')
        m_max = (data_range - 1) // 2
        offset = m_max
        
        horz_plot = 2
        vert_plot = 1
        for m in range(1, m_max + 1):
            yrang = [0, max(np.amax(abs_data[:, offset - m]), np.amax(abs_data[:, offset + m]),
                np.amax(abs_ref[:, offset - m]), np.amax(abs_ref[:, offset + m]))]
            plt.figure(figsize = (9.6, 4.8))
            ax = plt.subplot(vert_plot, horz_plot, 1)
            plt.plot(s, abs_ref[:, offset - m], 'r--')
            plt.plot(s, abs_data[:, offset - m])
            ax.ticklabel_format(style = 'sci', scilimits = (-3, 4))
            plt.ylim(yrang)
            plt.title('$m = {}$'.format(-m))
            plt.ylabel(r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx')
            plt.xlabel(r'$s$')
            index = [i for (i, val) in enumerate(m_resonant) if abs(val) == m]
            if len(index) == 1:
                ax.axvline(s_resonant[index], color = 'b', alpha = 0.5)
            ax = plt.subplot(vert_plot, horz_plot, 2)
            plt.plot(s, abs_ref[:, offset + m], 'r--')
            plt.plot(s, abs_data[:, offset + m])
            ax.ticklabel_format(style = 'sci', scilimits = (-3, 4))
            plt.ylim(yrang)
            plt.title('$m = {}$'.format(m))
            plt.xlabel(r'$s$')
            plt.ylabel(r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx')
            index = [i for (i, val) in enumerate(m_resonant) if abs(val) == m]
            if len(index) == 1:
                ax.axvline(s_resonant[index], color = 'b', alpha = 0.5)
            plt.tight_layout()
            plt.savefig(os.path.splitext(self.datafile)[0] + '_{}.pdf'.format(m))
            plt.close()
        plt.figure(figsize = (9.6, 4.8))
        ax = plt.subplot(vert_plot, horz_plot, 1)
        plt.plot(s, abs_ref[:, offset], 'r--')
        plt.plot(s, abs_data[:, offset])
        ax.ticklabel_format(style = 'sci', scilimits = (-3, 4))
        plt.title('$m = 0$')
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\left\vert \sqrt{g} B_{mn}^{\psi} \right\vert$ / Mx')
        ax = plt.subplot(vert_plot, horz_plot, 2)
        plt.plot(s, np.abs(q), label = 'kinetic')
        plt.plot(self.s_mhd, self.q_mhd, 'r--', label = 'MHD')
        ax.legend()
        plt.xlabel(r'$s$')
        plt.ylabel(r'$q$')
        plt.tight_layout()
        plt.savefig(os.path.splitext(self.datafile)[0] + '_0.pdf')
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
        fluxvar = np.loadtxt(os.path.join(self.datadir, self.config['fluxvar_file']))
        self.rho   = fluxvar[:,0]
        self.psi   = fluxvar[:,1]
        self.q     = fluxvar[:,2]
        self.dens  = fluxvar[:,3]
        self.temp  = fluxvar[:,4]
        self.pres0 = fluxvar[:,5]
        self.s = (self.psi - self.psi[0]) / (self.psi[-1] - self.psi[0])

    def load_mesh(self):
        with open(self.meshfile, 'r') as f:
            data = f.readlines()
        print('meshfile header: ', data[0])
        [NN, NT, NE] = np.loadtxt(StringIO(data[0]), dtype = int)
        node = np.loadtxt(StringIO(''.join(data[1:NN+1])), dtype = float)
        self.node = node[:, 0:2]
        tri = np.loadtxt(StringIO(''.join(data[NN+1:NN+NT+1])), dtype = int)
        self.tri = tri[:, 0:3]
        edge = np.loadtxt(StringIO(''.join(data[NN+NT+1:])), dtype = int)
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

    def generate_2d_triplots(self, datafile, start_column, infix_list, filename_decorator):
        print('reading contents of ', datafile)
        contents = np.loadtxt(os.path.join(self.datadir, datafile))
        for column, infix in enumerate(infix_list, start_column):
            self.plots.append(magdif_2d_triplot(
                node = self.node, tri = self.tri,
                data = contents[:, column],
                filename = filename_decorator(datafile, infix)
                ))

    def generate_default_plots(self):
        self.plots.append(magdif_1d_cutplot(self.rho, r'$\rho$ / cm', self.psi,
                r'$\psi$ / Mx', 'disc poloidal flux', 'plot_psi.pdf'
            ))
        self.plots.append(magdif_1d_cutplot(self.rho, r'$\rho$ / cm', self.q,
                r'$q$', 'safety factor', 'plot_q.pdf'
            ))
        self.plots.append(magdif_1d_cutplot(self.rho, r'$\rho$ / cm', self.dens,
                r'$n$ / cm\textsuperscript{-3}', 'particle density', 'plot_dens.pdf'
            ))
        self.plots.append(magdif_1d_cutplot(self.rho, r'$\rho$ / cm', self.temp,
                r'$T$ / eV', 'temperature', 'plot_temp.pdf'
            ))
        self.plots.append(magdif_1d_cutplot(self.rho, r'$\rho$ / cm', self.pres0,
                r'$p_{0}$ / dyn cm\textsuperscript{-2}', 'pressure', 'plot_pres0.pdf'
            ))
        self.generate_2d_triplots('j0phi.dat', 3, [''],
            self.__class__.decorate_default)
        
        Bn_datafiles = ['plot_Bn.dat', 'plot_Bn_001.dat', 'plot_Bn_vac.dat']
        currn_datafiles = ['plot_currn.dat', 'plot_currn_001.dat']
        presn_datafiles = ['presn.dat', 'presn_001.dat']
        vector_infix = ['_R_Re', '_R_Im', '_Z_Re', '_Z_Im', '_phi_Re', '_phi_Im',
                        '_contradenspsi_Re', '_contradenspsi_Im', '_cotheta_Re', '_cotheta_Im']
        scalar_infix = ['_Re', '_Im']

        for datafile in Bn_datafiles + currn_datafiles:
            self.generate_2d_triplots(datafile, 2, vector_infix,
                self.__class__.decorate_filename_vectorplot)
        for datafile in presn_datafiles:
            self.generate_2d_triplots(datafile, 0, scalar_infix,
                self.__class__.decorate_filename_presnplot)
        
        self.plots.append(magdif_conv_plot(self.datadir, 'conv_sum.dat', 'conv_diff.dat'))
        self.plots.append(magdif_poloidal_modes(self.config['n'], self.s, self.q,
            self.datadir, 'Bmn_psi.dat', 'Bmn_vac_psi.dat'))

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

#==============================================================================
#     testcase = magdif(
#         '/temp/lainer_p/NEO-EQ/coarse_grid/test_res_precon_lowdens',
#         'test_res_precon_lowdens.in',
#         '../PRELOAD/inputformaxwell.msh'
#     )
#     testcase.load_mesh()
#     vacuum = np.loadtxt('/temp/lainer_p/NEO-EQ/coarse_grid/test_res_precon_highdens/plot_Bn_vac.dat')
#     precon = np.loadtxt('/temp/lainer_p/NEO-EQ/coarse_grid/test_res_precon_lowdens/plot_Bn.dat')
#     direct = np.loadtxt('/temp/lainer_p/NEO-EQ/coarse_grid/test_res_direct_lowdens/plot_Bn.dat')
#     highdens = np.loadtxt('/temp/lainer_p/NEO-EQ/coarse_grid/test_res_precon_highdens/plot_Bn.dat')
#     vacuum = vacuum[:,4] + 1j * vacuum[:,5]
#     precon = precon[:,4] + 1j * precon[:,5]
#     direct = direct[:,4] + 1j * direct[:,5]
#     highdens = highdens[:,4] + 1j * highdens[:,5]
#     testcase.plots.append(magdif_2d_triplot(
#         node = testcase.node, tri = testcase.tri,
#         data = 1e-4 * np.real(vacuum),
#         title = r'$\mathrm{Re} \: \delta B_{\mathrm{v}}^{Z}$ / T',
#         filename = '/temp/lainer_p/NEO-EQ/coarse_grid/Bn_vac_Z_Re.pdf'
#     ))
#     testcase.plots.append(magdif_2d_triplot(
#         node = testcase.node, tri = testcase.tri,
#         data = 1e-4 * np.real(precon),
#         title = r'$\mathrm{Re} \: \delta B^{Z}$ / T',
#         filename = '/temp/lainer_p/NEO-EQ/coarse_grid/Bn_Z_Re.pdf'
#     ))
#     testcase.plots.append(magdif_2d_triplot(
#         node = testcase.node, tri = testcase.tri,
#         data = 1e-4 * np.real(precon - direct),
#         title = r'$\mathrm{Re} \: \Delta \delta B^{Z}$ / T',
#         filename = '/temp/lainer_p/NEO-EQ/coarse_grid/Bn_diff_Z_Re.pdf'
#     ))
#     testcase.plots.append(magdif_2d_triplot(
#         node = testcase.node, tri = testcase.tri,
#         data = 1e-4 * np.real(highdens),
#         title = r'$\mathrm{Re} \: \delta B^{Z}$ / T',
#         filename = '/temp/lainer_p/NEO-EQ/coarse_grid/Bn_highdens_Z_Re.pdf'
#     ))
#     testcase.dump_plots()
#==============================================================================
