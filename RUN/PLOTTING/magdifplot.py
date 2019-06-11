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
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from multiprocessing import Pool

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.use('Qt5Agg')
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
    def __init__(self, node, tri, data, filename):
        self.node = node
        self.tri = tri
        self.data = data
        #self.title = title
        self.filename = filename

    def dump_plot(self):
        print('plotting ', self.filename)
        plt.figure(figsize = (5.2, 4.8))
        plt.tripcolor(self.node[:, 0], self.node[:, 1], self.tri - 1,
            self.data, cmap = 'RdBu')
        plt.xlim([min(self.node[:,0]) * 0.9, max(self.node[:,0]) * 1.1])
        plt.axis('equal')
        plt.colorbar()
        plt.clim([-max(abs(self.data)), max(abs(self.data))])
        plt.xlabel(r'$R$ / cm')
        plt.ylabel(r'$Z$ / cm')
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
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.savefig(self.filename)
        plt.close()

class magdif_poloidal_modes:
    def __init__(self, datadir, datafile, reffile):
        self.datadir = datadir
        self.datafile = datafile
        self.reffile = reffile

    def dump_plot(self):
        print('plotting poloidal modes from', self.datafile)
        data = np.loadtxt(os.path.join(self.datadir, self.datafile))
        ref = np.loadtxt(os.path.join(self.datadir, self.reffile))
        
        psi_n = (np.max(data[:,0]) - data[:,0]) / np.max(data[:,0])
        
        M_data = (np.shape(data)[1] - 2) // 2
        abs_data = np.hypot(data[:, 2:(2 + M_data)], data[:, (2 + M_data):])
        M_ref = (np.shape(ref)[1] - 2) // 2
        abs_ref = np.hypot(ref[:, 2:(2 + M_ref)], ref[:, (2 + M_ref):])
        
        horz_plot = 3
        vert_plot = 2
        k2_max = horz_plot * vert_plot
        k1_max = M_data // 2 // k2_max
        for k1 in range(k1_max):
            yrang = [0, np.amax([
                np.amax(abs_data[:, (k1 * k2_max):((k1 + 1) * k2_max - 1)]),
                np.amax(abs_ref[:, (k1 * k2_max):((k1 + 1) * k2_max - 1)])])]
            plt.figure(figsize = (7.2, 3.6))
            for k2 in range(k2_max):
                m = k1 * k2_max + k2
                plt.subplot(vert_plot, horz_plot, k2 + 1)
                plt.plot(psi_n, abs_data[:, m])
                plt.plot(psi_n, abs_ref[:, m], 'r--')
                plt.title('$m = {}$'.format(m))
                plt.plot(np.array([0.608, 0.608])**2, yrang, 'b' ,alpha=.5)#1.5
                plt.plot(np.array([0.760, 0.760])**2, yrang, 'b' ,alpha=.5)#2
                plt.plot(np.array([0.823, 0.823])**2, yrang, 'b' ,alpha=.5)#2.5
                plt.plot(np.array([0.861, 0.861])**2, yrang, 'b' ,alpha=.5)#3
                plt.plot(np.array([0.891, 0.891])**2, yrang, 'b' ,alpha=.5)#3.5
                plt.plot(np.array([0.918, 0.918])**2, yrang, 'b' ,alpha=.5)#4
                plt.xlabel(r'$\psi / \psi_{\max}$')
                plt.ylabel(r'$\left\vert B_{mn}^{\psi} \right\vert$')
                plt.ylim(yrang)
                plt.tight_layout()
            plt.savefig(os.path.splitext(self.datafile)[0] + 
                '_{}.pdf'.format(k1+1))
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
        nml = p.read(os.path.join(self.datadir, self.configfile))
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
        self.s = self.rho / np.amax(self.rho)
        self.s_2 = np.square(self.s)

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
        vector_infix = ['_ReR', '_ImR', '_ReZ', '_ImZ', '_Rephi', '_Imphi',
                        '_Reproj', '_Improj']
        scalar_infix = ['_Re', '_Im']

        for datafile in Bn_datafiles + currn_datafiles:
            self.generate_2d_triplots(datafile, 2, vector_infix,
                self.__class__.decorate_filename_vectorplot)
        for datafile in presn_datafiles:
            self.generate_2d_triplots(datafile, 0, scalar_infix,
                self.__class__.decorate_filename_presnplot)
        
        self.plots.append(magdif_poloidal_modes(self.datadir,
            'Bmn_psi.dat', 'Bmn_vac_psi.dat'))

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
