#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:39:22 2019

@author: lainer_p
"""

from sys import argv
from os import path, cpu_count
from enum import Enum
import netCDF4
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
#     'pgf.preamble': fr"\input{{{latex_preamble}}}"
# }
# from matplotlib.backends.backend_pgf import FigureCanvasPgf
# matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
# matplotlib.rcParams.update(pgf_config)
# =============================================================================

c_cgs = 2.9979246e+10
cm_to_m = 1.0e-02
c1_statA_to_A = 1.0e+01
c1_statA_per_cm2_to_A_per_m2 = c1_statA_to_A / cm_to_m ** 2
statA_to_A = c1_statA_to_A / c_cgs
statA_per_cm2_to_A_per_m2 = c1_statA_per_cm2_to_A_per_m2 / c_cgs

def scifmt():
    fmt = matplotlib.ticker.ScalarFormatter()
    fmt.set_powerlimits((-3, 4))
    return fmt

class magdif_2d_triplot:
    def __init__(self, triangulation, data, label, filename, title=None,
                 clim_scale=None):
        self.triangulation = triangulation
        self.data = data
        self.label = label
        self.title = title
        self.filename = filename
        if clim_scale is not None:
            self.clim_scale = clim_scale
        else:
            self.clim_scale = (1.0, 1.0)

    def dump_plot(self):
        print(f"plotting {self.filename}")
        plt.figure(figsize=(3.3, 4.4))
        plt.tripcolor(self.triangulation, self.data, cmap=colorcet.cm.coolwarm)
        plt.gca().set_aspect('equal')
        cbar = plt.colorbar(format=scifmt())
        cbar.set_label(self.label, rotation=90)
        plt.clim([-max(abs(self.data)) * self.clim_scale[0],
                  max(abs(self.data)) * self.clim_scale[1]])
        plt.xlabel(r'$R$ / cm')
        plt.ylabel(r'$Z$ / cm')
        if self.title is not None:
            plt.title(self.title)
        plt.savefig(self.filename, dpi=300)
        plt.close()


class magdif_2d_rectplots:
    def __init__(self, R, Z, data, label, filename, title=None, clim_scale=None):
        self.R = R
        self.Z = Z
        self.data = data
        self.label = label
        self.title = title
        self.filename = filename
        if clim_scale is None:
            self.clim_scale = (1.0, 1.0)
        else:
            self.clim_scale = clim_scale

    def dump_plot(self):
        print(f"plotting {self.filename}")
        xlim = (min(R[0] for R in self.R), max(R[-1] for R in self.R))
        ylim = (min(Z[0] for Z in self.Z), max(Z[-1] for Z in self.Z))
        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8.0, 5.0))
        images = []
        for k in range(2):
            images.append(axs[k].imshow(self.data[k], cmap=colorcet.cm.coolwarm, interpolation='gaussian',
                                        extent=[self.R[k][0], self.R[k][-1], self.Z[k][0], self.Z[k][-1]]))
            axs[k].set_aspect('equal')
            axs[k].set_xlabel(r'$R$ / cm')
            axs[k].set_ylabel(r'$Z$ / cm')
            axs[k].set_xlim(xlim)
            axs[k].set_ylim(ylim)
            if self.title is not None:
                axs[k].set_title(self.title[k])
        axs[1].yaxis.set_tick_params(labelleft=True)
        clim = (-max(np.amax(np.abs(image.get_array())) for image in images) * self.clim_scale[0],
                max(np.amax(np.abs(image.get_array())) for image in images) * self.clim_scale[1])
        norm = matplotlib.colors.Normalize(vmin=clim[0], vmax=clim[1])
        for im in images:
            im.set_norm(norm)
        cbar = fig.colorbar(images[0], ax=axs, format=scifmt())
        cbar.set_label(self.label, rotation=90)
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
        print(f"plotting {self.filename}")
        sup_eigval = self.config['ritz_threshold']
        L2int_Bnvac = self.data['/iter/L2int_Bnvac'][()]
        L2int_Bn_diff = self.data['/iter/L2int_Bn_diff'][()]
        kiter = np.arange(0, len(L2int_Bn_diff))
        plt.figure(figsize=(6.6, 3.6))
        plt.semilogy(
                kiter, L2int_Bnvac * sup_eigval ** kiter, 'r-',
                label=r'$\lvert \lambda_{\text{max}} \rvert^{k}'
                + r' \lVert \mathbf{B}_{n}^{(0)} \rVert_{2}$'
        )
        plt.semilogy(
                kiter, L2int_Bn_diff, 'xk',
                label=r'$\lVert \delta \mathbf{B}_{n}^{(k)} \rVert_{2}$'
        )
        if self.xlim is not None:
            plt.xlim(self.xlim)
        if self.ylim is not None:
            plt.ylim(self.ylim)
        plt.gca().legend(loc='upper right')
        plt.xticks(kiter)
        plt.xlabel('iteration step $k$')
        plt.ylabel(r'$\lVert R \mathbf{B}^{\text{pol}} \rVert_{2}$ / \si{\maxwell}')
        if self.title is not None:
            plt.title(self.title)
        else:
            plt.title('estimation of convergence')

        plt.tight_layout()
        plt.savefig(self.filename)
        plt.close()


class fslabel(Enum):
    psi_norm = r'normalized poloidal flux $\hat{\psi}$'
    r = r'minor radius $r$ / \si{\centi\meter}'

class polmodes:
    def __init__(self, label, fmt):
        self.label = label
        self.fmt = fmt
        self.var = {}
        self.rho = {}

    def read_magdif(self, data, rad_coord, var_name='/postprocess/Bmn/coeff_rad'):
        self.type = 'MEPHIT'
        self.rad_coord = rad_coord
        if self.rad_coord == fslabel.psi_norm:
            rho = data['/cache/fs_half/psi'][()]
            # normalize psi
            rho = (rho - rho[0]) / (rho[-1] - rho[0])
        else:
            rho = data['/cache/fs_half/rad'][()]
        self.m_max = (data[var_name].shape[1] - 1) // 2
        for m in range(-self.m_max, self.m_max + 1):
            self.rho[m] = rho
            self.var[m] = np.array(data[var_name][:, m + self.m_max], dtype='D')

    def read_fouriermodes(self, data):
        self.type = 'amn.dat'
        self.rad_coord = fslabel.psi_norm
        rho = data['/debug_fouriermodes/psi_n'][()]
        var_name = '/debug_fouriermodes/comp_psi_contravar_dens'
        self.m_max = (data[var_name].shape[1] - 1) // 2
        for m in range(-self.m_max, self.m_max + 1):
            self.rho[m] = rho
            self.var[m] = np.array(data[var_name][:, m + self.m_max], dtype='D')

    def read_KiLCA(self, datafile, var_name='Br'):
        self.type = 'KiLCA'
        self.rad_coord = fslabel.r
        self.m_max = 0
        data = h5py.File(datafile, 'r')
        sgn_q = int(np.sign(data['/output/background/profiles/q_i'][0, -1]))
        for name, grp in data['/output'].items():
            if 'postprocessor' not in name:
                continue
            m = int(grp['mode'][0, 0]) * sgn_q
            self.m_max = max(self.m_max, abs(m))
            self.rho[m] = np.array(grp['r'][0, :], dtype='d')
            self.var[m] = np.zeros(self.rho[m].shape, dtype='D')
            self.var[m].real = grp[var_name][0, :]
            if grp[var_name].shape[0] == 2:
                self.var[m].imag = grp[var_name][1, :]
        data.close()

    def read_GPEC(self, datafile, var_name='Jbgradpsi'):
        self.type = 'GPEC'
        self.rad_coord = fslabel.psi_norm
        self.m_max = 0
        rootgrp = netCDF4.Dataset(datafile, 'r')
        rho = np.array(rootgrp.variables['psi_n'])
        sgn = int(rootgrp.getncattr('helicity'))
        for k, m_out in enumerate(rootgrp.variables['m_out'][:]):
            m = m_out * sgn
            self.m_max = max(self.m_max, abs(m))
            self.rho[m] = rho
            self.var[m] = np.empty(rho.shape, dtype='D')
            self.var[m].real = rootgrp.variables[var_name][0, k, :]
            self.var[m].imag = rootgrp.variables[var_name][1, k, :]
            # convert weber to maxwell
            self.var[m] *= 1e8
        rootgrp.close()


class magdif_poloidal_plots:
    def __init__(self, datadir, filename, config, data, rad_coord, ylabel,
                 comp, *poldata):
        self.datadir = datadir
        self.filename = filename
        self.config = config
        self.data = data
        self.rad_coord = rad_coord
        self.xlabel = self.rad_coord.value
        self.ylabel = ylabel
        self.comp = comp
        self.poldata = poldata

        psi_norm = self.data['/cache/fs/psi']
        psi_norm = (psi_norm - psi_norm[0]) / (psi_norm[-1] - psi_norm[0])
        rad = self.data['/cache/fs/rad']
        self.psi2rad = interpolate.interp1d(psi_norm, rad, kind='cubic')
        self.rad2psi = interpolate.interp1d(rad, psi_norm, kind='cubic')

    def interp_rho(self, data, m):
        if self.rad_coord is fslabel.psi_norm and data.rad_coord is fslabel.r:
            return self.rad2psi(data.rho[m])
        elif self.rad_coord is fslabel.r and data.rad_coord is fslabel.psi_norm:
            return self.psi2rad(data.rho[m])
        else:
            return data.rho[m]

    def dump_plot(self):
        sgn_m_res = -np.sign(self.data['/cache/fs/q'][-1])
        m_res = range(self.data['/mesh/m_res_min'][()],
                      self.data['/mesh/m_res_max'][()] + 1) * sgn_m_res
        if self.rad_coord is fslabel.psi_norm:
            rho = self.data['/cache/fs/psi']
            resonance = dict(zip(m_res, (self.data['/mesh/psi_res'] - rho[0]) /
                                 (rho[-1] - rho[0])))
            # normalize psi
            rho = (rho - rho[0]) / (rho[-1] - rho[0])
        else:
            rho = self.data['/cache/fs/rad'][()]
            resonance = dict(zip(m_res, self.data['/mesh/rad_norm_res'] *
                                 rho[-1]))

        fmt = path.basename(path.splitext(self.filename)[0] + '_{}' +
                            path.splitext(self.filename)[1])
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        # plot non-symmetric modes
        m_max = min(map(lambda d: d.m_max, self.poldata))
        for m_abs in range(1, m_max):
            filename = fmt.format(m_abs)
            print(f"plotting {filename}")
            fig, axs = plt.subplots(vert_plot, horz_plot, sharex=True,
                                    sharey=True, figsize=two_squares)
            for k in range(horz_plot):
                m = (2 * k - 1) * m_abs
                axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
                if m in resonance:
                    axs[k].axvline(resonance[m], color='b', alpha=0.5,
                                   label='resonance position', lw=0.5)
                for data in self.poldata:
                    if m in data.var:
                        axs[k].plot(self.interp_rho(data, m),
                                    self.comp(data.var[m]),
                                    data.fmt, label=data.label, lw=0.5)
                axs[k].legend(loc='upper left', fontsize='x-small')
                axs[k].ticklabel_format(style='sci', scilimits=(-3, 4))
                axs[k].set_title(('resonant ' if m in resonance else
                                  'non-resonant ') + fr"$m = {m}$")
                axs[k].set_xlabel(self.xlabel)
                axs[k].set_ylabel(self.ylabel)
            axs[1].yaxis.set_tick_params(labelleft=True)
            plt.tight_layout()
            plt.savefig(path.join(self.datadir, filename))
            plt.close()
        # plot symmetric mode and safety factor
        m = 0
        filename = fmt.format(m)
        print(f"plotting {filename}")
        plt.figure(figsize=two_squares)
        ax = plt.subplot(vert_plot, horz_plot, 1)
        for data in self.poldata:
            if m in data.var:
                ax.plot(self.interp_rho(data, m), self.comp(data.var[m]),
                        data.fmt, label=data.label, lw=0.5)
        ax.legend(loc='upper left', fontsize='x-small')
        ax.ticklabel_format(style='sci', scilimits=(-3, 4))
        ax.set_title(f"$m = {m}$")
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax = plt.subplot(vert_plot, horz_plot, 2)
        for res in resonance.values():
            ax.axvline(np.abs(res), color='b', alpha=0.5, lw=0.5)
        q = self.data['/cache/fs/q'][()]
        ax.plot(rho, q, 'k-')
        if (np.any(q > 0.0)):
            ax.set_ylim(bottom=0.0)
        else:
            ax.set_ylim(top=0.0)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(r'$q$')
        plt.tight_layout()
        plt.savefig(path.join(self.datadir, filename))
        plt.close()


class parcurr:
    def __init__(self):
        pass

    def process_magdif(self, datafile_func, m_range, nrad, rres, symfluxcoord=True):
        self.rad = {}
        self.jnpar = {}
        self.rres = rres
        self.part_int = {}
        self.bndry = {}
        if symfluxcoord:
            self.psi = {}
            self.Ichar = {}
            self.Delta = {}
        for m in m_range:
            data = np.loadtxt(datafile_func(m))
            self.rad[m] = data[:, 1].copy()
            self.jnpar[m] = np.empty((nrad), dtype='D')
            self.jnpar[m].real = data[:, 2].copy()
            self.jnpar[m].imag = data[:, 3].copy()
            self.part_int[m] = np.empty((nrad), dtype='D')
            self.part_int[m].real = data[:, 6].copy()
            self.part_int[m].imag = data[:, 7].copy()
            self.bndry[m] = np.empty((nrad), dtype='D')
            self.bndry[m].real = data[:, 10].copy()
            self.bndry[m].imag = data[:, 11].copy()
            if symfluxcoord:
                self.psi[m] = data[:, 0].copy()
                self.Ichar[m] = data[:, 14].copy()
                self.Delta[m] = np.empty((nrad), dtype='D')
                self.Delta[m].real = data[:, 15].copy()
                self.Delta[m].imag = data[:, 16].copy()
        const = 2.0 * np.pi * statA_to_A
        self.width = {}
        self.intJ = {}
        self.intB = {}
        self.bndryB = {}
        if symfluxcoord:
            self.GPEC = {}
        for m in m_range:
            mid = np.searchsorted(self.rad[m], self.rres[m])
            half = max(0, min(nrad - mid - 1, mid - 1))
            self.width[m] = np.empty(half)
            self.intJ[m] = np.empty(half, dtype='D')
            self.intB[m] = np.empty(half, dtype='D')
            self.bndryB[m] = np.empty(half, dtype='D')
            if symfluxcoord:
                self.GPEC[m] = np.empty(half, dtype='D')
                integrand = self.jnpar[m]
                diff = self.psi[m]
            else:
                integrand = self.jnpar[m] * self.rad[m]
                diff = self.rad[m]
            for w in range(0, half):
                lo = mid - 1 - w
                hi = mid + w
                self.width[m][w] = self.rad[m][hi] - self.rad[m][lo]
                self.intJ[m][w] = np.trapz(integrand[lo:hi], diff[lo:hi]) * const
                self.intB[m][w] = np.trapz(self.part_int[m][lo:hi], diff[lo:hi]) * const
                self.bndryB[m][w] = (self.bndry[m][hi] - self.bndry[m][lo]) * const
                if symfluxcoord:
                    self.GPEC[m][w] = -1j / m * self.Ichar[m][mid] * (
                            self.Delta[m][hi] - self.Delta[m][lo]) * statA_to_A


    def process_KiLCA(self, datafile, nrad):
        self.rad = {}
        self.jnpar = {}
        self.rres = {}
        self.d = {}
        self.jnpar_int = {}
        self.Ipar = {}
        data = h5py.File(datafile, 'r')
        self.m_max = 0
        ### sgn_q = int(np.sign(data['/output/background/profiles/q_i'][0, -1]))
        for name, grp in data['/output'].items():
            if 'postprocessor' not in name:
                continue
            m = int(grp['mode'][0, 0]) # * sgn_q
            self.m_max = max(self.m_max, abs(m))
            self.rad[m] = np.array(grp['r'][0, :], dtype='d')
            self.jnpar[m] = np.zeros(self.rad[m].shape, dtype='D')
            self.jnpar[m].real = grp['Jpar'][0, :]
            if grp['Jpar'].shape[0] == 2:
                self.jnpar[m].imag = grp['Jpar'][1, :]
            self.jnpar_int[m] = interpolate.interp1d(
                self.rad[m], self.jnpar[m], kind='cubic',
                fill_value='extrapolate', assume_sorted=True
            )
            self.rres[m] = grp['rres'][0, 0].copy()
            self.d[m] = grp['d'][0, 0].copy()
            self.Ipar[m] = grp['Ipar'][0, 0].copy()
        self.hz = (data['/output/background/b0z'][0, :] /
                   data['/output/background/b0'][0, :])
        self.rad[0] = data['/output/background/R'][0, :].copy()
        self.hz_int = interpolate.interp1d(
            self.rad[0], self.hz, kind='cubic',
            fill_value='extrapolate', assume_sorted=True
        )
        data.close()
        const = 2.0 * np.pi * c1_statA_to_A
        mid = nrad // 2
        half = nrad // 2 - 1
        self.width = {}
        self.intJ = {}
        self.Imnpar = {}
        for m in self.jnpar.keys():
            interval = np.linspace(self.rres[m] - 1.5 * self.d[m],
                                   self.rres[m] + 1.5 * self.d[m], nrad)
            self.intJ[m] = np.empty(half, dtype='D')
            # kilca_intB = np.empty(half, dtype='D')
            # kilca_bndryB = np.empty(half, dtype='D')
            self.width[m] = np.empty(half)
            for w in range(0, half):
                lo = mid - 1 - w
                hi = mid + w
                self.width[m][w] = interval[hi] - interval[lo]
                self.intJ[m][w] = np.trapz(self.jnpar_int[m](interval[lo:hi]) *
                                           interval[lo:hi] * self.hz_int(interval[lo:hi]),
                                           interval[lo:hi]) * const
            interval = np.linspace(self.rres[m] - 0.5 * self.d[m],
                                   self.rres[m] + 0.5 * self.d[m], nrad // 2)
            self.Imnpar[m] = np.trapz(self.jnpar_int[m](interval) * interval *
                                      self.hz_int(interval), interval) * const

    def process_GPEC(self, datafile):
        rootgrp = netCDF4.Dataset(datafile, 'r')
        n = int(rootgrp.getncattr('n'))
        q = np.array(rootgrp.variables['q_rational'][:])
        m_min = int(np.rint(n * np.amin(np.abs(q))))
        m_max = int(np.rint(n * np.amax(np.abs(q))))
        self.I_res = {}
        self.w_isl = {}
        for k, m in enumerate(range(m_min, m_max + 1)):
            self.I_res[m] = np.complex(rootgrp.variables['I_res'][0, k],
                                       rootgrp.variables['I_res'][1, k])
            self.w_isl[m] = rootgrp.variables['w_isl'][k]


unit_J = r'\statampere\per\centi\meter\squared'
unit_p = r'\dyne\per\centi\meter\squared'
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
        self.triangulation = matplotlib.tri.Triangulation(
            self.data['/mesh/node_R'], self.data['/mesh/node_Z'],
            self.data['/mesh/tri_node'][()] - 1)

    def generate_RT0_triplots(self, grp, label, unit, filename):
        for dataset, decorator in RT0_comp.items():
            nameparts = path.splitext(filename)
            self.plots.append(magdif_2d_triplot(
                triangulation=self.triangulation,
                data=self.data[grp + dataset][()].real,
                label=fr"$\Real {decorator['math'](label)}$ / \si{{{unit}}}",
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_Re' + nameparts[1])))
            self.plots.append(magdif_2d_triplot(
                triangulation=self.triangulation,
                data=self.data[grp + dataset][()].imag,
                label=fr"$\Imag {decorator['math'](label)}$ / \si{{{unit}}}",
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_Im' + nameparts[1])))
            self.plots.append(magdif_2d_triplot(
                triangulation=self.triangulation, clim_scale = (0.0, 1.0),
                data=np.abs(self.data[grp + dataset][()]),
                label=fr"$\lvert {decorator['math'](label)} \rvert$ / \si{{{unit}}}",
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_abs' + nameparts[1])))
            self.plots.append(magdif_2d_triplot(
                triangulation=self.triangulation,
                data=np.angle(self.data[grp + dataset][()]),
                label=fr"$\arg {decorator['math'](label)}$",
                filename=path.join(self.datadir, nameparts[0] +
                                   decorator['file'] + '_arg' + nameparts[1])))

    def generate_L1_triplots(self, grp, label, unit, filename):
        nameparts = path.splitext(filename)
        self.plots.append(magdif_2d_triplot(
            triangulation=self.triangulation,
            data=self.data[grp + '/L1_DOF'][()].real,
            label=fr"$\Real {label}$ / \si{{{unit}}}",
            filename=path.join(self.datadir, nameparts[0] +
                               '_Re' + nameparts[1])))
        self.plots.append(magdif_2d_triplot(
            triangulation=self.triangulation,
            data=self.data[grp + '/L1_DOF'][()].imag,
            label=fr"$\Imag {label}$ / \si{{{unit}}}",
            filename=path.join(self.datadir, nameparts[0] +
                               '_Im' + nameparts[1])))
        self.plots.append(magdif_2d_triplot(
            triangulation=self.triangulation, clim_scale = (0.0, 1.0),
            data=np.absolute(self.data[grp + '/L1_DOF'][()]),
            label=fr"$\lvert {label} \rvert$ / \si{{{unit}}}",
            filename=path.join(self.datadir, nameparts[0] +
                               '_abs' + nameparts[1])))
        self.plots.append(magdif_2d_triplot(
            triangulation=self.triangulation,
            data=np.angle(self.data[grp + '/L1_DOF'][()]),
            label=fr"$\arg {label}$ / \si{{{unit}}}",
            filename=path.join(self.datadir, nameparts[0] +
                               '_arg' + nameparts[1])))

    def generate_default_plots(self):
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'][()], r'$r$ / \si{\centi\meter}',
                self.data['/cache/fs/psi'][()], r'$\psi$ / \si{\maxwell}',
                'disc poloidal flux', path.join(self.datadir, 'plot_psi.pdf')
        ))
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'][()], r'$r$ / \si{\centi\meter}',
                self.data['/cache/fs/q'][()], r'$q$',
                'safety factor', path.join(self.datadir, 'plot_q.pdf')
        ))
        self.plots.append(magdif_1d_cutplot(
                self.data['/cache/fs/rad'][()], r'$r$ / \si{\centi\meter}',
                self.data['/cache/fs/p'][()],
                r'$p_{0}$ / \si{\dyne\per\centi\meter\squared}',
                'pressure', path.join(self.datadir, 'plot_p0.pdf')
        ))
        # TODO: j0phi edge plot
        self.generate_RT0_triplots('/Bnvac', r'B_{n}', r'\gauss',
                                   'plot_Bnvac.png')
        self.generate_RT0_triplots('/iter/Bnplas', r'B_{n}', r'\gauss',
                                   'plot_Bnplas.png')
        self.generate_RT0_triplots('/iter/Bn', r'B_{n}', r'\gauss',
                                   'plot_Bn.png')
        self.generate_RT0_triplots('/iter/Bn_000', r'B_{n}', r'\gauss',
                                   'plot_Bn_000.png')
        self.generate_RT0_triplots('/iter/Bn_001', r'B_{n}', r'\gauss',
                                   'plot_Bn_001.png')
        self.generate_RT0_triplots('/iter/jn', r'J_{n}', unit_J,
                                   'plot_Jn.png')
        self.generate_RT0_triplots('/iter/jn_000', r'J_{n}', unit_J,
                                   'plot_Jn_000.png')
        self.generate_RT0_triplots('/iter/jn_001', r'J_{n}', unit_J,
                                   'plot_Jn_001.png')
        self.generate_L1_triplots('/iter/pn', r'p_{n}', unit_p, 'plot_pn.pdf')
        self.generate_L1_triplots('/iter/pn_000', r'p_{n}', unit_p,
                                  'plot_pn_000.png')
        self.generate_L1_triplots('/iter/pn_001', r'p_{n}', unit_p,
                                  'plot_pn_001.png')

        self.plots.append(magdif_conv_plot(
            self.config, self.data, path.join(self.datadir, 'convergence.pdf'))
        )

        if 'kilca_scale_factor' in self.config:
            kilca_scale_factor = self.config['kilca_scale_factor']
        else:
            kilca_scale_factor = 0
        grp = '/postprocess'
        if kilca_scale_factor == 0:
            pert = polmodes('full perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.psi_norm, grp + '/Bmn/coeff_rad')
            vac = polmodes('vacuum perturbation', 'r--')
            vac.read_magdif(self.data, fslabel.psi_norm, grp + '/Bmn_vac/coeff_rad')
            self.plots.append(magdif_poloidal_plots(self.datadir,
                    'Bmn_psi.pdf', self.config, self.data, fslabel.psi_norm,
                    r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\maxwell}',
                    np.abs, pert, vac
            ))
            pert = polmodes('initial perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.psi_norm, grp + '/jmn_000/coeff_pol')
            self.plots.append(magdif_poloidal_plots(self.datadir,
                    'jmn_000_theta.pdf', self.config, self.data, fslabel.psi_norm,
                    r'$\lvert J_{mn \theta}^{(0)} \rvert$'
                    + r' / \si{\statampere\per\centi\meter}', np.abs, pert
            ))

        else:
            pert = polmodes('full perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.r, grp + '/Bmn/coeff_rad')
            vac = polmodes('vacuum perturbation', 'r--')
            vac.read_magdif(self.data, fslabel.r, grp + '/Bmn_vac/coeff_rad')
            self.plots.append(magdif_poloidal_plots(self.datadir, 'Bmn_r.pdf',
                    self.config, self.data, fslabel.r,
                    r'$\lvert B_{mn}^{r} \rvert$ / \si{\gauss}',
                    np.abs, pert, vac
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
