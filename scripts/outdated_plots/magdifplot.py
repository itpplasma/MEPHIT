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
from cycler import cycler
from matplotlib import rcParams
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.tri import Triangulation
from matplotlib.colors import Normalize
import colorcet
import numpy as np
from scipy import interpolate
from multiprocessing import Pool

# supply path to scripts directory
scripts_dir = path.dirname(path.realpath(__file__))
run_dir = path.realpath(scripts_dir + '/../run')

# complex values are stored as compound types in libneo/hdf5tools
h5py_hack().complex_names = ('real', 'imag')

# matplotlib default configuration
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
latex_preamble = path.join(scripts_dir, 'outdated_plots/magdifplot.tex')
rcParams['text.latex.preamble'] = fr"\input{{{latex_preamble}}}"

c_cgs = 2.9979246e+10
cm_to_m = 1.0e-02
G_to_T = 1.0e-04
Mx_to_Wb = 1.0e-08
c1_statA_to_A = 1.0e+01
c1_statA_per_cm2_to_A_per_m2 = c1_statA_to_A / cm_to_m ** 2
statA_to_A = c1_statA_to_A / c_cgs
statA_per_cm2_to_A_per_m2 = c1_statA_per_cm2_to_A_per_m2 / c_cgs


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
        fig = Figure(figsize=(3.3, 4.4))
        ax = fig.subplots()
        im = ax.tripcolor(self.triangulation, self.data, cmap=colorcet.cm.coolwarm)
        ax.set_aspect('equal')
        cbar = fig.colorbar(im)
        cbar.set_label(self.label, rotation=90)
        im.set_clim([-max(abs(self.data)) * self.clim_scale[0], max(abs(self.data)) * self.clim_scale[1]])
        ax.set_xlabel(r'$R$ / cm')
        ax.set_ylabel(r'$Z$ / cm')
        if self.title is not None:
            ax.set_title(self.title)
        canvas = FigureCanvas(fig)
        fig.savefig(self.filename, dpi=300)


class magdif_2d_rectplots:
    def __init__(self, R, Z, data, label, filename, title=None, clim_scale=None, centered=True):
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
        self.centered = centered

    def dump_plot(self):
        print(f"plotting {self.filename}")
        xlim = (min(R[0] for R in self.R), max(R[-1] for R in self.R))
        ylim = (min(Z[0] for Z in self.Z), max(Z[-1] for Z in self.Z))
        fig = Figure(figsize=(8.0, 5.0))
        axs = fig.subplots(1, 2, sharex='all', sharey='all')
        images = []
        for k in range(2):
            images.append(axs[k].imshow(self.data[k], cmap=colorcet.cm.coolwarm, interpolation='gaussian',
                                        extent=[self.R[k][0], self.R[k][-1], self.Z[k][0], self.Z[k][-1]],
                                        origin='lower'))
            axs[k].set_aspect('equal')
            axs[k].set_xlabel(r'$R$ / m')
            axs[k].set_ylabel(r'$Z$ / m')
            axs[k].set_xlim(xlim)
            axs[k].set_ylim(ylim)
            if self.title is not None:
                axs[k].set_title(self.title[k])
        axs[1].yaxis.set_tick_params(labelleft=True)
        axs[1].yaxis.offsetText.set_visible(True)
        if self.centered:
            clim = (-max(np.amax(np.abs(image.get_array())) for image in images) * self.clim_scale[0],
                    max(np.amax(np.abs(image.get_array())) for image in images) * self.clim_scale[1])
        else:
            clim = (min(np.amin(image.get_array()) for image in images) * self.clim_scale[0],
                    max(np.amax(image.get_array()) for image in images) * self.clim_scale[1])
        norm = Normalize(vmin=clim[0], vmax=clim[1])
        for im in images:
            im.set_norm(norm)
        for k in range(2):
            axs[k].contour(self.R[k], self.Z[k], self.data[k], np.linspace(clim[0], clim[1], 10),
                           colors='k', linewidths=0.1)
        cbar = fig.colorbar(images[0], ax=axs)
        cbar.set_label(self.label, rotation=90)
        canvas = FigureCanvas(fig)
        fig.savefig(self.filename, dpi=300)


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
        fig = Figure(figsize=(6.6, 3.6))
        ax = fig.subplots()
        ax.plot(self.x, self.y, '-k')
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        canvas = FigureCanvas(fig)
        fig.savefig(self.filename)


class magdif_conv_plot:
    def __init__(self, data, filename, xlim=None, ylim=None, title=None):
        self.data = data
        self.filename = filename
        self.xlim = xlim
        self.ylim = ylim
        self.title = title

    def dump_plot(self):
        print(f"plotting {self.filename}")
        sup_eigval = self.data['/config/ritz_threshold']
        L2int_Bnvac = self.data['/iter/L2int_Bnvac'][()]
        L2int_Bn_diff = self.data['/iter/L2int_Bn_diff'][()]
        kiter = np.arange(0, len(L2int_Bn_diff))
        fig = Figure(figsize=(6.6, 3.6))
        ax = fig.subplots()
        ax.semilogy(kiter, L2int_Bnvac * sup_eigval ** kiter, 'r-',
                    label=r'$\lvert \lambda_{\text{max}} \rvert^{k} \lVert \mathbf{B}_{n}^{(0)} \rVert_{2}$')
        ax.semilogy(kiter, L2int_Bn_diff, 'xk', label=r'$\lVert \delta \mathbf{B}_{n}^{(k)} \rVert_{2}$')
        if self.xlim is not None:
            ax.set_xlim(self.xlim)
        if self.ylim is not None:
            ax.set_ylim(self.ylim)
        ax.legend(loc='upper right')
        ax.set_xticks(kiter)
        ax.set_xlabel('iteration step $k$')
        ax.set_ylabel(r'$\lVert R \mathbf{B}^{\text{pol}} \rVert_{2}$ / \si{\maxwell}')
        if self.title is not None:
            ax.set_title(self.title)
        else:
            ax.set_title('estimation of convergence')
        canvas = FigureCanvas(fig)
        fig.savefig(self.filename)


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
        if self.rad_coord is fslabel.psi_norm:
            psi_axis = data['/cache/fs/psi'][0]
            psi_edge = data['/cache/fs/psi'][-1]
            rho = (data['/cache/fs_half/psi'][()] - psi_axis) / (psi_edge - psi_axis)
        else:
            rho = data['/cache/fs_half/rad'][()]
        self.m_max = (data[var_name].shape[1] - 1) // 2
        for m in range(-self.m_max, self.m_max + 1):
            self.rho[m] = rho
            self.var[m] = np.array(data[var_name][:, m + self.m_max], dtype='D') * Mx_to_Wb

    def read_fouriermodes(self, data):
        self.type = 'amn.dat'
        self.rad_coord = fslabel.psi_norm
        rho = data['/debug_fouriermodes/psi_n'][()]
        var_name = '/debug_fouriermodes/comp_psi_contravar_dens'
        self.m_max = (data[var_name].shape[1] - 1) // 2
        for m in range(-self.m_max, self.m_max + 1):
            self.rho[m] = rho
            self.var[m] = np.array(data[var_name][:, m + self.m_max], dtype='D') * Mx_to_Wb

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
            self.var[m].real = grp[var_name][0, :] * Mx_to_Wb
            if grp[var_name].shape[0] == 2:
                self.var[m].imag = grp[var_name][1, :] * Mx_to_Wb
        data.close()

    def read_GPEC(self, datafile, sgn_dpsi, var_name='Jbgradpsi'):
        self.type = 'GPEC'
        self.rad_coord = fslabel.psi_norm
        self.m_max = 0
        rootgrp = netCDF4.Dataset(datafile, 'r')
        rho = np.array(rootgrp.variables['psi_n'])
        helicity = int(rootgrp.getncattr('helicity'))
        for k, m_out in enumerate(rootgrp.variables['m_out'][:]):
            m = m_out * helicity
            self.m_max = max(self.m_max, abs(m))
            self.rho[m] = rho
            self.var[m] = np.empty(rho.shape, dtype='D')
            # GPEC always uses normal vectors pointing outwards
            # and includes the factor 2 for Fourier series of a real function in the coefficient
            self.var[m].real = rootgrp.variables[var_name][0, k, :] * 0.5 * sgn_dpsi
            # GPEC uses clockwise toroidal angle for positive helicity
            # and expands Fourier series in negative toroidal angle
            self.var[m].imag = rootgrp.variables[var_name][1, k, :] * 0.5 * sgn_dpsi * helicity
        rootgrp.close()


class magdif_poloidal_plots:
    def __init__(self, datadir, filename, data, rad_coord, ylabel, comp, *poldata):
        self.datadir = datadir
        self.filename = filename
        self.data = data
        self.rad_coord = rad_coord
        self.xlabel = self.rad_coord.value
        self.ylabel = ylabel
        self.comp = comp
        self.omit_res = False  # omit resonance position indicator when zooming in
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
        m_res = np.arange(self.data['/mesh/m_res_min'][()],
                          self.data['/mesh/m_res_max'][()] + 1) * sgn_m_res
        if self.rad_coord is fslabel.psi_norm:
            rho = self.data['/cache/fs/psi']
            resonance = dict(zip(m_res, (self.data['/mesh/psi_res'] - rho[0]) / (rho[-1] - rho[0])))
            # normalize psi
            rho = (rho - rho[0]) / (rho[-1] - rho[0])
        else:
            rho = self.data['/cache/fs/rad'][()]
            resonance = dict(zip(m_res, self.data['/mesh/rad_norm_res'] * rho[-1]))

        fmt = path.basename(path.splitext(self.filename)[0] + '_{}' + path.splitext(self.filename)[1])
        horz_plot = 2
        vert_plot = 1
        two_squares = (6.6, 3.3)
        # plot non-symmetric modes
        m_max = min(map(lambda d: d.m_max, self.poldata))
        for m_abs in range(1, m_max + 1):
            filename = fmt.format(m_abs)
            print(f"plotting {filename}")
            fig = Figure(figsize=two_squares)
            axs = fig.subplots(vert_plot, horz_plot, sharex='all', sharey='all')
            for k in range(horz_plot):
                m = (2 * k - 1) * m_abs
                axs[k].axhline(0.0, color='k', alpha=0.5, lw=0.5)
                if m in resonance and not self.omit_res:
                    axs[k].axvline(resonance[m], color='b', alpha=0.5, lw=0.5, label='resonance position')
                for data in self.poldata:
                    if m in data.var:
                        axs[k].plot(self.interp_rho(data, m), self.comp(data.var[m]), data.fmt, label=data.label)
                axs[k].legend(fontsize='x-small')
                axs[k].set_title(('resonant ' if m in resonance else 'non-resonant ') + fr"$m = {m}$")
                axs[k].set_xlabel(self.xlabel)
                axs[k].set_ylabel(self.ylabel)
            axs[1].yaxis.set_tick_params(labelleft=True)
            canvas = FigureCanvas(fig)
            fig.savefig(path.join(self.datadir, filename))
        # plot symmetric mode and safety factor
        m = 0
        filename = fmt.format(m)
        print(f"plotting {filename}")
        fig = Figure(figsize=two_squares)
        axs = fig.subplots(vert_plot, horz_plot)
        for data in self.poldata:
            if m in data.var:
                axs[0].plot(self.interp_rho(data, m), self.comp(data.var[m]), data.fmt, label=data.label)
        axs[0].legend(fontsize='x-small')
        axs[0].set_title(f"$m = {m}$")
        axs[0].set_xlabel(self.xlabel)
        axs[0].set_ylabel(self.ylabel)
        for res in resonance.values():
            axs[1].axvline(np.abs(res), color='b', alpha=0.5, lw=0.5)
        q = self.data['/cache/fs/q'][()]
        axs[1].plot(rho, q, 'k-')
        if np.any(q > 0.0):
            axs[1].set_ylim(bottom=0.0)
        else:
            axs[1].set_ylim(top=0.0)
        axs[1].set_xlabel(self.xlabel)
        axs[1].set_ylabel(r'$q$')
        canvas = FigureCanvas(fig)
        fig.savefig(path.join(self.datadir, filename))


class parcurr:
    def __init__(self):
        pass

    def process_magdif(self, data_func, m_range, nrad, rres, symfluxcoord=True):
        self.rad = {}
        self.jnpar = {}
        self.rres = rres
        self.part_int = {}
        self.bndry = {}
        if symfluxcoord:
            self.psi = {}
            self.Ichar = {}
            self.Delta = {}
            dataset_func = lambda m: f"/postprocess/Imn_par_{m}"
        else:
            dataset_func = lambda m: '/postprocess/Imn_par_KiLCA'
        for m in m_range:
            self.rad[m] = data_func(m)[dataset_func(m) + '/rad'][()]
            self.jnpar[m] = data_func(m)[dataset_func(m) + '/jmn_par_neg'][()]
            self.part_int[m] = data_func(m)[dataset_func(m) + '/part_int_neg'][()]
            self.bndry[m] = data_func(m)[dataset_func(m) + '/bndry_neg'][()]
            if symfluxcoord:
                self.psi[m] = data_func(m)[dataset_func(m) + '/psi'][()]
                self.Ichar[m] = data_func(m)[dataset_func(m) + '/I_char'][()]
                self.Delta[m] = data_func(m)[dataset_func(m) + '/Delta_mn_neg'][()]
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
            self.I_res[m] = complex(rootgrp.variables['I_res'][0, k],
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

    def __init__(self, datadir, datafile='mephit.h5'):
        self.plots = []
        self.datadir = datadir
        self.datafile = datafile

    def read_datafile(self):
        print(f"reading contents of {self.datafile}")
        self.data = h5py.File(path.join(self.datadir, self.datafile), 'r')
        self.triangulation = Triangulation(self.data['/mesh/node_R'], self.data['/mesh/node_Z'],
                                           self.data['/mesh/tri_node'][()] - 1)

    def read_convexwall(self, file=scripts_dir+'/../data/convexwall_asdex.dat'):
        print(f"reading contents of {file}")
        self.convexwall = np.loadtxt(file)

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
        self.generate_RT0_triplots('/vac/Bn', r'B_{n}', r'\gauss',
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
            self.data, path.join(self.datadir, 'convergence.pdf'))
        )

        grp = '/postprocess'
        kilca_scale_factor = self.data['/config/kilca_scale_factor'][()]
        fsl = fslabel.psi_norm if (kilca_scale_factor == 0) else fslabel.r
        if kilca_scale_factor == 0:
            pert = polmodes('full perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.psi_norm, grp + '/Bmn/coeff_rad')
            vac = polmodes('vacuum perturbation', 'r--')
            vac.read_magdif(self.data, fslabel.psi_norm, grp + '/Bmn_vac/coeff_rad')
            self.plots.append(magdif_poloidal_plots(self.datadir,
                    'Bmn_psi.pdf', self.data, fslabel.psi_norm,
                    r'$\lvert \sqrt{g} B_{mn}^{\psi} \rvert$ / \si{\maxwell}',
                    np.abs, pert, vac
            ))
            pert = polmodes('initial perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.psi_norm, grp + '/jmn_000/coeff_pol')
            self.plots.append(magdif_poloidal_plots(self.datadir,
                    'jmn_000_theta.pdf', self.data, fslabel.psi_norm,
                    r'$\lvert J_{mn \theta}^{(0)} \rvert$'
                    + r' / \si{\statampere\per\centi\meter}', np.abs, pert
            ))

        else:
            pert = polmodes('full perturbation', 'k-')
            pert.read_magdif(self.data, fslabel.r, grp + '/Bmn/coeff_rad')
            vac = polmodes('vacuum perturbation', 'r--')
            vac.read_magdif(self.data, fslabel.r, grp + '/Bmn_vac/coeff_rad')
            self.plots.append(magdif_poloidal_plots(self.datadir, 'Bmn_r.pdf',
                    self.data, fslabel.r,
                    r'$\lvert B_{mn}^{r} \rvert$ / \si{\gauss}',
                    np.abs, pert, vac
            ))
        pmn_initial = polmodes('initial iteration', 'k--')
        pmn_initial.read_magdif(self.data, fsl, grp + '/pmn_000/coeff')
        for k in pmn_initial.var.keys():
            pmn_initial.rho[k] = (self.data['/cache/fs/psi'][()] - self.data['/cache/fs/psi'][0]) /\
                                 (self.data['/cache/fs/psi'][-1] - self.data['/cache/fs/psi'][0])
            pmn_initial.var[k] /= Mx_to_Wb
        pmn_final = polmodes('final iteration', 'r--')
        pmn_final.read_magdif(self.data, fsl.psi_norm, grp + '/pmn/coeff')
        for k in pmn_final.var.keys():
            pmn_final.rho[k] = (self.data['/cache/fs/psi'][()] - self.data['/cache/fs/psi'][0]) /\
                               (self.data['/cache/fs/psi'][-1] - self.data['/cache/fs/psi'][0])
            pmn_final.var[k] /= Mx_to_Wb
        self.plots.append(magdif_poloidal_plots(self.datadir, 'pmn.pdf', self.data, fsl,
                                                r'$\lvert p_{mn} \rvert$ / \si{\dyne\per\centi\meter\squared}',
                                                np.abs, pmn_initial, pmn_final))

    def dump_plots(self):
        for p in self.plots:
            p.dump_plot()

    def wrapper(self, index):
        self.plots[index].dump_plot()

    def dump_plots_parallel(self):
        with Pool(max(1, cpu_count() - 1)) as p:
            p.map(self.wrapper, range(len(self.plots)))


if __name__ == '__main__':
    testcase = magdif(argv[1], argv[2])
    testcase.read_datafile()
    testcase.generate_default_plots()
    testcase.dump_plots_parallel()
