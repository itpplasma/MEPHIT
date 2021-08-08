#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 15:56:57 2020

@author: patrick
"""

from os import path
from h5py import File
import magdifplot
from matplotlib import ticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from scipy import interpolate

def compare(datafile, outdir, prefix, subtitle, m):
    data = File(datafile, 'r')
    kf = data['/mesh/res_ind'][m - data['/mesh/m_res_min'][()]]
    ktri_min = data['/mesh/kt_low'][kf]
    ktri_max = data['/mesh/kt_low'][kf] + data['/mesh/kt_max'][kf] - 1
    psi_n = (data['/cache/fs_half/psi'][kf] - data['/cache/fs/psi'][0]) / \
            (data['/cache/fs/psi'][-1] - data['/cache/fs/psi'][0])
    th = data['/cache/sample_polmodes/theta'][ktri_min:ktri_max]
    BR = data['/debug_RT0/Bn_R'][ktri_min:ktri_max]
    Bph = data['/debug_RT0/Bn_phi'][ktri_min:ktri_max]
    BZ = data['/debug_RT0/Bn_Z'][ktri_min:ktri_max]
    BR_int = data['/debug_RT0/Bn_R_RT0'][ktri_min:ktri_max]
    Bph_int = data['/debug_RT0/Bn_phi_RT0'][ktri_min:ktri_max]
    BZ_int = data['/debug_RT0/Bn_Z_RT0'][ktri_min:ktri_max]
    BR_spl = interpolate.interp1d(th, BR_int, kind='linear', copy=False, fill_value='extrapolate', assume_sorted=True)
    Bph_spl = interpolate.interp1d(th, Bph_int, kind='linear', copy=False, fill_value='extrapolate', assume_sorted=True)
    BZ_spl = interpolate.interp1d(th, BZ_int, kind='linear', copy=False, fill_value='extrapolate', assume_sorted=True)
    data.close()

    title = f"Vacuum perturbation field at fixed $\\hat{{\\psi}} = {psi_n:.3}$,\n{subtitle}"
    fname = path.join(outdir, prefix)

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.real(BR), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.real(BR), '--k', label='exact')
    ax.plot(np.rad2deg(th), np.real(BR_spl(th)) - np.real(BR), '-b', label='difference')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Real B_{n}^{R}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_R_Re.pdf')

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.imag(BR), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.imag(BR), '--k', label='exact')
    ax.plot(np.rad2deg(th), np.imag(BR_spl(th)) - np.imag(BR), '-b', label='difference')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Imag B_{n}^{R}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_R_Im.pdf')

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.real(BZ), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.real(BZ), '--k', label='exact')
    ax.plot(np.rad2deg(th), np.real(BZ_spl(th)) - np.real(BZ), '-b', label='difference')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Real B_{n}^{Z}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_Z_Re.pdf')

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.imag(BZ), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.imag(BZ), '--k', label='exact')
    ax.plot(np.rad2deg(th), np.imag(BZ_spl(th)) - np.imag(BZ), '-b', label='difference')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Imag B_{n}^{Z}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_Z_Im.pdf')

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.real(Bph_int), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.real(Bph_spl(th)) - np.real(Bph), '-b', label='difference')
    ax.plot(np.rad2deg(th), np.real(Bph), '--k', label='exact')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Real R B_{n}^{\varphi}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_phi_Re.pdf')

    fig = Figure()
    ax = fig.subplots()
    ax.plot(np.rad2deg(th), np.imag(Bph_int), '-r', label='interpolated')
    ax.plot(np.rad2deg(th), np.imag(Bph_spl(th)) - np.imag(Bph), '-b', label='difference')
    ax.plot(np.rad2deg(th), np.imag(Bph), '--k', label='exact')
    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    ax.legend()
    ax.set_xlabel(r'$\vartheta$ / \si{\degree}')
    ax.set_ylabel(r'$\Imag R B_{n}^{\varphi}$ / \si{\gauss}')
    ax.set_title(title)
    canvas = FigureCanvas(fig)
    fig.savefig(fname + '_Bn_phi_Im.pdf')

workdir = '/home/patrick/git/NEO-EQ/run/33353_2325'
compare(path.join(workdir, 'magdif_fix.h5'), workdir, 'res_fix', 'with fixed poloidal resolution', 6)
compare(path.join(workdir, 'magdif.h5'), workdir, 'res_opt', 'with optimized poloidal resolution', 6)
