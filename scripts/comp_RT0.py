#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 15:56:57 2020

@author: patrick
"""

from os import path
from matplotlib import use, rcParams, ticker, pyplot as plt
import numpy as np
from scipy import interpolate

rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'
use('Agg')
scifmt = ticker.ScalarFormatter()
scifmt.set_powerlimits((-3, 4))
canvas = (6.6, 3.6)
res = 300
thin = 0.5

def compare(datafile, outdir, prefix, subtitle, m, avg):
    try:
        data = np.loadtxt(datafile)
    except Exception as err:
        print('Error: {}'.format(err))
        exit
    fname = path.join(outdir, prefix)
    title = "Vacuum perturbation field at fixed $r$,\n{}".format(subtitle)
    
    sh = (data.shape[0])
    # r = data[:, 0]
    th = data[:, 1]
    BR = np.empty(sh, dtype=complex)
    BR.real = data[:, 2]
    BR.imag = data[:, 3]
    Bph = np.empty(sh, dtype=complex)
    Bph.real = data[:, 4]
    Bph.imag = data[:, 5]
    BZ = np.empty(sh, dtype=complex)
    BZ.real = data[:, 6]
    BZ.imag = data[:, 7]
    BR_int = np.empty(sh, dtype=complex)
    BR_int.real = data[:, 8]
    BR_int.imag = data[:, 9]
    Bph_int = np.empty(sh, dtype=complex)
    Bph_int.real = data[:, 10]
    Bph_int.imag = data[:, 11]
    BZ_int = np.empty(sh, dtype=complex)
    BZ_int.real = data[:, 12]
    BZ_int.imag = data[:, 13]
    
    if avg:
        BR_int = 0.5 * (BR_int[0::2] + BR_int[1::2])
        Bph_int = 0.5 * (Bph_int[0::2] + Bph_int[1::2])
        BZ_int = 0.5 * (BZ_int[0::2] + BZ_int[1::2])
        th_int = 0.5 * (th[0::2] + th[1::2])
    else:
        th_int = th
    
    BR_spl = interpolate.interp1d(th_int, BR_int, kind='cubic', copy=False,
                                  fill_value='extrapolate', assume_sorted=True)
    Bph_spl = interpolate.interp1d(th_int, Bph_int, kind='cubic', copy=False,
                                  fill_value='extrapolate', assume_sorted=True)
    BZ_spl = interpolate.interp1d(th_int, BZ_int, kind='cubic', copy=False,
                                  fill_value='extrapolate', assume_sorted=True)
    
    plt.figure(figsize=canvas)
    plt.plot(np.rad2deg(th_int), np.imag(BR_int), '-r', lw=thin,
             label='interpolated')
    plt.plot(np.rad2deg(th), np.imag(BR), '-k', lw=thin,
             label='exact')
    plt.plot(np.rad2deg(th), np.imag(BR_spl(th)) - np.imag(BR), '-b', lw=thin,
             label='difference')
    plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend(loc='lower right')
    plt.xlabel(r'$\theta$ / °')
    plt.ylabel(r'$\mathrm{Im} \, B_{n}^{R}$ / G')
    plt.title(title)
    plt.savefig(fname + '_BR.png', dpi=res)
    plt.close()
    
    plt.figure(figsize=canvas)
    plt.plot(np.rad2deg(th_int), np.imag(BZ_int), '-r', lw=thin,
             label='interpolated')
    plt.plot(np.rad2deg(th), np.imag(BZ), '-k', lw=thin,
             label='exact')
    plt.plot(np.rad2deg(th), np.imag(BZ_spl(th)) - np.imag(BZ), '-b', lw=thin,
             label='difference')
    plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend(loc='lower right')
    plt.xlabel(r'$\theta$ / °')
    plt.ylabel(r'$\mathrm{Im} \, B_{n}^{Z}$ / G')
    plt.title(title)
    plt.savefig(fname + '_BZ.png', dpi=res)
    plt.close()
    
    plt.figure(figsize=canvas)
    plt.plot(np.rad2deg(th_int), np.real(Bph_int), '-r', lw=thin,
             label='interpolated')
    plt.plot(np.rad2deg(th), np.real(Bph_spl(th)) - np.real(Bph), '-b', lw=thin,
             label='difference')
    plt.plot(np.rad2deg(th), np.real(Bph), '-k', lw=thin,
             label='exact')
    plt.gca().get_xaxis().set_major_locator(ticker.MultipleLocator(45))
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend(loc='lower right')
    plt.xlabel(r'$\theta$ / °')
    plt.ylabel(r'$\mathrm{Re} \, R B_{n}^{\varphi}$ / G')
    plt.title(title)
    plt.savefig(fname + '_Bph.png', dpi=res)
    plt.close()
    
    Br = BZ * np.sin(th) + BR * np.cos(th)
    Bth = BZ * np.cos(th) - BR * np.sin(th)
    Br_int = BZ_int * np.sin(th_int) + BR_int * np.cos(th_int)
    Bth_int = BZ_int * np.cos(th_int) - BR_int * np.sin(th_int)
    
    return {'Bmn_r': np.mean(Br * np.exp(-1j * m * th)),
            'Bmn_theta': np.mean(Bth * np.exp(-1j * m * th)),
            'Bmn_z': np.mean(Bph * np.exp(-1j * m * th)),
            'int_Bmn_r': np.mean(Br_int * np.exp(-1j * m * th_int)),
            'int_Bmn_theta': np.mean(Bth_int * np.exp(-1j * m * th_int)),
            'int_Bmn_z': np.mean(Bph_int * np.exp(-1j * m * th_int))}


coeff_RT0_GL2 = compare('cmp_RT0.dat', '.', 'cmp_RT0_GL2_cut',
        'RT0 elements, GL quadrature order 2, w/o averaging', 9, False)
coeff_RT0_GL2_avg = compare('cmp_RT0.dat', '.', 'cmp_RT0_GL2_cut_avg',
        'RT0 elements, GL quadrature order 2, w/ averaging', 9, True)
print("Bmn_r = {:.15e}".format(coeff_RT0_GL2['Bmn_r']))
print("Bmn_r_RT0_GL2 = {:.15e}".format(coeff_RT0_GL2['int_Bmn_r']))
print("Bmn_r_RT0_GL2_avg = {:.15e}".format(coeff_RT0_GL2_avg['int_Bmn_r']))
print("Bmn_z = {:.15e}".format(coeff_RT0_GL2['Bmn_z']))
print("Bmn_z_RT0_GL2 = {:.15e}".format(coeff_RT0_GL2['int_Bmn_z']))
print("Bmn_z_RT0_GL2_avg = {:.15e}".format(coeff_RT0_GL2_avg['int_Bmn_z']))

coeff_BDM1 = compare('cmp_RT0_ff.dat', '.', 'cmp_BDM1_cut',
        'BDM1 elements, w/o averaging', 9, False)
coeff_BDM1_avg = compare('cmp_RT0_ff.dat', '.', 'cmp_BDM1_cut_avg',
        'BDM1 elements, w/ averaging', 9, True)
print("Bmn_r = {:.15e}".format(coeff_BDM1['Bmn_r']))
print("Bmn_r_BDM1 = {:.15e}".format(coeff_BDM1['int_Bmn_r']))
print("Bmn_r_BDM1_avg = {:.15e}".format(coeff_BDM1_avg['int_Bmn_r']))
print("Bmn_z = {:.15e}".format(coeff_BDM1['Bmn_z']))
print("Bmn_z_BDM1 = {:.15e}".format(coeff_BDM1['int_Bmn_z']))
print("Bmn_z_BDM1_avg = {:.15e}".format(coeff_BDM1_avg['int_Bmn_z']))

coeff_optpolres = compare('cmp_polres.dat', '.', 'cmp_optpolres_cut',
        'BDM1 elements, optimized poloidal resolution', 9, False)
print("Bmn_r = {:.15e}".format(coeff_optpolres['Bmn_r']))
print("Bmn_r_optpolres = {:.15e}".format(coeff_optpolres['int_Bmn_r']))
print("Bmn_z = {:.15e}".format(coeff_optpolres['Bmn_z']))
print("Bmn_z_optpolres = {:.15e}".format(coeff_optpolres['int_Bmn_z']))
