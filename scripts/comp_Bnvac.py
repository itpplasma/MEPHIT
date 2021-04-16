#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:33:27 2020

@author: patrick
"""

import numpy as np
from scipy import interpolate as interp
import matplotlib
from matplotlib import pyplot as plt
from magdifplot import fslabel, magdif_mnDat as mnDat
import h5py
from os import path

scifmt = matplotlib.ticker.ScalarFormatter()
scifmt.set_powerlimits((-3, 4))
canvas = (6.6, 3.6)
res = 300
thin = 0.5

def get_vac_datapoints(kilca_file):
    m_min = 3
    m_max = 9
    kilca = h5py.File(kilca_file, 'r')
    for m in range(m_min, m_max+1):
        k = m - m_min + 1
        loc = '/output/postprocessor{}'.format(k)
        kilca_r = kilca[loc + '/r'][0, :]
        if kilca[loc + '/Br'].shape[0] == 2:
            kilca_Br = np.empty(kilca[loc + '/Br'].shape[1:], dtype=complex)
            kilca_Br.real = kilca[loc + '/Br'][0, ...]
            kilca_Br.imag = kilca[loc + '/Br'][1, ...]
        else:
            kilca_Br = kilca[loc + '/Br'][0, :]
        if kilca[loc + '/Bth'].shape[0] == 2:
            kilca_Bth = np.empty(kilca[loc + '/Bth'].shape[1:], dtype=complex)
            kilca_Bth.real = kilca[loc + '/Bth'][0, ...]
            kilca_Bth.imag = kilca[loc + '/Bth'][1, ...]
        else:
            kilca_Bth = kilca[loc + '/Bth'][0, :]
        if kilca[loc + '/Bz'].shape[0] == 2:
            kilca_Bz = np.empty(kilca[loc + '/Bz'].shape[1:], dtype=complex)
            kilca_Bz.real = kilca[loc + '/Bz'][0, ...]
            kilca_Bz.imag = kilca[loc + '/Bz'][1, ...]
        else:
            kilca_Bz = kilca[loc + '/Bz'][0, :]
        print("kilca_vac_r({}) = {:.15e}".format(m, kilca_r[-1]))
        print("kilca_vac_Bz({0}) = ({1.real:.15e}, {1.imag:.15e})".format(m, 
              kilca_Bz[-1]))
        plt.figure(figsize=canvas)
        plt.plot(kilca_r, np.imag(kilca_Br), '-k', lw=thin,
                 label=r'$\mathrm{Im} \, B_{mn}^{r}$')
        plt.plot(kilca_r, np.real(kilca_Bth), '-r', lw=thin,
                 label=r'$\mathrm{Re} \, r B_{mn}^{\theta}$')
        plt.plot(kilca_r, np.real(kilca_Bz), '-b', lw=thin,
                 label=r'$\mathrm{Re} \, B_{mn}^{z}$')
        plt.gca().legend(loc='upper left')
        plt.xlabel('$r$ / cm')
        plt.ylabel('$B_{mn}$ / G')
        plt.title('Vacuum perturbation field for n = 2, m = {}'.format(m))
        plt.savefig('KiLCA_B{}n.png'.format(m))
        plt.close()

def abserr(dev_x, dev_y, ref_x, ref_y):
    resampled = interp.interp1d(ref_x, ref_y, kind='cubic', copy=False,
                               fill_value='extrapolate', assume_sorted=True)
    return dev_y - resampled(dev_x)

def compare(workdir, m, gamma, nflux, ref):
    title = ("Vacuum perturbation field, error in postprocessing\n" +
        "$m$ = {}, ".format(m) + "scaling $\\gamma$ = {}, ".format(gamma) +
        "{} flux surfaces ".format(nflux) +
        ("(with refinement)" if ref else "(without refinement)"))
    
    pre_file = path.join(workdir, 'cmp_vac.dat')
    try:
        pre = np.loadtxt(pre_file)
    except Exception as err:
        print('Error: {}'.format(err))
        exit
    pre_mask = pre[:, 0] >= 5.0
    pre_r = pre[pre_mask, 0]
    pre_Br = np.empty(pre_r.shape, dtype=complex)
    pre_Br.real = pre[pre_mask, 7]
    pre_Br.imag = pre[pre_mask, 8]
    pre_Bth = np.empty(pre_r.shape, dtype=complex)
    pre_Bth.real = pre[pre_mask, 9]
    pre_Bth.imag = pre[pre_mask, 10]
    pre_Bz = np.empty(pre_r.shape, dtype=complex)
    pre_Bz.real = pre[pre_mask, 11]
    pre_Bz.imag = pre[pre_mask, 12]
    
    post_file = 'Bmn_vac_r.dat'
    post = mnDat(workdir, post_file, gamma, fslabel.r, 'magdif')
    post.process()
    post_mask = post.rho >= 5.0
    post_r = post.rho[post_mask]
    post_Br = post.abs[post_mask, post.offset + m]
    post_Br_arg = post.arg[post_mask, post.offset + m]
    post_file = 'Bmn_vac_theta.dat'
    post = mnDat(workdir, post_file, gamma, fslabel.r, 'magdif')
    post.process()
    post_mask = post.rho >= 5.0
    post_r = post.rho[post_mask]
    post_Bth = post.abs[post_mask, post.offset + m]
    post_Bth_arg = post.arg[post_mask, post.offset + m]
    post_file = 'Bmn_vac_z.dat'
    post = mnDat(workdir, post_file, gamma, fslabel.r, 'magdif')
    post.process()
    post_mask = post.rho >= 5.0
    post_r = post.rho[post_mask]
    post_Bz = post.abs[post_mask, post.offset + m]
    post_Bz_arg = post.arg[post_mask, post.offset + m]
    
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.abs(pre_Br), '-k', lw=thin, label='exact')
    plt.plot(post_r, np.abs(post_Br), '-r', lw=thin, label='interpolated')
    plt.plot(post_r, abserr(post_r, post_Br, pre_r, np.abs(pre_Br)),
             '-b', lw=thin, label='difference')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$|B_{mn}^{r}|$ / G')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_r.png', dpi=res)
    plt.close()
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.rad2deg(np.angle(pre_Br)), '-k', lw=thin,
             label='exact')
    plt.plot(post_r, np.rad2deg(post_Br_arg), '-r', lw=thin,
             label='interpolated')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$\arg B_{mn}^{r}$ / °')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_r_arg.png', dpi=res)
    plt.close()
    
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.abs(pre_Bth), '-k', lw=thin, label='exact')
    plt.plot(post_r, np.abs(post_Bth), '-r', lw=thin, label='interpolated')
    plt.plot(post_r, abserr(post_r, post_Bth, pre_r, np.abs(pre_Bth)),
             '-b', lw=thin, label='difference')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$|r B_{mn}^{\theta}|$')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_theta.png', dpi=res)
    plt.close()
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.rad2deg(np.angle(pre_Bth)), '-k', lw=thin,
             label='exact')
    plt.plot(post_r, np.rad2deg(post_Bth_arg), '-r', lw=thin,
             label='interpolated')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$\arg r B_{mn}^{\theta}$ / °')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_theta_arg.png', dpi=res)
    plt.close()
    
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.abs(pre_Bz), '-k', lw=thin, label='exact')
    plt.plot(post_r, np.abs(post_Bz), '-r', lw=thin, label='interpolated')
    plt.plot(post_r, abserr(post_r, post_Bz, pre_r, np.abs(pre_Bz)),
             '-b', lw=thin, label='difference')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$|B_{mn}^{z}|$ / G')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_z.png', dpi=res)
    plt.close()
    plt.figure(figsize=canvas)
    plt.plot(pre_r, np.rad2deg(np.angle(pre_Bz)), '-k', lw=thin,
             label='exact')
    plt.plot(post_r, np.rad2deg(post_Bz_arg), '-r', lw=thin,
             label='interpolated')
    plt.gca().get_yaxis().set_major_formatter(scifmt)
    plt.gca().legend()
    plt.xlabel(post.rad_coord.value)
    plt.ylabel(r'$\arg B_{mn}^{z}$ / °')
    plt.title(title)
    plt.savefig(path.basename(workdir) + '_cmp_Bmn_vac_z_arg.png', dpi=res)
    plt.close()
    
    return {'pre_r': pre_r, 'pre_Br': pre_Br, 'pre_Bth': pre_Bth,
            'pre_Bz': pre_Bz, 'post_r': post_r, 'post_Br': post_Br,
            'post_Br_arg': post_Br_arg, 'post_Bth': post_Bth,
            'post_Bth_arg': post_Bth_arg, 'post_Bz': post_Bz,
            'post_Bz_arg': post_Bz_arg}


kilca_file = '/home/patrick/git/NEO-EQ/run/geomint_TCFP/TCFP_vac_hip.hdf5'
get_vac_datapoints(kilca_file)
test_dir = '/home/patrick/itp-temp/NEO-EQ/run'
compare(path.join(test_dir, 'test_hires_hiscale_ref'), 9, 1000, 110, True)
compare(path.join(test_dir, 'test_hires_hiscale_unref'), 9, 1000, 100, False)
compare(path.join(test_dir, 'test_hires_loscale_ref'), 9, 100, 110, True)
compare(path.join(test_dir, 'test_hires_loscale_unref'), 9, 100, 100, False)
compare(path.join(test_dir, 'test_lores_hiscale_ref'), 9, 1000, 85, True)
compare(path.join(test_dir, 'test_lores_hiscale_unref'), 9, 1000, 75, False)
compare(path.join(test_dir, 'test_lores_loscale_ref'), 9, 100, 85, True)
compare(path.join(test_dir, 'test_lores_loscale_unref'), 9, 100, 75, False)
compare('/home/patrick/git/NEO-EQ/run/geomint_TCFP', 9, 1000, 110, True)

test2_hhf = compare(path.join(test_dir, 'test2_hires_half_flux'),
        9, 1000, 110, True)
test2_hhg = compare(path.join(test_dir, 'test2_hires_half_geom'),
        9, 1000, 110, True)
test2_hff = compare(path.join(test_dir, 'test2_hires_full_flux'),
        9, 1000, 110, True)
test2_hfg = compare(path.join(test_dir, 'test2_hires_full_geom'),
        9, 1000, 110, True)
test2_lhf = compare(path.join(test_dir, 'test2_lores_half_flux'),
        9, 1000, 110, True)
test2_lhg = compare(path.join(test_dir, 'test2_lores_half_geom'),
        9, 1000, 110, True)
test2_lff = compare(path.join(test_dir, 'test2_lores_full_flux'),
        9, 1000, 110, True)
test2_lfg = compare(path.join(test_dir, 'test2_lores_full_geom'),
        9, 1000, 110, True)

plt.figure(figsize=canvas)
plt.plot(test2_lhf['pre_r'], np.abs(test2_lhf['pre_Bz']), '-k', lw=thin,
         label='exact')
plt.plot(test2_lhf['post_r'], test2_lhf['post_Bz'], '-r', lw=thin,
         label='flux')
plt.plot(test2_lhg['post_r'], test2_lhg['post_Bz'], '--b', lw=thin,
         label='geom.')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend(loc='upper left')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$|B_{mn}^{z}|$ / G')
plt.title('Error in vacuum perturbation field\n'
          'flux symmetry angle vs. geometrical angle')
plt.savefig('test2_cmp_geom_flux_Bmn_vac_z.png', dpi=res)
plt.close()
plt.figure(figsize=canvas)
plt.plot(test2_lhf['pre_r'], np.angle(test2_lhf['pre_Bz'], True), '-k',
         lw=thin, label='exact')
plt.plot(test2_lhf['post_r'], np.rad2deg(test2_lhf['post_Bz_arg']), '-r',
         lw=thin, label='flux')
plt.plot(test2_lhg['post_r'], np.rad2deg(test2_lhg['post_Bz_arg']), '--b',
         lw=thin, label='geom.')
plt.gca().legend(loc='upper right')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$\arg B_{mn}^{z}$ / °')
plt.title('Error in vacuum perturbation field\n'
          'flux symmetry angle vs. geometrical angle')
plt.savefig('test2_cmp_geom_flux_Bmn_vac_z_arg.png', dpi=res)
plt.close()

plt.figure(figsize=canvas)
plt.plot(test2_lhf['pre_r'], np.abs(test2_lhf['pre_Bz']), '-k', lw=thin,
         label='exact')
plt.plot(test2_lhf['post_r'], test2_lhf['post_Bz'], '-r', lw=thin,
         label='offset: 1/2 triangle')
plt.plot(test2_lff['post_r'], test2_lff['post_Bz'], '--b', lw=thin,
         label='offset: 1/4 triangle')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend(loc='upper left')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$|B_{mn}^{z}|$ / G')
plt.title('Error in vacuum perturbation field\n'
          'effects of poloidal sample point offset')
plt.savefig('test2_cmp_half_full_Bmn_vac_z.png', dpi=res)
plt.close()
plt.figure(figsize=canvas)
plt.plot(test2_lhf['pre_r'], np.angle(test2_lhf['pre_Bz'], True), '-k',
         lw=thin, label='exact')
plt.plot(test2_lhf['post_r'], np.rad2deg(test2_lhf['post_Bz_arg']), '-r',
         lw=thin, label='offset: 1/2 triangle')
plt.plot(test2_lff['post_r'], np.rad2deg(test2_lff['post_Bz_arg']), '-b',
         lw=thin, label='offset: 1/4 triangle')
plt.gca().legend(loc='upper right')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$\arg B_{mn}^{z}$ / °')
plt.title('Error in vacuum perturbation field\n'
          'effects of poloidal sample point offset')
plt.savefig('test2_cmp_half_full_Bmn_vac_z_arg.png', dpi=res)
plt.close()

plt.figure(figsize=canvas)
plt.plot(test2_hhf['pre_r'], np.abs(test2_hhf['pre_Bz']), '-k', lw=thin,
         label='exact')
plt.plot(test2_lhf['post_r'], test2_lhf['post_Bz'], '-r', lw=thin,
         label='$n_{pol} = 300$')
plt.plot(test2_hhf['post_r'], test2_hhf['post_Bz'], '--b', lw=thin,
         label='$n_{pol} = 500$')
plt.gca().get_yaxis().set_major_formatter(scifmt)
plt.gca().legend(loc='upper left')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$|B_{mn}^{z}|$ / G')
plt.title('Error in vacuum perturbation field\neffects of poloidal resolution')
plt.savefig('test2_cmp_polres_Bmn_vac_z.png', dpi=res)
plt.close()
plt.figure(figsize=canvas)
plt.plot(test2_hhf['pre_r'], np.angle(test2_hhf['pre_Bz'], True), '-k',
         lw=thin, label='exact')
plt.plot(test2_lhf['post_r'], np.rad2deg(test2_lhf['post_Bz_arg']), '-r',
         lw=thin, label='$n_{pol} = 300$')
plt.plot(test2_hhf['post_r'], np.rad2deg(test2_hhf['post_Bz_arg']), '-b',
         lw=thin, label='$n_{pol} = 500$')
plt.gca().legend(loc='upper right')
plt.xlabel(fslabel.r.value)
plt.ylabel(r'$\arg B_{mn}^{z}$ / °')
plt.title('Error in vacuum perturbation field\neffects of poloidal resolution')
plt.savefig('test2_cmp_polres_Bmn_vac_z_arg.png', dpi=res)
plt.close()
