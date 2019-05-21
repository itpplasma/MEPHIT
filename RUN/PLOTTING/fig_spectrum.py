# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:43:45 2016

@author: Christopher Albert
"""

import numpy as np
import matplotlib.pyplot as plt
from fig_common import *
#from exportfig import exportfig

data = np.loadtxt('/temp/ert/KINEQ/21_arnoldi_test_fine_30k/RUN/hpsi_vac.spectrum')
r = np.sqrt((np.max(data[:,0])-data[:,0])/np.max(data[:,0]))
q = data[:,1]
Nm = data[:,2:].shape[1]/2

m = np.arange(-Nm/2+1,Nm/2+1)
n = 2
[M,Q] = np.meshgrid(m,q)

spec_vac = data[:,2:(Nm+2)] + 1j* data[:,(Nm+2):]


data = np.loadtxt('/temp/ert/KINEQ/21_arnoldi_test_fine_30k/RUN/hpsi_relaxed_71.spectrum')
spec_plas = data[:,2:(Nm+2)] + 1j* data[:,(Nm+2):]

def plotq():
    plt.plot(3, r[np.argmin((q+1.5)**2)], 'w+')
    plt.plot(4, r[np.argmin((q+2)**2)], 'w+')
    plt.plot(5, r[np.argmin((q+2.5)**2)], 'w+')
    plt.plot(6, r[np.argmin((q+3)**2)], 'w+')
    plt.plot(7, r[np.argmin((q+3.5)**2)], 'w+')
    plt.plot(8, r[np.argmin((q+4)**2)], 'w+')
    plt.plot(9, r[np.argmin((q+4.5)**2)], 'w+')
    plt.plot(10, r[np.argmin((q+5)**2)], 'w+')
    plt.plot(11, r[np.argmin((q+5.5)**2)], 'w+')
    plt.plot(12, r[np.argmin((q+6)**2)], 'w+')

plt.figure()
plt.subplot(1,2,1)
levels = np.linspace(0,225,226)
plt.contourf(m,r,np.abs(spec_vac),levels=levels)
plotq()
#plt.colorbar()
plt.xlim([-15,15])
plt.ylim([0,1])
plt.gca().tick_params(axis='x', direction='out')
plt.gca().xaxis.tick_bottom()
plt.gca().tick_params(axis='y', direction='out')
plt.gca().yaxis.tick_left()
plt.xlabel(r'$m$')
plt.ylabel(r'$\rho_{\rm pol}$')
plt.subplot(1,2,2)
plt.contourf(m,r,np.abs(spec_vac+spec_plas),levels=levels)
plt.xlabel(r'$m$')

#plt.colorbar()
plotq()

plt.xlim([-15,15])
plt.ylim([0,1])
plt.gca().get_yaxis().set_ticks([])
plt.gca().tick_params(axis='x', direction='out')
plt.gca().xaxis.tick_bottom()

plt.tight_layout()

plt.savefig('fig_spectrum.jpg', dpi=800)