# -*- coding: utf-8 -*-
"""
Created on Sat Aug 06 08:48:24 2016

@author: chris
"""
import numpy as np
import matplotlib.pyplot as plt
from loadmesh import loadmesh
from exportfig import exportfig

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
#plt.rc('font',**{'family':'sans-serif','size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[3.2, 2.13])
plt.rc('lines', markersize=3)
colors = plt.rcParams['axes.color_cycle']


[node,tri,edge] = loadmesh('inputformaxwell')
data = np.genfromtxt('hpsi_vac.dat')
datap = np.genfromtxt('hpsi_relaxed_37.dat')
#datap = np.genfromtxt('/temp/ert/KINEQ/41_arnoldi_test_fine_30k_MPI/RUN/PLOTTING/hpsi_relaxed.dat')

#tri_ind = np.array(data[:,0],dtype=int)-1

plt.figure()
xlims=[100,220]
ylims=[-150,150]
clims=[-500,500]

plt.subplot(1,2,1)
plt.tripcolor(node[:, 0], node[:, 1], tri-1, data[:,0], cmap='seismic')
plt.gca().text(170,-190,r'$R$')
plt.gca().text(20,-5,r'$Z$')
#plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
#plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.xlim(xlims)
plt.ylim(ylims)
plt.clim(clims)
#plt.colorbar()
#plt.clim([float(sys.argv[4]),float(sys.argv[5])])
plt.axis('equal')
plt.subplot(1,2,2)
plt.tripcolor(node[:, 0], node[:, 1], tri-1, datap[:,0], cmap='seismic')
plt.gca().text(170,-190,r'$R$')
plt.gca().text(20,-5,r'$Z$')
plt.xlim(xlims)
plt.ylim(ylims)
plt.clim(clims)
plt.axis('equal')
#plt.colorbar()
plt.tight_layout()

#plt.savefig('fig_bpsi_rev.jpg', dpi=800)
    
