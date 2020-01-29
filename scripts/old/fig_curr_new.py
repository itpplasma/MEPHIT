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

#tri_ind = np.array(data[:,0],dtype=int)-1

#plt.figure()
#
#plt.subplot(1,2,1)
#plt.tripcolor(node[:, 0], node[:, 1], tri-1, data1[:,0], cmap='seismic')
#plt.gca().text(170,-190,r'$R$')
#plt.gca().text(20,-5,r'$Z$')
##plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
##plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
#plt.xlim(xlims)
#plt.ylim(ylims)
#plt.clim(clims)
##plt.colorbar()
##plt.clim([float(sys.argv[4]),float(sys.argv[5])])
#plt.axis('equal')
#plt.subplot(1,2,2)
#plt.tripcolor(node[:, 0], node[:, 1], tri-1, data2[:,0], cmap='seismic')
#plt.gca().text(170,-190,r'$R$')
#plt.gca().text(20,-5,r'$Z$')
#plt.xlim(xlims)
#plt.ylim(ylims)
#plt.clim(clims)
#plt.axis('equal')
##plt.colorbar()
#plt.tight_layout()

#plt.savefig('fig_bpsi_rev.jpg', dpi=800)
 
if True:   
    xtri = np.mean(node[tri-1,0],1)
    ytri = np.mean(node[tri-1,1],1)
    xlims=[201,206]
    
    data1 = np.genfromtxt('currp1.dat')
    data2 = np.genfromtxt('currp2.dat')
    plt.figure()
    ylims=[3,13]
    clims=[0,5e10]
    plt.subplot(1,2,1)
    jmod = np.sqrt(data1[:,0]**2 + data1[:,2]**2)
    plt.quiver(node[:,0], node[:,1], data1[:,0], data1[:,2], jmod, 
               scale=5./2.*1.5e11, cmap='Blues')
    plt.xlabel('$R$')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.clim(clims)
    cb = plt.colorbar()
    cb.set_ticks(np.array([0,1,2,3,4,5])*1e10)
    #plt.gca().text(205,0,r'$R$')
    plt.gca().text(199.9,8,r'$Z$')
    plt.tight_layout()
    
    plt.subplot(1,2,2)
    clims=[0,2e9]
    jmod = np.sqrt(data2[:,0]**2 + data2[:,2]**2)
    plt.quiver(node[:,0], node[:,1], data2[:,0], data2[:,2], jmod, 
               scale=1.5e10, cmap='Reds')
    plt.xlabel('$R$')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.clim(clims)
    #plt.gca().text(205,0,r'$R$')
    plt.gca().text(199.9,8,r'$Z$')
    cb = plt.colorbar()
    cb.set_ticks(np.array([0,.4,.8,1.2,1.6,2])*1e9)
    plt.tight_layout()

#data1 = np.genfromtxt('curr1.dat')
#data2 = np.genfromtxt('curr2.dat')
#plt.figure()
#xlims=[201,206]
#ylims=[3,13]
##clims=[-1e10,5e10]
#plt.subplot(1,2,1)
#clims=[0,5e10]
#jmod = np.sqrt(data1[:,0]**2 + data1[:,2]**2)
#plt.quiver(xtri, ytri, data1[:,0], data1[:,2], jmod, scale=5e11)
#plt.xlim(xlims)
#plt.ylim(ylims)
#plt.clim(clims)
#plt.colorbar()
#plt.subplot(1,2,2)
#clims=[0,2e9]
#jmod = np.sqrt(data2[:,0]**2 + data2[:,2]**2)
#plt.quiver(xtri, ytri, data2[:,0], data2[:,2], jmod, scale=2e10)
#plt.xlim(xlims)
#plt.ylim(ylims)
#plt.clim(clims)
#plt.colorbar()
##plt.axis('equal')
plt.savefig('fig_curr_new.jpg', dpi=800)