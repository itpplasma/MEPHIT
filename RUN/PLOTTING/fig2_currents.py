# -*- coding: utf-8 -*-
"""
Created on Sat Aug 06 08:48:24 2016

@author: chris
"""
import numpy as np
import matplotlib.pyplot as plt
from loadmesh import loadmesh

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
#plt.rc('font',**{'family':'sans-serif','size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[3.0, 2.13])
plt.rc('lines', markersize=3)
colors = plt.rcParams['axes.color_cycle']


[node,tri,edge] = loadmesh('inputformaxwell')
data = np.genfromtxt('curr_toplot_45.dat')

#tri_ind = np.array(data[:,0],dtype=int)-1

plt.figure()
xlims=[100,220]
ylims=[-150,150]
clims=[-1e10,1e10]

#j1 = np.sqrt(data[:,0]**2 + data[:,2]**2 + data[:,4]**2 + data[:,6]**2)
#j2 = np.sqrt(data[:,1]**2 + data[:,3]**2 + data[:,5]**2 + data[:,7]**2)

j1 = data[:,0]
j2 = data[:,1]

plt.subplot(1,2,1)
plt.tripcolor(node[:, 0], node[:, 1], tri-1, j1, cmap='seismic')
plt.gca().text(170,-190,r'$R$')
plt.gca().text(20,-5,r'$Z$')
#plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
#plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.xlim(xlims)
plt.ylim(ylims)
plt.clim(clims)
plt.axis('equal')
#cb = plt.colorbar()
#cb.ax.yaxis.set_offset_position('right') 
#plt.clim([float(sys.argv[4]),float(sys.argv[5])])
plt.subplot(1,2,2)
#clims=[-2e9,2e9]
plt.tripcolor(node[:, 0], node[:, 1], tri-1, j2, cmap='seismic')
plt.gca().text(170,-190,r'$R$')
plt.gca().text(20,-5,r'$Z$')
#plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
plt.xlim(xlims)
plt.ylim(ylims)
plt.clim(clims)
#plt.xticks(np.arange(50,350,50))
plt.axis('equal')
#cb = plt.colorbar()
#cb.ax.yaxis.set_offset_position('right') 
plt.tight_layout()

plt.savefig('fig_currents_rev.jpg', dpi=800)
    