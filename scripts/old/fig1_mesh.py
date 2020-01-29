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
plt.subplot(1,2,1)
plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=0.05)
plt.plot([100,250],[0,0],'.w')
plt.axis('equal')
plt.xticks(np.arange(50,350,50))
plt.gca().text(170,-190,r'$R$')
plt.gca().text(20,-5,r'$Z$')

[node,tri,edge] = loadmesh('inputformaxwell_ext')
plt.subplot(1,2,2)
plt.triplot(node[:, 0], node[:, 1], tri-1, linewidth=0.05)
#plt.plot([0,404],[0,0],'.w')
plt.axis('equal')
plt.xticks(np.arange(0,500,100))
plt.ylim([-300,300])
plt.gca().text(200,-410,r'$R$')
plt.gca().text(-100,-5,r'$Z$')

plt.tight_layout()

#plt.savefig('fig_mesh.jpg', dpi=800)
