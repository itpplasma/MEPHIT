# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:19:54 2016

@author: Christopher Albert
"""

import matplotlib.pyplot as plt

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
#plt.rc('font',**{'family':'sans-serif','size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[3.2, 2.13])
plt.rc('lines', markersize=3)
colors = plt.rcParams['axes.color_cycle']