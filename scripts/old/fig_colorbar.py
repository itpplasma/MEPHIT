# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 13:56:15 2016

@author: Christopher Albert
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[0.60, 2.13])
plt.rc('lines', markersize=3)
colors = plt.rcParams['axes.color_cycle']

# Make a figure and axes with dimensions as desired.
fig = plt.figure()
ax1 = fig.add_axes([0.25, 0.15, 0.22, 0.75])

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
norm = mpl.colors.Normalize(vmin=-1e10, vmax=1e10)


# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1,
                                cmap='seismic',
                                norm=norm,
                                orientation='vertical')
                                

plt.savefig('fig_colorbar_m1e10_1e10_rev.jpg', dpi=800)