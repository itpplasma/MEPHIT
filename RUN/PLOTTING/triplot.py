#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:15:14 2016

@author: Christopher Albert
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from io import BytesIO
    
meshfile = sys.argv[1]
datafile = sys.argv[2]

with open(meshfile, 'r') as f:
    data = f.readlines()

print(data[0])
[NN,NT,NE] = np.loadtxt(StringIO(data[0]),dtype=int)

node = np.loadtxt(StringIO(''.join(data[1:NN+1])),dtype=float)
nlab = np.array(node[:,2],dtype=int)
node = node[:,0:2]

tri = np.loadtxt(StringIO(''.join(data[NN+1:NN+NT+1])),dtype=int)
tlab = tri[:,3]
tri = tri[:,0:3]

edge = np.loadtxt(StringIO(''.join(data[NN+NT+1:])),dtype=int)
elab = edge[:,2]
edge = edge[:,0:2]
                
data = np.genfromtxt(datafile)

print(sys.argv[3])

#tri_ind = np.array(data[:,0],dtype=int)-1
val = data[:,int(float(sys.argv[3]))]

plt.figure()
plt.tripcolor(node[:, 0], node[:, 1], tri-1, val, cmap='RdBu')
plt.colorbar()
plt.clim([-max(abs(val)), max(abs(val))])
plt.axis('equal')
if (len(sys.argv) > 4):
    try:
        plt.clim([float(sys.argv[4]), -float(sys.argv[4])])
        plt.show()
    except:
        plt.savefig(sys.argv[4], format = 'pdf')
else:
    plt.show()
