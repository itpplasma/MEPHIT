#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt
import sys

datafile = sys.argv[1]
data = np.genfromtxt(datafile)
rcoords = data[:, int(float(sys.argv[2]))]
zcoords = data[:, int(float(sys.argv[3]))]
rcomps = data[:, int(float(sys.argv[4]))]
zcomps = data[:, int(float(sys.argv[5]))]

plt.figure()
plt.quiver(rcoords, zcoords, rcomps, zcomps)
plt.axis('equal')
plt.show()
