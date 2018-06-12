# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 22:01:13 2017

@author: chral
"""

import numpy as np
import matplotlib.pyplot as plt

files = ['along_line_rmp.1e5','along_line_rmp.5e5','along_line_rmp.1e5.constsrc','along_line_rmp.dat']
data = []

for file in files:
    data.append(np.genfromtxt(file))
    
plt.figure()
#plt.plot(data[0][:,0],data[0][:,3], '--')
#plt.plot(data[1][:,0],data[1][:,3], '--')
#plt.plot(data[2][:,0],data[2][:,3])
plt.plot(data[3][:,0],data[3][:,3])

plt.figure()
#plt.plot(data[0][:,0],data[0][:,4], '--')
#plt.plot(data[1][:,0],data[1][:,4], '--')
#plt.plot(data[2][:,0],data[2][:,4])
plt.plot(data[3][:,0],data[3][:,4])

plt.figure()
#plt.plot(data[0][:,0],data[0][:,5], '--')
#plt.plot(data[1][:,0],data[1][:,5], '--')
#plt.plot(data[2][:,0],data[2][:,5])
plt.plot(data[3][:,0],data[3][:,5])

plt.figure()
#plt.plot(data[0][:,0],data[0][:,6], '--')
#plt.plot(data[1][:,0],data[1][:,6], '--')
#plt.plot(data[2][:,0],data[2][:,6])
plt.plot(data[3][:,0],data[3][:,6])
