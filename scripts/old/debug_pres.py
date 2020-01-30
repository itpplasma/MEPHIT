#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 19:16:51 2017

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('test_pressure_x.out')

n = data.shape[1]

#%%
plt.figure()
plt.plot(data[n-5:n-1,8],data[n-5:n-1,9], '-')
plt.figure()
plt.plot(data[1:5,8],data[1:5,9], '.-')

#%%
a = np.diff(data[:,8])
a = np.append(a,data[0,8]-data[-1,8])

b = np.diff(data[:,9])
b = np.append(b,data[0,9]-data[-1,9])
plt.figure()
plt.plot(np.arange(130), b[70:], '.-')
plt.plot(np.arange(70)+130, b[:70], '.-')

#%%
plt.figure(1)
plt.clf()
plt.plot(np.arange(70),np.sqrt(a[130:]**2+b[130:]**2))
plt.plot(np.arange(130)+70,np.sqrt(a[:130]**2+b[:130]**2))
