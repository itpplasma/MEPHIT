#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 16:50:11 2018

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fort.1000')
d = data[:,0] + 1j*data[:,1]

data = np.loadtxt('fort.1001')
du = data[:,0] + 1j*data[:,1]

data = np.loadtxt('fort.1002')
q = data[:,0] + 1j*data[:,1]

data = np.loadtxt('fort.1003')
a1 = data[:,0] + 1j*data[:,1]
data = np.loadtxt('fort.1004')
a2 = data[:,0] + 1j*data[:,1]
data = np.loadtxt('fort.1005')
b1 = data[:,0] + 1j*data[:,1]
data = np.loadtxt('fort.1006')
b2 = data[:,0] + 1j*data[:,1]
data = np.loadtxt('fort.1007')
x = data[:,0] + 1j*data[:,1]

A = np.diag(d)

plt.close()
plt.figure()
plt.subplot(2,2,1)
plt.plot(np.real(d),'.')
plt.plot(np.imag(d),'.')
plt.subplot(2,2,2)
plt.plot(np.real(du),'.')
plt.plot(np.imag(du),'.')
plt.subplot(2,2,3)
plt.plot(np.real(q),'.')
plt.plot(np.imag(q),'.')
plt.subplot(2,2,4)
plt.plot(np.real(x),'.')
plt.plot(np.imag(x),'.')