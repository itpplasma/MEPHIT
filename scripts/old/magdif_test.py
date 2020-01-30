#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 18:05:14 2017

@author: ert
"""

from ctypes import cdll, POINTER, c_double, c_int, byref
from _ctypes import dlclose
import numpy as np
import matplotlib.pyplot as plt
import os 
import subprocess
from pytchfort import fortran_module
#%%

try:
    dlclose(magdiflib._handle)
except:
    pass

#p = subprocess.Popen(['make'], stdout=subprocess.PIPE)
#print(p.stdout.read())

magdiflib = cdll.LoadLibrary('./libmagdif.so')

magdif = fortran_module(magdiflib, 'magdif')

magdif.init()
#magdif.test_pressure_profile()
#magdif.test_current_new()

#%%
nr_max = 74
psi = magdif.psi(c_double, nr_max)
q = magdif.q(c_double, nr_max)
dqdpsi = magdif.dqdpsi(c_double, nr_max)
psimin = magdif.psimin(c_double)

plt.plot(psi, q)