#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 11:57:44 2018

@author: ert
"""

import numpy as np
import sys

nkpol = 200 #TODO: read from config file

datafile = sys.argv[1]
data = np.genfromtxt(datafile)

dataavg = data
dataavg[nkpol:-1:2,:] = 0.5*(data[nkpol:-1:2,:]+data[nkpol+1::2,:])
dataavg[nkpol+1::2,:] = dataavg[nkpol:-1:2,:]

np.savetxt(datafile+'.avg', dataavg)
