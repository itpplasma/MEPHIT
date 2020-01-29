# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:13:45 2016

@author: Christopher Albert
"""

import numpy as np
from StringIO import StringIO

def loadmesh(infile):
    with open(infile + '.msh', 'r') as f:
        data = f.readlines()
        
        [NN,NT,NE] = np.genfromtxt(StringIO(data[0]),dtype=int)
        
        node = np.genfromtxt(StringIO(''.join(data[1:NN+1])),dtype=float)
        node = node[:,0:2]
        
        tri = np.genfromtxt(StringIO(''.join(data[NN+1:NN+NT+1])),dtype=int)
        tri = tri[:,0:3]
        
        edge = np.genfromtxt(StringIO(''.join(data[NN+NT+1:])),dtype=int)
        edge = edge[:,0:2]
        
    return [node,tri,edge]