#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 14:41:40 2018

@author: calbert
"""

import numpy as np
import codecs
from io import StringIO

# read msh file
def readmsh(infile):
    with codecs.open(infile, mode = 'r', encoding = 'utf-8') as f:
        data = f.readlines()
    
    # generate arrays for nodes, triangles, boundary edges    
    [NN,NT,NE] = np.genfromtxt(StringIO(data[0]),dtype=int)
    
    node = np.genfromtxt(StringIO(''.join(data[1:NN+1])),dtype=float)
    nlab = np.array(node[:,2],dtype=int)
    node = node[:,0:2]
    
    tri = np.genfromtxt(StringIO(''.join(data[NN+1:NN+NT+1])),dtype=int)
    tlab = tri[:,3]
    tri = tri[:,0:3]
    
    edge = np.genfromtxt(StringIO(''.join(data[NN+NT+1:])),dtype=int)
    elab = edge[:,2]
    edge = edge[:,0:2]
    
    return [node, tri, edge, nlab, tlab, elab]    
