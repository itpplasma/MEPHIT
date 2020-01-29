# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 11:39:03 2016

@author: Christopher Albert
"""

import numpy as np
import matplotlib.pyplot as plt

filename = 'FIELD/17151.3800.AUGD.EQI.00.efit'
#filename = 'fix_a100_k10_q2_bn010_prof1'
#filename = '15MA_Q10_T38.geqdsk.txt' # from GACODE github

with open(filename) as f:
    lines = f.readlines()

shotstr = lines[0][:25]
timestr = lines[0][25:51]
idum,nw,nh = [int(k) for k in lines[0][51:].split()]

def split(a):
    aflat = a.replace('\n','')
    alist = [aflat[16*k:16*k+16] for k in range(int(len(aflat)/16))]
    return np.array(alist, dtype=float)

rdim, zdim, rcentr, rleft, zmid = np.array([lines[1][16*k:16*k+16] for k in range(5)],dtype=float)
rmaxis, zmaxis, simag1, sibry1, bcentr = np.array([lines[2][16*k:16*k+16] for k in range(5)],dtype=float)
current, simag2, dummy, rmaxis, dummy = np.array([lines[3][16*k:16*k+16] for k in range(5)],dtype=float)
zmaxis, dummy, sibry2, dummy, dummy = np.array([lines[4][16*k:16*k+16] for k in range(5)],dtype=float)

offset = 5
fpol = split(''.join(lines[offset:offset+int(nw/5)+1]))
offset = offset+int(nw/5)+1
pres = split(''.join(lines[offset:offset+int(nw/5)+1]))
offset = offset+int(nw/5)+1
ffprim = split(''.join(lines[offset:offset+int(nw/5)+1]))
offset = offset+int(nw/5)+1
pprime = split(''.join(lines[offset:offset+int(nw/5)+1]))
offset = offset+int(nw/5)+1
psi2d = split(''.join(lines[offset:offset+int(nw*nh/5)+1])).reshape([nw,nh])
offset = offset+int(nw*nh/5)+1
q = split(''.join(lines[offset:offset+int(nw/5)+1]))

psi = np.linspace(0,sibry1-simag1,nw)

plt.figure(1)
plt.contour(psi2d,30)
plt.colorbar()
plt.axis('equal')

plt.figure(2)
plt.plot(psi,q)
