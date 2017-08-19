# -*- coding: utf-8 -*-
"""
\brief     Solve magnetic differential equations in FD and FV method
\author    Christopher Albert
\date      2017-03-23
\copyright EUPL v.1.1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mpltri

#%% connected chain magnetic DGL

R0 = 1.0e3
def r(R,Z): return np.sqrt((R-R0)**2+Z**2)
def th(R,Z): the = np.array(np.arctan2(Z,R-R0)); the[the<0] += 2*np.pi; return the
def absgradpsi(R,Z): return r(R,Z)/2.0 # psi = r**2/4
Bph = 1.0/R0
m = 1         # poloidal mode number
n = 2         # toroidal mode number
N = 128

th1 = np.linspace(0,2*np.pi,N+1)[:-1]
th2 = np.linspace(0,2*np.pi,N+1)[:-1]
th1[1:-1] = th1[1:-1] + 2*np.pi/(N+1)*0.5*(np.random.rand(len(th1)-2)-0.5)
th2[1:-1] = th1[1:-1] + 2*np.pi/(N+1)*0.2*(np.random.rand(len(th1)-2)-0.5)
R1 = R0 + (1.0-np.pi/(N+1))*np.cos(th1)
Z1 = (1.0-np.pi/(N+1))*np.sin(th1)
R2 = R0 + (1.0+np.pi/(N+1))*np.cos(th2)
Z2 = (1.0+np.pi/(N+1))*np.sin(th2)

R = np.c_[[R1],[R2]].flatten()
Z = np.c_[[Z1],[Z2]].flatten()

tri = []
neigh = [] # neighbours
for k in range(N-1):
    tri.append([N+k+1,k,N+k])
    neigh.append([2*k+1,2*k-1])
    tri.append([k,N+k+1,k+1])
    neigh.append([2*k,2*k+2])

tri.append([N,N-1,2*N-1])    
neigh.append([2*N-1,2*N-3])
tri.append([N-1,N,0])
neigh.append([2*N-2,0])
tri = np.array(tri)
neigh = np.array(neigh)


tria = mpltri.Triangulation(R,Z,tri)

plt.figure()
plt.plot(R,Z,'x')
plt.axis('equal')

plt.figure()
plt.triplot(tria)
plt.axis('equal')


def h0(R,Z): return absgradpsi(R,Z)*np.array([-np.sin(th(R,Z)), np.cos(th(R,Z))])

[U,V] = h0(R,Z)
plt.figure()
plt.quiver(R,Z,U,V)
plt.show()

RE = [] # edge midpoint positions R
ZE = [] # edge midpoint positions Z
REV = [] # edge vector R
ZEV = [] # edge vector Z
ST = [] # triangle surfaces
h0perp = [] # flux of poloidal h0 through edges
for k in range(2*N):
    RE.append(np.array([R[tri[k][1]]+R[tri[k][0]],
                 R[tri[k][2]]+R[tri[k][1]],
                 R[tri[k][0]]+R[tri[k][2]]])/2.0)
    ZE.append(np.array([Z[tri[k][1]]+Z[tri[k][0]],
                 Z[tri[k][2]]+Z[tri[k][1]],
                 Z[tri[k][0]]+Z[tri[k][2]]])/2.0)
    REV.append(np.array([R[tri[k][1]]-R[tri[k][0]],
                 R[tri[k][2]]-R[tri[k][1]],
                 R[tri[k][0]]-R[tri[k][2]]]))
    ZEV.append(np.array([Z[tri[k][1]]-Z[tri[k][0]],
                 Z[tri[k][2]]-Z[tri[k][1]],
                 Z[tri[k][0]]-Z[tri[k][2]]]))
    ST.append(np.linalg.det([R[tri[0]],Z[tri[0]],np.ones(3)])/2.0)
    h0perp.append(np.array([
            np.cross(h0(RE[-1][0],ZE[-1][0]).T,[REV[-1][0],ZEV[-1][0]]),
            np.cross(h0(RE[-1][1],ZE[-1][1]).T,[REV[-1][1],ZEV[-1][1]]),
            np.cross(h0(RE[-1][2],ZE[-1][2]).T,[REV[-1][2],ZEV[-1][2]])
            ]))

RE = np.array(RE)
ZE = np.array(ZE)
REV = np.array(REV)
ZEV = np.array(ZEV)
le = np.sqrt(REV**2+ZEV**2) # edge lengths

RT = np.sum(R[tri],1)/3.0
ZT = np.sum(Z[tri],1)/3.0
    
plt.figure()
plt.quiver(RE,ZE,REV,ZEV)

# construct matrix
A = np.zeros([2*N,2*N], complex)
for k in range(2*N):
    A[k,k] = 1.0j*n*0.66 # TODO: add hphi
    if(h0perp[k][0] >= 0):
        A[k,k] += h0perp[k][0]/ST[k]
    else:
        A[k,neigh[k][0]] = h0perp[k][0]/ST[neigh[k][1]]
    if(h0perp[k][1] >= 0):
        A[k,k] += h0perp[k][1]/ST[k]
    else:
        A[k,neigh[k][1]] = h0perp[k][1]/ST[neigh[k][1]]
        
# analytical comparison
qmn = 1.0+0.5j
q = RT/absgradpsi(RT,ZT) * qmn * np.exp(1.0j*m*th(RT,ZT))

Bth = 1.0/(2.0*R0) # psi = r**2/4
pmn = qmn/(1j*(m*Bth+n*Bph))
pn = pmn*np.exp(1j*m*th(RT,ZT))
        
p = np.linalg.solve(A,q)

plt.figure()
plt.plot(th(RT,ZT), np.real(p))
plt.plot(th(RT,ZT), np.imag(p), '--')
plt.show()
