# -*- coding: utf-8 -*-
"""
\brief     Solve magnetic differential equations in FD and FV method
\author    Christopher Albert
\date      2017-03-23
\copyright EUPL v.1.1
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

#%% connected chain magnetic DGL

R0 = 1.0e3
r = 1.0
absgradpsi = r/2.0 # psi = r**2/4
Bph = 1.0/R0
m = 7         # poloidal mode number
n = 2         # toroidal mode number
N = 100

th = np.linspace(0,2*np.pi,N+1)[:-1]
th[1:-1] = th[1:-1] + 2*np.pi/(N+1)*0.5*(np.random.rand(len(th)-2)-0.5)
R = R0 + r*np.cos(th)
Z = r*np.sin(th)

#plt.figure()
#plt.plot(R,Z)

# redundant, will be needed for more general case
th = np.arctan2(Z,R-R0) 
th[th<0] = th[th<0] + 2*np.pi

dR = np.roll(R,-1) - R 
dZ = np.roll(Z,-1) - Z
ds = np.sqrt(dR**2 + dZ**2)
dth = np.roll(th,-1) - th 
dth[-1] = dth[-1] + 2*np.pi


qmn = 1.0+0.5j
q = -1.0j*R/absgradpsi * qmn * np.exp(1.0j*m*th)


# analytical solution

#r = np.sqrt((R-R0)**2+Z**2)
Bth = 1.0/(2.0*R0) # psi = r**2/4
pmn = qmn/(1j*(m*Bth+n*Bph))
pn = pmn*np.exp(1j*m*th)

#%%
#plt.figure()
#plt.plot(th, np.real(pn))
#plt.plot(th, np.imag(pn), '--')

# forward difference scheme
A = sps.lil_matrix((N,N), dtype=complex)
for k in range(N):
    A[k,k] = n*R[k]*Bph/absgradpsi+1.0j/ds[k]
    if k != N-1:
        A[k,k+1] = -1.0j/ds[k]
A[N-1,0] = -1.0j/ds[N-1]

#q = (np.roll(q,-1) + q)/2.0
print(np.linalg.cond(A.toarray()))
p = spsl.spsolve(A.tocsr(),q)
#p = (np.roll(p,1) + p)/2.0

T = A.toarray()
Ti = np.linalg.inv(T)

#plt.figure()
#plt.plot(th, np.real(p))
#plt.plot(th, np.imag(p), '--')

plt.figure()
plt.plot(th, (np.real(pn)-np.real(p))/np.max(np.real(pn)))
plt.plot(th, (np.imag(pn)-np.imag(p))/np.max(np.imag(pn)), '--')

A = sps.lil_matrix((N,N), dtype=complex)

for k in range(N):
    A[k,k] = n*R[k]*Bph/absgradpsi
    if k != 0:
        A[k,k] = A[k,k] - 1.0j*(ds[k]-ds[k-1])/(ds[k]*ds[k-1])
        A[k,k-1] = 1.0j*ds[k]/(ds[k-1]*(ds[k-1]+ds[k]))
        
    if k != N-1:
        A[k,k+1] = -1.0j*ds[k-1]/(ds[k]*(ds[k-1]+ds[k]))
        
A[0,0] = A[0,0] - 1.0j*(ds[0]-ds[N-1])/(ds[0]*ds[N-1])
A[0,N-1] = 1.0j*ds[0]/(ds[N-1]*(ds[N-1]+ds[0]))
A[N-1,0] = -1.0j*ds[N-2]/(ds[N-1]*(ds[N-2]+ds[N-1]))

A = A.tocsr()
print(np.linalg.cond(A.toarray()))

p = spsl.spsolve(A,q)

#A2 = np.dot(A.toarray(),Ti)
#pp = np.linalg.solve(A2,q)
#p = np.linalg.solve(T,pp)

#A3 = np.dot(Ti,A.toarray())
#p = np.linalg.solve(A3,np.dot(Ti,q))


plt.figure()
plt.plot(th, np.real(p))
plt.plot(th, np.real(pn))
plt.plot(th, np.imag(p), '--')
plt.plot(th, np.imag(pn), '--')

plt.figure()
plt.plot(th, (np.real(pn)-np.real(p))/np.max(np.real(pn)))
plt.plot(th, (np.imag(pn)-np.imag(p))/np.max(np.imag(pn)), '--')
#%%
