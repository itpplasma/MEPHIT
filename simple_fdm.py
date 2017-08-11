import numpy as np
import matplotlib.pyplot as plt

ns = 100
ds = 2*np.pi/(ns+1)

s = np.linspace(0,2*np.pi-ds,ns)
n = 1
m = 1

qmn = 1.0
pmn = 1.0/(1.0j*(m+n))
qn = qmn*np.exp(1.0j*m*s)
pn = pmn*np.exp(1.0j*m*s)

dsinv = 1.0/ds
diag1 = np.ones(ns)*(0.5j*n-dsinv)
diag2 = np.ones(ns-1)*(0.5j*n+dsinv)
A = np.diag(diag1) + np.diag(diag2,1)
A[-1,0] = (0.5j*n+dsinv)

qnnum = 0.5*(np.roll(qn,-1)+qn)
pnnum = np.linalg.solve(A,qnnum)

plt.figure()
plt.plot(pn)
plt.plot(qn)
plt.plot(pnnum)
plt.show()


