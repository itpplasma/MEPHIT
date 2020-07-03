#%% Init
import numpy as np
import matplotlib.pyplot as plt
from libmc_fffi import libmc_efit, field_eq
from scipy.integrate import solve_ivp
from scipy.interpolate import SmoothBivariateSpline


r0 = 50.
p0 = 0.
z0 = 50.

Br = libmc_efit.new('double', 0.0)
Bp = libmc_efit.new('double', 0.0)
Bz = libmc_efit.new('double', 0.0)
dBrdR = libmc_efit.new('double', 0.0)
dBrdp = libmc_efit.new('double', 0.0)
dBrdZ = libmc_efit.new('double', 0.0)
dBpdR = libmc_efit.new('double', 0.0)
dBpdp = libmc_efit.new('double', 0.0)
dBpdZ = libmc_efit.new('double', 0.0)
dBzdR = libmc_efit.new('double', 0.0)
dBzdp = libmc_efit.new('double', 0.0)
dBzdZ = libmc_efit.new('double', 0.0)

tau = libmc_efit.new('double', 0.0)
vy = libmc_efit._ffi.new('double[5]')  # future: vz = np.zeros(5)
y = libmc_efit._ffi.new('double[5]')
y[0] = r0
y[1] = p0
y[2] = z0
y[3] = 0.7
y[4] = 0.3
vy[1] = 1e4
#
##%% Testing
#
#psi_pol = field_eq.psif - field_eq.psib
#
#print(f'Br before: {Br[0]}')
#print(f'psi_pol before: {psi_pol}')
#
#libmc_efit.field(r0, p0, z0, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
#
#psi_pol = field_eq.psif - field_eq.psib
#print(f'Br after: {Br[0]}')
#print(f'psi_pol after: {psi_pol}')

# check EFIT file with https://github.com/itpplasma/libneo/tree/master/matlab/EFIT

#%% spline psi

# grad A_(phi) = grad A_phi/R = grad (-psi/R)

nr = 20
nz=25
r0 = np.linspace(80,220,nr)
z0 = np.linspace(-100,100,nz)
R, Z = np.meshgrid(r0,z0)
aphi = np.zeros([nz,nr])
for k in range(len(r0)):
    for j in range(len(z0)):
        r = R[j,k]
        z = Z[j,k]
        libmc_efit.field(r, p0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
        psi_pol = field_eq.psif - field_eq.psib
        aphi[j,k] = psi_pol#-psi_pol/r

plt.pcolor(R, Z, aphi, cmap='plasma')
plt.colorbar()
plt.show()

spl = SmoothBivariateSpline(np.ravel(R),np.ravel(Z),np.ravel(aphi))


#psifun = spline(RR, ZZ, PSI)
# e.g. scipy.interpolate.SmoothBivariateSpline

# plot contours of psi






""" ############################## Constants ############################## """

c = 2.9979e10           # cm/s
qe = 4.8032e-10   # franklin ( = 3.336e-10C)
e_mass = 9.1094e-28     # g
p_mass = 1.6726e-24     # g 
ev = 1.6022e-12         # erg ( = 1e-7J)
am = 2                  # Atomic mass 2 of deuterium ions
Zb = 1                    # Atomic charge 1 of deuterium ions
tempi1 = 0.17549561306e4  # ion temperature
v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))        # Reference (thermal) velocity

""" ####################################################################### """


## Initial conditions
#R0 = 0.01  # cm
#phi0 = 0.0
#Z0 = 0.  # cm
#
#vR0 = 0.45*v0
#vph0 = 0.1*v0
#vZ0 = 0.45*v0
#
#libmc_efit.field(R0, phi0, Z0, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#psi_pol = field_eq.psif - field_eq.psib
#Aphi = -psi_pol/R0
#
m = am*p_mass
#pph0 = m*vph0 + qe/c*Aphi
#
#def ydot(t,y):
#    """ y = [R, Z, vR, vZ] """
#    
#    libmc_efit.field(y[0], y[1], y[2], Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,
#                     dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ) 
#    psi_pol = field_eq.psif - field_eq.psib
#    Aphi = -psi_pol/R0
#    dAphidr = spl.__call__(y[0], y[2], dx=1)
#    dAphidz = spl.__call__(y[0], y[2], dy=1)
#    
#    Rdot = y[2]
#    Zdot = y[3]
#    vRdot = -qe/c*y[3]*Bp[0] + 1/m*(pph0 - qe/c*Aphi)*(pph0/y[0] + qe/c*dAphidr)
#    vZdot = qe/c*y[2]*Bp[0] + 1/m*(pph0 - qe/c*Aphi)*qe/c*dAphidz
#    
#    return [Rdot, Zdot, vRdot, vZdot]
#
#y0 = np.array([R0, Z0, vR0, vZ0])
#integ = solve_ivp(ydot, (0,1e-3), y0, max_step=1e-5)

t=1000
tau=np.sqrt(2*tempi1/m)*t


libmc_efit.velo(tau,y,vy)










































