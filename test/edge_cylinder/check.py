"""Check the MEPHIT boundary sheet (edge_moved_current) in the KiLCA cylinder.

At a cut boundary r = a the equilibrium current density j0 drops to zero. A displacement
xi_r exp(i (m theta + n phi)) of the boundary moves the surface current K_mv = xi_r j0
past it; its surface divergence is closed by a field-aligned sheet kappa B0. In the
cylinder, with B0.grad = i k_par, k_par = m B_theta / a + n B_phi / R0, and ideal
B0.grad xi_r = B_r(a):
    xi_r = B_r(a) / (i k_par),
    kappa = -(m j_theta / a + n j_phi / R0) xi_r / k_par,
    K_theta = xi_r j_theta + kappa B_theta = xi_r n (j_theta B_phi - B_theta j_phi) / (R0 k_par),
    K_phi = -m R0 K_theta / (n a)   (zero surface divergence),
so the toroidal sheet current per radian of theta is I_tor = a K_phi. The field of that
sheet with B_r = 0 on the outer FreeFem++ boundary r_w is the Bessel solution of
test/sheet_cylinder/check.py.
"""
import sys

import h5py
import numpy as np
from scipy.special import ivp, kvp

C = 2.99792458e10


def cplx(d):
    return d["real"] + 1j * d["imag"]


work = sys.argv[1]
f = h5py.File(f"{work}/mephit.h5", "r")
e = f["sheet/edge"]
xi_mn = cplx(e["xi_mn"][()])
m_max = (xi_mn.size - 1) // 2
k0 = int(np.argmax(np.abs(xi_mn)))
m = k0 - m_max
n = int(f["mesh/n"][()])
R0, Z0 = float(f["mesh/R_O"][()]), float(f["mesh/Z_O"][()])
kp_low = f["mesh/kp_low"][()]
kp_max = f["mesh/kp_max"][()]
lo, nkp = int(kp_low[-1]), int(kp_max[-1])
R = f["mesh/node_R"][()][lo:lo + nkp]
Z = f["mesh/node_Z"][()][lo:lo + nkp]
a = np.hypot(R - R0, Z - Z0).mean()
# equilibrium on the boundary edge midpoints (poloidal edges kp_low + 1 .. kp_low + kp_max)
mf = f["cache/mid_fields"]
sl = slice(lo - 1, lo - 1 + nkp)
Rm = 0.5 * (R + np.roll(R, -1))
Zm = 0.5 * (Z + np.roll(Z, -1))
th = np.arctan2(Zm - Z0, Rm - R0)
B_th = np.mean(-np.sin(th) * mf["B0_R"][()][sl] + np.cos(th) * mf["B0_Z"][()][sl])
B_ph = np.mean(mf["B0_phi"][()][sl])
j_th = np.mean(-np.sin(th) * mf["j0_R"][()][sl] + np.cos(th) * mf["j0_Z"][()][sl])
j_ph = np.mean(mf["j0_phi"][()][sl])
k_par = m * B_th / a + n * B_ph / R0
# B_r(a) of the total field, extrapolated from the last two half-grid cells
coeff = cplx(f["iter/Bmn/coeff_rad"][()])
if coeff.shape[0] == 2 * m_max + 1:
    coeff = coeff.T
rad_half = f["cache/fs_half/rad"][()]
prof = coeff[-len(rad_half):, k0]
B_r_a = prof[-1] + (prof[-1] - prof[-2]) * (a - rad_half[-1]) / (rad_half[-1] - rad_half[-2])
xi_an = B_r_a / (1j * k_par)
K_th = xi_an * n * (j_th * B_ph - B_th * j_ph) / (R0 * k_par)
Itor_an = -a * m * R0 * K_th / (n * a)
Itor = cplx(e["Itor_mn"][()])[k0]

failures = []
print(f"m = {m}, n = {n}, a = {a:.3f} cm, R0 = {R0:.4g} cm, k_par = {k_par:.4e} G/cm")
print(f"j_theta B_phi - B_theta j_phi = {j_th * B_ph - B_th * j_ph:.4e} statA G cm^-2")
print(f"xi_r MEPHIT / analytic = {xi_mn[k0] / xi_an:.5f}")
print(f"I_tor MEPHIT / analytic = {Itor / Itor_an:.5f}")
if abs(xi_mn[k0] / xi_an - 1) > 0.02:
    failures.append("boundary displacement differs from B_r / (i k_par)")
if abs(Itor / Itor_an - 1) > 0.03:
    failures.append("toroidal boundary sheet current differs from the analytic sheet")
# field of the sheet alone against the Bessel solution with MEPHIT's own sheet current
outer = np.loadtxt(f"{work}/outer_boundary.dat", skiprows=1)
r_w = np.hypot(outer[:, 0] - R0, outer[:, 1] - Z0).max()
k = n / R0
# z = -phi in the right-handed (r, theta, z) of the Bessel solution, hence -m
D = 4 * np.pi * Itor / (1j * C * -m)


def B_r(r):
    lo_, hi_ = k * min(r, a), k * max(r, a)
    return k * D * k * a * ivp(m, lo_) * (kvp(m, hi_) - kvp(m, k * r_w) * ivp(m, hi_) / ivp(m, k * r_w))


Bmn = cplx(e["Bmn_rad"][()])
if Bmn.shape[0] != 2 * m_max + 1:
    Bmn = Bmn.T
profile = Bmn[k0]
inside = [j for j in range(len(rad_half)) if 2.0 < a - rad_half[j] < 8.0]
ratios = np.array([profile[j] / B_r(rad_half[j]) for j in inside])
print(f"B_r MEPHIT / analytic at {len(inside)} cells 2-8 cm inside the boundary: "
      f"{ratios.real.min():.4f} .. {ratios.real.max():.4f}, max |imag| {np.abs(ratios.imag).max():.1e}")
# The sheet sits on the interface to the unstructured vacuum triangulation, which is
# coarser than the core mesh and not rotationally symmetric: measured 1.8 % low and 1 deg
# in phase, uniformly inside (the resonant sheet deep in the core mesh agrees to 0.5 %).
if len(inside) == 0 or np.abs(ratios - 1).max() > 0.03:
    failures.append("boundary sheet field differs from the analytic field")
if failures:
    sys.exit("FAIL: " + "; ".join(failures))
print("boundary sheet matches the analytic cylinder")
