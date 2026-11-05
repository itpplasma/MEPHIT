"""Check the MEPHIT sheet response in the KiLCA cylinder against the analytic field.

A divergence-free surface current K exp(i (m theta + k z)) at r = r_s in vacuum, with
B_r = 0 on a conducting wall at r = r_w (the outer FreeFem++ boundary), has (Gaussian)
    B_r(r) = k D x_s I_m'(x_<) [K_m'(x_>) - K_m'(x_w) I_m'(x_>) / I_m'(x_w)],
D = 4 pi K_z r_s / (i c m), x = k r, x_< = k min(r, r_s), x_> = k max(r, r_s).
In KiLCA geometry MEPHIT's coeff_rad is the physical B_r harmonic, and Itor_mn / r_s
is the K_z harmonic. The MEPHIT sheet has a finite width (the hat functions of the
nodes bracketing the resonance), which lowers B_r inside the layer by about 1.5 %.
"""
import sys

import h5py
import numpy as np
from scipy.special import ivp, kvp


def cplx(d):
    return d["real"] + 1j * d["imag"]


work = sys.argv[1]
f = h5py.File(f"{work}/mephit.h5", "r")
g = f["sheet"]
m = 4
i = list(g["m"][()]).index(m)
assert g["active"][()][i] == 1
R0, Z0 = float(f["mesh/R_O"][()]), float(f["mesh/Z_O"][()])
k = int(f["mesh/n"][()]) / R0
kf = int(f["mesh/res_ind"][()][i])
rad = f["cache/fs_half/rad"][()]
r_s = rad[kf - 1]
outer = np.loadtxt(f"{work}/outer_boundary.dat", skiprows=1)
r_w = np.hypot(outer[:, 0] - R0, outer[:, 1] - Z0).max()
Itor = cplx(g["Itor_mn"][()])[i]
Itor_lift = cplx(g["Itor_mn_lift"][()])[i]
D = 4 * np.pi * Itor / (1j * 2.99792458e10 * m)
x_s, x_w = k * r_s, k * r_w


def B_r(r):
    lo, hi = k * min(r, r_s), k * max(r, r_s)
    return k * D * x_s * ivp(m, lo) * (kvp(m, hi) - kvp(m, x_w) * ivp(m, hi) / ivp(m, x_w))


failures = []
ratio_res = cplx(g["Phi_mn"][()])[i, i] / B_r(r_s)
print(f"r_s = {r_s:.3f} cm, r_w = {r_w:.3f} cm, k r_s = {x_s:.4f}")
print(f"Itor_mn / Itor_mn_lift = {Itor / Itor_lift:.6f}")
print(f"B_r(r_s) MEPHIT / analytic = {ratio_res:.5f}")
if abs(Itor / Itor_lift - 1) > 1e-3:
    failures.append("toroidal sheet current from zero divergence differs from the lift")
if abs(ratio_res - 1) > 0.03:
    failures.append("resonant flux of the sheet differs from the analytic field")
profile = cplx(g["Bmn_rad"][()])[i]
outside = [j for j in range(len(rad)) if 2.0 < abs(rad[j] - r_s) < 8.0]
ratios = np.array([profile[j] / B_r(rad[j]) for j in outside])
print(f"B_r MEPHIT / analytic at {len(outside)} cells 2-8 cm from the sheet: "
      f"{abs(ratios).min():.4f} .. {abs(ratios).max():.4f}")
if np.abs(ratios - 1).max() > 0.02:
    failures.append("sheet field outside the layer differs from the analytic field")
singcur = cplx(g["singcur"][()])[i]
if abs(abs(singcur) / (4 * np.pi * int(f["mesh/n"][()]) * abs(Itor) / m) - 1) > 1e-12:
    failures.append("|singcur| is not 4 pi n |Itor_mn| / m")
if failures:
    sys.exit("FAIL: " + "; ".join(failures))
print("sheet response matches the analytic cylinder")
