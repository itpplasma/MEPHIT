"""Check the bordered shielding solve in the KiLCA cylinder (runs jump, ideal and layer of run_test.sh).

- ideal: the resonant flux harmonic in the resonant cell vanishes (< 1e-5 of the vacuum flux there);
- ideal and jump (pressure-jump sheet) give the same field away from the resonance, since both shield;
- layer (resistive-inertial layer, slow rotation): the field partly penetrates and the sheet current is
  smaller than in the ideal case.
"""
import sys

import h5py
import numpy as np


def cplx(d):
    return d["real"] + 1j * d["imag"]


work = sys.argv[1]
f = {case: h5py.File(f"{work}/{case}/mephit.h5", "r") for case in ("jump", "ideal", "layer")}
g = f["ideal"]
m, n_res = 4, 0
i = list(range(int(g["mesh/m_res_min"][()]), int(g["mesh/m_res_max"][()]) + 1)).index(m)
kf = int(g["mesh/res_ind"][()][i])
m_max = int(g["iter/Bmn/m_max"][()])
col = -int(np.sign(g["cache/fs/q"][()][-1])) * m + m_max


def residual(case):
    tot = cplx(f[case]["iter/Bmn/coeff_rad"][()])[kf - 1, col]
    vac = cplx(f[case]["iter/Bmn_vac/coeff_rad"][()])[kf - 1, col]
    return abs(tot / vac)


failures = []
res_jump, res_ideal, res_layer = residual("jump"), residual("ideal"), residual("layer")
amp_ideal = abs(cplx(f["ideal"]["iter/sheet_amplitude"][()])[i])
amp_layer = abs(cplx(f["layer"]["iter/sheet_amplitude"][()])[i])
print(f"resonant flux / vacuum: jump {res_jump:.3e}, ideal {res_ideal:.3e}, layer {res_layer:.3e}")
print(f"sheet amplitude: ideal {amp_ideal:.4e}, layer {amp_layer:.4e}")
if res_ideal > 1e-5:
    failures.append("ideal bordered solve leaves a resonant flux")
if not 0.05 < res_layer < 0.5:
    failures.append("slowly rotating resistive layer does not let the field partly penetrate")
if not amp_layer < 0.5 * amp_ideal:
    failures.append("resistive layer does not reduce the sheet current")
a = cplx(f["jump"]["iter/Bmn/coeff_rad"][()])
b = cplx(f["ideal"]["iter/Bmn/coeff_rad"][()])
far = np.abs(np.arange(a.shape[0]) - (kf - 1)) > 6
diff = np.max(np.abs(a - b)[far]) / np.max(np.abs(a))
print(f"field away from the resonance, ideal vs jump: max relative difference {diff:.2e}")
if diff > 2e-3:
    failures.append("ideal bordered field differs from the pressure-jump field away from the resonance")
if failures:
    sys.exit("FAIL: " + "; ".join(failures))
print("bordered shielding in the cylinder: OK")
