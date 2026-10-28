"""Write the synthetic KiLCA vacuum file and flat profile files for the sheet test."""
import sys

import h5py
import numpy as np

out = sys.argv[1] if len(sys.argv) > 1 else "."
with h5py.File(f"{out}/kilca_vac.h5", "w") as f:
    g = f.create_group("output/postprocessor1")
    g["mode"] = np.array([[4.0], [2.0]])          # Fortran (1, 2): m, n
    r = np.linspace(0.0, 70.0, 8)
    g["r"] = r[np.newaxis, :]                      # Fortran (N, 1)
    g["Bz"] = np.vstack([np.full_like(r, 1.0), np.zeros_like(r)])  # Fortran (N, 2): Re, Im
x = np.linspace(0.0, 1.0, 101)
for name, y in (("n", 5e19 * (1.0 - 0.9 * x**2)), ("Te", 1e3 * (1.0 - 0.9 * x**2)),
                ("Ti", 1e3 * (1.0 - 0.9 * x**2)), ("Er", 0.0 * x)):
    np.savetxt(f"{out}/{name}.dat", np.column_stack([x, y]))
