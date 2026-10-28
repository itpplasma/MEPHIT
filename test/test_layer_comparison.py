"""Independent harmonic/current oracles for the offline comparison script."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy.io import netcdf_file

SPEC = importlib.util.spec_from_file_location(
    "comparison", Path(__file__).resolve().parents[1] / "scripts/compare_layer_integrated.py")
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


class LayerComparisonTest(unittest.TestCase):
    def gpec_case(self, root, current=1.0):
        path = root / "gpec_singfld_n3.out"
        path.write_text("q psi spot real(singflx) imag(singflx) ...\n"
                        f"1.333 0.5 0.001 3D-4 4D-4 {current} 0 6D-6 8D-6\n")
        with netcdf_file(root / "gpec_profile_output_n3.nc", "w") as f:
            f.createDimension("surface", 1)
            f.createVariable("q_rational", "d", ("surface",))[:] = [4 / 3]
            f.createVariable("area_rational", "d", ("surface",))[:] = [5.0]
        return path

    def test_gpec_total_flux_divided_by_surface_area(self):
        with tempfile.TemporaryDirectory() as folder:
            row = comparison.read_gpec(self.gpec_case(Path(folder)), 3)[4]
        self.assertAlmostEqual(abs(complex(*row["total_T"])), 2e-6)
        self.assertAlmostEqual(abs(complex(*row["sheet_T"])), 5e-4)

    def test_parallel_current_against_helical_cylinder_surface_integral(self):
        # Independent surface-current oracle: K = grad(kappa) x normal.
        # theta,z coordinates of a cylinder, q=m/n, toroidal circumference 2 pi R.
        m, n, radius, major, kappa_peak = 4, 3, 1.25, 5.0, 12.5
        theta = np.linspace(0, 2 * np.pi, 16384, endpoint=False)
        k_phi = 1j * m / radius * kappa_peak * np.exp(1j * m * theta)
        k_theta = 1j * n / major * kappa_peak * np.exp(1j * m * theta)
        pitch = n * radius / (m * major)
        k_parallel = (k_phi + pitch * k_theta) / np.sqrt(1 + pitch**2)
        integrated = np.mean(k_parallel * np.exp(-1j * m * theta)) * radius / 2
        # GPEC peak I_res = 2 pi n kappa_peak, from its potential definition.
        gpec_current = 2 * np.pi * n * kappa_peak
        surface = dict(m=m, q=m / n, h_cm=0.01, refinement_width_cm=0.1,
                       geometry=np.sqrt(1 + pitch**2),
                       parallel_current_A=comparison.pair(integrated), total_T=[1e-5, 0])
        base = dict(n=n, damp=True, layer_model=1, layer_factor=1, surfaces=[surface])
        fine = dict(surface, h_cm=0.005)
        samples = [dict(base, name="coarse"), dict(base, name="fine", surfaces=[fine])]
        with tempfile.TemporaryDirectory() as folder:
            output = comparison.compare(samples, self.gpec_case(Path(folder), gpec_current),
                                        "coarse", "fine")["comparison"][0]
        self.assertAlmostEqual(output["current_ratio"], 1.0, places=13)
        self.assertTrue(output["locally_refined"])

    def test_collect_harmonic_units_fortran_indices_and_nonfinite_current(self):
        theta = np.linspace(0, 2 * np.pi, 1024, endpoint=False)
        # Flux-density harmonic of this real field, independent Fourier oracle.
        field = 4 * np.cos(4 * theta) - 4 * np.sin(4 * theta)
        coefficient = np.mean(field * np.exp(-4j * theta))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "manufactured.h5"
            with h5py.File(path, "w") as f:
                values = {
                    "mesh/n": 3, "mesh/res_modes": [4], "mesh/res_ind": [2],
                    "cache/fs/rad": [0, 0.9, 1, 1.2], "cache/fs_half/psi": [3, 2, 1],
                    "mesh/psi_res": [2], "sheet/m": [4], "sheet/area_norm": [10],
                    "sheet/Phi_tot": [coefficient], "sheet/Phi_vac": [coefficient],
                    "sheet/singflx": [2 * coefficient / 10],
                    "iter/Ires": [complex(float("nan"), 0)],
                    "equil/cocos/sgn_q": -1, "mesh/kt_low": [0, 2, 4],
                    "mesh/kt_max": [2, 2, 2], "mesh/Delta_psi_res_curr": [0.3],
                    "config/Delta_rad_res": [0.1], "config/add_fine": [4],
                    "config/damp": 1, "config/res_layer_model": 1,
                    "config/res_layer_factor": 1, "config/iter_rel_err": 1e-9,
                    "iter/rel_err": 1e-10, "iter/niter": 5,
                }
                for key, value in values.items():
                    f[key] = value
                f["iter/Ires"].attrs["lbounds"] = [4]
                profile = np.zeros((3, 13), dtype=complex)
                profile[:, 10] = coefficient
                f["iter/Bmn/coeff_rad"] = profile
                f["iter/Bmn/coeff_rad"].attrs["lbounds"] = [-6, 1]
                for key, value in (("B0_R", 0), ("B0_phi", -3), ("B0_Z", 4)):
                    f["cache/sample_polmodes_half/" + key] = np.full(6, value)
            row = comparison.collect_run("manufactured", path)["surfaces"][0]
        np.testing.assert_allclose(row["total_T"], [4e-5, 4e-5], rtol=1e-13)
        np.testing.assert_allclose(row["exact_surface_total_T"], row["total_T"], rtol=1e-13)
        self.assertAlmostEqual(row["h_cm"], 0.1)
        self.assertAlmostEqual(row["geometry"], 5 / 3)
        self.assertIsNone(row["parallel_current_A"])

    def test_phase_change_cannot_be_called_amplitude_convergence(self):
        change = comparison.relative_changes([1, 0], [0, 1])
        self.assertEqual(change["amplitude"], 0)
        self.assertAlmostEqual(change["complex"], np.sqrt(2))
        self.assertIsNone(comparison.relative_changes(None, [1, 0])["amplitude"])
        self.assertIsNone(comparison.relative_changes([0, 0], [1, 0])["amplitude"])


if __name__ == "__main__":
    unittest.main()
