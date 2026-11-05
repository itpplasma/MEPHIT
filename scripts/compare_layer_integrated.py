#!/usr/bin/env python3
"""Compare damped MEPHIT layer currents and total resonant flux with GPEC.

Collect small observables where the large HDF5 files reside::

    python compare_layer_integrated.py collect samples.json NAME=RUN/mephit.h5 ...

Then compare against an explicitly selected GPEC case (no case-name inference)::

    python compare_layer_integrated.py compare samples.json GPEC/gpec_singfld_n3.out \
        comparison.json --coarse NAME --fine NAME

Fluxes are peak, area-normalized amplitudes in T. Ires is a one-sided parallel
current harmonic integrated over MEPHIT's configured physical window, in A.
For a GPEC helical sheet, K_phi dl/dtheta = i m kappa, K_parallel/K_phi = B/B_phi,
and peak I_res = 4 pi n kappa. Thus the corresponding one-sided parallel-current
magnitude is |I_res| m <B/|B_phi|> / (4 pi n). This compares a finite integration
window with an ideal sheet; it does not reconstruct omitted layer tails.

GPEC singflx is the sheet-only field, NOT its total field. Total MEPHIT flux /
GPEC singflx is a penetration scale, not a like-for-like field ratio. In GPEC's
gpout_singfld, singbwp is the unnormalized total flux (despite its NetCDF unit
label); divide by area_rational for the like-for-like total-field comparison.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.io import netcdf_file

STATA_PER_A = 2.99792458e9


def complex_data(dataset):
    value = dataset[()]
    if value.dtype.names:
        return value["real"] + 1j * value["imag"]
    return value


def pair(value):
    return [float(value.real), float(value.imag)] if np.isfinite(value) else None


def relative_changes(coarse, fine):
    if coarse is None or fine is None or abs(complex(*coarse)) == 0:
        return {"amplitude": None, "complex": None}
    a, b = complex(*coarse), complex(*fine)
    return {"amplitude": abs(abs(b) / abs(a) - 1), "complex": abs(b / a - 1)}


def collect_run(name, path):
    with h5py.File(path, "r") as f:
        n = int(f["mesh/n"][()])
        modes = f["mesh/res_modes"][()]
        indices = f["mesh/res_ind"][()]
        radial = f["cache/fs/rad"][()]
        psi = f["cache/fs_half/psi"][()]
        psi_res = f["mesh/psi_res"][()]
        sheet_modes = list(f["sheet/m"][()])
        area = f["sheet/area_norm"][()]
        total = complex_data(f["sheet/Phi_tot"])
        vacuum = complex_data(f["sheet/Phi_vac"])
        sheet = complex_data(f["sheet/singflx"])
        current = complex_data(f["iter/Ires"])
        current_low = int(f["iter/Ires"].attrs["lbounds"][0])
        bmn = f["iter/Bmn/coeff_rad"]
        coefficients = complex_data(bmn)
        m_low = int(bmn.attrs["lbounds"][0])
        sign_q = int(f["equil/cocos/sgn_q"][()])
        result = dict(name=name, path=str(path), n=n,
                      damp=bool(f["config/damp"][()]),
                      layer_model=int(f["config/res_layer_model"][()]),
                      layer_factor=float(f["config/res_layer_factor"][()]),
                      niter=int(f["iter/niter"][()]),
                      rel_err=float(f["iter/rel_err"][()]),
                      target_rel_err=float(f["config/iter_rel_err"][()]),
                      surfaces=[])
        for j, (m, kf) in enumerate(zip(modes, indices)):
            m, kf = int(m), int(kf)
            i = sheet_modes.index(m)
            start = int(f["mesh/kt_low"][kf - 1])
            stop = start + int(f["mesh/kt_max"][kf - 1])
            sample = f["cache/sample_polmodes_half"]
            br, bp, bz = [sample[k][start:stop] for k in ("B0_R", "B0_phi", "B0_Z")]
            geometry = float(np.mean(np.sqrt(br**2 + bp**2 + bz**2) / abs(bp)))
            # Interpolate the complex harmonic, rather than its magnitude.
            profile = coefficients[:, -sign_q * m - m_low]
            order = np.argsort(psi)
            exact = np.interp(psi_res[j], psi[order], profile[order])
            norm = 2e-4 / area[i]  # one-sided G cm^2 -> peak T / area in cm^2
            result["surfaces"].append(dict(
                m=m, q=m / n, h_cm=float(radial[kf] - radial[kf - 1]),
                refinement_width_cm=float(f["config/Delta_rad_res"][j]),
                add_fine=int(f["config/add_fine"][j]),
                layer_width_psi=float(f["mesh/Delta_psi_res_curr"][j]),
                area_norm_cm2=float(area[i]), geometry=geometry,
                total_T=pair(total[i] * norm), exact_surface_total_T=pair(exact * norm),
                vacuum_T=pair(vacuum[i] * norm), sheet_T=pair(sheet[i] * 1e-4),
                parallel_current_A=pair(current[m - current_low] / STATA_PER_A)))
    return result


def read_gpec(path, n):
    rows = {}
    for line in Path(path).read_text().splitlines():
        try:
            values = [float(x.replace("D", "E")) for x in line.split()]
        except ValueError:
            continue
        if len(values) < 9:
            continue
        q = values[0]
        m = round(n * q)
        if abs(m / n - q) > 0.000501:
            raise ValueError(f"GPEC q={q} is not a resonance of n={n}")
        if m in rows:
            raise ValueError(f"Duplicate GPEC surface m={m}")
        rows[m] = dict(q=q, sheet_T=pair(complex(*values[3:5])),
                       current_A=pair(complex(*values[5:7])),
                       total_flux_Wb=pair(complex(*values[7:9])))
    if not rows:
        raise ValueError(f"No resonant data in {path}")
    # gpout_singfld leaves singbwp unnormalized; use its own surface areas.
    profile = Path(path).with_name(f"gpec_profile_output_n{n}.nc")
    with netcdf_file(profile, "r", mmap=False) as f:
        qs = f.variables["q_rational"][:].copy()
        areas = f.variables["area_rational"][:].copy()
    for m, row in rows.items():
        i = int(np.argmin(abs(qs - m / n)))
        if abs(qs[i] - m / n) > 1e-6 or areas[i] <= 0:
            raise ValueError(f"Missing/invalid GPEC surface area for m={m}")
        row["area_m2"] = float(areas[i])
        row["total_T"] = pair(complex(*row["total_flux_Wb"]) / areas[i])
    return rows


def magnitude(value):
    return abs(complex(*value)) if value is not None else None


def compare(samples, gpec_path, coarse_name, fine_name):
    runs = {r["name"]: r for r in samples}
    coarse, fine = runs[coarse_name], runs[fine_name]
    if (coarse["n"], coarse["damp"], coarse["layer_model"], coarse["layer_factor"]) != (
            fine["n"], fine["damp"], fine["layer_model"], fine["layer_factor"]):
        raise ValueError("Convergence pair must use the same n and layer model")
    gpec = read_gpec(gpec_path, fine["n"])
    previous = {r["m"]: r for r in coarse["surfaces"]}
    rows = []
    for row in fine["surfaces"]:
        m = row["m"]
        ref, old = gpec[m], previous[m]
        gpar = magnitude(ref["current_A"]) * m * row["geometry"] / (4 * np.pi * fine["n"])
        current = magnitude(row["parallel_current_A"])
        total = magnitude(row["total_T"])
        sheet_flux = magnitude(ref["sheet_T"])
        flux_change = relative_changes(old["total_T"], row["total_T"])
        rows.append(dict(
            m=m, q=row["q"], h_cm=row["h_cm"],
            parallel_current_A=current, gpec_parallel_sheet_A=gpar,
            current_ratio=current / gpar if current is not None and gpar > 0 else None,
            total_T=total, gpec_total_T=magnitude(ref["total_T"]),
            penetration_scale=total / sheet_flux if sheet_flux > 0 else None,
            current_change=relative_changes(old["parallel_current_A"], row["parallel_current_A"]),
            flux_change=flux_change,
            locally_refined=row["h_cm"] < old["h_cm"] * 0.99,
            same_refinement_width=row["refinement_width_cm"] == old["refinement_width_cm"],
            flux_amplitude_within_20pct=(flux_change["amplitude"] <= 0.2
                                        if flux_change["amplitude"] is not None else None)))
    return dict(gpec=str(gpec_path), coarse=coarse_name, fine=fine_name,
                gpec_surfaces=gpec, runs=samples, comparison=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    collect = sub.add_parser("collect")
    collect.add_argument("output", type=Path)
    collect.add_argument("runs", nargs="+", metavar="NAME=PATH")
    analyze = sub.add_parser("compare")
    analyze.add_argument("samples", type=Path)
    analyze.add_argument("gpec", type=Path)
    analyze.add_argument("output", type=Path)
    analyze.add_argument("--coarse", required=True)
    analyze.add_argument("--fine", required=True)
    args = parser.parse_args()
    if args.command == "collect":
        result = [collect_run(*spec.split("=", 1)) for spec in args.runs]
    else:
        result = compare(json.loads(args.samples.read_text()), args.gpec,
                         args.coarse, args.fine)
        print("m   Ipar/A    GPEC_Ipar/A   I_ratio   total/T    total/GPEC_sheet"
              "  d|flux|    dflux_complex    d|Ipar|")
        for r in result["comparison"]:
            def fmt(x):
                return f"{x:.4g}" if x is not None else "undefined"
            print(r["m"], *(fmt(x) for x in (
                r["parallel_current_A"], r["gpec_parallel_sheet_A"], r["current_ratio"],
                r["total_T"], r["penetration_scale"], r["flux_change"]["amplitude"],
                r["flux_change"]["complex"], r["current_change"]["amplitude"])))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
