# MEPHIT

MEPHIT computes the linear plasma response of an axisymmetric tokamak
equilibrium to external non-axisymmetric magnetic perturbations, such as the
fields of ELM-control (RMP) coils. It is used for ASDEX Upgrade, MAST Upgrade,
ITER and large-aspect-ratio (KiLCA) test cases.

Status: research code under active development; interfaces and input options
may change.

## Method and features

- Equilibrium from a g-file (EQDSK); a triangular mesh of the poloidal plane
  aligned with flux surfaces, optionally refined around the resonant surfaces
  of the selected toroidal harmonic.
- Vacuum perturbation from precomputed Fourier modes of coil fields (via
  libneo's `vacfield.x`), from GPEC, or from other sources.
- Plasma response currents from ideal MHD or from the KiLCA kinetic model.
- Self-consistent field by fixed-point iteration between response current and
  Maxwell solver (FreeFem++, optionally MFEM), accelerated by a preconditioner
  built from dominant eigenvalues found with the Arnoldi method.
- Optional post-processing of the resonant sheet currents and their fields.
- Results in HDF5 (`mephit.h5`), with Python plotting routines (`mephit_plot`).

## Building

Requirements: GNU/Linux, a recent GNU Fortran/C/C++ toolchain, CMake, Ninja,
BLAS/LAPACK, SuiteSparse, Boost and FreeFem++. libneo, fortnum and Triangle
are fetched and built automatically. Package lists for Debian/Ubuntu and Arch
Linux are in [setup/](setup/).

```bash
make
make test
```

The binaries and the driver script are placed in `build/`. Details, including
how to select the libneo version and enable MFEM, are in
[docs/building.md](docs/building.md).

## Usage

```bash
export MEPHIT_DIR="$(git rev-parse --show-toplevel)/build"
export MEPHIT_RUN_DIR="$(git rev-parse --show-toplevel)/run"

$MEPHIT_DIR/scripts/mephit.bash init -c data/mephit_g33353_2900_EQH.in \
    -g <gfile> -d asdex $MEPHIT_RUN_DIR/33353_2900_EQH
$MEPHIT_DIR/scripts/mephit.bash run $MEPHIT_RUN_DIR/33353_2900_EQH
```

`init` copies the configuration, the g-file and auxiliary files into a working
directory; `run` performs meshing, preconditioner construction and iterations.
Most configurations need a coil file generated once with libneo, and coil
currents and kinetic profiles. Configuration templates are in [data/](data/).
See [docs/usage.md](docs/usage.md).

## Documentation

- [docs/building.md](docs/building.md): dependencies, build options, tests,
  call graphs.
- [docs/usage.md](docs/usage.md): coil files, working directories,
  configuration, simulation phases, output files and plotting.
- [docs/contributing.md](docs/contributing.md): pull requests and merging.
- Configuration options are documented in [src/mephit_conf.f90](src/mephit_conf.f90).

## Related projects

- [libneo](https://github.com/itpplasma/libneo): field and utility library,
  including the coil-field tool `vacfield.x`.

## License

MIT License, see [LICENSE](LICENSE).
