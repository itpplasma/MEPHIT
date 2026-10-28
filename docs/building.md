# Building MEPHIT

## Dependencies

- current GNU/Linux environment (Bash, coreutils, getopt, ...)
- current Fortran, C and C++ compilers (tested with `gfortran` >= 13.1.0)
- [CMake](https://cmake.org/) and Ninja
- [LAPACK](https://www.netlib.org/lapack/) and BLAS
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)
- [Boost](https://www.boost.org/) (1.74 or newer)
- [FreeFem++](https://github.com/FreeFem/FreeFem-sources)
- [zlib](https://zlib.net/) for Fortio's compressed NetCDF-4/HDF5 subset
- Python 3 with the packages in [requirements.txt](../requirements.txt) for
  plotting
- optional: [MFEM](https://mfem.org/)
- optional: [Doxygen](https://doxygen.nl/) and
  [TeX Live](https://www.tug.org/texlive/) for call graphs

The following are fetched and built by CMake:

- [libneo](https://github.com/itpplasma/libneo), which in turn fetches Fortio,
- [fortnum](https://github.com/lazy-fortran/fortnum), the numerical library
  (special functions, quadrature, FFT) that replaces GSL and FFTW, pinned by
  `FORTNUM_REF` in `CMakeLists.txt`,
- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html), unless
  `TRIANGLE_DIR` points to an existing build.

The scripts in [setup/](../setup/) install the system packages on Debian/Ubuntu
(`debian.sh`) and Arch Linux (`archlinux.sh`), and load modules on the MPCDF
Viper cluster (`mpcdf_viper.sh`).

## Build

```bash
make            # configure with the CMake preset "default" (Ninja) and build in build/
make test       # run ctest in build/
make install
make clean      # remove build/
```

`build/` then contains the executables and the driver script
`build/scripts/mephit.bash`.

### libneo version

The build pins a tested libneo commit (see `cmake/SetupCODE.cmake`). To fetch
a different branch, tag or commit, pass `-DLIBNEO_REF=<ref>` to CMake or
`make LIBNEO_REF=<ref>`. To use a local libneo source tree instead, pass
`-DLIBNEO_PATH=<dir>` or `make LIBNEO_PATH=<dir>`. The `Makefile` forwards
these only when given on the command line, not from the environment.

libneo's tools are built as part of MEPHIT; in particular `vacfield.x`, used
to precompute coil fields, is at `build/libneo/vacfield.x`.

### MFEM

MFEM support is off by default. Enable it with `-DWITH_MFEM=ON` when
configuring with CMake.

## Tests

`make test` runs the CTest suite: unit tests for mesh refinement and
quadrature, helical and resonant currents and layer widths, a test of the
FreeFem++ failure handling, and, where available, cylinder test cases for
sheet and edge currents (`test/sheet_cylinder`, `test/edge_cylinder`).

The driver option `mephit.bash run --test` runs MEPHIT's internal tests on a
working directory. These are expensive.

## Call graphs

The Doxygen documentation of the source code is out of date, but the
`EXTRACT_ALL = YES` option of Doxygen is useful for generating call graphs:

```bash
doxygen doxygen.conf
```
