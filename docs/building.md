# Building MEPHIT

## Dependencies

- current GNU/Linux environment (Bash, coreutils, getopt, ...)
- current Fortran, C and C++ compilers (tested with `gfortran` >= 13.1.0)
- [CMake](https://cmake.org/) >= 3.22 and Ninja
- [LAPACK](https://www.netlib.org/lapack/) and BLAS
- [AMD](https://github.com/DrTimothyAldenDavis/SuiteSparse/tree/dev/AMD) ordering headers/library (BSD); no external sparse factors
- [GSL](https://www.gnu.org/software/gsl/)
- [FFTW3](http://fftw.org/)
- [Boost](https://www.boost.org/) (1.74 or newer)
- [FreeFem++](https://github.com/FreeFem/FreeFem-sources), including its plugin SDK (`ff-c++` and headers)
- [zlib](https://zlib.net/) for Fortio's compressed NetCDF-4/HDF5 subset
- Python 3 with the packages in [requirements.txt](../requirements.txt) for
  plotting
- optional: [MFEM](https://mfem.org/)
- optional: [Doxygen](https://doxygen.nl/) and
  [TeX Live](https://www.tug.org/texlive/) for call graphs

The following are fetched and built by CMake:

- [libneo](https://github.com/itpplasma/libneo), which in turn fetches Fortio,
- [FortSparse](https://github.com/lazy-fortran/fortsparse) at a pinned revision,
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

### Sparse solvers (FortSparse)

Sparse solves use FortSparse's native general LU exclusively in the Fortran
FDM solver, the optional MFEM path and the real FreeFem Maxwell daemon.
The build pins FortSparse commit
`d5ab98a8b5d12550eb9aacf795eba97d782f1735`; `FORTSPARSE_PATH` selects a local
source tree and `FORTSPARSE_REF` overrides the fetched revision. SuperLU,
UMFPACK helpers and MUMPS are disabled, with no factorization fallback. Dense
BLAS/LAPACK remain dependencies. MFEM is built serially without external sparse
factorization libraries; prebuilt MFEM configurations with those libraries are
rejected. More than one MPI process is unsupported by the MFEM adapter.

The existing Fortran `sparse_solve` APIs and factor/solve/free lifecycle remain.
Method 3 (default) disables refinement; method 2 enables original-matrix
refinement. Historical `sparse_solve_suitesparse` names remain compatibility
aliases to native code. Both MFEM solves enable refinement, as their previous
solver did. AMD ordering replaces the previous MFEM METIS ordering;
this migration makes no performance claim. Matrix assembly, RHS definitions and
normal-operation solve conventions are unchanged. Factorization or solve errors
fail explicitly. The legacy Fortran interface still owns one retained factor at
a time and is not reentrant.

FortSparse owns the optional `fortsparse_native` FreeFem module, built from
the installed SDK compile definitions, while
linking FortSparse directly without the SDK tool's external solver link flags.
Unused UMFPACK/CHOLMOD SDK adapters are excluded from the plugin headers;
SDK ABI definitions are preserved. AMD supplies only the ordering permutation.
The launcher uses the exact configured serial `FREEFEM_EXECUTABLE`, loads the
native plugin and has no
MUMPS/UMFPACK fallback. Its retained real factors support normal/transpose
solves; the existing three Maxwell Richardson corrections are unchanged. The
installed FreeFem executable may itself contain unused external solver code;
linkage checks apply to the MEPHIT library and native plugin. The daemon must
successfully load this plugin. For standalone `electrostatics.edp` use the build's
`scripts` directory in `FF_LOADPATH` as the launcher does.

For full runtime independence, provide a FreeFem runtime/SDK built without
external sparse solvers, then configure both `FREEFEM_EXECUTABLE` and
`FREEFEM_CXX_EXECUTABLE` to that prefix's `bin/FreeFem++` and `bin/ff-c++`.
A stock interpreter with unused vendor libraries remains valid for native factor
execution. The verified minimal-runtime recipe and migration measurement scope
are described in [test/native_migration](../test/native_migration/README.md).

Focused checks after building are `ctest --test-dir build -R '^native_sparse_'`
plus `ctest --test-dir build -R '^native_freefem_'`; with `WITH_MFEM=ON`,
run `ctest --test-dir build -R '^native_mfem_adapter$'`.
These cover independent known solutions, retained factors, input/singular
failures, and actual shared-library linkage. They do not replace scientific
cylinder regression or qualify other configurations. The matched migration
measurement recipe is in [test/native_migration](../test/native_migration/README.md).

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
