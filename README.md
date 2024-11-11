# MEPHIT

## Installation

Prerequisites from external sources for running MEPHIT are as follows.

- current GNU/Linux environment (Bash, coreutils, getopt, ...)
- current Fortran compiler (tested with `gfortran` >= 13.1.0)
- [CMake](https://cmake.org/)
- [libneo](https://github.com/itpplasma/libneo)
- [LAPACK](https://www.netlib.org/lapack/)
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)
- [SuperLU](https://github.com/xiaoyeli/superlu)
- [GSL](https://www.gnu.org/software/gsl/) and [FGSL](https://github.com/reinh-bader/fgsl)
- [FFTW3](http://fftw.org/)
- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)
- [Boost](https://www.boost.org/)
- [FreeFem++](https://github.com/FreeFem/FreeFem-sources)
- [MFEM](https://mfem.org/)
- [HDF5](https://www.hdfgroup.org/downloads/hdf5)
- [NetCDF](https://github.com/Unidata/netcdf-fortran)
- Python 3 including packages listed in [`requirements.txt`](requirements.txt) for plotting
- [Doxygen](https://doxygen.nl/) and [TeX Live](https://www.tug.org/texlive/) for call graphs

### Initial build

In the following sections, it is assumed that the environment variable `MEPHIT_DIR` points to the desired build directory and simulations are saved in `$MEPHIT_DIR/run`. If libneo, MFEM & FGSL are not in their default location (adjacent to MEPHIT), the environment variables `LIBNEO_DIR`, `MFEM_DIR`, and `FGSL_DIR` need to be set to the corresponding build directories as well. At ITPcp, you can refer to the `.gitlab-ci.yml` in [CODE](https://gitlab.tugraz.at/plasma/code).

To build MEPHIT, run:

```bash
mkdir $MEPHIT_DIR
cmake -B $MEPHIT_DIR
cd $MEPHIT_DIR
make -j
```

### Coil geometry

In order to use the default configuration in `mephit.in` (see below) for pre-computed Fourier modes for the vacuum field, i.e., `config%vac_src = 2`, you need to generate a coil file once to be used for `config%coil_file`. Create a namelist file for `coil_field`, like [libneo](https://github.com/itpplasma/libneo)'s `tools/vacfield_AUG.in`, and run `$LIBNEO_DIR/vacfield.x`. To generate the coil file for ASDEX Upgrade at ITPcp, run:

```bash
$LIBNEO_DIR/vacfield.x AUG 16 /proj/plasma/DATA/AUG/COILS/B{u,l}{1..8}n.asc Fourier vacfield_AUG.in $MEPHIT_DIR/run/AUG_B_coils.h5
```

This should take about 20 minutes and the resulting file uses about 2 GB of disk space. Note that the coils are consecutively numbered in the order given on the command line, starting from 1.

You can also use [GPEC](https://princetonuniversity.github.io/GPEC/developers.html) coil data. Assuming the environment variable `GPECHOME` is set to the GPEC directory:

```bash
$LIBNEO_DIR/vacfield.x GPEC 2 $GPECHOME/coil/aug_b{u,l}.dat Fourier vacfield_AUG.in $MEPHIT_DIR/run/AUG_B_coils.h5
```

## Setting up working directories

The general syntax is:

```text
$MEPHIT_DIR/scripts/mephit.bash init { -c | --config } <config> { -g | --g-eqdsk } <gfile> { -d | --device } { asdex | kilca | mastu } <working_directory> ...
```

This copies the `<config>`, `<gfile>`, and other necessary files to each given `<working_directory>`, specifically:

- `convexwall_$device.dat` containing the convex hull, where `$device` is one of those given above,
- [`field_divB0.inp`](data/field_divB0.inp), which is modified to point to `<gfile>` and `convexwall_$device.dat`, and
- `preload_for_SYNCH.inp`, which contains the configuration for the [field line integration](src/preload_for_SYNCH.f90) and doesn't need to be changed.

The `<config>` file and some sample gfiles can be taken from a list of templates `mephit_<gfile>.in` in the `data` directory, e.g. [`data/mephit_g33353_2900_EQH.in`](data/mephit_g33353_2900.in). Currently, only geometry files for ASDEX Upgrade, MAST Upgrade, and KiLCA (large aspect ratio) are available.

At ITPcp, you can reproduce the “standard” test case via:

```bash
$MEPHIT_DIR/scripts/mephit.bash init -c data/mephit_g33353_2900_EQH.in -g /proj/plasma/DATA/BALANCE/EQUI/33353/g33353.2900_EQH_MARKL -d asdex $MEPHIT_DIR/run/33353_2900_EQH
```

## Simulations

The config file `mephit.in` needs do be adapated for each simulation:

- In the `arrays` namelist, array indices must be within the range of `m_res_min` and `m_res_max`.
- For ASDEX Upgrade and MAST Upgrade, the data files containing coil currents and kinetic profiles must be supplied. See [the provided example](data/mephit_g33353_2900_EQH.in).
- For KiLCA, the requested poloidal mode must be set in `config%kilca_pol_mode`. For `config%kilca_scale_factor`, a value of `1000` yields reasonable results. Last but not least, the HDF5 output of the KiLCA vacuum run has to provied via `config%kilca_vac_output`. See [the provided example](data/mephit_g000001.0001_TCFP_hip.in).
- Some more options, mostly for debugging, are given in [`mephit_conf`](src/mephit_conf.f90).

The general syntax is:

```text
$MEPHIT_DIR/scripts/mephit.bash run [-m | --meshing] [-p | --preconditioner] [-i | --iterations] [--debug | --memcheck | --test] [<working_directory_or_file> ...]
```

Simulations consist of three phases:

1. meshing (including calculation of the vacuum field),
2. construction of preconditioner,
3. preconditioned iterations.

Each phase can be run separately by specifying the corresponding command line switch; if none are given, all phases are run by default. Use `--debug` to start a [GDB](https://www.gnu.org/software/gdb) session, `--memcheck` to run [Valgrind's Memcheck](https://valgrind.org/info/tools.html#memcheck), or `--test` to run internal tests (expensive, ignores phases). If no working directory is given, the current directory is used. If a file is given, it is used as an input namelist; if a directory is given, all files with pattern `mephit*.in` in it are used consecutively. Non-existing directories or files are skipped.

### Plots

Make the plotting routines available via:

```bash
python3 -m pip install -e .
```

To generate Jupyter notebooks for plotting, run:

```bash
jupytext -s $MEPHIT_DIR/scripts/{arnoldi,kinetic,magf,parcurr,tri}_plots.py
```

To update the files actually under version control, run:

```bash
jupytext -s $MEPHIT_DIR/scripts/*.ipynb
```

## Call graphs

While the Doxygen documentation of the source code is horribly out of date, the `EXTRACT_ALL = YES` option of Doxygen is useful for generating call graphs. To generate call graphs, run:

```bash
doxygen doxygen.conf
```
