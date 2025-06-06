# MEPHIT

## Installation

Prerequisites from external sources for running MEPHIT are as follows.

- current GNU/Linux environment (Bash, coreutils, getopt, ...)
- current Fortran compiler (tested with `gfortran` >= 13.1.0)
- [CMake](https://cmake.org/)
- [libneo](https://github.com/itpplasma/libneo)
- [LAPACK](https://www.netlib.org/lapack/)
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)
- [GSL](https://www.gnu.org/software/gsl/)
- [FFTW3](http://fftw.org/)
- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)
- [Boost](https://www.boost.org/)
- [FreeFem++](https://github.com/FreeFem/FreeFem-sources)
- [MFEM](https://mfem.org/) is optional
- [HDF5](https://www.hdfgroup.org/downloads/hdf5)
- [NetCDF](https://github.com/Unidata/netcdf-fortran)
- Python 3 including packages listed in [`requirements.txt`](requirements.txt) for plotting
- [Doxygen](https://doxygen.nl/) and [TeX Live](https://www.tug.org/texlive/) for call graphs

### Initial build

In the following sections, it is assumed that the environment variable `MEPHIT_DIR` points to the absolute path of the `build` directory containing the binaries and `MEPHIT_RUN_DIR` points to the absolute path of the `run` directory containing the simulations. If libneo is not in its default location (adjacent to MEPHIT), the environment variable `LIBNEO_DIR` needs to be set to the corresponding build directory as well. At ITPcp, you can refer to the `.gitlab-ci.yml` in [CODE](https://gitlab.tugraz.at/plasma/code).

To build MEPHIT, run:

```bash
make
export MEPHIT_DIR="$(git rev-parse --show-toplevel)/build"
export MEPHIT_RUN_DIR="$(git rev-parse --show-toplevel)/run"
```

### Coil geometry

In order to use the default configuration in `mephit.in` (see below) for pre-computed Fourier modes for the vacuum field, i.e., `config%vac_src = 2`, you need to generate a coil file once to be used for `config%coil_file`. **Create a namelist file for `coil_field`, like [libneo](https://github.com/itpplasma/libneo)'s `tools/vacfield_AUG.in`**, and run `$LIBNEO_DIR/vacfield.x`. To generate the coil file for ASDEX Upgrade at ITPcp, run:

```bash
$LIBNEO_DIR/vacfield.x AUG 16 /proj/plasma/DATA/AUG/COILS/B{u,l}{1..8}n.asc Fourier vacfield_AUG.in $MEPHIT_RUN_DIR/AUG_B_coils.h5
```

This should take about 20 minutes and the resulting file uses about 2 GB of disk space. The coils are consecutively numbered in the order given on the command line, starting from 1, and the currents in `config%currents_file` are expected in this order. The `config%Biot_savart_prefactor` needs to be set such that the result is given in statampere-turns with the speed of light set to 1. Note that this factor is also given in `vacfield_AUG.in`, but it is only used when running with `sum` instead of `Fourier`, together with the currents file. This corresponds to the setting `config%vac_src = 0` in `mephit.in`, which then ignores its `config%currents_file`.

You can also use [GPEC](https://princetonuniversity.github.io/GPEC/developers.html) coil data. Assuming the environment variable `GPECHOME` is set to the GPEC directory, the command above changes to:

```bash
$LIBNEO_DIR/vacfield.x GPEC 2 $GPECHOME/coil/aug_b{u,l}.dat Fourier vacfield_AUG.in $MEPHIT_RUN_DIR/AUG_B_coils.h5
```

## Setting up working directories

The general syntax is:

```text
$MEPHIT_DIR/scripts/mephit.bash init { -c | --config } <config> { -g | --g-eqdsk } <gfile>
    { -d | --device } { asdex | kilca | mastu } <working_directory> ...
```

This copies the `<config>`, `<gfile>`, and other necessary files to each given `<working_directory>`, specifically:

- `convexwall_$device.dat` containing the convex hull, where `$device` is one of those given above,
- [`field_divB0.inp`](data/field_divB0.inp), which is modified to point to `<gfile>` and `convexwall_$device.dat`, and
- `preload_for_SYNCH.inp`, which contains the configuration for the [field line integration](src/preload_for_SYNCH.f90) and doesn't need to be changed.

The `<config>` file and some sample gfiles can be taken from a list of templates `mephit_<gfile>.in` in the `data` directory, e.g. [`data/mephit_g33353_2900_EQH.in`](data/mephit_g33353_2900.in). Currently, only geometry files for ASDEX Upgrade, MAST Upgrade, and KiLCA (large aspect ratio) are available.

At ITPcp, you can reproduce the “standard” test case via:

```bash
$MEPHIT_DIR/scripts/mephit.bash init -c data/mephit_g33353_2900_EQH.in -g /proj/plasma/DATA/BALANCE/EQUI/33353/g33353.2900_EQH_MARKL -d asdex $MEPHIT_RUN_DIR/33353_2900_EQH
```

## Running simulations

The config file `mephit.in` needs do be adapated for each simulation:

- In the `arrays` namelist, array indices must be within the range of `m_res_min` and `m_res_max`, the minimal and maximal poloidal mode numbers for which resonances occur.
- For ASDEX Upgrade and MAST Upgrade, the data files containing coil currents and kinetic profiles must be supplied. See [the provided example](data/mephit_g33353_2900_EQH.in).
- For KiLCA, the requested poloidal mode must be set in `config%kilca_pol_mode`. For `config%kilca_scale_factor`, a value of `1000` yields reasonable results. Last but not least, the HDF5 output of the KiLCA vacuum run has to provied via `config%kilca_vac_output`. See [the provided example](data/mephit_g000001.0001_TCFP_hip.in).
- Some more options, mostly for debugging, are given in [`mephit_conf`](src/mephit_conf.f90).

The general syntax is:

```text
$MEPHIT_DIR/scripts/mephit.bash run [-m | --meshing] [-p | --preconditioner] [-i | --iterations]
    [--debug | --memcheck | --test] [<working_directory_or_file> ...]
```

Simulations consist of three phases:

1. meshing (including calculation of the vacuum field),
2. construction of preconditioner,
3. preconditioned iterations.

Each phase can be run separately by specifying the corresponding command line switch; if none are given, all phases are run by default. For testing and debugging, there are additionally three mutually exclusive flags:

- `--debug` starts a [GDB](https://www.gnu.org/software/gdb) session, using [a script](scripts/mephit.gdb) to navigate to the entry point of the main Fortran subroutine,
- `--memcheck` runs [Valgrind's Memcheck](https://valgrind.org/info/tools.html#memcheck),
- `--test` runs internal tests (expensive, ignores phases).

If no working directory is given, the current directory is used. If a file is given, it is used as an input namelist; if a directory is given, all files with pattern `mephit*.in` in it are used consecutively. Non-existing directories or files are skipped.

The following files are generated for `mephit*.in`, with `*` a possible suffix:

- `mephit*.h5` contains all relevant numerical results,
- `mephit*.log` contains the text output also displayed on the screen,
- `core_plasma*.msh`, `outer*.msh`, and `maxwell*.msh` contain the meshes for the core plasma, the surrounding elliptical boundary for imposition of boundary conditions, and the union of the two, for use with FreeFem++,
- `edgemap*.dat` contains the mapping of edge degrees of freedom between MEPHIT and FreeFem++,
- `fglut_dump*` contains graphics generated by FreeFem++ which can be viewed with the command `ffglut`.

### Plotting results

Make the plotting routines available via:

```bash
python3 -m pip install --no-build-isolation -e .
```

To generate Jupyter notebooks for plotting, run:

```bash
jupytext -s scripts/{arnoldi,kinetic,magf,parcurr,tri}_plots.py
```

The input files for the plots are usually set near the top of the Jupyter notebook.
To update the files actually under version control, run:

```bash
jupytext -s scripts/*.ipynb
```

*Note that the above commands are relative to the repository root, not to `$MEPHIT_DIR`.* The plotting routines themselves, however, usually read from the `$MEPHIT_RUN_DIR/*/mephit.h5` files.

## Call graphs

While the Doxygen documentation of the source code is horribly out of date, the `EXTRACT_ALL = YES` option of Doxygen is useful for generating call graphs. To generate call graphs, run:

```bash
doxygen doxygen.conf
```
