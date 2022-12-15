# How to run

Prerequisites for running MEPHIT are as follows.

- current GNU/Linux environment (Bash, coreutils, getopt, ...)
- current Fortran compiler (tested with `gfortran` >= 9.2.0 and `ifort` 18.0.1)
- [Doxygen](https://doxygen.nl/) and [TeX Live](https://www.tug.org/texlive/) for documentation
- [CMake](https://cmake.org/)
- [LAPACK](https://www.netlib.org/lapack/)
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)
- [SuperLU](https://github.com/xiaoyeli/superlu)
- [GSL](https://www.gnu.org/software/gsl/) and [FGSL](https://github.com/reinh-bader/fgsl)
- [FFTW3](http://fftw.org/)
- [FreeFem++](https://github.com/FreeFem/FreeFem-sources)
- [HDF5](https://www.hdfgroup.org/downloads/hdf5)
- [NetCDF](https://github.com/Unidata/netcdf-fortran)
- Python 3 including packages listed in [`requirements.txt`](requirements.txt)

## Initial build

In the following sections, it is assumed that the environment variable `MEPHIT_DIR` points to the desired build directory. If libneo and FGSL are not in their default location (adjacent to MEPHIT and in the system default, respectively), the environment variables `LIBNEO_DIR` and `FGSL_PATH` need to be set as well. At ITPcp, you can put the following into your `~/.bashrc`:

```bash
export LIBNEO_DIR=/temp/AG-plasma/codes/libneo/build-master
export FGSL_PATH=/temp/AG-plasma/codes/contrib/fgsl-1.5.0/LIB
```

To build MEPHIT, run

```bash
mkdir $MEPHIT_DIR
cmake -B $MEPHIT_DIR
cd $MEPHIT_DIR
make
```

## Coil geometry conversions

```bash
$MEPHIT_DIR/scripts/mephit.bash convert <input_type> <output_type> <source_directory> [<target_directory>]
```

Converts coil geometry of `<input_type>` to `<output_type>`, reading appropriately named files from `<source_directory>` and writing to `<target_directory>`. If `<target_directory>` does not exist, it is created, and if it is omitted, it is set to the same value as `<source_directory>`. Possible values for `<input_type>` are `AUG`, `GPEC`, and `Nemov`. Possible values for `<output_type>` are `Fourier`, `GPEC`, and `Nemov`.

In order to use the default configuration `vac_src = 2` (pre-computed Fourier modes), run

```bash
$MEPHIT_DIR/scripts/mephit.bash convert AUG Fourier $MEPHIT_DIR/data
```

once to generate `$MEPHIT_DIR/data/AUG_B_coils.h5` which should take about 20 minutes and uses about 2 GB of disk space.

## Setting up working directories

```bash
$MEPHIT_DIR/scripts/mephit.bash init { -c | --config } <config> { -g | --g-eqdsk } <gfile> { -t | --type } { asdex | kilca } <working_directory> ...
```

This copies the `<config>`, `<gfile>`, and other necessary files to each given `<working_directory>`. The `<config>` file and some sample gfiles can be taken from a list of templates `mephit_<gfile>.in` in the `data` directory, e.g. [`data/mephit_g33353.2335.in`](data/mephit_g33353.2335.in). Currently, only geometry files for ASDEX Upgrade and KiLCA (large aspect ratio) are available.

The config file `mephit.in` needs do be adapated for each simulation:

- In the `arrays` namelist, array indices must be within the range of `m_res_min` and `m_res_max`.
- For ASDEX Upgrade, the coil currents must be set in `config%Ic`.
- For KiLCA, the requested poloidal mode must be set in `config%kilca_pol_mode`. For `config%kilca_scale_factor`, a value of `1000` yields reasonable results. Last but not least, the HDF5 output of the KiLCA vacuum run has to provied via `config%kilca_vac_output`.

## Simulations

Simulations consist of three phases:

1. meshing (includes calculation of the vacuum field)
2. iterations
3. analysis (poloidal modes, parallel currents)

Each phase can be run separately by specifying the corresponding command line switch; if none are given, all phases are run by default. If no working directory is given, the current directory is used; non-existing directories are skipped.

```bash
$MEPHIT_DIR/scripts/mephit.bash run [-m | --meshing] [-i | --iterations] [-a | --analysis] [<working_directory> ...]
```

## Tests

To run some default testing/debugging routines after the meshing phase, use

```bash
$MEPHIT_DIR/scripts/mephit.bash test [<working_directory> ...]
```

## Plots

Default plots are generated by running

```bash
$MEPHIT_DIR/scripts/mephit.bash plot [<working_directory> ...]
```

## Cleanup

```bash
$MEPHIT_DIR/scripts/mephit.bash clean [<working_directory> ...]
```

This removes all files generated by the `run` and `plot` commands in each `<working_directory>`.

# Documentation

To generate documentation, run the following commands.

```bash
cd DOC
latexmk magdif.tex
doxygen doxygen.conf
```

While both the LaTeX documentation of the method and the Doxygen documentation of the source code are horribly out of date, the `EXTRACT_ALL = YES` option of Doxygen is useful for generating call graphs.
