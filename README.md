# How to run

Prerequisites for running NEO-EQ are as follows.

-   current GNU/Linux environment (Bash, coreutils, getopt, ...)
-   CMake
-   current Fortran compiler (tested with `gfortran` and `ifort`)
-   LAPACK, SuiteSparse and SuperLU
-   FFTW3
-   GSL and FGSL
-   Python 3 including matplotlib
-   FreeFem++

## Initial build

    mkdir build
    cd build
    cmake ..
    make

## Coil geometry conversions

    build/bin/magdif convert <input_type> <output_type> <source_directory> [<target_directory>]

Converts coil geometry of `<input_type>` to `<output_type>`, reading appropriately named files from `<source_directory>` and writing to `<target_directory>`. If `<target_directory>` does not exist, it is created, and if it is omitted, it is set to the same value as `<source_directory>`. Possible values for `<input_type>` are `AUG`, `GPEC`, and `Nemov`. Possible values for `<output_type>` are `Fourier`, `GPEC`, and `Nemov`.

In order to use the configuration `vac_src = 2` (pre-computed Fourier modes), a symbolic link to the HDF5 file generated by the following command should be placed in the `<working_directory>`:

    build/bin/magdif convert AUG Fourier data <target_directory>
    ln -s <target_directory>/AUG_B_coils.h5 <working_directory>

## Setting up working directories

    build/bin/magdif init -c <config> -g <gfile> [-v <vacuum_field>] -w <convex_wall> <working_directory> ...

This copies the `<config>`, `<gfile>`, `<vacuum_field>`, and
`<convex_wall>` and other necessary files to each given
`<working_directory>`. The `<config>` file can be taken from a list of
templates in the `data` directory. `<convex_wall>` is usually either
`data/convexwall.asdex` or `data/convexwall.kilca`. For a KiLCA test
case, the `<vacuum_field>` should be omitted.

## Run

To generate the mesh and vacuum field data, and run the actual calculations including poloidal mode postprocessing, use

    build/bin/magdif run <working_directory> ...

## Plots

Default plots are generated by running

    build/bin/magdif plot <working_directory> ...

## Cleanup

    build/bin/magdif clean <working_directory> ...

This removes all files generated by the `run` and `plot` commands in each `<working_directory>`.
