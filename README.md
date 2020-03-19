# How to run

Prerequisites for running NEO-EQ are as follows.

-   current GNU/Linux environment (Bash, coreutils, getopt, ...)
-   CMake
-   current Fortran compiler (tested with `gfortran` and `ifort`)
-   LAPACK, SuiteSparse and SuperLU
-   GSL and FGSL
-   Python 3 including matplotlib
-   FreeFem++

## Initial build

    mkdir build
    cd build
    cmake ..
    make

## Setting up working directories

    build/bin/magdif init -c <config> -g <gfile> [-v <vacuum_field>] -w <convex_wall> <working_directory> ...

This copies the `<config>`, `<gfile>`, `<vacuum_field>`, and
`<convex_wall>` and other necessary files to each given
`<working_directory>`. The `<config>` file can be taken from a list of
templates in the `data` directory. `<convex_wall>` is usually either
`data/convexwall.asdex` or `data/convexwall.kilca`. For a KiLCA test
case, the `<vacuum_field>` can be omitted.

## Preprocessing

When the configuration has been adapted as needed, the mesh and vacuum
field data are generated via

    build/bin/magdif prepare -c <config> <working_directory> ...

in each `<working_directory>`. `<config>` is assumed to be a path
relative to each respective `<working_directory>`.

## Run

To run the actual calculations including poloidal mode postprocessing,
use

    build/bin/magdif run -c <config> <working_directory> ...

Again, `<config>` is assumed to be a path relative to each respective
`<working_directory>`.

## Plots

Default plots are generated by running

    build/bin/magdif plot -c <config> <working_directory> ...

Again, `<config>` is assumed to be a path relative to each respective
`<working_directory>`.

## Cleanup

    build/bin/magdif clean <working_directory> ...

This removes all files generated by the `prepare`, `run`, and `plot`
commands in each `<working_directory>`.
