# Using MEPHIT

The examples below assume that `MEPHIT_DIR` is the absolute path of the
`build` directory with the binaries and `MEPHIT_RUN_DIR` the absolute path of
a directory for simulations:

```bash
export MEPHIT_DIR="$(git rev-parse --show-toplevel)/build"
export MEPHIT_RUN_DIR="$(git rev-parse --show-toplevel)/run"
```

At ITPcp, the `.gitlab-ci.yml` in [CODE](https://gitlab.tugraz.at/plasma/code)
shows a complete setup.

## Coil geometry

The default configuration (`config%vac_src = 2`) uses precomputed Fourier
modes of the vacuum field, read from the file given in `config%coil_file`.
Generate this file once with libneo's `vacfield.x`, which is built at
`$MEPHIT_DIR/libneo/vacfield.x`. It needs a namelist file for `coil_field`,
such as libneo's `tools/vacfield_AUG.in`. With
`LIBNEO_DIR=$MEPHIT_DIR/libneo`, the ASDEX Upgrade coil file is generated at
ITPcp with

```bash
$LIBNEO_DIR/vacfield.x AUG 16 /proj/plasma/DATA/AUG/COILS/B{u,l}{1..8}n.asc Fourier vacfield_AUG.in $MEPHIT_RUN_DIR/AUG_B_coils.h5
```

This takes about 20 minutes and produces a file of about 2 GB. The coils are
numbered consecutively in the order given on the command line, starting from
1, and the currents in `config%currents_file` are expected in this order.
`config%Biot_savart_prefactor` must be set such that the result is in
statampere-turns with the speed of light set to 1. The same factor in
`vacfield_AUG.in` is only used with `sum` instead of `Fourier`, together with
a currents file; this corresponds to `config%vac_src = 0`, in which case
`config%currents_file` is ignored.

[GPEC](https://princetonuniversity.github.io/GPEC/developers.html) coil data
can be used as well. With `GPECHOME` pointing to the GPEC directory:

```bash
$LIBNEO_DIR/vacfield.x GPEC 2 $GPECHOME/coil/aug_b{u,l}.dat Fourier vacfield_AUG.in $MEPHIT_RUN_DIR/AUG_B_coils.h5
```

## Setting up working directories

```text
$MEPHIT_DIR/scripts/mephit.bash init { -c | --config } <config> { -g | --g-eqdsk } <gfile>
    { -d | --device } { asdex | iter | kilca | mastu } <working_directory> ...
```

This copies `<config>`, `<gfile>` and the following files into each
`<working_directory>`:

- `convexwall_$device.dat`, the convex hull of the first wall,
- [field_divB0.inp](../data/field_divB0.inp), modified to point to `<gfile>`
  and `convexwall_$device.dat`,
- `preload_for_SYNCH.inp`, the settings for the
  [field line integration](../src/preload_for_SYNCH.f90), which normally need
  no change.

Configuration templates `mephit_<gfile>.in` and sample g-files are in
[data/](../data/). Geometry files are available for ASDEX Upgrade, ITER, MAST
Upgrade and KiLCA (large aspect ratio). For ITER, the convex wall is the hull
of the first-wall limiter of the ITER TC24 g-file, used as a proxy for the
first wall. [data/mephit_iter_tc24.in](../data/mephit_iter_tc24.in) and
[data/vacfield_ITER_BL2.in](../data/vacfield_ITER_BL2.in) are the configuration
and coil-grid templates for the 27 `bl2` ELM coils (GPEC
`coil/iter_bl2{u,m,l}.dat`, one winding, currents in ampere-turns, hence
`Biot_Savart_prefactor = 0.1`).

At ITPcp, the standard ASDEX Upgrade test case is set up with

```bash
$MEPHIT_DIR/scripts/mephit.bash init -c data/mephit_g33353_2900_EQH.in -g /proj/plasma/DATA/BALANCE/EQUI/33353/g33353.2900_EQH_MARKL -d asdex $MEPHIT_RUN_DIR/33353_2900_EQH
```

## Configuration

Adapt `mephit.in` for each simulation:

- In the `arrays` namelist, array indices must lie between `m_res_min` and
  `m_res_max`, the minimal and maximal poloidal mode numbers with resonances.
- For ASDEX Upgrade and MAST Upgrade, supply the files with coil currents and
  kinetic profiles; see
  [data/mephit_g33353_2900_EQH.in](../data/mephit_g33353_2900_EQH.in).
- For KiLCA, set the poloidal mode in `config%kilca_pol_mode`;
  `config%kilca_scale_factor = 1000` gives reasonable results. Provide the HDF5
  output of the KiLCA vacuum run in `config%kilca_vac_output`; see
  [data/mephit_g000001.0001_TCFP_hip.in](../data/mephit_g000001.0001_TCFP_hip.in).
- For an equilibrium cut at an inner flux surface without an X-point there
  (e.g. the ITER TC24 g-file cut at s_pol = 0.997), set
  `config%truncate_at_sibry = .true.`. Otherwise the flux surfaces extend to
  the outermost closed surface inside the g-file box, not to `SIBRY`.
- The profile files `config%dens_file`, `temp_e_file`, `temp_i_file` and
  `E_r_file` are read in every run, also for ideal MHD (`currn_model = 0`),
  where they are not used.

All options, including debugging options, are documented in
[src/mephit_conf.f90](../src/mephit_conf.f90).

## Running simulations

```text
$MEPHIT_DIR/scripts/mephit.bash run [-m | --meshing] [-p | --preconditioner] [-i | --iterations]
    [-s | --sheet] [--debug | --memcheck | --test] [<working_directory_or_file> ...]
```

A simulation has three phases:

1. meshing, including the calculation of the vacuum field,
2. construction of the preconditioner,
3. preconditioned iterations.

Each phase can be run separately with its switch; without switches all three
run.

The optional phase `--sheet` (not run by default) post-processes a finished
run. It rebuilds each resonant sheet (shielding) current from the stored
pressure perturbation, solves for its field alone and writes the group `sheet`
of `mephit.h5`. Its `singflx` and `singcur` follow GPEC's `Phi_res` and
`I_res` (peak amplitude, flux normalised by the surface area) in gauss and
statampere; see `sheet_response` in [src/mephit_iter.F90](../src/mephit_iter.F90).
It appends to `mephit.h5`, so run it on a copy of a run directory that must
stay unchanged.

Three mutually exclusive flags help with testing and debugging:

- `--debug` starts a [GDB](https://www.gnu.org/software/gdb) session, using
  [scripts/mephit.gdb](../scripts/mephit.gdb) to go to the entry point of the
  main Fortran subroutine,
- `--memcheck` runs [Valgrind's Memcheck](https://valgrind.org/info/tools.html#memcheck),
- `--test` runs internal tests (expensive, ignores phases).

Without a working directory, the current directory is used. A file argument is
used as input namelist; for a directory, all files matching `mephit*.in` in it
are run in turn. Non-existing directories or files are skipped.

## Output files

For each `mephit*.in`, with `*` a possible suffix:

- `mephit*.h5`: all relevant numerical results,
- `mephit*.log`: the text output also shown on screen,
- `core_plasma*.msh`, `outer*.msh`, `maxwell*.msh`: meshes for the core
  plasma, for the surrounding elliptical boundary where boundary conditions are
  imposed, and for the union of both, for use with FreeFem++,
- `edgemap*.dat`: mapping of edge degrees of freedom between MEPHIT and
  FreeFem++,
- `fglut_dump*`: FreeFem++ graphics, viewable with `ffglut`.

## Plotting

Install the plotting package `mephit_plot` from the repository root:

```bash
python3 -m pip install --no-build-isolation -e .
```

Generate Jupyter notebooks from the plotting scripts with

```bash
jupytext -s scripts/{arnoldi,kinetic,magf,parcurr,tri}_plots.py
```

The input files for the plots are usually set near the top of each notebook.
To update the scripts under version control from edited notebooks, run

```bash
jupytext -s scripts/*.ipynb
```

These commands are relative to the repository root, not to `$MEPHIT_DIR`. The
plotting routines usually read `$MEPHIT_RUN_DIR/*/mephit.h5`.
