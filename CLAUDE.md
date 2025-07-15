# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MEPHIT is a magnetohydrodynamic (MHD) stability code for tokamak plasmas written in Fortran, C, and C++. It simulates plasma instabilities and magnetic field perturbations in fusion devices like ASDEX Upgrade, MAST Upgrade, and KiLCA.

## Build System and Development Commands

### Build Commands
```bash
# Full build (uses CMake with presets)
make

# Build only
make ninja

# Run tests
make test

# Install
make install

# Clean build directory
make clean
```

### Environment Variables Required
```bash
export MEPHIT_DIR="$(git rev-parse --show-toplevel)/build"
export MEPHIT_RUN_DIR="$(git rev-parse --show-toplevel)/run"
# Optional: export LIBNEO_DIR="/path/to/libneo/build" if not adjacent to MEPHIT
```

### Testing
- C tests: `cd build && ctest`
- Python tests: `pytest test/` (requires Python plotting dependencies)
- Internal tests: `$MEPHIT_DIR/scripts/mephit.bash run --test`

### Python Development
```bash
# Install plotting package in development mode
python3 -m pip install --no-build-isolation -e .

# Generate Jupyter notebooks from Python scripts
jupytext -s scripts/{arnoldi,kinetic,magf,parcurr,tri}_plots.py

# Sync notebooks back to Python scripts
jupytext -s scripts/*.ipynb
```

## Code Architecture

### Core Components
- **Configuration**: `src/mephit_conf.f90` - Main configuration module with parameters and types
- **Main executable**: `src/mephit_run.c` - C wrapper that manages child processes (FreeFem++ and Fortran)
- **Finite Element**: `src/mephit_fem.{c,cpp}` - Finite element mesh and matrix operations
- **Iteration**: `src/mephit_iter.F90` - Core MHD iteration algorithms
- **Mesh**: `src/mephit_mesh.F90` - Mesh generation and manipulation
- **Perturbation**: `src/mephit_pert.f90` - Magnetic perturbation calculations
- **Utilities**: `src/mephit_util.{c,f90}` - Utility functions and data structures

### Simulation Workflow
1. **Meshing**: Generate computational mesh and calculate vacuum field
2. **Preconditioner**: Construct preconditioner for iterative solver
3. **Iterations**: Run preconditioned iterations to solve MHD equations

### Key Data Structures
- `config_t`: Main configuration type containing all simulation parameters
- Various profile types for pressure, current, and safety factor
- Mesh data structures for finite element calculations

## Development Scripts

### Main Script
`scripts/mephit.bash` - Master script for initialization and running:
```bash
# Initialize working directory
$MEPHIT_DIR/scripts/mephit.bash init -c <config> -g <gfile> -d <device> <workdir>

# Run simulation (all phases)
$MEPHIT_DIR/scripts/mephit.bash run [<workdir>]

# Run specific phases
$MEPHIT_DIR/scripts/mephit.bash run -m  # meshing only
$MEPHIT_DIR/scripts/mephit.bash run -p  # preconditioner only
$MEPHIT_DIR/scripts/mephit.bash run -i  # iterations only

# Debug modes
$MEPHIT_DIR/scripts/mephit.bash run --debug     # GDB session
$MEPHIT_DIR/scripts/mephit.bash run --memcheck  # Valgrind
```

### Plotting Scripts
- `scripts/arnoldi_plots.py` - Eigenvalue analysis plots
- `scripts/kinetic_plots.py` - Kinetic effects visualization
- `scripts/magf_plots.py` - Magnetic field plots
- `scripts/parcurr_plots.py` - Parallel current plots
- `scripts/tri_plots.py` - Triangular mesh plots

## Configuration Files

### Input Files
- `mephit*.in` - Main namelist configuration files
- `field_divB0.inp` - Magnetic field configuration
- `convexwall_*.dat` - Device geometry (asdex, kilca, mastu)
- `preload_for_SYNCH.inp` - Field line integration parameters

### Output Files
- `mephit*.h5` - HDF5 files with numerical results
- `mephit*.log` - Text output logs
- `core_plasma*.msh`, `outer*.msh`, `maxwell*.msh` - FreeFem++ meshes
- `edgemap*.dat` - Edge mapping data
- `fglut_dump*` - Graphics for visualization

## Dependencies

### External Libraries
- libneo (adjacent directory or set LIBNEO_DIR)
- LAPACK/BLAS
- SuiteSparse (UMFPACK)
- GSL (GNU Scientific Library)
- FFTW3
- HDF5 (C and Fortran)
- NetCDF (Fortran)
- Triangle (mesh generation)
- Boost (C++)
- FreeFem++ (finite element)
- Optional: MFEM (configure with -DWITH_MFEM=ON)

### Python Dependencies
See `requirements.txt` for plotting dependencies including matplotlib, numpy, scipy, h5py, etc.

## Compiler Settings

### Fortran Flags
- gfortran: `-O2 -march=native -g -fno-realloc-lhs -fbacktrace -fcheck=all`
- ifort: `-fpp -g -assume norealloc_lhs -traceback -check all`
- Warning flags: `-Wall -Wextra -pedantic -fmax-errors=1`

### C/C++ Flags
- `-O2 -g -march=native -Wconversion -Wfloat-equal -Wshadow`

## Device Support

Supported tokamak devices:
- **asdex**: ASDEX Upgrade
- **kilca**: KiLCA (large aspect ratio)
- **mastu**: MAST Upgrade

Each device requires specific geometry files and configuration parameters.