# KIM-MEPHIT Interface Design

## Summary

Integrate KIM (from KAMEL) as a third response current model in MEPHIT, alongside iMHD and KiLCA/FLR2. KIM pre-computes integral kernels externally; MEPHIT reads them and solves the electrostatic Poisson problem within its iteration loop.

Branch name (both repos): `feature/kim-mephit-interface`

## Data Flow (Per Iteration, Per Resonant Mode m)

1. MEPHIT extracts poloidal Fourier mode B_mn^psi from the 2D magnetic field
2. Transform to KIM's radial coordinate: B^r = (dr/dpsi_pol) * B^psi, interpolate onto KIM's xl_grid
3. Solve Poisson equation on KIM's grid: (L + 4pi K_rho_phi) Phi = -4pi K_rho_B * B_r with zero-misalignment-E BC
4. Compute parallel current: j_par = M^{-1} (K_j_phi * Phi + K_j_B * B_r)
5. Interpolate Phi_m and j_par/B0 back to MEPHIT's flux surface grid

## KAMEL Changes

### Write B_ref to KIM Output

`btor` from `KIM_SETUP` namelist is the reference magnetic field used in the effective radius definition r = sqrt(2 psi_tor / B_ref). Currently not written to output.

Changes:
- Text output: write `btor` to `setup/B_ref.dat` (single value)
- HDF5 output: write `btor` as scalar dataset `/setup/B_ref`

Files to modify:
- KIM I/O routines (the write/output section)

## MEPHIT Changes

### 1. Configuration (mephit_conf.f90)

New constant:
```fortran
integer, parameter :: currn_model_kim = 2
```

New namelist parameters in `&scalars`:
```fortran
character(len = 1024) :: kim_data_dir = ''
```

`kim_data_dir` points to base directory containing per-mode subdirectories (`kim_m6/`, `kim_m7/`, ...) or HDF5 files (`kim_m6.h5`, `kim_m7.h5`, ...).

### 2. New Module: mephit_kim.f90

#### Data Structures

```fortran
type :: kim_mode_t
    integer :: mpol                             ! poloidal mode number
    integer :: npts                             ! KIM xl_grid size
    integer :: npts_rg                          ! KIM rg_grid size
    real(dp), allocatable :: xl_grid(:)         ! KIM radial grid [cm]
    real(dp), allocatable :: rg_grid(:)         ! KIM profile grid [cm]
    real(dp), allocatable :: mass_matrix(:,:)   ! FEM mass matrix (real)
    real(dp), allocatable :: laplace(:,:)       ! Laplace operator in FEM (real)
    complex(dp), allocatable :: K_rho_phi(:,:)  ! charge-potential kernel
    complex(dp), allocatable :: K_rho_B(:,:)    ! charge-field kernel
    complex(dp), allocatable :: K_j_phi(:,:)    ! current-potential kernel
    complex(dp), allocatable :: K_j_B(:,:)      ! current-field kernel
    real(dp), allocatable :: Er0(:)             ! equilibrium E_r on rg_grid [statV/cm]
    real(dp), allocatable :: B0(:)              ! equilibrium B0 on rg_grid [G]
    real(dp), allocatable :: kp(:)              ! parallel wavenumber on rg_grid [1/cm]
    real(dp) :: B_ref                           ! reference field at axis [G]
    ! Pre-computed for solver:
    complex(dp), allocatable :: A_sparse(:)     ! CSC values of (L + 4pi K_rho_phi) with BC
    integer, allocatable :: A_row(:), A_col(:)  ! CSC indices
    integer :: A_nnz                            ! number of nonzeros
    real(dp), allocatable :: mass_matrix_inv(:,:) ! M^{-1}
    real(dp), allocatable :: Er0_xl(:)          ! E_r interpolated to xl_grid
    real(dp), allocatable :: B0_xl(:)           ! B0 interpolated to xl_grid
    real(dp), allocatable :: kp_xl(:)           ! k_par interpolated to xl_grid
    ! Coordinate mapping (MEPHIT flux surfaces -> KIM grid):
    real(dp), allocatable :: r_eff_mephit(:)    ! effective radius on MEPHIT grid
    real(dp), allocatable :: dr_dpsi_pol(:)     ! Jacobian on MEPHIT grid
end type

type :: kim_t
    type(kim_mode_t), allocatable :: modes(:)   ! one per resonant mode
    integer :: m_min, m_max                     ! resonant mode range
end type
```

#### File I/O

Two code paths:

**Text mode** (directory per mode):
```
kim_data_dir/
  kim_m6/
    kernel/K_rho_phi_re.dat, K_rho_phi_im.dat
    kernel/K_rho_B_re.dat, K_rho_B_im.dat
    kernel/K_j_phi_re.dat, K_j_phi_im.dat
    kernel/K_j_B_re.dat, K_j_B_im.dat
    kernel/mass_matrix.dat
    kernel/Laplace_in_FEM.dat
    grid/xl_xb.dat
    grid/rg_xb.dat
    backs/Er.dat, B0.dat, kp.dat  (on rg_grid)
    setup/B_ref.dat
  kim_m7/
    ...
```

Matrix files: space-separated, one row per line, Fortran scientific notation.
Complex kernels: separate real/imaginary files.
Profile files: 4-column format (x, Re, Im, |y|) or 2-column (x, value).

**HDF5 mode** (one file per mode):
```
kim_data_dir/
  kim_m6.h5
    /kernel/K_rho_phi, K_rho_B, K_j_phi, K_j_B  (complex 2D)
    /kernel/mass_matrix, Laplace_in_FEM           (real 2D)
    /grid/xl_xb, rg_xb                            (real 1D)
    /backs/Er, B0, kp                              (real 1D)
    /setup/B_ref                                   (real scalar)
  kim_m7.h5
    ...
```

Detection: check if `kim_m<m>/` directory exists (text mode) or `kim_m<m>.h5` file exists (HDF5 mode).

#### Initialization (called once)

```
subroutine kim_init(kim, mesh, equil)
```

For each resonant mode m (from mesh%m_res_min to mesh%m_res_max):
1. Read KIM data (kernels, grid, mass matrix, Laplace, profiles, B_ref)
2. Interpolate profiles (Er0, B0, kp) from rg_grid to xl_grid via linear interpolation
3. Invert mass matrix: M_inv = M^{-1}
4. Assemble system matrix: A = L + 4*pi * K_rho_phi
5. Apply BC structure to A (boundary rows set to identity)
6. Convert A to sparse CSC format (threshold 1e-14)
7. Build coordinate mapping: compute r_eff and dr/dpsi_pol on MEPHIT's flux surfaces
   - r_eff = sqrt(2 * psi_tor / B_ref)
   - dr/dpsi_pol = q / (B_ref * r_eff)

#### Poisson Solve (called each iteration, per mode)

```
subroutine kim_solve_poisson(mode, Br_kim, Phi_kim)
```

1. Assemble RHS: b = -4*pi * K_rho_B @ Br_kim
2. Apply BC: b(1) = -Phi_aligned(1), b(n) = -Phi_aligned(n)
   where Phi_aligned = c * Er0 * Br / (B0 * kp)
3. Solve A * Phi = b using UMFPACK (umf4zsym, umf4znum, umf4zsolr)
4. Return Phi_kim

#### Current Density (called each iteration, per mode)

```
subroutine kim_compute_current(mode, Phi_kim, Br_kim, jpar_over_B0)
```

1. Compute j_par = M_inv @ (K_j_phi @ Phi_kim + K_j_B @ Br_kim)
2. Divide by B0 (interpolated to xl_grid): jpar_over_B0 = j_par / B0

#### Top-Level Response Current (called from mephit_iter.F90)

```
subroutine compute_kim_current(Bn, jmnpar_over_Bmod, kim, Phi_mn, Phi_aligned_mn)
```

For each resonant mode m:
1. Extract B_mn^psi from Bn via RT0_poloidal_modes (same as KiLCA path)
2. Compute B^r = (dr/dpsi_pol) * B^psi on MEPHIT's flux surfaces
3. Recover B^r by multiplying with B0^phi
4. Interpolate B^r onto KIM's xl_grid
5. Call kim_solve_poisson -> Phi on KIM grid
6. Call kim_compute_current -> j_par/B0 on KIM grid
7. Interpolate Phi and j_par/B0 back to MEPHIT's flux surfaces
8. Store in jmnpar_over_Bmod and Phi_mn arrays

### 3. Integration in mephit_iter.F90

In `compute_currn` (~line 1092), add:
```fortran
case (currn_model_kim)
  call compute_kim_current(perteq%Bn, resonant_jmnpar_over_Bmod, &
    kim, perteq%Phi_mn, perteq%Phi_aligned_mn)
```

In initialization (~line 262), allocate Phi_mn and Phi_aligned_mn for currn_model_kim (same as KiLCA case).

### 4. CMakeLists.txt

Add `src/mephit_kim.f90` to the source file list.

## Coordinate Transformation Details

KIM's effective radius: r = sqrt(2 * psi_tor / B_ref)

Contravariant magnetic field transformation:
  B^r = (dr/dpsi_pol) * B^psi
  dr/dpsi_pol = q / (B_ref * r)

This uses MEPHIT's equilibrium data (psi_tor, q on flux surfaces) combined with KIM's B_ref.

Verified against KiLCA-FLR2 standalone (radial_variable == 1 mode):
- avR2nabpsi2 = R_tor^2
- bcovar_phi = B_tor * R_tor
- bpsi_over_bphi rescaled by 1 / (B_tor / R_tor)

## Boundary Conditions

Zero misalignment electric field at KIM grid boundaries:
  Phi(boundary) = -Phi_aligned(boundary)
  Phi_aligned = c * E_r * B_r / (B0 * k_parallel)

E_r, B0, k_parallel come from KIM's backs/ profiles (interpolated to xl_grid).
B_r at boundaries comes from interpolated MEPHIT field.

## Solver

- System matrix A = L + 4*pi * K_rho_phi is constant across iterations (pre-assembled)
- Uses complex UMFPACK (already available in MEPHIT via umf4_f77zwrapper.c)
- Symbolic factorization done once at init; numeric factorization + solve each iteration
- Sparsity threshold: 1e-14 (same as KIM)
- Current: j_par = M^{-1} @ (K_j_phi @ Phi + K_j_B @ B_r), using pre-inverted mass matrix

## Units

Both KIM and MEPHIT use CGS Gaussian units. No conversion needed.
- Magnetic field: Gauss
- Electric field: statV/cm
- Potential: statV
- Current density: statA/cm^2
- Lengths: cm
- Wavenumbers: 1/cm

KIM writes temperatures in eV to output files, but temperature is not directly used in MEPHIT's Poisson solve (it's encoded in the kernels).

## Files Modified

### KAMEL (branch: feature/kim-mephit-interface)
- KIM I/O module: add B_ref (= btor) to text and HDF5 output

### MEPHIT (branch: feature/kim-mephit-interface)
- `src/mephit_conf.f90`: add currn_model_kim, kim_data_dir
- `src/mephit_kim.f90`: new module (I/O, solver, coordinate mapping)
- `src/mephit_iter.F90`: add currn_model_kim case in compute_currn, allocations
- `CMakeLists.txt`: add mephit_kim.f90
