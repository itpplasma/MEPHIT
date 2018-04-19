module magdif
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       mesh_element_rmp, bphicovar, knot
  use for_macrostep, only : t_min, d_min
  use sparse_mod, only: remap_rc, sparse_solve, sparse_matmul

  implicit none

  integer log_level
  logical :: log_err, log_warn, log_info, log_debug ! specify log levels
  logical :: nonres = .false.  !< use non-resonant test case

  character(len=1024) :: point_file   !< input data file for mesh points
  character(len=1024) :: tri_file     !< input data file for triangles and edges
  character(len=1024) :: Bnflux_file  !< input data file for magnetic field perturbation
  character(len=1024) :: hpsi_file    !< input data file for \f$ h_{n}^{\psi} \f$
  character(len=1024) :: config_file  !< input config file for namelist settings
  character(len=1024) :: presn_file   !< output data file for pressure perturbation

  integer  :: n               !< harmonic index of perturbation
  integer  :: nkpol           !< number of knots per poloidal loop
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density

  namelist / settings / log_level, nonres, point_file, tri_file, Bnflux_file, hpsi_file, presn_file, &
       n, nkpol, nflux, ti0, di0  !< namelist for input parameters

  integer, parameter :: logfile = 6             !< log to stdout, TODO: make this configurable

  real(dp), allocatable :: pres0(:)             !< unperturbed pressure \f$ p_{0} \f$ in dyn cm^-1
  real(dp), allocatable :: dpres0_dpsi(:)       !< derivative of unperturbed pressure w.r.t. flux surface label, \f$ p_{0}'(\psi) \f$
  real(dp), allocatable :: dens(:)              !< density \f$ \frac{N}{V} \f$ on flux surface in cm^-3
  real(dp), allocatable :: temp(:)              !< temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV
  real(dp), allocatable :: psi(:)               !< flux surface label \f$ \psi \f$
  complex(dp), allocatable :: presn(:)          !< pressure perturbation \f$ p_{n} \f$ in each mesh point
  complex(dp), allocatable :: currn(:,:)        !< edge currents \f$ R \vec{j}_{n} \cdot \vec{n} \f$ weighted by \f$ R \f$
  complex(dp), allocatable :: Bnflux(:,:)       !< edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$ weighted by \f$ R \f$
  complex(dp), allocatable :: Bnphi(:)          !< physical toroidal component of magnetic perturbation, \f$ B_{n \phi} \f$

  real(dp) :: psimin  !< minimum flux surface label, located at the separatrix
  real(dp) :: psimax  !< maximum flux surface label, located at the magnetic axis
  real(dp), parameter :: R0 = 172.74467899999999d0  !< distance of magnetic axis from center, \f$ R_{0} \f$
  real(dp), parameter :: ideal_gas_factor = 1.6021766208d-12  !< unit conversion factor in ideal gas law \f$ p = \frac{N}{V} k_{\mathrm{B}} T \f$
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< imaginary unit in double precision

contains

  !> Initialize magdif module
  subroutine magdif_init
    call read_config(config_file)
    call read_mesh
    call read_bnflux
    call read_hpsi ! TODO: get rid of this due to redundancy with bnflux
    call init_flux_variables
    !call init_safety_factor ! TODO
    if (log_info) write(logfile, *) 'magdif initialized'
  end subroutine magdif_init

  !> Final cleanup of magdif module
  subroutine magdif_cleanup
    deallocate(pres0)
    deallocate(dpres0_dpsi)
    deallocate(dens)
    deallocate(temp)
    deallocate(psi)
    deallocate(presn)
    deallocate(currn)
    deallocate(Bnflux)
    deallocate(Bnphi)
    deallocate(mesh_point)
    deallocate(mesh_element)
    deallocate(mesh_element_rmp)
    if (log_info) write(logfile, *) 'magdif cleanup finished'
  end subroutine magdif_cleanup

  !> Read configuration file for magdif
  !! @param config_filename file name of config file
  subroutine read_config(config_filename)
    character(len = *) :: config_filename

    open(1, file = config_filename)
    read(1, nml = settings)
    close(1)

    log_err = .false.
    log_warn = .false.
    log_info = .false.
    log_debug = .false.
    if (log_level > 0) log_err = .true.
    if (log_level > 1) log_warn = .true.
    if (log_level > 2) log_info = .true.
    if (log_level > 3) log_debug = .true.
  end subroutine read_config

  !> Read mesh points and triangles
  subroutine read_mesh
    open(1, file = point_file, form = 'unformatted')
    read(1) npoint
    if (log_info) write(logfile, *) 'npoint = ', npoint
    allocate(mesh_point(npoint))
    allocate(presn(npoint))
    read(1) mesh_point
    close(1)

    open(1, file = tri_file, form = 'unformatted')
    read(1) ntri
    if (log_info) write(logfile, *) 'ntri   = ', ntri
    allocate(mesh_element(ntri))
    read(1) mesh_element
    read(1) bphicovar
    close(1)

    allocate(mesh_element_rmp(ntri))
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    allocate(currn(ntri, 3))

    Bnflux = 0d0
    Bnphi = 0d0
    currn = 0d0
  end subroutine read_mesh

  !> Read fluxes of perturbation field
  subroutine read_bnflux
    integer :: k
    real(dp) :: dummy8(8)

    open(1, file = Bnflux_file)
    do k = 1, ntri
       read(1, *) dummy8
       Bnflux(k,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(k,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(k,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(k) = cmplx(dummy8(7), dummy8(8), dp)

       if (abs((sum(Bnflux(k,:)) + imun * n * Bnphi(k) * mesh_element(k)%det_3 * 0.5d0) &
            / sum(Bnflux(k,:))) > 1d-10) then
          if (log_err) write(logfile, *) 'vacuum B_n not divergence free'
          stop 'vacuum B_n not divergence free'
       end if
    end do
    close(1)
  end subroutine read_bnflux

  !> TODO: get rid of this due to redundancy with bnflux
  subroutine read_hpsi
    integer :: k
    real(dp) :: dummy_re, dummy_im

    open(1, file = hpsi_file)
    do k = 1, ntri
       read (1, *) dummy_re, dummy_im
       if (.not. nonres) then
          mesh_element_rmp(k)%bnorm_vac = cmplx(dummy_re, dummy_im, dp)
       else
       !Test case: completely non-resonant perturbation
          mesh_element_rmp(k)%bnorm_vac = 3.d0 * R0 * abs(bphicovar) &
               / sum(mesh_point(mesh_element(k)%i_knot(:))%rcoord ** 2 &
               * mesh_point(mesh_element(k)%i_knot(:))%b_mod)
       endif
    end do
    close(1)
  end subroutine read_hpsi

  !>Initialize poloidal psi and unperturbed pressure p0
  subroutine init_flux_variables
    integer :: k
    real(dp) :: ddens_dpsi, dtemp_dpsi

    psimin = minval(mesh_point%psi_pol)
    psimax = maxval(mesh_point%psi_pol)

    if (log_info) write(logfile,*) 'psimin = ', psimin, '  psimax = ', psimax

    ddens_dpsi = di0 / psimax
    dtemp_dpsi = ti0 / psimax

    allocate(pres0(0:nflux+1))
    allocate(dpres0_dpsi(0:nflux+1))
    allocate(dens(0:nflux+1))
    allocate(temp(0:nflux+1))
    allocate(psi(0:nflux+1))

    psi(0) = mesh_point(1)%psi_pol  ! magnetic axis at k == 0 is not counted as flux surface
    do k = 1, nflux+1
       ! average over the loop to smooth out numerical errors
       psi(k) = sum(mesh_point((1 + (k-1) * nkpol + 1):(1 + k * nkpol))%psi_pol) / nkpol
    end do
    dens = (psi - psimin) / psimax * di0 + d_min
    temp = (psi - psimin) / psimax * ti0 + t_min
    pres0 = dens * temp * ideal_gas_factor
    dpres0_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ideal_gas_factor
  end subroutine init_flux_variables

  !>TODO
  subroutine assemble_system_first_order(nrow, a, b, d, du)
    integer, intent(in) :: nrow                    !< number of system rows
    complex(dp), intent(in), dimension(nrow) :: a  !< system coefficient \f$ a_{k} \f$
    complex(dp), intent(in), dimension(nrow) :: b  !< system coefficient \f$ b_{k} \f$
    complex(dp), intent(out) :: d(nrow)            !< diagonal of stiffness matrix \f$ A \f$
    complex(dp), intent(out) :: du(nrow)           !< superdiagonal of stiffness matrix \f$ A \f$ and \f$ A_{n, 1} \f$

    d = -a + b * 0.5d0
    du = a + b * 0.5d0
  end subroutine assemble_system_first_order

  !>TODO
  subroutine assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    integer, intent(in)  :: nrow                          !< number of system rows
    complex(dp), intent(in)  :: d(nrow)                   !< diagnonal of stiffness matrix \f$ A \f$
    complex(dp), intent(in)  :: du(nrow)                  !< superdiagonal of stiffness matrix \f$ A \f$ and \f$ A_{n, 1} \f$
    integer, intent(out) :: nz                            !< number of non-zero entries
    integer, intent(out) :: irow(2*nrow), icol(2*nrow)    !< matrix index representation
    complex(dp), intent(out) :: aval(2*nrow)              !< matrix values

    integer :: k

    irow(1) = 1
    icol(1) = 1
    aval(1) = d(1)

    irow(2) = nrow
    icol(2) = 1
    aval(2) = du(nrow)

    do k = 2,nrow
       ! off-diagonal
       irow(2*k-1) = k-1
       icol(2*k-1) = k
       aval(2*k-1) = du(k-1)

       ! diagonal
       irow(2*k) = k
       icol(2*k) = k
       aval(2*k) = d(k)
    end do

    nz = 2*nrow
  end subroutine assemble_sparse

  !>TODO
  subroutine common_triangles(knot1, knot2, common_tri)
    type(knot), intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)

    integer :: k, l, kcom
    kcom = 0
    do k = 1, knot1%n_owners
       do l = 1, knot2%n_owners
          if (knot1%i_owner_tri(k) == knot2%i_owner_tri(l)) then
             kcom = kcom+1
             if (kcom > 2) stop "Error: more than two common triangles for knots"
             common_tri(kcom) = knot1%i_owner_tri(k)
          end if
       end do
    end do
  end subroutine common_triangles

  !>Compute pressure perturbation. TODO: documentation
  subroutine compute_presn
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    real(dp) :: rold, zold, Brold, Bpold, Bzold ! previous values in loop
    real(dp) :: lr, lz  ! edge vector components
    complex(dp) :: Bnpsi
    complex(dp), dimension(nkpol) :: a, b, x, d, du
    integer :: kp, kl, kpold
    integer :: nz
    integer, dimension(2*nkpol) :: irow, icol
    complex(dp), dimension(2*nkpol) :: aval
    type(knot) :: oldknot, curknot
    integer :: common_tri(2)
    real(dp) :: perps(2)

    open(1, file = presn_file, recl = 1024)
    do kl = 1, nflux ! loop through flux surfaces

       curknot = mesh_point(1 + kl*nkpol)
       r = curknot%rcoord; p = 0d0; z = curknot%zcoord
       call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
            dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

       do kp = 1, nkpol ! loop through poloidal points along ring
          rold = r; zold = z
          Brold = Br; Bpold = Bp; Bzold = Bz
          kpold = kp - 1
          if (kpold == 0) kpold = nkpol
          oldknot = curknot
          curknot = mesh_point(1 + (kl-1)*nkpol + kp)
          r = curknot%rcoord; p = 0d0; z = curknot%zcoord

          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

          call common_triangles(oldknot, curknot, common_tri)

          lr = r - rold
          lz = z - zold
          perps = mesh_element(common_tri(:))%det_3 / hypot(lr, lz)
          if (nonres) then
             Bnflux(minval(common_tri), mod(mesh_element(minval(common_tri))%knot_h, 3) + 1) = &
                  R0 / (r + rold) * 2 * abs(bphicovar) * hypot(lr, lz) * sum(perps) / (psi(kl+1) - psi(kl-1))
          endif
          Bnpsi = Bnflux(minval(common_tri), mod(mesh_element(minval(common_tri))%knot_h, 3) + 1) / &
               (r + rold) * 2 / hypot(lr, lz) * (psi(kl+1) - psi(kl-1)) / sum(perps)

          x(kpold) = -dpres0_dpsi(kl) * Bnpsi

          a(kpold) = ((Br + Brold) * 0.5d0 * lr + (Bz + Bzold) * 0.5d0 * lz) / &
               (lr ** 2 + lz ** 2)

          b(kpold) = imun * n * (Bp / r + Bpold / rold) * 0.5d0
       end do ! kp

       ! solve linear system
       call assemble_system_first_order(nkpol, a, b, d, du)
       call assemble_sparse(nkpol, d, du, nz, irow, icol, aval)
       call sparse_solve(nkpol, nkpol, nz, irow, icol, aval, x)

       if (kl == 1) then ! first point on axis before actual output
          write(1, *) psi(0), dens(0), temp(0), pres0(0), &
               real(sum(x) / size(x)), aimag(sum(x) / size(x))
       end if
       do kp = 1, nkpol
          presn((kl - 1) * nkpol + 1 + kp) = x(kp)
          write(1, *) psi(kl), dens(kl), temp(kl), pres0(kl), &
               real(x(kp)), aimag(x(kp))
       end do
    end do ! kl
    do kp = 1, (npoint - nkpol * nflux - 1) ! write zeroes in remaining points until end
       write(1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine compute_presn

end module magdif
