module magdif_config
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90

  implicit none

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer, parameter :: logfile = 6             !< log to stdout, TODO: make this configurable

  integer :: log_level
  integer :: runmode 
  logical :: log_err, log_warn, log_info, log_debug ! specify log levels
  logical :: nonres = .false.  !< use non-resonant test case

  character(len=1024) :: point_file       !< input data file for mesh points
  character(len=1024) :: tri_file         !< input data file for triangles and edges
  character(len=1024) :: Bnflux_file      !< input data file for magnetic field perturbation
  character(len=1024) :: Bnflux_vac_file  !< input data file for magnetic field perturbation
  character(len=1024) :: Bn_sum_file      !< output data file for accumulated magnetic field perturbation
  character(len=1024) :: hpsi_file        !< input data file for \f$ h_{n}^{\psi} \f$
  character(len=1024) :: config_file      !< input config file for namelist settings
  character(len=1024) :: presn_file       !< output data file for pressure perturbation
  character(len=1024) :: currn_file       !< output data file for current perturbation

  integer  :: niter = 20
  integer  :: n               !< harmonic index of perturbation
  integer  :: nkpol           !< number of knots per poloidal loop
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density

  namelist / settings / log_level, runmode, nonres, point_file, tri_file, Bnflux_file, &
       Bnflux_vac_file, Bn_sum_file, hpsi_file, presn_file, currn_file, niter, n, nkpol, &
       nflux, ti0, di0
       !< namelist for input parameters

  contains

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
end module magdif_config

module magdif
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi, ev2erg                               ! PRELOAD/SRC/orbit_mod.f90
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       mesh_element_rmp, bphicovar, knot, triangle
  use for_macrostep, only : t_min, d_min
  use sparse_mod, only: remap_rc, sparse_solve, sparse_matmul
  use magdif_config

  implicit none

  real(dp), allocatable :: pres0(:)             !< unperturbed pressure \f$ p_{0} \f$ in dyn cm^-1
  real(dp), allocatable :: dpres0_dpsi(:)       !< derivative of unperturbed pressure w.r.t. flux surface label, \f$ p_{0}'(\psi) \f$
  real(dp), allocatable :: dens(:)              !< density \f$ \frac{N}{V} \f$ on flux surface in cm^-3
  real(dp), allocatable :: temp(:)              !< temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV
  real(dp), allocatable :: psi(:)               !< flux surface label \f$ \psi \f$
  real(dp), allocatable :: B2avg(:)             !< flux surface average \f$ \langle B_{0}^{2} \rangle \f$
  complex(dp), allocatable :: presn(:)          !< pressure perturbation \f$ p_{n} \f$ in each mesh point
  complex(dp), allocatable :: Bnflux(:,:)       !< edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$
  complex(dp), allocatable :: Bnphi(:)          !< physical toroidal component of magnetic perturbation, \f$ B_{n (\phi)} \f$
  complex(dp), allocatable :: Bnflux_vac(:,:)   !< vacuum edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$
  complex(dp), allocatable :: Bnphi_vac(:)      !< vacuum physical toroidal component of magnetic perturbation, \f$ B_{n (\phi)} \f$
  complex(dp), allocatable :: jnflux(:,:)       !< edge currents \f$ R \vec{j}_{n} \cdot \vec{n} \f$
  complex(dp), allocatable :: jnphi(:)          !< physical toroidal component of current perturbation, \f$ j_{n (\phi)} \f$
  real(dp), allocatable :: j0phi(:)             !< physical toroidal component of equilibrium current, \f$ j_{0 (\phi)} \f$

  real(dp) :: psimin  !< minimum flux surface label, located at the separatrix
  real(dp) :: psimax  !< maximum flux surface label, located at the magnetic axis
  real(dp), parameter :: R0 = 172.74467899999999d0  !< distance of magnetic axis from center, \f$ R_{0} \f$
  real(dp), parameter :: clight = 2.99792458d10  !< speed of light in cm/s
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< imaginary unit in double precision

  real(dp), allocatable :: q(:), dqdpsi(:) !< safety factor and derivative over psi

  interface
     subroutine sub_assemble_flux_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, l, &
          edge_index, edge_name)
       import :: dp, triangle
       complex(dp), dimension(:), intent(inout) :: x, d, du
       integer, intent(in) :: kl, kp, k_low
       real(dp), intent(in) :: Deltapsi
       type(triangle), intent(in) :: elem
       integer, dimension(2), intent(in) :: l
       integer, intent(in) :: edge_index
       character, intent(in) :: edge_name
     end subroutine sub_assemble_flux_coeff
  end interface

  interface
     subroutine sub_assign_flux(x, kl, kp, kp_max, k_low, ei, eo, ef)
       import :: dp
       complex(dp), dimension(:), intent(in) :: x
       integer, intent(in) :: kl, kp, kp_max, k_low
       integer, intent(in) :: ei, eo, ef
     end subroutine sub_assign_flux
  end interface

contains

  !> Initialize magdif module
  subroutine magdif_init
    call read_mesh
    call read_hpsi ! TODO: get rid of this due to redundancy with bnflux
    call init_flux_variables
    call init_safety_factor
    call compute_j0phi
    if (nonres) then
       call compute_Bn_nonres
    else
       call read_bnflux(Bnflux_vac_file)
    end if
    Bnflux_vac = Bnflux
    Bnphi_vac = Bnphi
    if (log_info) write(logfile, *) 'magdif initialized'
  end subroutine magdif_init

  !> Final cleanup of magdif module
  subroutine magdif_cleanup
    if (allocated(q)) deallocate(q)
    if (allocated(dqdpsi)) deallocate(dqdpsi)
    if (allocated(pres0)) deallocate(pres0)
    if (allocated(dpres0_dpsi)) deallocate(dpres0_dpsi)
    if (allocated(dens)) deallocate(dens)
    if (allocated(temp)) deallocate(temp)
    if (allocated(psi)) deallocate(psi)
    if (allocated(B2avg)) deallocate(B2avg)
    if (allocated(presn)) deallocate(presn)
    if (allocated(jnflux)) deallocate(jnflux)
    if (allocated(Bnflux)) deallocate(Bnflux)
    if (allocated(Bnphi)) deallocate(Bnphi)
    if (allocated(Bnflux_vac)) deallocate(Bnflux_vac)
    if (allocated(Bnphi_vac)) deallocate(Bnphi_vac)
    if (allocated(jnphi)) deallocate(jnphi)
    if (allocated(j0phi)) deallocate(j0phi)
    if (allocated(mesh_point)) deallocate(mesh_point)
    if (allocated(mesh_element)) deallocate(mesh_element)
    if (allocated(mesh_element_rmp)) deallocate(mesh_element_rmp)
    if (log_info) write(logfile, *) 'magdif cleanup finished'
  end subroutine magdif_cleanup

  subroutine magdif_single
    integer :: stat

    call compute_presn
    call compute_currn
    call compute_Bn(stat)
    call read_bnflux(Bnflux_file)
  end subroutine magdif_single

  subroutine magdif_direct
    integer :: kiter, kt
    integer :: stat
    complex(dp) :: Bnflux_sum(ntri, 3)
    complex(dp) :: Bnphi_sum(ntri)

    Bnflux_sum = 0.0d0
    Bnphi_sum = 0.0d0
    do kiter = 1, niter
       if (log_info) write(logfile, *) 'Iteration ', kiter, ' of ', niter
       call compute_presn             ! compute pressure based on previous perturbation field
       call compute_currn
       call compute_Bn(stat)          ! use field code to generate new field from currents
       call read_bnflux(Bnflux_file)  ! read new bnflux from field code
       Bnflux_sum = Bnflux_sum + Bnflux
       Bnphi_sum = Bnphi_sum + Bnphi
    end do

    open(1, file = Bn_sum_file)
    do kt = 1, ntri
       write(1, *) &
            real(Bnflux_sum(kt, 1)), aimag(Bnflux_sum(kt, 1)), &
            real(Bnflux_sum(kt, 2)), aimag(Bnflux_sum(kt, 2)), &
            real(Bnflux_sum(kt, 3)), aimag(Bnflux_sum(kt, 3)), &
            real(Bnphi_sum(kt) * mesh_element(kt)%det_3 * 0.5d0), &
            aimag(Bnphi_sum(kt) * mesh_element(kt)%det_3 * 0.5d0)
    end do
    close(1)
  end subroutine magdif_direct

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
    allocate(Bnflux_vac(ntri, 3))
    allocate(Bnphi_vac(ntri))
    allocate(jnphi(ntri))
    allocate(j0phi(ntri))
    allocate(jnflux(ntri, 3))

    Bnflux = 0d0
    Bnphi = 0d0
    Bnflux_vac = 0d0
    Bnphi_vac = 0d0
    jnphi = 0d0
    j0phi = 0d0
    jnflux = 0d0
  end subroutine read_mesh

  subroutine compute_Bn(stat)
    integer, intent(out) :: stat

    call execute_command_line (&
         "PATH=/temp/ert/local/bin:$PATH /temp/ert/local/bin/FreeFem++ " // &
         "../FEM/maxwell.edp ../PRELOAD/inputformaxwell_ext.msh " // currn_file // &
         "2 > /tmp/freefem.out 2>&1 && cd ..", &
         exitstat = stat)

  end subroutine compute_Bn

  !> Read fluxes of perturbation field
  subroutine read_bnflux(filename)
    character(len = 1024) :: filename
    integer :: k
    real(dp) :: dummy8(8)

    open(1, file = filename)
    do k = 1, ntri
       read(1, *) dummy8
       Bnflux(k,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(k,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(k,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(k) = cmplx(dummy8(7), dummy8(8), dp) / (mesh_element(k)%det_3 * 0.5d0)

       if (abs((sum(Bnflux(k,:)) + imun * n * Bnphi(k) * mesh_element(k)%det_3 * 0.5d0) &
            / sum(Bnflux(k,:))) > 1d-3) then
          if (log_err) write(logfile, *) 'B_n not divergence free: ',&
               abs((sum(Bnflux(k,:)) + imun * n * Bnphi(k) * mesh_element(k)%det_3 * 0.5d0) &
            / sum(Bnflux(k,:)))
          stop 'B_n not divergence free'
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

    subroutine init_safety_factor
    integer :: kl, kp, kp_max, k_low
    type(triangle) :: elem
    real(dp) :: r, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    allocate(q(nflux+1))
    allocate(dqdpsi(0:nflux))
    q = 0d0
    dqdpsi = 0d0

    do kl = 1, nflux+1
       select case(kl)
       case (1)
          kp_max = nkpol
          k_low = 0
       case default
          kp_max = 2 * nkpol
          k_low = (kl-1) * kp_max - nkpol
       end select
       do kp = 1, kp_max
          elem = mesh_element(k_low + kp)
          call ring_centered_avg_coord(elem, r, z)
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          q(kl) = q(kl) + Bp * elem%det_3 * 0.5d0
       end do
       q(kl) = -q(kl) * 0.5d0 / pi / (psi(kl) - psi(kl-1))  ! check sign
    end do

    do kl = 1, nflux
       ! q on triangle strip, psi on edge ring
       dqdpsi(kl) = (q(kl+1) - q(kl)) / (psi(kl+1) - psi(kl-1)) * 0.5d0
    end do

    if (log_debug) then
       open(1, file = 'qsafety.out')
       write(1,*) psi(0), 0.0d0, dqdpsi(0)
       do kl = 1, nflux
          write(1,*) psi(kl), q(kl), dqdpsi(kl)
       end do
       close(1)
    end if
  end subroutine init_safety_factor

  subroutine ring_centered_avg_coord(elem, r, z)
    type(triangle), intent(in) :: elem
    real(dp), intent(out) :: r, z
    type(knot), dimension(3) :: knots

    knots = mesh_point(elem%i_knot)
    r = (sum(knots%rcoord) + knots(elem%knot_h)%rcoord) * 0.25d0
    z = (sum(knots%zcoord) + knots(elem%knot_h)%zcoord) * 0.25d0
  end subroutine ring_centered_avg_coord

  !>Initialize poloidal psi and unperturbed pressure p0
  subroutine init_flux_variables
    integer :: kl
    real(dp) :: ddens_dpsi, dtemp_dpsi
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kp, kp_max, k_low

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
    allocate(B2avg(nflux))

    psi(0) = mesh_point(1)%psi_pol  ! magnetic axis at k == 0 is not counted as flux surface
    do kl = 1, nflux+1
       ! average over the loop to smooth out numerical errors
       psi(kl) = sum(mesh_point((1 + (kl-1) * nkpol + 1):(1 + kl * nkpol))%psi_pol) / nkpol

       select case(kl)
       case (1)
          kp_max = nkpol
          k_low = 0
       case default
          kp_max = 2 * nkpol
          k_low = (kl-1) * kp_max - nkpol
       end select
       if (kl > nflux) cycle
       B2avg(kl) = 0d0
       do kp = 1, kp_max
          call ring_centered_avg_coord(mesh_element(k_low + kp), r, z)
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          B2avg(kl) = B2avg(kl) + Br**2 + Bp**2 + Bz**2
       end do !kp
       B2avg(kl) = B2avg(kl) / kp_max
    end do
    dens = (psi - psimin) / psimax * di0 + d_min
    temp = (psi - psimin) / psimax * ti0 + t_min
    pres0 = dens * temp * ev2erg
    dpres0_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
  end subroutine init_flux_variables

  subroutine compute_j0phi
    integer :: kl, kp, kp_max, k_low
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    real(dp) :: Bpol

    do kl = 1, nflux
       select case(kl)
       case (1)
          kp_max = nkpol
          k_low = 0
       case default
          kp_max = 2 * nkpol
          k_low = (kl-1) * kp_max - nkpol
       end select
       do kp = 1, kp_max
          call ring_centered_avg_coord(mesh_element(k_low + kp), r, z)
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          Bpol = hypot(Br, Bz)
          j0phi(k_low + kp) = clight * (pres0(kl) - pres0(kl-1)) / (psi(kl) - psi(kl-1)) * &
               (Bp ** 2 / B2avg(kl) + (Bpol ** 2 - Bp ** 2) / (Bpol ** 2 + Bp ** 2))
       end do
    end do
  end subroutine compute_j0phi

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

  subroutine unshared_knots(common_tri, outer_knot, inner_knot)
    integer, intent(in) :: common_tri(2)
    type(knot), intent(out) :: outer_knot, inner_knot
    integer, dimension(2, 3) :: set
    integer, dimension(2) :: i_unshared
    integer :: tri, k, n_unshared
    character(len = *), parameter :: errmsg = &
         'not exactly one unshared knot between triangles'

    set(1, :) = mesh_element(common_tri(1))%i_knot(:)
    set(2, :) = mesh_element(common_tri(2))%i_knot(:)
    do tri = 1, 2
       n_unshared = 0
       do k = 1, 3
          if (.not. any(set(tri, k) == set(3-tri, :))) then
             n_unshared = n_unshared + 1
             i_unshared(tri) = set(tri, k)
          end if
       end do
       if (n_unshared /= 1) then
          if (log_debug) write(logfile, *) errmsg
          stop errmsg
       end if
    end do
    if (i_unshared(1) > i_unshared(2)) then
       outer_knot = mesh_point(i_unshared(1))
       inner_knot = mesh_point(i_unshared(2))
    else
       outer_knot = mesh_point(i_unshared(2))
       inner_knot = mesh_point(i_unshared(1))
    end if
  end subroutine unshared_knots

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

       do kp = 1, nkpol ! loop through poloidal points along loop
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
          Bnpsi = Bnflux(minval(common_tri), mod(mesh_element(minval(common_tri))%knot_h, 3) + 1) / &
               (r + rold) * 2 / hypot(lr, lz) * (psi(kl+1) - psi(kl-1)) / sum(perps)

          x(kp) = -dpres0_dpsi(kl) * Bnpsi

          a(kp) = ((Br + Brold) * 0.5d0 * lr + (Bz + Bzold) * 0.5d0 * lz) / &
               (lr ** 2 + lz ** 2)

          b(kp) = imun * n * (Bp / r + Bpold / rold) * 0.5d0
       end do ! kp

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
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

  subroutine get_labeled_edges(elem, li, lo, lf, ei, eo, ef)
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    integer, dimension(3) :: i_knot_diff
    integer :: knot_i, knot_o, knot_f
    integer :: i1, i2
    logical :: closing_loop
    character(len = *), parameter :: errmsg = 'cannot find correct label for triangle edges'

    i1 = 0
    i2 = 0
    select case (elem%knot_h)
    case (1)
       i1 = 2
       i2 = 3
    case (2)
       i1 = 3
       i2 = 1
    case (3)
       i1 = 1
       i2 = 2
    end select
    if (elem%i_knot(i1) == elem%i_knot(i2)) then
       if (log_debug) write(logfile, *) errmsg
       stop errmsg
    end if
    closing_loop = abs(elem%i_knot(i1) - elem%i_knot(i2)) /= 1  ! indices not next to each other
    if ((elem%i_knot(i1) > elem%i_knot(i2)) .neqv. closing_loop) then  ! i1 is next counter-clockwise
       knot_i = i1
       knot_o = i2
    else  ! i2 is next counter-clockwise
       knot_i = i2
       knot_o = i1
    end if
    knot_f = elem%knot_h
    i_knot_diff = elem%i_knot - elem%i_knot(knot_f)
    if (all(i_knot_diff >= 0)) then  ! knot_f lies on inner surface
       ei = knot_f
       eo = knot_i
       ef = knot_o
       li = (/ elem%i_knot(knot_f), elem%i_knot(knot_o) /)
       lo = (/ elem%i_knot(knot_i), elem%i_knot(knot_f) /)
       lf = (/ elem%i_knot(knot_o), elem%i_knot(knot_i) /)
    else if (all(i_knot_diff <= 0)) then ! elem%i_knot(knot_f) lies on outer surface
       ei = knot_o
       eo = knot_f
       ef = knot_i
       li = (/ elem%i_knot(knot_o), elem%i_knot(knot_f) /)
       lo = (/ elem%i_knot(knot_f), elem%i_knot(knot_i) /)
       lf = (/ elem%i_knot(knot_i), elem%i_knot(knot_o) /)
    else
       if (log_debug) write(logfile, *) errmsg
       stop errmsg
    end if
  end subroutine get_labeled_edges

  subroutine compute_triangle_flux(assemble_flux_coeff, assign_flux, pol_flux, tor_flux, &
       outfile)
    procedure(sub_assemble_flux_coeff) :: assemble_flux_coeff
    procedure(sub_assign_flux) :: assign_flux
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_flux(:)
    character(len = 1024), intent(in), optional :: outfile

    complex(dp), dimension(2 * nkpol) :: x, d, du
    integer :: kl, kp, kp_max, k_low
    integer :: nz
    integer, dimension(4 * nkpol) :: irow, icol
    complex(dp), dimension(4 * nkpol) :: aval
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    real(dp) :: Deltapsi
    integer :: ktri

    ktri = 0

    if (present(outfile)) open(1, file = outfile, recl = 1024)
    do kl = 1, nflux ! loop through flux surfaces
       select case(kl)
       case (1)
          kp_max = nkpol
          k_low = 0
       case (2)
          kp_max = 2 * nkpol
          k_low = (kl-1) * kp_max - nkpol
       case default
          kp_max = 2 * nkpol
          k_low = (kl-1) * kp_max - nkpol
       end select
       Deltapsi = psi(kl) - psi(kl-1)
       do kp = 1, kp_max
          elem = mesh_element(k_low + kp)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef)
          call assemble_flux_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, lf, ef, 'f')
          call assemble_flux_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, li, ei, 'i')
          call assemble_flux_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, lo, eo, 'o')
       end do
       call assemble_sparse(kp_max, d(:kp_max), du(:kp_max), nz, irow(:2*kp_max),&
            icol(:2*kp_max), aval(:2*kp_max))
       call sparse_solve(kp_max, kp_max, nz, irow(:nz), icol(:nz), aval(:nz), x(:kp_max))
       do kp = 1, kp_max
          elem = mesh_element(k_low + kp)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef)
          call assign_flux(x, kl, kp, kp_max, k_low, ei, eo, ef)
          if (present(outfile)) write(1, *) &
               real(pol_flux(k_low + kp, 1)), aimag(pol_flux(k_low + kp, 1)), &
               real(pol_flux(k_low + kp, 2)), aimag(pol_flux(k_low + kp, 2)), &
               real(pol_flux(k_low + kp, 3)), aimag(pol_flux(k_low + kp, 3)), &
               real(tor_flux(k_low + kp)), aimag(tor_flux(k_low + kp))
          ktri = ktri+1
       end do
    end do

    if (present(outfile)) then
       do while (ktri < ntri)
          write(1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
          ktri = ktri+1
       end do
       close(1)
    end if
  end subroutine compute_triangle_flux

  subroutine assemble_currn_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, l, &
       edge_index, edge_name)

    complex(dp), dimension(:), intent(inout) :: x, d, du
    integer, intent(in) :: kl, kp, k_low
    real(dp), intent(in) :: Deltapsi
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(in) :: l
    integer, intent(in) :: edge_index
    character, intent(in) :: edge_name

    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    type(knot) :: base, tip

    p = 0d0
    ! use midpoint of edge
    base = mesh_point(l(1))
    tip = mesh_point(l(2))
    r = (base%rcoord + tip%rcoord) * 0.5d0
    z = (base%zcoord + tip%zcoord) * 0.5d0
    call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
    select case (edge_name)
    case ('f')
       ! first term on source side: flux through edge f
       x(kp) = -clight * r / Bp * (presn(l(2)) - presn(l(1))) - &
            j0phi(k_low + kp) / Bp / r * Bnflux(k_low + kp, edge_index)
       jnflux(k_low + kp, edge_index) = x(kp)
    case ('i')
       ! diagonal matrix element
       d(kp) = -1d0 - imun * n * elem%det_3 * 0.25d0 * Bp / Deltapsi
       ! additional term from edge i on source side
       x(kp) = x(kp) - imun * n * elem%det_3 * 0.25d0 * (clight * r / (-Deltapsi) * &
            (presn(l(2)) - presn(l(1)) - Bnphi(k_low + kp) / r / Bp * &
            (pres0(kl) - pres0(kl-1))) + j0phi(k_low + kp) * (Bnphi(k_low + kp) / Bp + &
            Bnflux(k_low + kp, edge_index) / r / (-Deltapsi)))
    case ('o')
       ! superdiagonal matrix element
       du(kp) = 1d0 - imun * n * elem%det_3 * 0.25d0 * Bp / Deltapsi
       ! additional term from edge o on source side
       x(kp) = x(kp) - imun * n * elem%det_3 * 0.25d0 * (clight * r / Deltapsi * &
            (presn(l(2)) - presn(l(1)) - Bnphi(k_low + kp) / r / Bp * &
            (pres0(kl-1) - pres0(kl))) + j0phi(k_low + kp) * (Bnphi(k_low + kp) / Bp + &
            Bnflux(k_low + kp, edge_index) / r / Deltapsi))
    end select
  end subroutine assemble_currn_coeff

  subroutine assign_currn(x, kl, kp, kp_max, k_low, ei, eo, ef)
    complex(dp), dimension(:), intent(in) :: x
    integer, intent(in) :: kl, kp, kp_max, k_low
    integer, intent(in) :: ei, eo, ef
    jnflux(k_low + kp, ei) = -x(kp)
    jnflux(k_low + kp, eo) = x(mod(kp, kp_max) + 1)
    jnphi(k_low + kp) = sum(jnflux(k_low + kp, :)) * imun / n
  end subroutine assign_currn

  subroutine compute_currn
    call compute_triangle_flux(assemble_currn_coeff, assign_currn, jnflux, jnphi, currn_file)
  end subroutine compute_currn

  subroutine assemble_Bn_nonres_coeff(x, d, du, kl, kp, k_low, Deltapsi, elem, l, &
       edge_index, edge_name)

    complex(dp), dimension(:), intent(inout) :: x, d, du
    integer, intent(in) :: kl, kp, k_low
    real(dp), intent(in) :: Deltapsi
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(in) :: l
    integer, intent(in) :: edge_index
    character, intent(in) :: edge_name

    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    type(knot) :: base, tip, outer, inner
    integer :: common_tri(2)
    real(dp) :: dpsi_dr, dpsi_dz, nr, nz

    select case (edge_name)
    case ('f')
       p = 0d0
       ! use midpoint of edge
       base = mesh_point(l(1))
       tip = mesh_point(l(2))
       call common_triangles(base, tip, common_tri)
       call unshared_knots(common_tri, outer, inner)
       r = (base%rcoord + tip%rcoord) * 0.5d0
       z = (base%zcoord + tip%zcoord) * 0.5d0
       call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
            dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
       dpsi_dr = (psi(kl+1) - psi(kl-1)) / (outer%rcoord - inner%rcoord)
       dpsi_dz = (psi(kl+1) - psi(kl-1)) / (outer%zcoord - inner%zcoord)
       nr = (tip%zcoord - base%zcoord)
       nz = -(tip%rcoord - base%rcoord)  ! normal pointing out of the triangle
       Bnflux(k_low + kp, edge_index) = Bp * R0 / r * (nr * dpsi_dr + nz * dpsi_dz) / &
            (dpsi_dr ** 2 + dpsi_dz ** 2)
       ! use weighted triangle midpoint
       call ring_centered_avg_coord(elem, r, z)
       call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
            dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
       Bnphi(k_low + kp) = Bp * R0 / r ** 2 * imun / n / q(kl) * &
            (dqdpsi(kl) + dqdpsi(kl-1)) * 0.5d0
       x(kp) = -Bnflux(k_low + kp, edge_index) - imun * n * elem%det_3 * 0.5d0 * r * &
            Bnphi(k_low + kp)
    case ('i')
       d(kp) = -1d0
    case ('o')
       du(kp) = 1d0
    end select
  end subroutine assemble_Bn_nonres_coeff

  subroutine assign_Bn_nonres(x, kl, kp, kp_max, k_low, ei, eo, ef)
    complex(dp), dimension(:), intent(in) :: x
    integer, intent(in) :: kl, kp, kp_max, k_low
    integer, intent(in) :: ei, eo, ef
    Bnflux(k_low + kp, ei) = -x(kp)
    Bnflux(k_low + kp, eo) = x(mod(kp, kp_max) + 1)
  end subroutine assign_Bn_nonres

  subroutine compute_Bn_nonres
    integer :: k
    call compute_triangle_flux(assemble_Bn_nonres_coeff, assign_Bn_nonres, Bnflux, Bnphi)
    do k = 1, ntri
       if (abs((sum(Bnflux(k,:)) + imun * n * Bnphi(k) * mesh_element(k)%det_3 * 0.5d0) &
            / sum(Bnflux(k,:))) > 1d-10) then
          if (log_err) write(logfile, *) 'nonresonant B_n not divergence free'
          stop 'nonresonant B_n not divergence free'
       end if
    end do

  end subroutine compute_Bn_nonres

end module magdif
