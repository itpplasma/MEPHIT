module magdif
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi, ev2erg                               ! PRELOAD/SRC/orbit_mod.f90
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       mesh_element_rmp, bphicovar, knot, triangle
  use for_macrostep, only : t_min, d_min
  use sparse_mod, only: remap_rc, sparse_solve, sparse_matmul
  use magdif_config

  implicit none

  !> Unperturbed pressure \f$ p_{0} \f$ in dyn cm^-1.
  !>
  !> Values are taken on flux surfaces with indices running from 0 to #magdif_config::nflux
  !> +1, i.e. from the magnetic axis to just outside the last closed flux surface.
  real(dp), allocatable :: pres0(:)

  !> Derivative of unperturbed pressure w.r.t. flux surface label, \f$ p_{0}'(\psi) \f$.
  !>
  !> Values are taken on flux surfaces with indices running from 0 to #magdif_config::nflux
  !> +1, i.e. from the magnetic axis to just outside the last closed flux surface.
  real(dp), allocatable :: dpres0_dpsi(:)

  !> Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
  !>
  !> Values are taken on flux surfaces with indices running from 0 to #magdif_config::nflux
  !> +1, i.e. from the magnetic axis to just outside the last closed flux surface.
  real(dp), allocatable :: dens(:)

  !> Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
  !>
  !> Values are taken on flux surfaces with indices running from 0 to #magdif_config::nflux
  !> +1, i.e. from the magnetic axis to just outside the last closed flux surface.
  real(dp), allocatable :: temp(:)

  !> Magnetic flux surface label \f$ \psi \f$ in G cm^2.
  !>
  !> \f$ \psi \f$ is the ribbon poloidal flux divided by \f$ 2 \pi \f$. Its sign is
  !> positive and its magnitude is growing in the radially inward direction. The indices
  !> are running from 0 for the magnetic axis to #magdif_config::nflux +1 for the first
  !> non-closed flux surface.
  real(dp), allocatable :: psi(:)

  real(dp) :: psimin  !< Minimum flux surface label, located at the LCFS
  real(dp) :: psimax  !< Maximum flux surface label, located at the magnetic axis

  !> Safety factor \f$ q \f$ (dimensionless).
  !>
  !> Values are taken between two flux surfaces with indices running from 1 to
  !> #magdif_config::nflux, i.e. from the triangle strip surrounding the magnetic axis
  !> to the triangle strip just inside the last closed flux surface.
  real(dp), allocatable :: q(:)

  !> Flux surface average \f$ \langle B_{0}^{2} \rangle \f$.
  !>
  !> Values are taken between two flux surfaces with indices running from 1 to
  !> #magdif_config::nflux, i.e. from the triangle strip surrounding the magnetic axis to
  !> the triangle strip contained by the last closed flux surface.
  real(dp), allocatable :: B2avg(:)

  !> Pressure perturbation \f$ p_{n} \f$ in dyn cm^-1.
  !>
  !> Values are taken at each mesh point and the indexing scheme is the same as for
  !> #mesh_mod::mesh_point.
  complex(dp), allocatable :: presn(:)

  !> Edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$ in G cm^2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: Bnflux(:,:)

  !> Physical toroidal component of magnetic perturbation \f$ B_{n (\phi)} \f$ in G.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: Bnphi(:)

  !> Vacuum perturbation edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$ in G cm^2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: Bnflux_vac(:,:)

  !> Physical toroidal component of vacuum magnetic perturbation \f$ B_{n(\phi)} \f$ in G.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: Bnphi_vac(:)

  !> Edge currents \f$ R \vec{j}_{n} \cdot \vec{n} \f$ in statampere.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: jnflux(:,:)

  !> Physical toroidal component of current perturbation \f$ j_{n (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: jnphi(:)

  !> Physical toroidal component of equilibrium current \f$ j_{0 (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  real(dp), allocatable :: j0phi(:)

  integer, allocatable :: kp_max(:)
  integer, allocatable :: kt_max(:)
  integer, allocatable :: kp_low(:)
  integer, allocatable :: kt_low(:)

  !> Distance of magnetic axis from center \f$ R_{0} \f$ in cm.
  real(dp), parameter :: R0 = 172.74467899999999d0

  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

  interface
     subroutine sub_assemble_flux_coeff(x, d, du, kf, kt, elem, lk, ek, edge_name)
       import :: dp, triangle
       complex(dp), dimension(:), intent(inout) :: x, d, du
       integer, intent(in) :: kf, kt
       type(triangle), intent(in) :: elem
       integer, dimension(2), intent(in) :: lk
       integer, intent(in) :: ek
       character, intent(in) :: edge_name
     end subroutine sub_assemble_flux_coeff
  end interface

  interface
     subroutine sub_assign_flux(x, kf, kt, ei, eo, ef)
       import :: dp
       complex(dp), dimension(:), intent(in) :: x
       integer, intent(in) :: kf, kt
       integer, intent(in) :: ei, eo, ef
     end subroutine sub_assign_flux
  end interface

contains

  !> Initialize magdif module
  subroutine magdif_init
    call init_indices
    call read_mesh
    call init_flux_variables
    call init_safety_factor
    call compute_j0phi
    if (nonres) then
       call compute_Bn_nonres
       if (log_debug) call dump_triangle_flux(Bnflux, Bnphi, 'plot_Bn_nonres.dat')
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
    if (allocated(kp_max)) deallocate(kp_max)
    if (allocated(kt_max)) deallocate(kt_max)
    if (allocated(kp_low)) deallocate(kp_low)
    if (allocated(kt_low)) deallocate(kt_low)
    if (log_info) write(logfile, *) 'magdif cleanup finished'
  end subroutine magdif_cleanup

  subroutine magdif_single
    call compute_presn
    call compute_currn
    if (log_debug) call dump_triangle_flux(jnflux, jnphi, 'plot_jn_nonres.dat')
    call compute_Bn
    call read_bnflux(Bnflux_file)
  end subroutine magdif_single

  subroutine magdif_direct
    integer :: kiter, kt
    complex(dp) :: Bnflux_sum(ntri, 3)
    complex(dp) :: Bnphi_sum(ntri)

    Bnflux_sum = 0.0d0
    Bnphi_sum = 0.0d0
    do kiter = 1, niter
       if (log_info) write(logfile, *) 'Iteration ', kiter, ' of ', niter
       call compute_presn             ! compute pressure based on previous perturbation field
       call compute_currn
       call compute_Bn                ! use field code to generate new field from currents
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

  subroutine init_indices
    integer :: kf

    allocate (kp_max(nflux+2))
    allocate (kt_max(nflux+1))
    allocate (kp_low(nflux+2))
    allocate (kt_low(nflux+1))

    kp_max = nkpol
    kt_max = 2 * nkpol
    kt_max(1) = nkpol

    kp_low(1) = 1
    do kf = 2, nflux+2
       kp_low(kf) = kp_low(kf-1) + kp_max(kf-1)
    end do
    kt_low(1) = 0
    do kf = 2, nflux+1
       kt_low(kf) = kt_low(kf-1) + kt_max(kf-1)
    end do
  end subroutine init_indices

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

  subroutine compute_Bn
    integer :: stat, dummy
    call execute_command_line("./maxwell.sh " // currn_file, exitstat = stat, &
         cmdstat = dummy)
    if (stat /= 0) then
       if (log_err) write(logfile, *) 'FreeFem++ failed with exit code ', stat
       stop 'FreeFem++ failed'
    end if
  end subroutine compute_Bn


  !> Check if divergence-freeness of the given vector field is fulfilled on each
  !> triangle, otherwise halt the program.
  !>
  !> @param pol_flux poloidal flux components, e.g. #jnflux - first index refers to the
  !> triangle, second index refers to the edge on that triangle, as per
  !> #mesh_mod::mesh_element
  !> @param tor_comp toroidal field components, e.g. #jnphi - index refers to the triangle
  !> @param abs_err absolute error threshold, i.e. maximum acceptable value for divergence
  !> @param field_name name given to the vector field in the error message
  !>
  !> This subroutine calculates the divergence via the divergence theorem, i.e. by adding
  !> up the fluxes of the vector field through triangle edges. If this sum is higher than
  !> \p abs_err on any triangle, it halts the program.

  subroutine check_div_free(pol_flux, tor_comp, abs_err, field_name)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    real(dp), intent(in) :: abs_err
    character(len = *), intent(in) :: field_name

    integer :: kt
    real(dp) :: div
    character(len = len_trim(field_name) + 20) :: err_msg
    err_msg = trim(field_name) // ' not divergence-free'

    do kt = 1, ntri
       div = abs((sum(pol_flux(kt,:)) + imun * n * tor_comp(kt) * &
            mesh_element(kt)%det_3 * 0.5d0) / sum(pol_flux(kt,:)))
       if (div > abs_err) then
          if (log_err) write(logfile, *) err_msg, ': ', div
          stop err_msg
       end if
    end do
  end subroutine check_div_free

  !> Read fluxes of perturbation field
  subroutine read_bnflux(filename)
    character(len = 1024) :: filename
    integer :: kt
    real(dp) :: dummy8(8)

    open(1, file = filename)
    do kt = 1, ntri
       read(1, *) dummy8
       Bnflux(kt,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(kt,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(kt,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(kt) = cmplx(dummy8(7), dummy8(8), dp) / mesh_element(kt)%det_3 * 2d0
    end do
    close(1)

    call check_div_free(Bnflux, Bnphi, 1d-3, 'B_n')
  end subroutine read_bnflux

  subroutine init_safety_factor
    integer :: kf, kt
    type(triangle) :: elem
    real(dp) :: r, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    allocate(q(nflux))
    q = 0d0

    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call ring_centered_avg_coord(elem, r, z)
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          q(kf) = q(kf) + Bp * elem%det_3 * 0.5d0
       end do
       q(kf) = -q(kf) * 0.5d0 / pi / (psi(kf) - psi(kf-1))  ! check sign
    end do

    if (log_debug) then
       open(1, file = 'qsafety.out')
       write(1,*) psi(0), 0.0d0
       do kf = 1, nflux
          write(1,*) psi(kf), q(kf)
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
    integer :: kf, kt
    real(dp) :: ddens_dpsi, dtemp_dpsi
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

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

    ! magnetic axis at k == 0 is not counted as flux surface
    psi(0) = mesh_point(1)%psi_pol

    do kf = 1, nflux+1
       ! average over the loop to smooth out numerical errors
       psi(kf) = sum(mesh_point((kp_low(kf) + 1):kp_low(kf+1))%psi_pol) / kp_max(kf)

       if (kf > nflux) cycle
       B2avg(kf) = 0d0
       do kt = 1, kt_max(kf)
          call ring_centered_avg_coord(mesh_element(kt_low(kf) + kt), r, z)
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          B2avg(kf) = B2avg(kf) + Br**2 + Bp**2 + Bz**2
       end do
       B2avg(kf) = B2avg(kf) / kt_max(kf)
    end do
    dens = (psi - psimin) / psimax * di0 + d_min
    temp = (psi - psimin) / psimax * ti0 + t_min
    pres0 = dens * temp * ev2erg
    dpres0_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
  end subroutine init_flux_variables

  subroutine compute_j0phi
    integer :: kf, kt
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    real(dp) :: Bpol

    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          call ring_centered_avg_coord(mesh_element(kt_low(kf) + kt), r, z)
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          Bpol = hypot(Br, Bz)
          j0phi(kt_low(kf) + kt) = &
               clight * (pres0(kf) - pres0(kf-1)) / (psi(kf) - psi(kf-1)) * &
               (Bp ** 2 / B2avg(kf) + (Bpol ** 2 - Bp ** 2) / (Bpol ** 2 + Bp ** 2)) * r
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

  !>Compute pressure perturbation. TODO: documentation
  subroutine compute_presn
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    real(dp) :: r_prev, z_prev, Br_prev, Bp_prev, Bz_prev ! previous values in loop
    real(dp) :: lr, lz  ! edge vector components
    complex(dp) :: Bnpsi
    complex(dp), dimension(nkpol) :: a, b, x, d, du
    integer :: kf, kp, kp_prev
    integer :: nz
    integer, dimension(2*nkpol) :: irow, icol
    complex(dp), dimension(2*nkpol) :: aval
    type(knot) :: cur_knot, prev_knot
    integer :: common_tri(2)
    real(dp) :: perps(2)

    open(1, file = presn_file, recl = 1024)
    do kf = 1, nflux ! loop through flux surfaces

       cur_knot = mesh_point(kp_low(kf+1))
       r = cur_knot%rcoord; p = 0d0; z = cur_knot%zcoord
       call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
            dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

       do kp = 1, kp_max(kf) ! loop through poloidal points along loop
          r_prev = r; z_prev = z
          Br_prev = Br; Bp_prev = Bp; Bz_prev = Bz
          kp_prev = kp - 1
          if (kp_prev == 0) kp_prev = kp_max(kf)  ! index "wraps around"
          prev_knot = cur_knot
          cur_knot = mesh_point(kp_low(kf) + kp)
          r = cur_knot%rcoord; p = 0d0; z = cur_knot%zcoord

          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

          call common_triangles(prev_knot, cur_knot, common_tri)

          lr = r - r_prev
          lz = z - z_prev
          perps = mesh_element(common_tri(:))%det_3 / hypot(lr, lz) * 0.5d0
          Bnpsi = Bnflux(minval(common_tri), mod(mesh_element(minval(common_tri))% &
               knot_h, 3) + 1) / (r + r_prev) * 2 / hypot(lr, lz) * &
               (psi(kf+1) - psi(kf-1)) / sum(perps)

          x(kp) = -dpres0_dpsi(kf) * Bnpsi

          a(kp) = ((Br + Br_prev) * 0.5d0 * lr + (Bz + Bz_prev) * 0.5d0 * lz) / &
               (lr ** 2 + lz ** 2)

          b(kp) = imun * n * (Bp / r + Bp_prev / r_prev) * 0.5d0
       end do ! kp

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       call assemble_sparse(nkpol, d, du, nz, irow, icol, aval)
       call sparse_solve(nkpol, nkpol, nz, irow, icol, aval, x)

       if (kf == 1) then ! first point on axis before actual output
          presn(1) = sum(x) / size(x)
          write(1, *) psi(0), dens(0), temp(0), pres0(0), &
               real(sum(x) / size(x)), aimag(sum(x) / size(x))
       end if
       do kp = 1, kp_max(kf)
          presn(kp_low(kf) + kp) = x(kp)
          write(1, *) psi(kf), dens(kf), temp(kf), pres0(kf), &
               real(x(kp)), aimag(x(kp))
       end do
    end do ! kf
    do kp = kp_low(nflux+1) + 1, npoint ! write zeroes in remaining points until end
       write(1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine compute_presn


  !> Map edge symbols to integer indices.
  !>
  !> @param elem the triangle for which indices are to be obtained
  !> @param li knot indices for base (1) and tip (2) of edge i
  !> @param lo knot indices for base (1) and tip (2) of edge o
  !> @param lf knot indices for base (1) and tip (2) of edge f
  !> @param ei index of edge i, e.g. for #bnflux
  !> @param eo index of edge o, e.g. for #bnflux
  !> @param ef index of edge f, e.g. for #bnflux
  !> @param orient true if edge f lies on the outer flux surface, false otherwise
  !>
  !> It is assumed that the knots are numbered in ascending order starting from the
  !> magnetic axis and going counter-clockwise around each flux surface. Furthermore it
  !> is assumed that edge 1 goes from knot 1 to knot 2, edge 2 from knot 2 to knot 3 and
  !> edge 3 from knot 3 to knot 1.

  subroutine get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    logical, intent(out) :: orient
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
    ! last triangle in strip if indices not next to each other
    closing_loop = abs(elem%i_knot(i1) - elem%i_knot(i2)) /= 1
    if ((elem%i_knot(i1) > elem%i_knot(i2)) .neqv. closing_loop) then
       ! i1 is next counter-clockwise
       knot_i = i1
       knot_o = i2
    else
       ! i2 is next counter-clockwise
       knot_i = i2
       knot_o = i1
    end if
    knot_f = elem%knot_h
    i_knot_diff = elem%i_knot - elem%i_knot(knot_f)
    if (all(i_knot_diff >= 0)) then
       ! knot_f lies on inner surface
       orient = .true.
       ei = knot_f
       eo = knot_i
       ef = knot_o
       li = (/ elem%i_knot(knot_f), elem%i_knot(knot_o) /)
       lo = (/ elem%i_knot(knot_i), elem%i_knot(knot_f) /)
       lf = (/ elem%i_knot(knot_o), elem%i_knot(knot_i) /)
    else if (all(i_knot_diff <= 0)) then
       ! knot_f lies on outer surface
       orient = .false.
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


  !> Compute the fluxes of a complex vector field through each triangle, separately for
  !> each flux surface.
  !>
  !> @param assemble_flux_coeff subroutine that assembles the system of linear equations
  !> edge by edge, e.g. assemble_currn_coeff()
  !> @param assign_flux subroutine that assigns values to the fluxes through each triangle
  !> after the system of linear equations is solved, e.g. assign_currn()
  !> @param pol_flux poloidal flux components, e.g. #jnflux - first index refers to the
  !> triangle, second index refers to the edge on that triangle, as per
  !> #mesh_mod::mesh_element
  !> @param tor_comp toroidal field components, e.g. #jnphi - index refers to the triangle
  !> @param outfile if present, write \p pol_flux and area-weighted \p tor_comp to the
  !> given filename
  !>
  !> For each flux surface, this subroutine loops through every triangle and calls \p
  !> assemble_flux_coeff for edges f, i, o (in that order). Then, it assembles and solves
  !> the resulting sparse system of linear equations. With the known solution, it loops
  !> over the triangle strip again and assigns values to \p pol_flux and \p tor_comp via
  !> \p assign_flux. Note that \p assemble_flux_coeff can also assign flux components that
  !> are already known, usually the one through edge f.

  subroutine compute_triangle_flux(assemble_flux_coeff, assign_flux, pol_flux, tor_comp, &
       outfile)
    procedure(sub_assemble_flux_coeff) :: assemble_flux_coeff
    procedure(sub_assign_flux) :: assign_flux
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)
    character(len = 1024), intent(in), optional :: outfile

    complex(dp), dimension(2 * nkpol) :: x, d, du
    integer :: kf, kt
    integer :: nz
    integer, dimension(4 * nkpol) :: irow, icol
    complex(dp), dimension(4 * nkpol) :: aval
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    complex(dp) :: tor_flux

    if (present(outfile)) open(1, file = outfile, recl = 1024)
    do kf = 1, nflux ! loop through flux surfaces
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          call assemble_flux_coeff(x, d, du, kf, kt, elem, lf, ef, 'f')
          call assemble_flux_coeff(x, d, du, kf, kt, elem, li, ei, 'i')
          call assemble_flux_coeff(x, d, du, kf, kt, elem, lo, eo, 'o')
       end do
       call assemble_sparse(kt_max(kf), d(:kt_max(kf)), du(:kt_max(kf)), nz, &
            irow(:(2*kt_max(kf))), icol(:(2*kt_max(kf))), aval(:(2*kt_max(kf))))
       call sparse_solve(kt_max(kf), kt_max(kf), nz, irow(:nz), icol(:nz), aval(:nz), &
            x(:kt_max(kf)))
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          call assign_flux(x, kf, kt, ei, eo, ef)
          if (present(outfile)) then
             tor_flux = tor_comp(kt_low(kf) + kt) * elem%det_3 * 0.5d0
             write(1, *) &
               real(pol_flux(kt_low(kf) + kt, 1)), aimag(pol_flux(kt_low(kf) + kt, 1)), &
               real(pol_flux(kt_low(kf) + kt, 2)), aimag(pol_flux(kt_low(kf) + kt, 2)), &
               real(pol_flux(kt_low(kf) + kt, 3)), aimag(pol_flux(kt_low(kf) + kt, 3)), &
               real(tor_flux), aimag(tor_flux)
          end if
       end do
    end do

    if (present(outfile)) then
       do kt = kt_low(nflux+1) + 1, ntri
          write(1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
       end do
       close(1)
    end if
  end subroutine compute_triangle_flux

  subroutine assemble_currn_coeff(x, d, du, kf, kt, elem, lk, ek, edge_name)

    complex(dp), dimension(:), intent(inout) :: x, d, du
    integer, intent(in) :: kf, kt
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(in) :: lk
    integer, intent(in) :: ek
    character, intent(in) :: edge_name

    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    real(dp) :: Deltapsi
    type(knot) :: base, tip

    p = 0d0
    ! use midpoint of edge
    base = mesh_point(lk(1))
    tip = mesh_point(lk(2))
    r = (base%rcoord + tip%rcoord) * 0.5d0
    z = (base%zcoord + tip%zcoord) * 0.5d0
    call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
    Deltapsi = psi(kf) - psi(kf-1)
    select case (edge_name)
    case ('f')
       ! first term on source side: flux through edge f
       jnflux(kt_low(kf) + kt, ek) = clight * r / Bp * (presn(lk(2)) - presn(lk(1))) + &
            j0phi(kt_low(kf) + kt) / Bp * Bnflux(kt_low(kf) + kt, ek)
       x(kt) = -jnflux(kt_low(kf) + kt, ek)
    case ('i')
       ! diagonal matrix element
       d(kt) = -1d0 - imun * n * elem%det_3 * 0.25d0 * Bp / Deltapsi
       ! additional term from edge i on source side
       x(kt) = x(kt) - imun * n * elem%det_3 * 0.25d0 * (clight * r / (-Deltapsi) * &
            (presn(lk(2)) - presn(lk(1)) - Bnphi(kt_low(kf) + kt) / Bp * &
            (pres0(kf) - pres0(kf-1))) + j0phi(kt_low(kf) + kt) * &
            (Bnphi(kt_low(kf) + kt) / Bp + Bnflux(kt_low(kf) + kt, ek) / r / (-Deltapsi)))
    case ('o')
       ! superdiagonal matrix element
       du(kt) = 1d0 - imun * n * elem%det_3 * 0.25d0 * Bp / Deltapsi
       ! additional term from edge o on source side
       x(kt) = x(kt) - imun * n * elem%det_3 * 0.25d0 * (clight * r / Deltapsi * &
            (presn(lk(2)) - presn(lk(1)) - Bnphi(kt_low(kf) + kt) / Bp * &
            (pres0(kf-1) - pres0(kf))) + j0phi(kt_low(kf) + kt) * &
            (Bnphi(kt_low(kf) + kt) / Bp + Bnflux(kt_low(kf) + kt, ek) / r / Deltapsi))
    end select
  end subroutine assemble_currn_coeff

  subroutine assign_currn(x, kf, kt, ei, eo, ef)
    complex(dp), dimension(:), intent(in) :: x
    integer, intent(in) :: kf, kt
    integer, intent(in) :: ei, eo, ef
    jnflux(kt_low(kf) + kt, ei) = -x(kt)
    jnflux(kt_low(kf) + kt, eo) = x(mod(kt, kt_max(kf) + 1))
    jnphi(kt_low(kf) + kt) = sum(jnflux(kt_low(kf) + kt, :)) * imun / n / &
         mesh_element(kt_low(kf) + kt)%det_3 * 2d0
  end subroutine assign_currn

  subroutine compute_currn
    call compute_triangle_flux(assemble_currn_coeff, assign_currn, jnflux, jnphi, currn_file)
  end subroutine compute_currn

  subroutine compute_Bn_nonres
    integer :: kf, kt
    type(triangle) :: elem
    type(knot) :: base, tip
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    integer :: common_tri(2)
    real(dp) :: lr, lz, Deltapsi, perps(2)
    complex(dp) :: Bnpsi, Bnflux_avg, tor_flux_avg, tor_flux_diff
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    tor_flux_avg = (0d0, 0d0); tor_flux_diff = (0d0, 0d0)  ! to suppress -Wuninitialized
    p = 0d0
    open(1, file = 'Bn_nonres.dat', recl = 1024)
    do kf = 1, nflux ! loop through flux surfaces
       Bnflux_avg = (0d0, 0d0)
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          if (orient) then
             Deltapsi = psi(kf+1) - psi(kf-1)
          else
             Deltapsi = psi(kf-2) - psi(kf)
          end if
          ! use midpoint of edge f
          base = mesh_point(lf(1))
          tip = mesh_point(lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          z = (base%zcoord + tip%zcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord
          call common_triangles(base, tip, common_tri)
          perps = mesh_element(common_tri(:))%det_3 / hypot(lr, lz) * 0.5d0
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          Bnpsi = Bp / r
          Bnflux(kt_low(kf) + kt, ef) = Bnpsi * r * hypot(lr, lz) * sum(perps) / Deltapsi
          Bnflux_avg = Bnflux_avg + Bnflux(kt_low(kf) + kt, ef)
          Bnphi(kt_low(kf) + kt) = imun / n * Bnflux(kt_low(kf) + kt, ef) &
               / elem%det_3 * 2d0
       end do
       Bnflux_avg = Bnflux_avg / kt_max(kf)
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          Bnflux(kt_low(kf) + kt, ei) = -2d0 * abs(Bnflux_avg)
          Bnflux(kt_low(kf) + kt, eo) =  2d0 * abs(Bnflux_avg)
          if (quad_avg .and. (kf > 1)) then
             select case (mod(kt, 2))
             case (1)
                tor_flux_avg = 0.25d0 * (Bnphi(kt_low(kf) + kt) * elem%det_3 + &
                     Bnphi(kt_low(kf) + kt + 1) * mesh_element(kt_low(kf) + kt + 1)%det_3)
                tor_flux_diff = 0.25d0 * (-Bnphi(kt_low(kf) + kt) * elem%det_3 + &
                     Bnphi(kt_low(kf) + kt + 1) * mesh_element(kt_low(kf) + kt + 1)%det_3)
                Bnphi(kt_low(kf) + kt) = tor_flux_avg / elem%det_3 * 2d0
                Bnflux(kt_low(kf) + kt, eo) = Bnflux(kt_low(kf) + kt, eo) - imun * n * &
                     tor_flux_diff
             case (0)
                Bnphi(kt_low(kf) + kt) = tor_flux_avg / elem%det_3 * 2d0
                Bnflux(kt_low(kf) + kt, ei) = Bnflux(kt_low(kf) + kt, ei) + imun * n * &
                     tor_flux_diff
             end select
          end if
          write (1, *) &
               real(Bnflux(kt_low(kf) + kt, 1)), real(Bnflux(kt_low(kf) + kt, 2)), &
               real(Bnflux(kt_low(kf) + kt, 3)), aimag(Bnphi(kt_low(kf) + kt))
       end do
    end do
    close(1)

    call check_div_free(Bnflux, Bnphi, 1d-09, 'non-resonant B_n')
  end subroutine compute_Bn_nonres

  subroutine dump_triangle_flux(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: ktri
    type(triangle) :: elem
    type(knot) :: tri(3)
    real(dp) :: length(3)
    complex(dp) :: pol_comp_r, pol_comp_z
    real(dp) :: r, z

    open(1, file = outfile, recl = 1024)
    do ktri = 1, kt_low(nflux+1)
       elem = mesh_element(ktri)
       tri = mesh_point(elem%i_knot(:))
       length(1) = hypot(tri(2)%rcoord - tri(1)%rcoord, tri(2)%zcoord - tri(1)%zcoord)
       length(2) = hypot(tri(3)%rcoord - tri(2)%rcoord, tri(3)%zcoord - tri(2)%zcoord)
       length(3) = hypot(tri(1)%rcoord - tri(3)%rcoord, tri(1)%zcoord - tri(3)%zcoord)
       call ring_centered_avg_coord(elem, r, z)
       pol_comp_r = 2d0 / elem%det_3 * ( &
            pol_flux(ktri, 1) * (r - tri(3)%rcoord) * length(1) + &
            pol_flux(ktri, 2) * (r - tri(1)%rcoord) * length(2) + &
            pol_flux(ktri, 3) * (r - tri(2)%rcoord) * length(3))
       pol_comp_z = 2d0 / elem%det_3 * ( &
            pol_flux(ktri, 1) * (z - tri(3)%zcoord) * length(1) + &
            pol_flux(ktri, 2) * (z - tri(1)%zcoord) * length(2) + &
            pol_flux(ktri, 3) * (z - tri(2)%zcoord) * length(3))
       write (1, *) r, z, real(pol_comp_r), aimag(pol_comp_r), real(tor_comp(ktri)), &
            aimag(tor_comp(ktri)), real(pol_comp_z), aimag(pol_comp_z)
    end do
    r = mesh_point(1)%rcoord
    z = mesh_point(1)%zcoord
    do ktri = kt_low(nflux+1) + 1, ntri
       write (1, *) r, z, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine dump_triangle_flux

end module magdif
