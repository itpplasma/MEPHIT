module magdif
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi, ev2erg                               ! PRELOAD/SRC/orbit_mod.f90
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       mesh_element_rmp, bphicovar, knot, triangle
  use for_macrostep, only: t_min, d_min
  use sparse_mod, only: sparse_solve, sparse_matmul
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

  !> \f$ R \$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0r(:,:)

  !> \f$ \phi \$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0phi(:,:)

  !> \f$ Z \$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0z(:,:)

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
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: j0phi(:,:)

  !> Number of knots on the flux surface given by the array index.
  !>
  !> The array index ranges from 1 for the innermost flux surface to
  !> #magdif_config::nflux+2 for consistency with #kp_low.
  integer, allocatable :: kp_max(:)

  !> Number of triangles inside the flux surface given by the array index.
  !>
  !> The array index ranges from 1 for the innermost flux surface to
  !> #magdif_config::nflux+1 for consistency with #kt_low.
  integer, allocatable :: kt_max(:)

  !> Global index of the last knot of the previous flux surface given by the array index.
  !>
  !> The global index of knots in #mesh_mod::mesh_point on the flux surface kf runs from
  !> #kp_low (kf)+1 to #kp_low (kf)+#kp_max (kf), so #kp_low is determined by cumulatively
  !> adding consecutive values of #kp_max. The array index ranges from 1, giving the
  !> global index of the knot on the magnetic axis (which has to be 1), to
  !> #magdif_config::nflux+2, giving the last knot on the flux surface just outside the
  !> last closed flux surface. The latter is needed for some interpolations.
  integer, allocatable :: kp_low(:)

  !> Global index of the last triangle of the previous flux surface given by the array
  !> index.
  !>
  !> The global index of triangles in #mesh_mod::mesh_element inside the flux surface kf
  !> runs from #kt_low (kf)+1 to #kt_low (kf)+#kt_max (kf), so #kt_low is determined by
  !> cumulatively adding consecutive values of #kt_max. The array index ranges from 1,
  !> giving the global index of the non-existent triangle on the magnetic axis (which is
  !> therefore 0), to #magdif_config::nflux+1, giving the last triangle just outside the
  !> last closed flux surface. The latter is needed for some interpolations.
  integer, allocatable :: kt_low(:)

  !> Distance of magnetic axis from center \f$ R_{0} \f$ in cm.
  real(dp), parameter :: R0 = 172.74467899999999d0

  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

contains

  !> Initialize magdif module
  subroutine magdif_init
    integer :: kf, kp

    call init_indices
    call read_mesh
    call init_flux_variables
    call cache_equilibrium_field
    call compute_safety_factor
    call compute_j0phi

    open(1, file = fluxvar_file, recl = longlines)
    write (1, *) psi(0), 0d0, dens(0), temp(0), pres0(0)
    do kf = 1, nflux
       do kp = 1, kp_max(kf)
          write (1, *) psi(kf), q(kf), dens(kf), temp(kf), pres0(kf)
       end do
    end do
    do kp = kp_low(nflux+1) + 1, npoint
       write (1, *) 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)

    if (nonres) then
       call compute_Bn_nonres
    else
       call read_Bn(Bn_vac_file)
    end if
    Bnflux_vac = Bnflux
    Bnphi_vac = Bnphi
    if (log_debug) then
       call write_triangle_flux(Bnflux_vac, Bnphi_vac, Bn_nonres_file, .false.)
       call write_triangle_flux(Bnflux_vac, Bnphi_vac, decorate_filename(Bn_nonres_file, &
            'plot_'), .true.)
    end if
    if (log_info) write(logfile, *) 'magdif initialized'
  end subroutine magdif_init

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    if (allocated(q)) deallocate(q)
    if (allocated(pres0)) deallocate(pres0)
    if (allocated(dpres0_dpsi)) deallocate(dpres0_dpsi)
    if (allocated(dens)) deallocate(dens)
    if (allocated(temp)) deallocate(temp)
    if (allocated(psi)) deallocate(psi)
    if (allocated(B0r)) deallocate(B0r)
    if (allocated(B0phi)) deallocate(B0phi)
    if (allocated(B0z)) deallocate(B0z)
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
    call compute_Bn
    call read_Bn(Bn_file)
  end subroutine magdif_single

  subroutine magdif_direct
    integer :: kiter
    complex(dp) :: Bnflux_sum(ntri, 3)
    complex(dp) :: Bnphi_sum(ntri)
    complex(dp) :: Bnflux_prev(ntri, 3)
    complex(dp) :: Bnphi_prev(ntri)
    complex(dp) :: Bnflux_diff(ntri, 3)
    complex(dp) :: Bnphi_diff(ntri)

    Bnflux_sum = 0.0d0
    Bnphi_sum = 0.0d0
    do kiter = 1, niter
       if (log_info) write(logfile, *) 'Iteration ', kiter, ' of ', niter
       call compute_presn     ! compute pressure based on previous perturbation field
       call compute_currn     ! compute currents based on previous perturbation field
       Bnflux_prev = Bnflux
       Bnphi_prev = Bnphi
       call compute_Bn        ! use field code to generate new field from currents
       call read_Bn(Bn_file)  ! read new bnflux from field code
       Bnflux_diff = Bnflux - Bnflux_prev
       Bnphi_diff = Bnphi - Bnphi_prev
       call write_triangle_flux(Bnflux_diff, Bnphi_diff, &
            decorate_filename(Bn_abs_err_file, '', kiter), .false.)
       call write_triangle_flux(Bnflux_diff, Bnphi_diff, &
            decorate_filename(Bn_abs_err_file, 'plot_', kiter), .true.)
       call write_triangle_flux(Bnflux_diff / Bnflux_prev, Bnphi_diff / Bnphi_prev, &
            decorate_filename(Bn_rel_err_file, '', kiter), .false.)
       call write_triangle_flux(Bnflux, Bnphi, &
            decorate_filename(Bn_file, 'plot_', kiter), .true.)
       call write_triangle_flux(jnflux, jnphi, &
            decorate_filename(currn_file, 'plot_', kiter), .true.)
       call write_nodal_dof(presn, decorate_filename(presn_file, '', kiter))
       Bnflux_sum = Bnflux_sum + Bnflux
       Bnphi_sum = Bnphi_sum + Bnphi
    end do
    call write_triangle_flux(Bnflux_sum, Bnphi_sum, Bn_sum_file, .false.)
    call write_triangle_flux(Bnflux_sum, Bnphi_sum, decorate_filename(Bn_sum_file, &
         'plot_'), .true.)

    Bnflux = Bnflux_sum
    Bnphi = Bnphi_sum
    call compute_presn
    call write_nodal_dof(presn, presn_file)
    call compute_currn
    call write_triangle_flux(jnflux, jnphi, decorate_filename(currn_file, 'plot_'), .true.)
  end subroutine magdif_direct

  !> Allocates and initializes #kp_low, #kp_max, #kt_low and #kt_max based on the values
  !> of #magdif_config::nflux and #magdif_config::nkpol. Deallocation is done in
  !> magdif_cleanup().
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

    if (log_info) then
       write (logfile, *) "Number of points up to LCFS: ", kp_low(nflux+1), &
            "Number of triangles up to LCFS: ", kt_low(nflux+1)
    end if
  end subroutine init_indices

  !> Reads mesh points and triangles.
  !>
  !> #mesh_mod::npoint and #mesh_mod::mesh_point are read directly from an unformatted
  !> #magdif_config::point_file, while #mesh_mod::ntri and #mesh_mod::mesh_element are
  !> read directly from an unformatted #magdif_config::tri_file. #presn, #bnflux, #bnphi,
  !> #bnflux_vac, #bnphi_vac, #j0phi, #jnphi and #jnflux are allocated and initialized to
  !> zero. Deallocation is done in magdif_cleanup().
  subroutine read_mesh
    open(1, file = point_file, form = 'unformatted')
    read(1) npoint
    if (log_info) write(logfile, *) 'npoint = ', npoint
    allocate(mesh_point(npoint))
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
    allocate(B0r(ntri, 3))
    allocate(B0phi(ntri, 3))
    allocate(B0z(ntri, 3))
    allocate(presn(npoint))
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    allocate(Bnflux_vac(ntri, 3))
    allocate(Bnphi_vac(ntri))
    allocate(jnphi(ntri))
    allocate(j0phi(ntri, 3))
    allocate(jnflux(ntri, 3))

    B0r = 0d0
    B0phi = 0d0
    B0z = 0d0
    presn = 0d0
    Bnflux = 0d0
    Bnphi = 0d0
    Bnflux_vac = 0d0
    Bnphi_vac = 0d0
    jnphi = 0d0
    j0phi = 0d0
    jnflux = 0d0
  end subroutine read_mesh

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi via an external program. No data
  !> is read yet; this is done by read_bn().
  !>
  !> Currently, this subroutine Calls FreeFem++ via shell script maxwell.sh in the current
  !> directory and handles the exit code. For further information see maxwell.sh and the
  !> script called therein.
  subroutine compute_Bn
    integer :: stat = 0, dummy = 0
    call execute_command_line("./maxwell.sh", exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       if (log_err) write(logfile, *) 'FreeFem++ failed with exit code ', stat
       stop 'FreeFem++ failed'
    end if
  end subroutine compute_Bn

  !> Checks if divergence-freeness of the given vector field is fulfilled on each
  !> triangle, otherwise halts the program.
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

  subroutine check_redundant_edges(pol_quant, comp_factor, name)
    complex(dp), intent(in) :: pol_quant(:,:)
    real(dp), intent(in) :: comp_factor
    character(len = *), intent(in) :: name
    integer :: ktri, ktri_adj, ke, ke_adj
    logical :: checked(ntri, 3)
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    character(len = len_trim(name) + 30) :: err_msg
    err_msg = trim(name) // ': inconsistent redundant edges'

    checked = .false.
    do ktri = 1, kt_low(nflux+1)
       elem = mesh_element(ktri)
       call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
       do ke = 1, 3
          if (ktri > kt_low(nflux) .and. ke == ef) cycle
          if (.not. checked(ktri, ke)) then
             checked(ktri, ke) = .true.
             ktri_adj = mesh_element(ktri)%neighbour(ke)
             ke_adj = mesh_element(ktri)%neighbour_edge(ke)
             checked(ktri_adj, ke_adj) = .true.
             if (pol_quant(ktri, ke) /= comp_factor * pol_quant(ktri_adj, ke_adj)) then
                if (log_err) write (logfile, *) err_msg, ' - ', &
                     name, '(', ktri, ',', ke, ') = ', pol_quant(ktri, ke), &
                     name, '(', ktri_adj, ',', ke_adj, ') = ', pol_quant(ktri_adj, ke_adj)
                stop err_msg
             end if
          end if
       end do
    end do
  end subroutine check_redundant_edges

  !> Reads fluxes of perturbation field and checks divergence-freeness.
  !>
  !> @param filename name of the formatted file containing the data
  !>
  !> Line numbers in the given file correspond to the global triangle index of
  !> #mesh_mod::mesh_element. The content of each line is read into #bnflux and #bnphi
  !> with numbering of edges in #bnflux as in #mesh_mod::mesh_element and the imaginary
  !> part of each value immediately following the real part.
  subroutine read_Bn(filename)
    character(len = 1024) :: filename
    integer :: kt
    real(dp) :: dummy8(8)

    open(1, file = filename, recl = longlines)
    do kt = 1, ntri
       read(1, *) dummy8
       Bnflux(kt,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(kt,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(kt,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(kt) = cmplx(dummy8(7), dummy8(8), dp) / mesh_element(kt)%det_3 * 2d0
    end do
    close(1)

    call check_div_free(Bnflux, Bnphi, 1d-3, 'B_n')
    call check_redundant_edges(Bnflux, -1d0, 'B_n')
  end subroutine read_Bn

  !> Allocates and computes the safety factor #q. Deallocation is done in
  !> magdif_cleanup().
  subroutine compute_safety_factor
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
  end subroutine compute_safety_factor

  !> Computes the "weighted" centroid for a triangle so that it is approximately
  !> equidistant between the enclosing flux surfaces, independent of triangle orientation.
  !>
  !> @param elem the triangle for which the centroid is to be computed
  !> @param r radial cylindrical coordinate of the centroid
  !> @param z axial cylindrical coordinate of the centroid
  !>
  !> Depending on the orientation of the triangle (see also \p orient of
  !> get_labeled_edges()), two knots lie on the inner flux surface and one on the outer
  !> one, or vice versa. A simple arithmetic mean of the three knots' coordinates would
  !> place the centroid closer to the inner flux surface for one orientation and closer
  !> to the outer one for the other. To counteract this, the "lonely" knot is counted
  !> twice in the averaging procedure, i.e. with double weighting.
  subroutine ring_centered_avg_coord(elem, r, z)
    type(triangle), intent(in) :: elem
    real(dp), intent(out) :: r, z
    type(knot), dimension(3) :: knots

    knots = mesh_point(elem%i_knot)
    r = (sum(knots%rcoord) + knots(elem%knot_h)%rcoord) * 0.25d0
    z = (sum(knots%zcoord) + knots(elem%knot_h)%zcoord) * 0.25d0
  end subroutine ring_centered_avg_coord

  !> Initializes quantities that are constant on each flux surface. The safety factor #q
  !> is initialized separately in compute_safety_factor().
  !>
  !> #psimin, #psimax and #psi are initialized from #mesh_mod::mesh_point::psi_pol.
  !> #pres0, #dpres0_dpsi, #dens and #temp are given an assumed profile based on #psimin,
  !> #psimax, #psi, #magdif_conf::di0, #magdif_conf::d_min, #magdif_conf::ti0 and
  !> #magdif_conf::t_min.
  subroutine init_flux_variables
    integer :: kf
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

    ! magnetic axis at k == 0 is not counted as flux surface
    psi(0) = mesh_point(1)%psi_pol

    do kf = 1, nflux+1
       ! average over the loop to smooth out numerical errors
       psi(kf) = sum(mesh_point((kp_low(kf) + 1):kp_low(kf+1))%psi_pol) / kp_max(kf)
    end do
    dens = (psi - psimin) / psimax * di0 + d_min
    temp = (psi - psimin) / psimax * ti0 + t_min
    pres0 = dens * temp * ev2erg
    dpres0_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
  end subroutine init_flux_variables

  subroutine cache_equilibrium_field
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kf, kt, ke
    type(triangle) :: elem
    type(knot) :: base, tip

    p = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          do ke = 1, 3
             base = mesh_point(elem%i_knot(ke))
             tip = mesh_point(elem%i_knot(mod(ke, 3) + 1))
             r = (base%rcoord + tip%rcoord) * 0.5d0
             z = (base%zcoord + tip%zcoord) * 0.5d0
             call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
                  dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
             B0r(kt_low(kf) + kt, ke) = Br
             B0phi(kt_low(kf) + kt, ke) = Bp
             B0z(kt_low(kf) + kt, ke) = Bz
          end do
       end do
    end do
  end subroutine cache_equilibrium_field

  !> Computes equilibrium current density #j0phi from given equilibrium magnetic field and
  !> assumed equilibrium pressure #pres0.
  !>
  !> This step is necessary because equilibrium pressure is not given experimentally as is
  !> \f$ \vec{B}_{0} \f$; arbitrary values are assumed. Consistency of MHD equilibrium is
  !> necessary in the derivation, while Ampere's equation is not used.
  subroutine compute_j0phi
    integer :: kf, kt
    real(dp) :: Bpol2, Btor2
    real(dp) :: B2avg(nflux)
    real(dp) :: B2avg_half(nflux)
    real(dp) :: r
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient

    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0r(kt_low(kf) + kt, ef) ** 2 + &
               B0phi(kt_low(kf) + kt, ef) ** 2 + B0z(kt_low(kf) + kt, ef) ** 2
          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0r(kt_low(kf) + kt, ei) ** 2 + &
               B0phi(kt_low(kf) + kt, ei) ** 2 + B0z(kt_low(kf) + kt, ei) ** 2
       end do
       B2avg(kf) = B2avg(kf) / kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / kt_max(kf)

       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          Bpol2 = B0r(kt_low(kf) + kt, ef) ** 2 + B0z(kt_low(kf) + kt, ef) ** 2
          Btor2 = B0phi(kt_low(kf) + kt, ef) ** 2
          if (kf > 1 .and. .not. orient) then
             j0phi(kt_low(kf) + kt, ef) = clight * dpres0_dpsi(kf-1) * (Bpol2 / &
                  B2avg(kf-1) + (Bpol2 - Btor2) / (Bpol2 + Btor2)) * r
          else
             j0phi(kt_low(kf) + kt, ef) = clight * dpres0_dpsi(kf) * (Bpol2 / &
                  B2avg(kf) + (Bpol2 - Btor2) / (Bpol2 + Btor2)) * r
          end if
          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          Bpol2 = B0r(kt_low(kf) + kt, ei) ** 2 + B0z(kt_low(kf) + kt, ei) ** 2
          Btor2 = B0phi(kt_low(kf) + kt, ei) ** 2
          j0phi(kt_low(kf) + kt, ei) = clight * (pres0(kf) - pres0(kf-1)) / &
               (psi(kf) - psi(kf-1)) * (Btor2 / B2avg_half(kf) + &
               (Bpol2 - Btor2) / (Bpol2 + Btor2)) * r
          r = sum(mesh_point(lo(:))%rcoord) * 0.5d0
          Bpol2 = B0r(kt_low(kf) + kt, eo) ** 2 + B0z(kt_low(kf) + kt, eo) ** 2
          Btor2 = B0phi(kt_low(kf) + kt, eo) ** 2
          j0phi(kt_low(kf) + kt, eo) = clight * (pres0(kf) - pres0(kf-1)) / &
               (psi(kf) - psi(kf-1)) * (Btor2 / B2avg_half(kf) + &
               (Bpol2 - Btor2) / (Bpol2 + Btor2)) * r
       end do
    end do

    call check_redundant_edges(cmplx(j0phi, 0d0, dp), 1d0, 'j0phi')
  end subroutine compute_j0phi

  !> Assembles a sparse matrix in coordinate list (COO) representation for use with
  !> sparse_mod::sparse_solve().
  !>
  !> @param nrow number \f$ n \f$ of rows in the matrix
  !> @param d \f$ n \f$ diagnonal elements of stiffness matrix \f$ A \f$
  !> @param du \f$ n-1 \f$ superdiagonal elements of stiffness matrix \f$ A \f$ and
  !> \f$ A_{n, 1} \f$ (lower left corner)
  !> @param nz number of non-zero entries (2*nrow)
  !> @param irow nz row indices of non-zero entries
  !> @param icol nz column indices of non-zero entries
  !> @param aval nz values of non-zero entries
  !>
  !> The input is a stiffness matrix \f$ K \f$ with non-zero entries on the main diagonal,
  !> the upper diagonal and, due to periodicity, in the lower left corner. This shape
  !> results from the problems in compute_presn() and compute_currn().
  subroutine assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    integer, intent(in)  :: nrow
    complex(dp), intent(in)  :: d(nrow)
    complex(dp), intent(in)  :: du(nrow)
    integer, intent(out) :: nz
    integer, intent(out) :: irow(2*nrow), icol(2*nrow)
    complex(dp), intent(out) :: aval(2*nrow)

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

  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  !>
  !> #psi, #dens, #temp, #pres0 and real and imaginary part of #presn are written, in that
  !> order, to #magdif_conf::presn_file, where line number corresponds to the knot index
  !> in #mesh_mod::mesh_point.
  subroutine compute_presn
    real(dp) :: r
    real(dp) :: lr, lz  ! edge vector components
    complex(dp) :: Bnpsi
    complex(dp), dimension(nkpol) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, kp
    integer :: nz
    integer, dimension(2*nkpol) :: irow, icol
    complex(dp), dimension(2*nkpol) :: aval
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    type(knot) :: base, tip
    integer :: common_tri(2)
    real(dp) :: perps(2)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux
inner: do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          if (.not. orient) cycle inner
          common_tri = (/ kt_low(kf) + kt, elem%neighbour(ef) /)
          ! use midpoint of edge f
          base = mesh_point(lf(1))
          tip = mesh_point(lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord

          perps = mesh_element(common_tri(:))%det_3 / hypot(lr, lz) * 0.5d0
          Bnpsi = Bnflux(kt_low(kf) + kt, ef) / r / hypot(lr, lz) * &
               (psi(kf+1) - psi(kf-1)) / sum(perps)

          if (kf == 1) then
             kp = kt
          else
             kp = kt / 2 + mod(kt, 2)
          end if

          x(kp) = -dpres0_dpsi(kf) * Bnpsi

          a(kp) = (B0r(kt_low(kf) + kt, ef) * lr + B0z(kt_low(kf) + kt, ef) * lz) / &
               (lr ** 2 + lz ** 2)

          b(kp) = imun * (n + imun * damp) * B0phi(kt_low(kf) + kt, ef) / r
       end do inner

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       call assemble_sparse(nkpol, d, du, nz, irow, icol, aval)
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(nkpol, nkpol, nz, irow, icol, aval, x)
       call sparse_matmul(nkpol, nkpol, irow, icol, aval, x, resid)
       resid = resid - inhom
       rel_err = abs(resid) / abs(inhom)
       max_rel_err = max(max_rel_err, maxval(rel_err))
       avg_rel_err = avg_rel_err + sum(rel_err)

       if (kf == 1) then ! first point on axis before actual output
          presn(1) = sum(x) / size(x)
       end if
       do kp = 1, kp_max(kf)
          presn(kp_low(kf) + kp) = x(kp)
       end do
    end do

    avg_rel_err = avg_rel_err / sum(kp_max(1:nflux))
    if (log_debug) write (logfile, *) 'compute_presn: diagonalization' // &
         ' max_rel_err = ', max_rel_err, ' avg_rel_err = ', avg_rel_err
    if (allocated(resid)) deallocate(resid)

    call write_nodal_dof(presn, presn_file)
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


  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #magdif_conf::currn_file.
  subroutine compute_currn
    complex(dp), dimension(2 * nkpol) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(2 * nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt
    integer :: nz
    integer, dimension(4 * nkpol) :: irow, icol
    complex(dp), dimension(4 * nkpol) :: aval
    type(triangle) :: elem
    real(dp) :: area
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    real(dp) :: Deltapsi
    real(dp) :: r

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux ! loop through flux surfaces
       Deltapsi = psi(kf) - psi(kf-1)
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          area = elem%det_3 * 0.5d0
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          ! first term on source side: flux through edge f
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          jnflux(kt_low(kf) + kt, ef) = clight * r / B0phi(kt_low(kf) + kt, ef) * &
          (presn(lf(2)) - presn(lf(1))) + j0phi(kt_low(kf) + kt, ef) / &
          B0phi(kt_low(kf) + kt, ef) * Bnflux(kt_low(kf) + kt, ef)
          x(kt) = -jnflux(kt_low(kf) + kt, ef)
          ! diagonal matrix element - edge i
          d(kt) = -1d0 - imun * (n + imun * damp) * area * 0.5d0 * &
               B0phi(kt_low(kf) + kt, ei) / Deltapsi
          ! additional term from edge i on source side
          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          x(kt) = x(kt) - imun * n * area * 0.5d0 * (clight * r / (-Deltapsi) * &
               (presn(li(2)) - presn(li(1)) - Bnphi(kt_low(kf) + kt) / &
               B0phi(kt_low(kf) + kt, ei) * (pres0(kf) - pres0(kf-1))) + &
               j0phi(kt_low(kf) + kt, ei) * (Bnphi(kt_low(kf) + kt) / &
               B0phi(kt_low(kf) + kt, ei) + Bnflux(kt_low(kf) + kt, ei) / r / (-Deltapsi)))
          ! superdiagonal matrix element - edge o
          du(kt) = 1d0 - imun * (n + imun * damp) * area * 0.5d0 * &
               B0phi(kt_low(kf) + kt, eo) / Deltapsi
          ! additional term from edge o on source side
          r = sum(mesh_point(lo(:))%rcoord) * 0.5d0
          x(kt) = x(kt) - imun * n * area * 0.5d0 * (clight * r / Deltapsi * &
               (presn(lo(2)) - presn(lo(1)) - Bnphi(kt_low(kf) + kt) / &
               B0phi(kt_low(kf) + kt, ei) * (pres0(kf-1) - pres0(kf))) + &
               j0phi(kt_low(kf) + kt, eo) * (Bnphi(kt_low(kf) + kt) / &
               B0phi(kt_low(kf) + kt, ei) + Bnflux(kt_low(kf) + kt, eo) / r / Deltapsi))
       end do
       call assemble_sparse(kt_max(kf), d(:kt_max(kf)), du(:kt_max(kf)), nz, &
            irow(:(2*kt_max(kf))), icol(:(2*kt_max(kf))), aval(:(2*kt_max(kf))))
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(kt_max(kf), kt_max(kf), nz, irow(:nz), icol(:nz), aval(:nz), &
            x(:kt_max(kf)))
       call sparse_matmul(kt_max(kf), kt_max(kf), irow(:nz), icol(:nz), aval(:nz), &
            x(:kt_max(kf)), resid)
       resid = resid - inhom(:kt_max(kf))
       rel_err(:kt_max(kf)) = abs(resid(:kt_max(kf))) / abs(inhom(:kt_max(kf)))
       max_rel_err = max(max_rel_err, maxval(rel_err(:kt_max(kf))))
       avg_rel_err = avg_rel_err + sum(rel_err(:kt_max(kf)))
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          area = elem%det_3 * 0.5d0
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          jnflux(kt_low(kf) + kt, ei) = -x(kt)
          jnflux(kt_low(kf) + kt, eo) = x(mod(kt, kt_max(kf)) + 1)
          jnphi(kt_low(kf) + kt) = sum(jnflux(kt_low(kf) + kt, :)) * imun / n / area
       end do
    end do
    avg_rel_err = avg_rel_err / sum(kt_max(1:nflux))
    if (log_debug) write (logfile, *) 'compute_triangle_flux: diagonalization' // &
         ' max_rel_err = ', max_rel_err, ' avg_rel_err = ', avg_rel_err
    if (allocated(resid)) deallocate(resid)

    call check_redundant_edges(jnflux, -1d0, 'jnflux')

    call write_triangle_flux(jnflux, jnphi, currn_file, .false.)
    call write_triangle_flux(jnflux, jnphi, decorate_filename(currn_file, 'plot_'), &
         .true.)
  end subroutine compute_currn

  subroutine avg_flux_on_quad(pol_flux, tor_comp)
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)

    integer :: kf, kt
    complex(dp) :: tor_flux_avg, tor_flux_diff
    type(triangle) :: elem1, elem2
    integer, dimension(2) :: li, lo, lf
    integer :: ei1, eo1, ef1, ei2, eo2, ef2
    logical :: orient

    do kf = 2, nflux
       do kt = 1, kt_max(kf), 2
          elem1 = mesh_element(kt_low(kf) + kt)
          elem2 = mesh_element(kt_low(kf) + kt + 1)
          call get_labeled_edges(elem1, li, lo, lf, ei1, eo1, ef1, orient)
          call get_labeled_edges(elem2, li, lo, lf, ei2, eo2, ef2, orient)
          tor_flux_avg = 0.25d0 * (tor_comp(kt_low(kf) + kt) * elem1%det_3 + &
               tor_comp(kt_low(kf) + kt + 1) * elem2%det_3)
          tor_flux_diff = 0.25d0 * (-tor_comp(kt_low(kf) + kt) * elem1%det_3 + &
               tor_comp(kt_low(kf) + kt + 1) * elem2%det_3)
          tor_comp(kt_low(kf) + kt) = tor_flux_avg / elem1%det_3 * 2d0
          pol_flux(kt_low(kf) + kt, eo1) = pol_flux(kt_low(kf) + kt, eo1) - imun * n * &
                  tor_flux_diff
          tor_comp(kt_low(kf) + kt + 1) = tor_flux_avg / elem2%det_3 * 2d0
          pol_flux(kt_low(kf) + kt + 1, ei2) = pol_flux(kt_low(kf) + kt + 1, ei2) + &
               imun * n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  subroutine compute_Bn_nonres
    integer :: kf, kt
    type(triangle) :: elem
    type(knot) :: base, tip
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    integer :: common_tri(2)
    real(dp) :: lr, lz, Deltapsi, perps(2)
    complex(dp) :: Bnpsi
    real(dp) :: r

    do kf = 1, nflux ! loop through flux surfaces
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
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord
          common_tri = (/ kt_low(kf) + kt, elem%neighbour(ef) /)
          perps = mesh_element(common_tri(:))%det_3 / hypot(lr, lz) * 0.5d0
          Bnpsi = R0 * B0phi(kt_low(kf) + kt, ef) / r
          Bnflux(kt_low(kf) + kt, ef) = Bnpsi * r * hypot(lr, lz) * sum(perps) / Deltapsi
          Bnphi(kt_low(kf) + kt) = imun / n * Bnflux(kt_low(kf) + kt, ef) &
               / elem%det_3 * 2d0
       end do
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          Bnflux(kt_low(kf) + kt, ei) = (0d0, 0d0)
          Bnflux(kt_low(kf) + kt, eo) = (0d0, 0d0)
       end do
    end do
    if (quad_avg) call avg_flux_on_quad(Bnflux, Bnphi)
    call check_div_free(Bnflux, Bnphi, 1d-09, 'non-resonant B_n')
    call check_redundant_edges(Bnflux, -1d0, 'non-resonant B_n')
  end subroutine compute_Bn_nonres

  subroutine write_triangle_flux(pol_flux, tor_comp, outfile, interpolate)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile
    logical, intent(in) :: interpolate

    integer :: ktri
    type(triangle) :: elem
    type(knot) :: tri(3)
    complex(dp) :: pol_comp_r, pol_comp_z
    real(dp) :: r, z

    open(1, file = outfile, recl = longlines)
    if (interpolate) then
       do ktri = 1, kt_low(nflux+1)
          elem = mesh_element(ktri)
          tri = mesh_point(elem%i_knot(:))
          call ring_centered_avg_coord(elem, r, z)
          pol_comp_r = 1d0 / elem%det_3 * ( &
               pol_flux(ktri, 1) * (r - tri(3)%rcoord) + &
               pol_flux(ktri, 2) * (r - tri(1)%rcoord) + &
               pol_flux(ktri, 3) * (r - tri(2)%rcoord))
          pol_comp_z = 1d0 / elem%det_3 * ( &
               pol_flux(ktri, 1) * (z - tri(3)%zcoord) + &
               pol_flux(ktri, 2) * (z - tri(1)%zcoord) + &
               pol_flux(ktri, 3) * (z - tri(2)%zcoord))
          write (1, *) r, z, real(pol_comp_r), aimag(pol_comp_r), real(pol_comp_z), &
               aimag(pol_comp_z), real(tor_comp(ktri)), aimag(tor_comp(ktri))
       end do
       r = mesh_point(1)%rcoord
       z = mesh_point(1)%zcoord
       do ktri = kt_low(nflux+1) + 1, ntri
          write (1, *) r, z, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
       end do
    else
       do ktri = 1, kt_low(nflux+1)
          write(1, *) &
               real(pol_flux(ktri, 1)), aimag(pol_flux(ktri, 1)), &
               real(pol_flux(ktri, 2)), aimag(pol_flux(ktri, 2)), &
               real(pol_flux(ktri, 3)), aimag(pol_flux(ktri, 3)), &
               real(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0), &
               aimag(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0)
       end do
       do ktri = kt_low(nflux+1) + 1, ntri
          write (1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
       end do
    end if
    close(1)
  end subroutine write_triangle_flux

  subroutine write_nodal_dof(nodal_dof, outfile)
    complex(dp), intent(in) :: nodal_dof(:)
    character(len = *), intent(in) :: outfile

    integer :: kpoint

    open(1, file = outfile, recl = longlines)
    do kpoint = 1, kp_low(nflux+1)
       write (1, *) real(nodal_dof(kpoint)), aimag(nodal_dof(kpoint))
    end do
    do kpoint = kp_low(nflux+1) + 1, npoint ! write zeroes in remaining points until end
       write (1, *) 0d0, 0d0
    end do
    close(1)
  end subroutine write_nodal_dof
end module magdif
