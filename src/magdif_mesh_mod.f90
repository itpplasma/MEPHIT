module magdif_mesh

  use iso_fortran_env, only: dp => real64
  use magdif_util, only: g_eqdsk, flux_func

  implicit none

  private

  public :: equil, fluxvar, flux_func_cache, fs, fs_half, mesh_t, mesh, &
       B0r, B0phi, B0z, B0r_Omega, B0phi_Omega, B0z_Omega, B0flux, j0phi, &
       flux_func_cache_init, flux_func_cache_check, flux_func_cache_destructor, generate_mesh, &
       refine_eqd_partition, refine_resonant_surfaces, write_kilca_convexfile, &
       create_mesh_points, init_indices, add_node_owner, common_triangles, &
       connect_mesh_points, get_labeled_edges, cache_mesh_data, write_mesh_data, &
       magdif_mesh_destructor, init_flux_variables, compute_pres_prof, &
       compute_safety_factor, check_safety_factor, cache_equilibrium_field, &
       compute_j0phi, check_curr0, write_fluxvar, point_location, point_in_triangle

  type(g_eqdsk) :: equil

    !> Structure containing flux functions evaluated at a specific flux surface, indicated
  !> by a common array index. For details see flux_func_cache_init().
  type, public :: flux_func_cache
     private
     integer :: nflux

     !> Magnetic flux surface label \f$ \psi \f$ in maxwell.
     !>
     !> \f$ \psi \f$ is the disc poloidal flux divided by \f$ 2 \pi \f$. Its sign is
     !> positive and its magnitude is growing in the radially inward direction.
     real(dp), dimension(:), allocatable, public :: psi

     !> Minor radius \f$ r \f$ in centimeter.
     real(dp), dimension(:), allocatable, public :: rad

     real(dp), dimension(:), allocatable, public :: F

     !> Equilibrium pressure \f$ p_{0} \f$ in barye.
     real(dp), dimension(:), allocatable, public :: p

     real(dp), dimension(:), allocatable, public :: FdF_dpsi

     !> Derivative of equilibrium pressure w.r.t. flux surface label,
     !> \f$ p_{0}'(\psi) \f$, in barye per maxwell.
     real(dp), dimension(:), allocatable, public :: dp_dpsi

     !> Safety factor \f$ q \f$ (dimensionless).
     real(dp), dimension(:), allocatable, public :: q
   contains
     procedure :: init => flux_func_cache_init
     final :: flux_func_cache_destructor
  end type flux_func_cache

  type(flux_func) :: fluxvar
  type(flux_func_cache) :: fs
  type(flux_func_cache) :: fs_half

  type :: mesh_t

     !> R coordinate of the O point in cm.
     real(dp) :: R_O

     !> Z coordinate of the O point in cm.
     real(dp) :: Z_O

     !> R coordinate of the X point in cm.
     real(dp) :: R_X

     !> Z coordinate of the X point in cm.
     real(dp) :: Z_X

     !> Minimal R value on computational grid in cm.
     real(dp) :: R_min

     !> Maximal Z value on computational grid in cm.
     real(dp) :: R_max

     !> Minimal R value on computational grid in cm.
     real(dp) :: Z_min

     !> Maximal Z value on computational grid in cm.
     real(dp) :: Z_max

     !> Number of flux surfaces. May differ from #magdif_conf::magdif_config::nflux due
     !> to refinement of flux surfaces.
     integer :: nflux

     !> Number of triangle edges. TODO: Add npoint and ntri along with mesh_element_rmp.
     integer :: nedge

     !> Toroidal mode number. May differ from #magdif_conf::magdif_config#n due to large
     !> aspect ratio scaling.
     integer :: n

     !> Minimal poloidal mode number in resonance.
     integer :: m_res_min

     !> Maximal poloidal mode number in resonance
     integer :: m_res_max

     !> Number of unrefined flux surfaces to be replaced by refined ones.
     integer, allocatable :: deletions(:)

     !> Number of refined flux surfaces.
     integer, allocatable :: additions(:)

     !> Relative size of most refined flux surface.
     real(dp), allocatable :: refinement(:)

     !> Poloidal mode number \f$ m \f$ (dimensionless) in resonance at given flux surface.
     !>
     !> Indexing is the same as for #q, on which the values depend. If no resonances are
     !> expected at a given index, #m_res is 0.
     integer, allocatable :: m_res(:)

     !> Indices of flux surfaces where resonance corresponding to a poloidal mode (given as
     !> array index) occurs.
     integer, allocatable :: res_ind(:)

     !> Number of knots on the flux surface given by the array index.
     !>
     !> The array index ranges from 1 for the innermost flux surface to
     !> #magdif_config::nflux for the last closed flux surface.
     integer, allocatable :: kp_max(:)

     !> Number of triangles inside the flux surface given by the array index.
     !>
     !> The array index ranges from 1 for the innermost flux surface to
     !> #magdif_config::nflux for the last closed flux surface.
     integer, allocatable :: kt_max(:)

     !> Global index of the last knot of the previous flux surface given by the array index.
     !>
     !> The global index of knots in #mesh_mod::mesh_point on the flux surface kf runs from
     !> #kp_low (kf)+1 to #kp_low (kf)+#kp_max (kf), so #kp_low is determined by cumulatively
     !> adding consecutive values of #kp_max. The array index ranges from 1, giving the
     !> global index of the knot on the magnetic axis (which has to be 1), to
     !> #magdif_config::nflux+1, effectively giving the last knot on the last closed
     !> flux surface.
     integer, allocatable :: kp_low(:)

     !> Global index of the last triangle of the previous flux surface given by the array
     !> index.
     !>
     !> The global index of triangles in #mesh_mod::mesh_element inside the flux surface kf
     !> runs from #kt_low (kf)+1 to #kt_low (kf)+#kt_max (kf), so #kt_low is determined by
     !> cumulatively adding consecutive values of #kt_max. The array index ranges from 1,
     !> giving the global index of the non-existent triangle on the magnetic axis (which is
     !> therefore 0), to #magdif_config::nflux+1, giving the last triangle inside the last
     !> closed flux surface.
     integer, allocatable :: kt_low(:)

     integer, allocatable :: edge_map2global(:, :)
     integer, allocatable :: edge_map2ktri(:, :)
     integer, allocatable :: edge_map2ke(:, :)

   contains
     final :: magdif_mesh_destructor
  end type mesh_t

  type(mesh_t) :: mesh

  !> \f$ R \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0r(:,:)

  !> \f$ \phi \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0phi(:,:)

  !> \f$ Z \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0z(:,:)

  !> \f$ \phi \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0r_Omega(:)

  !> \f$ Z \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0phi_Omega(:)

  !> \f$ R \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0z_Omega(:)

  real(dp), allocatable :: B0flux(:,:)

  !> Physical toroidal component of equilibrium current \f$ j_{0 (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: j0phi(:,:)

contains

  !> Set up arrays of cached values of flux functions.
  !>
  !> @nflux number of flux surfaces
  !> @half_step values are taken at flux surfaces (false) or between flux surfaces (true)
  !>
  !> For full-grid quantities, values are taken on flux surfaces with indices running
  !> from 0 to \p nflux, i.e. from the magnetic axis to the separatrix. An exception is
  !> made for \psi, where the index runs up to \p nflux +1. This value is extrapolated for
  !> finite differences in magdif::compute_presn() and magdif::compute_bn_nonres().
  !> For half-grid quantities, values are taken between two flux surfaces with indices
  !> running from 1 to \p nflux, i.e. from the triangle strip surrounding the magnetic
  !> axis to the triangle strip just inside the separatrix.
  subroutine flux_func_cache_init(this, nflux, half_step)
    class(flux_func_cache), intent(inout) :: this
    integer, intent(in) :: nflux
    logical, intent(in) :: half_step

    call flux_func_cache_destructor(this)
    if (half_step) then
       allocate(this%psi(nflux))
       allocate(this%rad(nflux))
       allocate(this%F(nflux))
       allocate(this%p(nflux))
       allocate(this%FdF_dpsi(nflux))
       allocate(this%dp_dpsi(nflux))
       allocate(this%q(nflux))
    else
       allocate(this%psi(0:nflux))
       allocate(this%rad(0:nflux))
       allocate(this%F(0:nflux))
       allocate(this%p(0:nflux))
       allocate(this%FdF_dpsi(0:nflux))
       allocate(this%dp_dpsi(0:nflux))
       allocate(this%q(0:nflux))
    end if
    this%psi = 0d0
    this%rad = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
  end subroutine flux_func_cache_init

  subroutine flux_func_cache_destructor(this)
    type(flux_func_cache), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%F)) deallocate(this%F)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%FdF_dpsi)) deallocate(this%FdF_dpsi)
    if (allocated(this%dp_dpsi)) deallocate(this%dp_dpsi)
    if (allocated(this%q)) deallocate(this%q)
  end subroutine flux_func_cache_destructor

  subroutine generate_mesh(unprocessed_geqdsk)
    use magdif_conf, only: conf, log
    use magdif_util, only: get_equil_filenames, initialize_globals

    character(len = *), intent(in) :: unprocessed_geqdsk
    character(len = 1024) :: gfile, convexfile

    call get_equil_filenames(gfile, convexfile)
    log%msg = 'attempting to read unprocessed G EQDSK file ' // trim(unprocessed_geqdsk)
    if (log%info) call log%write
    call equil%read(trim(unprocessed_geqdsk))
    call equil%classify
    call equil%standardise
    if (conf%kilca_scale_factor /= 0) then
       call equil%scale(conf%kilca_scale_factor)
    end if
    log%msg = 'attempting to write processed G EQDSK file ' // trim(gfile)
    if (log%info) call log%write
    call equil%write(trim(gfile))
    call initialize_globals(equil%rmaxis, equil%zmaxis)

    if (conf%kilca_scale_factor /= 0) then
       mesh%n = conf%n * conf%kilca_scale_factor
    else
       mesh%n = conf%n
    end if
    call create_mesh_points(convexfile)
    call connect_mesh_points
    call cache_mesh_data
    call write_mesh_data
  end subroutine generate_mesh

  subroutine refine_eqd_partition(nref, deletions, additions, refinement, res, refined, &
       ref_ind)
    use magdif_conf, only: conf, log
    use magdif_util, only: linspace
    integer, intent(in) :: nref
    integer, dimension(:), intent(in) :: deletions, additions
    real(dp), dimension(:), intent(in) :: refinement, res
    real(dp), dimension(:), allocatable, intent(out) :: refined
    integer, dimension(:), intent(out) :: ref_ind
    integer :: kref, k
    integer, dimension(nref) :: coarse_lo, coarse_hi, fine_lo, fine_hi
    real(dp) :: coarse_sep
    real(dp), dimension(nref) :: fine_sep, factor
    real(dp), dimension(maxval(additions) + 1, nref) :: geom_ser

    if (allocated(refined)) deallocate(refined)
    if (nref < 1) then
       mesh%nflux = conf%nflux_unref
       allocate(refined(0:mesh%nflux))
       refined = linspace(0d0, 1d0, mesh%nflux + 1, 0, 0)
       return
    end if
    if (nref /= size(deletions)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(deletions)', nref, &
            size(deletions))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(additions)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(additions)', nref, &
            size(additions))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(refinement)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(refinement)', nref, &
            size(refinement))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(res)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(res)', nref, size(res))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(ref_ind)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(ref_ind)', nref, &
            size(ref_ind))
       if (log%err) call log%write
       error stop
    end if
    mesh%nflux = conf%nflux_unref + 2 * sum(additions - deletions)
    allocate(refined(0:mesh%nflux))
    refined = 0d0
    coarse_lo = floor(res * (conf%nflux_unref + 1)) - deletions
    coarse_hi = ceiling(res * (conf%nflux_unref + 1)) + deletions
    ! compute upper and lower array indices of refined regions
    fine_lo = coarse_lo + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 0, nref - 1)]
    fine_hi = coarse_hi + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 1, nref)]
    ref_ind = fine_hi - additions
    ! compute separations between flux surfaces using geometric series
    coarse_sep = 1d0 / dble(conf%nflux_unref)
    fine_sep = coarse_sep * refinement
    factor = (coarse_sep / fine_sep) ** (1d0 / dble(additions + 1))
    geom_ser = reshape([(((factor(kref) ** k - factor(kref)) / (factor(kref) - 1d0), &
         k = 1, maxval(additions) + 1), kref = 1, nref)], [maxval(additions) + 1, nref])
    ! compute refined regions around resonant flux surfaces
    do kref = 1, nref
       refined(fine_lo(kref):fine_lo(kref)+additions(kref)) = res(kref) &
            - (geom_ser(additions(kref)+1:1:-1, kref) + 0.5d0) * fine_sep(kref)
       refined(fine_hi(kref)-additions(kref):fine_hi(kref)) = res(kref) &
            + (geom_ser(1:additions(kref)+1, kref) + 0.5d0) * fine_sep(kref)
    end do
    ! compute equidistant positions between refined regions
    refined(0:fine_lo(1)-1) = linspace(0d0, refined(fine_lo(1)), fine_lo(1), 0, 1)
    do kref = 2, nref
       refined(fine_hi(kref-1)+1:fine_lo(kref)-1) = linspace(refined(fine_hi(kref-1)), &
            refined(fine_lo(kref)), fine_lo(kref) - fine_hi(kref-1) - 1, 1, 1)
    end do
    refined(fine_hi(nref)+1:mesh%nflux) = linspace(refined(fine_hi(nref)), 1d0, &
         mesh%nflux - fine_hi(nref), 1, 0)
  end subroutine refine_eqd_partition

  subroutine refine_resonant_surfaces(psi_sample, q_sample, psi2rho_norm, rho_norm_ref)
    use magdif_conf, only: conf_arr, log
    use magdif_util, only: flux_func
    use netlib_mod, only: zeroin
    real(dp), dimension(:), intent(in) :: psi_sample
    real(dp), dimension(:), intent(in) :: q_sample
    interface
       function psi2rho_norm(psi)
         import :: dp
         real(dp), intent(in) :: psi
         real(dp) :: psi2rho_norm
       end function psi2rho_norm
    end interface
    real(dp), dimension(:), allocatable, intent(out) :: rho_norm_ref
    real(dp) :: psi_min, psi_max
    integer :: m, kref, m_lo, m_hi
    type(flux_func) :: psi_eqd
    real(dp), dimension(:), allocatable :: psi_res, rho_res
    integer, dimension(:), allocatable :: ref_ind
    logical, dimension(:), allocatable :: mask

    if (size(psi_sample) /= size(q_sample)) then
       call log%msg_arg_size('refine_resonant_surfaces', 'size(psi_sample)', &
            'size(q_sample)', size(psi_sample), size(q_sample))
       if (log%err) call log%write
       error stop
    end if
    psi_min = minval(psi_sample)
    psi_max = maxval(psi_sample)
    call psi_eqd%init(4, psi_sample)
    mesh%m_res_min = ceiling(minval(abs(q_sample)) * dble(mesh%n))
    mesh%m_res_max = floor(maxval(abs(q_sample)) * dble(mesh%n))
    if (conf_arr%m_min < mesh%m_res_min) then
       write (log%msg, '("Ignoring configuration values for ", i0, " <= m < ", i0, ".")') &
            conf_arr%m_min, mesh%m_res_min
       if (log%warn) call log%write
       m_lo = mesh%m_res_min
    else
       m_lo = conf_arr%m_min
    end if
    if (conf_arr%m_max > mesh%m_res_max) then
       write (log%msg, '("Ignoring configuration values for ", i0, " < m <= ", i0, ".")') &
            mesh%m_res_max, conf_arr%m_max
       if (log%warn) call log%write
       m_hi = mesh%m_res_max
    else
       m_hi = conf_arr%m_max
    end if
    if (allocated(mesh%refinement)) deallocate(mesh%refinement)
    if (allocated(mesh%deletions)) deallocate(mesh%deletions)
    if (allocated(mesh%additions)) deallocate(mesh%additions)
    allocate(mesh%refinement(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%deletions(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%additions(mesh%m_res_min:mesh%m_res_max))
    mesh%refinement = 0d0
    mesh%deletions = 0
    mesh%additions = 0
    mesh%refinement(m_lo:m_hi) = conf_arr%refinement(m_lo:m_hi)
    mesh%deletions(m_lo:m_hi) = conf_arr%deletions(m_lo:m_hi)
    mesh%additions(m_lo:m_hi) = conf_arr%additions(m_lo:m_hi)
    allocate(psi_res(mesh%m_res_min:mesh%m_res_max))
    allocate(rho_res(mesh%m_res_min:mesh%m_res_max))
    log%msg = 'resonance positions:'
    if (log%debug) call log%write
    do m = mesh%m_res_min, mesh%m_res_max
       psi_res(m) = zeroin(psi_min, psi_max, q_interp_resonant, 1d-9)
       rho_res(m) = psi2rho_norm(psi_res(m))
       write (log%msg, '("m = ", i2, ", psi_m = ", es24.16e3, ", rho_m = ", f19.16)') &
            m, psi_res(m), rho_res(m)
       if (log%debug) call log%write
    end do
    allocate(mask(mesh%m_res_min:mesh%m_res_max))
    mask = 0d0 < mesh%refinement .and. mesh%refinement < 1d0
    allocate(ref_ind(count(mask)))
    call refine_eqd_partition(count(mask), pack(mesh%deletions, mask), &
         pack(mesh%additions, mask), pack(mesh%refinement, mask), &
         pack(rho_res, mask), rho_norm_ref, ref_ind)

    log%msg = 'refinement positions:'
    if (log%debug) call log%write
    kref = 0
    do m = mesh%m_res_min, mesh%m_res_max
       if (.not. mask(m)) cycle
       kref = kref + 1
       write (log%msg, '("m = ", i0, ", kf = ", i0, ' // &
            '", rho: ", f19.16, 2(" < ", f19.16))') m, ref_ind(kref), &
            rho_norm_ref(ref_ind(kref) - 1), rho_res(m), rho_norm_ref(ref_ind(kref))
       if (log%debug) call log%write
    end do

    if (allocated(ref_ind)) deallocate(ref_ind)
    if (allocated(mask)) deallocate(mask)
    if (allocated(rho_res)) deallocate(rho_res)
    if (allocated(psi_res)) deallocate(psi_res)

  contains
    function q_interp_resonant(psi)
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant
      q_interp_resonant = psi_eqd%interp(abs(q_sample), psi) - dble(m) / dble(mesh%n)
    end function q_interp_resonant
  end subroutine refine_resonant_surfaces

  subroutine write_kilca_convexfile(rho_max, convexfile)
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90

    integer, parameter :: nrz = 96  ! at most 100 values are read in by field_divB0.f90
    real(dp), intent(in) :: rho_max
    character(len = *), intent(in) :: convexfile
    real(dp) :: theta
    integer :: fid, kp

    open(newunit = fid, file = convexfile, status = 'replace')
    do kp = 1, nrz
       theta = dble(kp) / dble(nrz) * 2d0 * pi
       write (fid, '(2(1x, es24.16e3))') equil%rmaxis + rho_max * cos(theta), &
            equil%zmaxis + rho_max * sin(theta)
    end do
    close(fid)
  end subroutine write_kilca_convexfile

  subroutine create_mesh_points(convexfile)
    use constants, only: pi
    use mesh_mod, only: npoint, mesh_point, ntri, mesh_element, mesh_element_rmp
    use magdif_conf, only: conf, log
    use magdif_util, only: interp_psi_pol, flux_func
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
         raxis, zaxis, load_magdata_in_symfluxcoord
    use field_line_integration_mod, only: circ_mesh_scale, o_point, x_point, &
         theta0_at_xpoint, theta_axis
    use points_2d, only: s_min, create_points_2d

    character(len = *), intent(in) :: convexfile
    integer :: kf, kpoi, fid
    integer, dimension(:), allocatable :: n_theta
    real(dp), dimension(:), allocatable :: rho_norm_eqd, rho_norm_ref
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi
    type(flux_func) :: psi_interpolator
    real(dp) :: psi_axis, rho_max

    theta0_at_xpoint = .true.
    circ_mesh_scale = conf%kilca_scale_factor
    if (conf%kilca_scale_factor /= 0) then
       ! calculate maximal extent from magnetic axis
       rho_max = min(equil%rmaxis - equil%rleft, &
            equil%rleft + equil%rdim - equil%rmaxis, &
            equil%zdim * 0.5d0 + equil%zmid - equil%zmaxis, &
            equil%zdim * 0.5d0 - equil%zmid + equil%zmaxis)
       call write_kilca_convexfile(rho_max, convexfile)
       o_point = [equil%rmaxis, equil%zmaxis]
       x_point = o_point + [rho_max, 0d0]
    end if
    ! calculates points on a fine grid in the core region by integrating along field lines
    call preload_for_SYNCH
    ! loads points that are calculated in preload_for_SYNCH into module variables
    ! and spline interpolates them for use with magdata_in_symfluxcoord_ext
    call load_magdata_in_symfluxcoord
    mesh%R_O = raxis
    mesh%Z_O = zaxis
    mesh%R_X = X_point(1)
    mesh%Z_X = X_point(2)
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    call psi_interpolator%init(4, psisurf(1:) * psipol_max + psi_axis)
    ! interpolate between psi and rho
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd = rbeg / hypot(theta_axis(1), theta_axis(2))

    call refine_resonant_surfaces(psisurf(1:) * psipol_max + psi_axis, qsaf, &
         psi2rho_norm, rho_norm_ref)
    call fs%init(mesh%nflux, .false.)
    call fs_half%init(mesh%nflux, .true.)
    fs%rad = rho_norm_ref
    fs%psi = [(interp_psi_pol(raxis + fs%rad(kf) * theta_axis(1), &
         zaxis + fs%rad(kf) * theta_axis(2)), kf = 0, mesh%nflux)]
    fs_half%rad = 0.5d0 * (fs%rad(0:mesh%nflux-1) + fs%rad(1:mesh%nflux))
    fs_half%psi = [(interp_psi_pol(raxis + fs_half%rad(kf) * theta_axis(1), &
         zaxis + fs_half%rad(kf) * theta_axis(2)), kf = 1, mesh%nflux)]
    fs%rad = fs%rad * hypot(theta_axis(1), theta_axis(2))
    fs_half%rad = fs_half%rad * hypot(theta_axis(1), theta_axis(2))
    call flux_func_cache_check
    ! dump presumedly optimal values for poloidal resolution
    open(newunit = fid, file = 'optpolres.dat', status = 'replace')
    do kf = 1, mesh%nflux - 1
       write (fid, '(2(1x, es24.16e3))') fs%rad(kf), 2d0 * pi * fs%rad(kf) / &
            (fs_half%rad(kf + 1) - fs_half%rad(kf))
    end do
    write (fid, '(2(1x, es24.16e3))') fs%rad(mesh%nflux), 2d0 * pi * fs%rad(mesh%nflux) / &
         (fs%rad(mesh%nflux) - fs%rad(mesh%nflux - 1))
    close(fid)

    call init_indices
    ntri = mesh%kt_low(mesh%nflux + 1)
    npoint = mesh%kp_low(mesh%nflux + 1)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))

    allocate(n_theta(mesh%nflux))
    allocate(points(3, npoint))
    allocate(points_s_theta_phi(3, npoint))
    n_theta = conf%nkpol
    s_min = 1d-16
    ! inp_label 2 to use poloidal psi with magdata_in_symfluxcoord_ext
    ! psi is not normalized by psipol_max, but shifted by -psi_axis
    call create_points_2d(2, n_theta, points, points_s_theta_phi, r_scaling_func = psi_ref)
    mesh_point(1)%rcoord = mesh%R_O
    mesh_point(1)%zcoord = mesh%Z_O
    do kpoi = 2, npoint
       mesh_point(kpoi)%rcoord = points(1, kpoi)
       mesh_point(kpoi)%zcoord = points(3, kpoi)
    end do
    mesh%R_min = minval(mesh_point%rcoord)
    mesh%R_max = maxval(mesh_point%rcoord)
    mesh%Z_min = minval(mesh_point%zcoord)
    mesh%Z_max = maxval(mesh_point%zcoord)

    if (allocated(n_theta)) deallocate(n_theta)
    if (allocated(points)) deallocate(points)
    if (allocated(points_s_theta_phi)) deallocate(points_s_theta_phi)
    if (allocated(rho_norm_ref)) deallocate(rho_norm_ref)
    if (allocated(rho_norm_eqd)) deallocate(rho_norm_eqd)

  contains
    function psi2rho_norm(psi) result(rho_norm)
      real(dp), intent(in) :: psi
      real(dp) :: rho_norm
      rho_norm = psi_interpolator%interp(rho_norm_eqd, psi)
    end function psi2rho_norm
    function psi_ref(psi_eqd)
      real(dp), dimension(:), intent(in) :: psi_eqd
      real(dp), dimension(size(psi_eqd)) :: psi_ref
      integer :: kf
      if (mesh%nflux /= size(psi_eqd)) then
         call log%msg_arg_size('psi_ref', 'nflux', 'size(psi_eqd)', &
              mesh%nflux, size(psi_eqd))
         if (log%err) call log%write
         error stop
      end if
      psi_ref = [(interp_psi_pol(raxis + rho_norm_ref(kf) * theta_axis(1), &
           zaxis + rho_norm_ref(kf) * theta_axis(2)) - psi_axis, kf = 1, mesh%nflux)]
    end function psi_ref
  end subroutine create_mesh_points

  !> Allocates and initializes #kp_low, #kp_max, #kt_low and #kt_max based on the values
  !> of #magdif_config::nflux and #magdif_config::nkpol. Deallocation is done in
  !> magdif_cleanup().
  subroutine init_indices
    use magdif_conf, only: conf, log
    integer :: kf

    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux + 1))
    allocate(mesh%kt_low(mesh%nflux + 1))

    mesh%kp_max = conf%nkpol
    mesh%kt_max = 2 * conf%nkpol
    mesh%kt_max(1) = conf%nkpol

    mesh%kp_low(1) = 1
    do kf = 2, mesh%nflux + 1
       mesh%kp_low(kf) = mesh%kp_low(kf-1) + mesh%kp_max(kf-1)
    end do
    mesh%kt_low(1) = 0
    do kf = 2, mesh%nflux + 1
       mesh%kt_low(kf) = mesh%kt_low(kf-1) + mesh%kt_max(kf-1)
    end do

    write (log%msg, '("Number of points up to LCFS: ", i0)') mesh%kp_low(mesh%nflux + 1)
    if (log%info) call log%write
    write (log%msg, '("Number of triangles up to LCFS: ", i0)') mesh%kt_low(mesh%nflux + 1)
    if (log%info) call log%write
  end subroutine init_indices

  elemental subroutine calculate_det_3(elem)
    use mesh_mod, only: triangle, knot, mesh_point
    type(triangle), intent(inout) :: elem
    real(dp) :: e1_r, e1_z, e2_r, e2_z
    e1_r = mesh_point(elem%i_knot(1))%rcoord - mesh_point(elem%i_knot(3))%rcoord
    e1_z = mesh_point(elem%i_knot(1))%zcoord - mesh_point(elem%i_knot(3))%zcoord
    e2_r = mesh_point(elem%i_knot(2))%rcoord - mesh_point(elem%i_knot(3))%rcoord
    e2_z = mesh_point(elem%i_knot(2))%zcoord - mesh_point(elem%i_knot(3))%zcoord
    elem%det_3 = abs(e1_r * e2_z - e1_z * e2_r)
  end subroutine calculate_det_3

  subroutine add_node_owner(kpoint, ktri)
    use mesh_mod, only: knot, mesh_point, n_owners_max
    use magdif_conf, only: log
    integer, intent(in) :: kpoint, ktri
    if (mesh_point(kpoint)%n_owners < n_owners_max) then
       mesh_point(kpoint)%i_owner_tri(mesh_point(kpoint)%n_owners + 1) = ktri
       mesh_point(kpoint)%n_owners = mesh_point(kpoint)%n_owners + 1
    else
       write (log%msg, '("Maximal number of owning triangles exceeded at point ", i0)') &
            kpoint
       if (log%warn) call log%write
    end if
  end subroutine add_node_owner

  !> Returns the indices of the two triangles sharing an edge.
  !>
  !> @param knot1 index of first knot of the edge
  !> @param knot2 index of second knot of the edge
  !> @param common_tri indices of the triangles sharing the given edge
  !>
  !> The program is halted if the input data is invalid, i.e. if more than two triangles
  !> appear to share the edge.
  subroutine common_triangles(knot1, knot2, common_tri)
    use mesh_mod, only: mesh_point
    use magdif_conf, only: log
    integer, intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)
    integer :: k, l, kcom

    common_tri = [0, 0]
    kcom = 0
    do k = 1, mesh_point(knot1)%n_owners
       do l = 1, mesh_point(knot2)%n_owners
          if (mesh_point(knot1)%i_owner_tri(k) == mesh_point(knot2)%i_owner_tri(l)) then
             kcom = kcom + 1
             if (kcom > 2) then
                write (log%msg, '("More than two common triangles for knots ", ' // &
                     'i0, " and ", i0)') knot1, knot2
                if (log%err) call log%write
                error stop
             end if
             common_tri(kcom) = mesh_point(knot1)%i_owner_tri(k)
          end if
       end do
    end do
  end subroutine common_triangles

  subroutine connect_mesh_points
    use mesh_mod, only: mesh_element
    integer :: kf, kp, kt, ktri, ktri_adj, common_tri(2)

    ktri = mesh%kt_low(1)
    ! define trianlges on innermost flux surface
    kf = 1
    do kp = 1, mesh%kp_max(kf)
       ktri = ktri + 1
       mesh_element(ktri)%i_knot(1) = mesh%kp_low(kf) + kp
       mesh_element(ktri)%i_knot(2) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
       mesh_element(ktri)%i_knot(3) = mesh%kp_low(kf)
       mesh_element(ktri)%knot_h = 3
       call add_node_owner(mesh_element(ktri)%i_knot(1), ktri)
       call add_node_owner(mesh_element(ktri)%i_knot(2), ktri)
       call calculate_det_3(mesh_element(ktri))
    end do
    ! define triangles on outer flux surfaces
    do kf = 2, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          ktri = ktri + 1
          mesh_element(ktri)%i_knot(1) = mesh%kp_low(kf) + kp
          mesh_element(ktri)%i_knot(2) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
          mesh_element(ktri)%i_knot(3) = mesh%kp_low(kf - 1) + kp
          mesh_element(ktri)%knot_h = 3
          call add_node_owner(mesh_element(ktri)%i_knot(1), ktri)
          call add_node_owner(mesh_element(ktri)%i_knot(2), ktri)
          call add_node_owner(mesh_element(ktri)%i_knot(3), ktri)
          call calculate_det_3(mesh_element(ktri))
          ktri = ktri + 1
          mesh_element(ktri)%i_knot(1) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
          mesh_element(ktri)%i_knot(2) = mesh%kp_low(kf - 1) + mod(kp, mesh%kp_max(kf - 1)) + 1
          mesh_element(ktri)%i_knot(3) = mesh%kp_low(kf - 1) + kp
          mesh_element(ktri)%knot_h = 1
          call add_node_owner(mesh_element(ktri)%i_knot(1), ktri)
          call add_node_owner(mesh_element(ktri)%i_knot(2), ktri)
          call add_node_owner(mesh_element(ktri)%i_knot(3), ktri)
          call calculate_det_3(mesh_element(ktri))
       end do
    end do
    ! set neighbours for edges i and o on innermost flux surface
    kf = 1
    do kt = 1, mesh%kt_max(kf)
       ktri = mesh%kt_low(kf) + kt
       ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
       mesh_element(ktri)%neighbour(2) = ktri_adj
       mesh_element(ktri)%neighbour_edge(2) = 3
       mesh_element(ktri_adj)%neighbour(3) = ktri
       mesh_element(ktri_adj)%neighbour_edge(3) = 2
    end do
    ! set neighbours for edges i and o on outer flux surfaces
    do kf = 2, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
          if (mod(kt, 2) == 1) then
             mesh_element(ktri)%neighbour(2) = ktri_adj
             mesh_element(ktri)%neighbour_edge(2) = 3
             mesh_element(ktri_adj)%neighbour(3) = ktri
             mesh_element(ktri_adj)%neighbour_edge(3) = 2
          else
             mesh_element(ktri)%neighbour(1) = ktri_adj
             mesh_element(ktri)%neighbour_edge(1) = 3
             mesh_element(ktri_adj)%neighbour(3) = ktri
             mesh_element(ktri_adj)%neighbour_edge(3) = 1
          end if
       end do
    end do
    ! set neighbours for edges f
    do kf = 1, mesh%nflux - 1
       do kp = 1, mesh%kp_max(kf)
          call common_triangles(mesh%kp_low(kf) + kp, &
               mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1, common_tri)
          ktri = minval(common_tri)
          ktri_adj = maxval(common_tri)
          mesh_element(ktri)%neighbour(1) = ktri_adj
          mesh_element(ktri)%neighbour_edge(1) = 2
          mesh_element(ktri_adj)%neighbour(2) = ktri
          mesh_element(ktri_adj)%neighbour_edge(2) = 1
       end do
    end do
    ! set dummy values for edges on boundary
    kf = mesh%nflux
    do kt = 1, mesh%kt_max(kf)
       if (mod(kt, 2) == 1) then
          ktri = mesh%kt_low(kf) + kt
          ktri_adj = ktri + mesh%kt_max(kf) + 1
          mesh_element(ktri)%neighbour(1) = ktri_adj
          mesh_element(ktri)%neighbour_edge(1) = 2
       end if
    end do
  end subroutine connect_mesh_points

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
  !> It is assumed that the knots are globally numbered in ascending order starting from
  !> the magnetic axis and going counter-clockwise around each flux surface. Furthermore
  !> it is assumed that edge 1 goes from knot 1 to knot 2, edge 2 from knot 2 to knot 3
  !> and edge 3 from knot 3 to knot 1. It is not assumed that knots are orderd
  !> counter-clockwise locally.

  subroutine get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
    use mesh_mod, only: triangle
    use magdif_conf, only: log
    type(triangle), intent(in) :: elem
    integer, dimension(:), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    logical, intent(out) :: orient
    integer, dimension(3) :: i_knot_diff
    integer :: knot_i, knot_o, knot_f
    integer :: i1, i2
    logical :: closing_loop

    if (2 /= size(li)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(li)', 'actual size(li)', &
            2, size(li))
       if (log%err) call log%write
       error stop
    end if
    if (2 /= size(lo)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(lo)', 'actual size(lo)', &
            2, size(lo))
       if (log%err) call log%write
       error stop
    end if
    if (2 /= size(lf)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(lf)', 'actual size(lf)', &
            2, size(lf))
       if (log%err) call log%write
       error stop
    end if
    log%msg = 'cannot find correct label for triangle edges'

    ! initialize to suppress compiler warnings
    i1 = 0
    i2 = 0
    knot_f = elem%knot_h
    select case (knot_f)
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
       if (log%err) call log%write
       error stop
    end if
    ! last triangle in strip if indices not next to each other
    closing_loop = abs(elem%i_knot(i1) - elem%i_knot(i2)) /= 1
    i_knot_diff = elem%i_knot - elem%i_knot(knot_f)
    if (all(i_knot_diff >= 0)) then
       ! knot_f lies on inner surface
       orient = .true.
       if ((elem%i_knot(i1) < elem%i_knot(i2)) .neqv. closing_loop) then
          ! i1 is next after knot_f counter-clockwise
          knot_o = i1
          knot_i = i2
       else
          ! i2 is next after knot_f counter-clockwise
          knot_o = i2
          knot_i = i1
       end if
       ei = knot_f
       eo = knot_i
       ef = knot_o
       li = [elem%i_knot(knot_f), elem%i_knot(knot_o)]
       lo = [elem%i_knot(knot_i), elem%i_knot(knot_f)]
       lf = [elem%i_knot(knot_o), elem%i_knot(knot_i)]
    else if (all(i_knot_diff <= 0)) then
       ! knot_f lies on outer surface
       orient = .false.
       if ((elem%i_knot(i1) > elem%i_knot(i2)) .neqv. closing_loop) then
          ! i1 is next after knot_f counter-clockwise
          knot_i = i1
          knot_o = i2
       else
          ! i2 is next after knot_f counter-clockwise
          knot_i = i2
          knot_o = i1
       end if
       ei = knot_o
       eo = knot_f
       ef = knot_i
       li = [elem%i_knot(knot_o), elem%i_knot(knot_f)]
       lo = [elem%i_knot(knot_f), elem%i_knot(knot_i)]
       lf = [elem%i_knot(knot_i), elem%i_knot(knot_o)]
    else
       if (log%err) call log%write
       error stop
    end if
  end subroutine get_labeled_edges

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
  pure subroutine ring_centered_avg_coord(elem, r, z)
    use mesh_mod, only: triangle, knot, mesh_point
    type(triangle), intent(in) :: elem
    real(dp), intent(out) :: r, z
    type(knot), dimension(3) :: knots

    knots = mesh_point(elem%i_knot)
    r = (sum(knots%rcoord) + knots(elem%knot_h)%rcoord) * 0.25d0
    z = (sum(knots%zcoord) + knots(elem%knot_h)%zcoord) * 0.25d0
  end subroutine ring_centered_avg_coord

  subroutine cache_mesh_data
    use mesh_mod, only: triangle, triangle_rmp, mesh_element, mesh_element_rmp
    integer :: ktri, kedge, ke, ke_adj, ktri_adj
    type(triangle) :: elem
    type(triangle_rmp) :: tri

    mesh%nedge = (3 * mesh%kt_low(mesh%nflux + 1) + mesh%kp_max(mesh%nflux)) / 2
    allocate(mesh%edge_map2global(mesh%kt_low(mesh%nflux + 1), 3))
    allocate(mesh%edge_map2ktri(mesh%nedge, 2))
    allocate(mesh%edge_map2ke(mesh%nedge, 2))
    mesh%edge_map2global = 0
    mesh%edge_map2ktri = 0
    mesh%edge_map2ke = 0
    kedge = 1

    do ktri = 1, mesh%kt_low(mesh%nflux + 1)
       elem = mesh_element(ktri)
       tri%area = 0.5d0 * elem%det_3
       call get_labeled_edges(elem, tri%li, tri%lo, tri%lf, tri%ei, tri%eo, tri%ef, &
            tri%orient)
       call ring_centered_avg_coord(elem, tri%R_Omega, tri%Z_Omega)
       mesh_element_rmp(ktri) = tri

       do ke = 1, 3
          if (mesh%edge_map2global(ktri, ke) == 0) then
             ktri_adj = elem%neighbour(ke)
             ke_adj = elem%neighbour_edge(ke)
             if (ktri_adj > mesh%kt_low(mesh%nflux + 1)) then
                mesh%edge_map2global(ktri, ke) = kedge
                mesh%edge_map2ktri(kedge, :) = [ktri, -1]
                mesh%edge_map2ke(kedge, :) = [ke, -1]
                kedge = kedge + 1
             else
                mesh%edge_map2global(ktri, ke) = kedge
                mesh%edge_map2global(ktri_adj, ke_adj) = kedge
                mesh%edge_map2ktri(kedge, :) = [ktri, ktri_adj]
                mesh%edge_map2ke(kedge, :) = [ke, ke_adj]
                kedge = kedge + 1
             end if
          end if
       end do
    end do
  end subroutine cache_mesh_data

  function point_location(r, z) result(location)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf
    use magdif_util, only: interp_psi_pol, binsearch
    use mesh_mod, only: mesh_point
    real(dp), intent(in) :: r, z
    integer :: location

    integer :: kf, kq, k, k_max, ktri, candidates(6)
    real(dp) :: psi, pol_frac, pol_offset, pol_interval(0:conf%nkpol), thickness

    location = -1
    if (R < mesh%R_min .or. R > mesh%R_max .or. Z < mesh%Z_min .or. Z > mesh%Z_max) return
    location = -2
    psi = interp_psi_pol(r, z)
    if (equil%cocos%sgn_dpsi == +1) then
       if (psi > fs%psi(mesh%nflux)) return
    else
       if (psi < fs%psi(mesh%nflux)) return
    end if
    call binsearch(fs%psi, lbound(fs%psi, 1), psi, kf)
    location = -3

    pol_interval = 0d0
    pol_interval(0:mesh%kp_max(kf)-1) = 0.5d0 / pi * atan2( &
         mesh_point((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))%zcoord - mesh%Z_O, &
         mesh_point((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))%rcoord - mesh%R_O)
    pol_offset = pol_interval(0)
    pol_interval = pol_interval - pol_offset
    pol_interval(mesh%kp_max(kf)) = 1d0
    where (pol_interval < 0d0) pol_interval = pol_interval + 1d0
    where (pol_interval > 1d0) pol_interval = pol_interval - 1d0
    pol_frac = 0.5d0 / pi * atan2(Z - mesh%Z_O, R - mesh%R_O) - pol_offset
    if (pol_frac < 0d0) pol_frac = pol_frac + 1d0
    if (pol_frac > 1d0) pol_frac = pol_frac - 1d0
    call binsearch(pol_interval, lbound(pol_interval, 1), pol_frac, kq)

    ! Triangle edges do not lie exactly on flux surfaces, so we include the two adjacent
    ! triangle strips in the search. The candidates are ordered by decreasing likelihood
    ! of being the correct guess, i.e., the current loop, the outer loop and the inner
    ! loop, or filler if any of these is not applicable.
    if (kf == 1) then
       k_max = 3
       candidates = [mesh%kt_low(kf) + kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            -1, -1, -1]
    elseif (kf == 2) then
       k_max = 5
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            mesh%kt_low(kf - 1) + kq, -1]
    elseif (kf == mesh%nflux) then
       k_max = 4
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf - 1) + 2 * kq - 1, mesh%kt_low(kf - 1) + 2 * kq, &
            -1, -1]
    else
       k_max = 6
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            mesh%kt_low(kf - 1) + 2 * kq - 1, mesh%kt_low(kf - 1) + 2 * kq]
    end if

    thickness = fs%rad(mesh%nflux) * sqrt(epsilon(1d0)) * 8d0
    do k = 1, k_max
       ktri = candidates(k)
       if (point_in_triangle(ktri, R, Z, thickness)) then
          location = ktri
          exit
       end if
    end do
  end function point_location

  ! based on http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
  function point_in_triangle(ktri, R, Z, thickness) result(probably)
    use mesh_mod, only: triangle, mesh_element, mesh_point
    integer, intent(in) :: ktri
    real(dp), intent(in) :: R, Z
    real(dp), intent(in), optional :: thickness
    logical :: probably
    real(dp), dimension(1:4) :: node_R, node_Z, dist_R, dist_Z
    real(dp), dimension(1:3) :: edge_R, edge_Z, edge_2, dist_2, dotprod
    type(triangle) :: elem

    probably = .false.
    if (ktri <= 0) return
    elem = mesh_element(ktri)
    node_R = mesh_point([elem%i_knot, elem%i_knot(1)])%rcoord
    node_Z = mesh_point([elem%i_knot, elem%i_knot(1)])%zcoord
    dist_R = R - node_R
    dist_Z = Z - node_Z
    edge_R = node_R(2:4) - node_R(1:3)
    edge_Z = node_Z(2:4) - node_Z(1:3)
    edge_2 = edge_R ** 2 + edge_Z ** 2
    ! perp_R = edge_Z, perp_Z = -edge_R
    dotprod = edge_Z * dist_R(1:3) - edge_R * dist_Z(1:3)
    probably = all(dotprod <= 0d0)
    if (probably .or. .not. present(thickness)) return
    ! reuse dotprod as parameter of edge vector in linear equation
    dotprod = edge_R * dist_R(1:3) + edge_Z * dist_Z(1:3)
    where (dotprod < 0d0)
       dist_2 = dist_R(1:3) ** 2 + dist_Z(1:3) ** 2
    elsewhere (dotprod > edge_2)
       dist_2 = dist_R(2:4) ** 2 + dist_Z(2:4) ** 2
    elsewhere
       dist_2 = dist_R(1:3) ** 2 + dist_Z(1:3) ** 2 - dotprod ** 2 / edge_2
    end where
    probably = any(dist_2 < thickness ** 2)
  end function point_in_triangle

  subroutine write_mesh_data
    use mesh_mod, only: npoint, ntri, knot, triangle, mesh_point, mesh_element, &
         mesh_element_rmp
    use magdif_conf, only: longlines
    use hdf5_tools, only: HID_T, h5_create, h5_define_group, h5_close_group, h5_add, &
         h5_close

    integer(HID_T) :: h5id_magdif, h5id_mesh, h5id_cache, h5id_fs, h5id_fs_half
    integer :: tri_node(3, ntri), adj_tri(3, ntri), adj_edge(3, ntri), &
         tri_li(2, ntri), tri_lo(2, ntri), tri_lf(2, ntri), tri_orient(ntri)
    integer :: fid, kpoi, ktri, kp, ke
    type(triangle) :: elem

    open(newunit = fid, file = 'mesh_new.asc', recl = longlines, status = 'replace')
    do ktri = 1, ntri
       elem = mesh_element(ktri)
       write (fid, '(3(1x, i5), 1x, i1, 3(1x, i5), 3(1x, i3))') &
            (elem%i_knot(ke), ke = 1, 3), elem%knot_h, &
            (elem%neighbour(ke), ke = 1, 3), (elem%neighbour_edge(ke), ke = 1, 3)
    end do
    close(fid)

    call flux_func_cache_check
    ! intermediate until mesh_mod is refactored into magdif_mesh
    do ktri = 1, ntri
       tri_node(:, ktri) = mesh_element(ktri)%i_knot
       adj_tri(:, ktri) = mesh_element(ktri)%neighbour
       adj_edge(:, ktri) = mesh_element(ktri)%neighbour_edge
       tri_li(:, ktri) = mesh_element_rmp(ktri)%li
       tri_lo(:, ktri) = mesh_element_rmp(ktri)%lo
       tri_lf(:, ktri) = mesh_element_rmp(ktri)%lf
    end do
    where (mesh_element_rmp%orient)
       tri_orient = 1
    elsewhere
       tri_orient = 0
    end where
    call h5_create('magdif.h5', h5id_magdif)
    call h5_define_group(h5id_magdif, 'mesh', h5id_mesh)
    call h5_add(h5id_mesh, 'R_O', mesh%R_O, &
         comment = 'R coordinate of O point', unit = 'cm')
    call h5_add(h5id_mesh, 'Z_O', mesh%Z_O, &
         comment = 'Z coordinate of O point', unit = 'cm')
    call h5_add(h5id_mesh, 'R_min', mesh%R_min, &
         comment = 'minimal R coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_mesh, 'Z_min', mesh%Z_min, &
         comment = 'minimal Z coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_mesh, 'R_max', mesh%R_max, &
         comment = 'maximal R coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_mesh, 'Z_max', mesh%Z_max, &
         comment = 'maximal Z coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_mesh, 'nflux', mesh%nflux, comment = 'number of flux surfaces')
    call h5_add(h5id_mesh, 'npoint', npoint, comment = 'number of points')
    call h5_add(h5id_mesh, 'ntri', ntri, comment = 'number of triangles')
    call h5_add(h5id_mesh, 'nedge', mesh%nedge, comment = 'number of edges')
    call h5_add(h5id_mesh, 'n', mesh%n, comment = 'toroidal mode number')
    call h5_add(h5id_mesh, 'm_res_min', mesh%m_res_min, &
         comment = 'minimal absolute poloidal mode number')
    call h5_add(h5id_mesh, 'm_res_max', mesh%m_res_max, &
         comment = 'maximal absolute poloidal mode number')
    call h5_add(h5id_mesh, 'kp_max', mesh%kp_max, lbound(mesh%kp_max), ubound(mesh%kp_max), &
         comment = 'number of nodes on flux surface')
    call h5_add(h5id_mesh, 'kp_low', mesh%kp_low, lbound(mesh%kp_low), ubound(mesh%kp_low), &
         comment = 'index of last node on previous flux surface')
    call h5_add(h5id_mesh, 'kt_max', mesh%kt_max, lbound(mesh%kt_max), ubound(mesh%kt_max), &
         comment = 'number of triangles on flux surface')
    call h5_add(h5id_mesh, 'kt_low', mesh%kt_low, lbound(mesh%kt_low), ubound(mesh%kt_low), &
         comment = 'index of last triangle on previous flux surface')
    call h5_add(h5id_mesh, 'node_R', mesh_point%rcoord, [1], [npoint], &
         comment = 'R coordinates of mesh points', unit = 'cm')
    call h5_add(h5id_mesh, 'node_Z', mesh_point%zcoord, [1], [npoint], &
         comment = 'Z coordinates of mesh points', unit = 'cm')
    call h5_add(h5id_mesh, 'tri_node', tri_node, lbound(tri_node), ubound(tri_node), &
         comment = 'triangle node indices')
    call h5_add(h5id_mesh, 'tri_node_F', mesh_element%knot_h, [1], [ntri], &
         comment = 'local node index of node F')
    call h5_add(h5id_mesh, 'tri_li', tri_li, lbound(tri_li), ubound(tri_li), &
         comment = 'local node indices of edge i, counter-clockwise')
    call h5_add(h5id_mesh, 'tri_lo', tri_lo, lbound(tri_lo), ubound(tri_lo), &
         comment = 'local node indices of edge o, counter-clockwise')
    call h5_add(h5id_mesh, 'tri_lf', tri_lf, lbound(tri_lf), ubound(tri_lf), &
         comment = 'local node indices of edge f, counter-clockwise')
    call h5_add(h5id_mesh, 'tri_ei', mesh_element_rmp%ei, &
         lbound(mesh_element_rmp), ubound(mesh_element_rmp), &
         comment = 'local edge index of edge i')
    call h5_add(h5id_mesh, 'tri_eo', mesh_element_rmp%eo, &
         lbound(mesh_element_rmp), ubound(mesh_element_rmp), &
         comment = 'local edge index of edge o')
    call h5_add(h5id_mesh, 'tri_ef', mesh_element_rmp%ef, &
         lbound(mesh_element_rmp), ubound(mesh_element_rmp), &
         comment = 'local edge index of edge f')
    call h5_add(h5id_mesh, 'tri_orient', tri_orient, [1], [ntri], &
         comment = 'triangle orientation: true (1) if edge f lies on outer flux surface')
    call h5_add(h5id_mesh, 'adj_tri', adj_tri, lbound(adj_tri), ubound(adj_tri), &
         comment = 'adjacent triangle indices')
    call h5_add(h5id_mesh, 'adj_edge', adj_edge, lbound(adj_edge), ubound(adj_edge), &
         comment = 'adjacent triangle edge indices')
    call h5_add(h5id_mesh, 'edge_map2kedge', mesh%edge_map2global, &
         lbound(mesh%edge_map2global), ubound(mesh%edge_map2global), &
         comment = 'mapping of triangle & local edge index to global edge index')
    call h5_add(h5id_mesh, 'edge_map2ktri', mesh%edge_map2ktri, &
         lbound(mesh%edge_map2ktri), ubound(mesh%edge_map2ktri), &
         comment = 'mapping of global edge index to triangle index')
    call h5_add(h5id_mesh, 'edge_map2ke', mesh%edge_map2ke, &
         lbound(mesh%edge_map2ke), ubound(mesh%edge_map2ke), &
         comment = 'mapping of global edge index to local edge index')
    call h5_add(h5id_mesh, 'tri_centr_R', mesh_element_rmp%R_Omega, &
         lbound(mesh_element_rmp), ubound(mesh_element_rmp), &
         'R coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_mesh, 'tri_centr_Z', mesh_element_rmp%Z_Omega, &
         lbound(mesh_element_rmp), ubound(mesh_element_rmp), &
         'Z coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_mesh, 'det', mesh_element%det_3, [1], [ntri], &
         comment = 'determinant of triangle node coordinates')
    call h5_close_group(h5id_mesh)
    call h5_define_group(h5id_magdif, 'cache', h5id_cache)
    call h5_define_group(h5id_cache, 'fs', h5id_fs)
    call h5_add(h5id_fs, 'psi', fs%psi, lbound(fs%psi), ubound(fs%psi), &
         comment = 'poloidal flux on flux surfaces', unit = 'Mx')
    call h5_add(h5id_fs, 'rad', fs%rad, lbound(fs%rad), ubound(fs%rad), &
         comment = 'radial position on OX line on flux surfaces', unit = 'cm')
    call h5_close_group(h5id_fs)
    call h5_define_group(h5id_cache, 'fs_half', h5id_fs_half)
    call h5_add(h5id_fs_half, 'psi', fs_half%psi, lbound(fs_half%psi), ubound(fs_half%psi), &
         comment = 'poloidal flux between flux surfaces', unit = 'Mx')
    call h5_add(h5id_fs_half, 'rad', fs_half%rad, lbound(fs_half%rad), ubound(fs_half%rad), &
         comment = 'radial position on OX line between flux surfaces', unit = 'cm')
    call h5_close_group(h5id_fs_half)
    call h5_close_group(h5id_cache)
    call h5_close(h5id_magdif)

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace')
    write (fid, '(3(1x, i0))') npoint, ntri, mesh%kp_max(mesh%nflux) - 1
    do kpoi = 1, mesh%kp_low(mesh%nflux + 1)
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 0
    end do
    do kpoi = mesh%kp_low(mesh%nflux + 1) + 1, npoint
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 1
    end do
    do ktri = 1, ntri
       write (fid, '(4(1x, i0))') mesh_element(ktri)%i_knot(:), 0
    end do
    do kp = 1, mesh%kp_max(mesh%nflux) - 1
       write (fid, '(4(1x, i0))') mesh%kp_low(mesh%nflux) + kp, mesh%kp_low(mesh%nflux) + kp + 1, 1
    end do
    close(fid)
  end subroutine write_mesh_data

  subroutine magdif_mesh_destructor(this)
    type(mesh_t), intent(inout) :: this
    if (allocated(this%deletions)) deallocate(this%deletions)
    if (allocated(this%additions)) deallocate(this%additions)
    if (allocated(this%refinement)) deallocate(this%refinement)
    if (allocated(this%m_res)) deallocate(this%m_res)
    if (allocated(this%res_ind)) deallocate(this%res_ind)
    if (allocated(this%kp_max)) deallocate(this%kp_max)
    if (allocated(this%kt_max)) deallocate(this%kt_max)
    if (allocated(this%kp_low)) deallocate(this%kp_low)
    if (allocated(this%kt_low)) deallocate(this%kt_low)
  end subroutine magdif_mesh_destructor

  subroutine init_flux_variables
    integer :: kf

    ! initialize fluxvar with equidistant psi values
    call fluxvar%init(4, equil%psi_eqd)

    call compute_pres_prof
    call compute_safety_factor
    fs%F = [(fluxvar%interp(equil%fpol, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%F = [(fluxvar%interp(equil%fpol, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    fs_half%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs_half%psi(kf)), &
         kf = 1, mesh%nflux)]
    call write_fluxvar
  end subroutine init_flux_variables

  subroutine compute_pres_prof
    use constants, only: ev2erg  ! orbit_mod.f90
    use magdif_conf, only: conf, pres_prof_eps, pres_prof_par, pres_prof_geqdsk, log
    integer :: kf
    real(dp) :: ddens_dpsi, dtemp_dpsi, psi_int, psi_ext

    !> Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
    real(dp) :: dens(0:mesh%nflux)

    !> Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
    real(dp) :: temp(0:mesh%nflux)

    dens = 0d0
    temp = 0d0
    if (equil%cocos%sgn_dpsi == -1) then
       psi_ext = minval(equil%psirz)
       psi_int = maxval(equil%psirz)
    else
       psi_ext = maxval(equil%psirz)
       psi_int = minval(equil%psirz)
    end if
    select case (conf%pres_prof)
    case (pres_prof_eps)
       ddens_dpsi = conf%dens_max / psi_int
       dtemp_dpsi = conf%temp_max / psi_int
       dens = (fs%psi - psi_ext) / psi_int * conf%dens_max + conf%dens_min
       temp = (fs%psi - psi_ext) / psi_int * conf%temp_max + conf%temp_min
       write (log%msg, '("temp@axis: ", es24.16e3, ", dens@axis: ", es24.16e3)') &
            temp(0), dens(0)
       if (log%info) call log%write
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psi_ext) / psi_int * conf%dens_max + conf%dens_min
       temp(1:) = (fs_half%psi - psi_ext) / psi_int * conf%temp_max + conf%temp_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_par)
       ddens_dpsi = (conf%dens_max - conf%dens_min) / (psi_int - psi_ext)
       dtemp_dpsi = (conf%temp_max - conf%temp_min) / (psi_int - psi_ext)
       dens = (fs%psi - psi_ext) / (psi_int - psi_ext) * (conf%dens_max - conf%dens_min) &
            + conf%dens_min
       temp = (fs%psi - psi_ext) / (psi_int - psi_ext) * (conf%temp_max - conf%temp_min) &
            + conf%temp_min
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
            (conf%dens_max - conf%dens_min) + conf%dens_min
       temp(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
            (conf%temp_max - conf%temp_min) + conf%temp_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_geqdsk)
       fs%p = [(fluxvar%interp(equil%pres, fs%psi(kf)), kf = 0, mesh%nflux)]
       fs%dp_dpsi = [(fluxvar%interp(equil%pprime, fs%psi(kf)), kf = 0, mesh%nflux)]
       fs_half%p = [(fluxvar%interp(equil%pres, fs_half%psi(kf)), kf = 1, mesh%nflux)]
       fs_half%dp_dpsi = [(fluxvar%interp(equil%pprime, fs_half%psi(kf)), &
            kf = 1, mesh%nflux)]
    case default
       write (log%msg, '("unknown pressure profile selection", i0)') conf%pres_prof
       if (log%err) call log%write
       error stop
    end select
  end subroutine compute_pres_prof

  !> Allocates and computes the safety factor #q and #m_res.
  !>
  !> Also allocates #magdif_config::sheet_current_factor, to be read in via
  !> magdif_config::read_delayed_config() in magdif_init(). All deallocation is done in
  !> magdif_cleanup().
  subroutine compute_safety_factor
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, conf_arr, q_prof_flux, q_prof_rot, q_prof_geqdsk, &
         log, cmplx_fmt
    use magdif_util, only: flux_func
    use magdata_in_symfluxcoor_mod, only: psipol_max, psisurf, qsaf
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    integer :: kf, kt, ktri, m, kf_res
    type(triangle_rmp) :: tri
    type(flux_func) :: psi_interpolator
    real(dp), dimension(mesh%nflux) :: abs_err

    select case (conf%q_prof)
    case (q_prof_flux)
       fs_half%q = 0d0
       do kf = 1, mesh%nflux
          do kt = 1, mesh%kt_max(kf)
             ktri = mesh%kt_low(kf) + kt
             tri = mesh_element_rmp(ktri)
             fs_half%q(kf) = fs_half%q(kf) + B0phi_Omega(ktri) * tri%area
          end do
          fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
       end do
       call psi_interpolator%init(4, fs_half%psi)
       ! Lagrange polynomial extrapolation for values at separatrix and magnetic axis
       fs%q = [(psi_interpolator%interp(fs_half%q, fs%psi(kf)), kf = 0, mesh%nflux)]
    case (q_prof_rot)
       ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
       ! load_magdata_in_symfluxcoord_ext divides by psipol_max
       call psi_interpolator%init(4, psisurf(1:) * psipol_max + fs%psi(0))
       ! Lagrange polynomial extrapolation for value at magnetic axis
       fs%q = [(psi_interpolator%interp(qsaf, fs%psi(kf)), kf = 0, mesh%nflux)]
       fs_half%q = [(psi_interpolator%interp(qsaf, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    case (q_prof_geqdsk)
       fs%q = [(fluxvar%interp(equil%qpsi, fs%psi(kf)), kf = 0, mesh%nflux)]
       fs_half%q = [(fluxvar%interp(equil%qpsi, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    case default
       write (log%msg, '("unknown q profile selection: ", i0)') conf%q_prof
       if (log%err) call log%write
       error stop
    end select

    allocate(mesh%m_res(mesh%nflux))
    mesh%m_res = 0
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    mesh%res_ind = 0
    log%msg = 'resonance positions:'
    if (log%debug) call log%write
    do m = mesh%m_res_max, mesh%m_res_min, -1
       abs_err = [(abs(abs(fs_half%q(kf)) - dble(m) / dble(mesh%n)), kf = 1, mesh%nflux)]
       kf_res = minloc(abs_err, 1)
       mesh%res_ind(m) = kf_res
       mesh%m_res(kf_res) = m
       write (log%msg, '("m = ", i0, ", kf_res = ", i0, ' // &
            '", rho: ", f19.16, 2(" < ", f19.16))') m, kf_res, &
            [fs%rad(kf_res - 1), fs_half%rad(kf_res), fs%rad(kf_res)] / fs%rad(mesh%nflux)
       if (log%debug) call log%write
    end do

    if (log%info) then
       log%msg = 'absolute poloidal mode number, sheet current factor'
       call log%write
       do m = conf%m_min, conf%m_max
          write (log%msg, '(i2, 1x, ' // cmplx_fmt // ')') m, &
               conf_arr%sheet_current_factor(m)
          call log%write
       end do
    end if
  end subroutine compute_safety_factor

  subroutine check_safety_factor
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: nlabel, psipol_max, psisurf, rbeg, qsaf
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    integer :: kf, kt, ktri, fid
    type(triangle_rmp) :: tri
    real(dp), allocatable :: q(:)

    allocate(q(mesh%nflux))
    q = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          q(kf) = q(kf) + B0phi_Omega(ktri) * tri%area
       end do
       q(kf) = q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    open(newunit = fid, file = 'check_q_step.dat', status = 'replace')
    do kf = 1, mesh%nflux
       write (fid, '(3(1x, es24.16e3))') (fs_half%psi(kf) - fs%psi(0)) / &
            (fs%psi(mesh%nflux) - fs%psi(0)), fs_half%rad(kf), q(kf)
    end do
    close(fid)
    deallocate(q)
    allocate(q(nlabel))
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    q = [(fluxvar%interp(equil%qpsi, psisurf(kf) * psipol_max + fs%psi(0)), &
         kf = 1, nlabel)]
    open(newunit = fid, file = 'check_q_cont.dat', status = 'replace')
    do kf = 1, nlabel
       write (fid, '(4(1x, es24.16e3))') psisurf(kf), rbeg(kf), qsaf(kf), q(kf)
    end do
    close(fid)
    deallocate(q)
  end subroutine check_safety_factor

  subroutine cache_equilibrium_field
    use mesh_mod, only: knot, triangle, mesh_point, mesh_element
    use magdif_conf, only: longlines
    real(dp) :: r, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kf, kt, ktri, ke, fid
    type(triangle) :: elem
    type(knot) :: base, tip
    real(dp) :: n_r, n_z

    open(newunit = fid, file = 'plot_B0.dat', recl = longlines, status = 'replace')
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          elem = mesh_element(ktri)
          do ke = 1, 3
             base = mesh_point(elem%i_knot(ke))
             tip = mesh_point(elem%i_knot(mod(ke, 3) + 1))
             r = (base%rcoord + tip%rcoord) * 0.5d0
             z = (base%zcoord + tip%zcoord) * 0.5d0
             call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
                  dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
             B0r(ktri, ke) = Br
             B0phi(ktri, ke) = Bp
             B0z(ktri, ke) = Bz
             n_r = tip%zcoord - base%zcoord
             n_z = base%rcoord - tip%rcoord
             B0flux(ktri, ke) = r * (Br * n_r + Bz * n_z)
          end do
          call ring_centered_avg_coord(elem, r, z)
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          B0r_Omega(ktri) = Br
          B0phi_Omega(ktri) = Bp
          B0z_Omega(ktri) = Bz
          write (fid, '(5(1x, es24.16e3))') r, z, Br, Bz, Bp
       end do
    end do
    close(fid)
  end subroutine cache_equilibrium_field

  !> Computes equilibrium current density #j0phi from given equilibrium magnetic field and
  !> assumed equilibrium pressure #pres0.
  !>
  !> This step is necessary because equilibrium pressure is not given experimentally as is
  !> \f$ \vec{B}_{0} \f$; arbitrary values are assumed. Consistency of MHD equilibrium is
  !> necessary in the derivation, while Ampere's equation is not used.
  subroutine compute_j0phi
    use constants, only: pi  ! orbit_mod.f90
    use mesh_mod, only: triangle_rmp, mesh_element_rmp, mesh_point
    use magdif_conf, only: conf, curr_prof_ps, curr_prof_rot, curr_prof_geqdsk, longlines
    use magdif_util, only: clight
    integer :: kf, kt, ktri, fid
    real(dp) :: r, z
    real(dp) :: Btor2
    real(dp), dimension(mesh%nflux) :: B2avg, B2avg_half
    real(dp) :: plot_j0phi
    type(triangle_rmp) :: tri

    open(newunit = fid, file = conf%j0phi_file, recl = longlines, status = 'replace')
    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0r(ktri, tri%ef) ** 2 + &
               B0phi(ktri, tri%ef) ** 2 + B0z(ktri, tri%ef) ** 2
          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0r(ktri, tri%ei) ** 2 + &
               B0phi(ktri, tri%ei) ** 2 + B0z(ktri, tri%ei) ** 2
       end do
       B2avg(kf) = B2avg(kf) / mesh%kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / mesh%kt_max(kf)

       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)

          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
             if (tri%orient) then
                j0phi(ktri, tri%ef) = clight * (fs%dp_dpsi(kf) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf) / r)
             else
                j0phi(ktri, tri%ef) = clight * (fs%dp_dpsi(kf-1) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf-1) / r)
             end if
          case (curr_prof_rot)
             z = sum(mesh_point(tri%lf(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%ef) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%ef) ** 2
             if (kf > 1 .and. .not. tri%orient) then
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf-1) * (1d0 - &
                     Btor2 / B2avg(kf-1))
             else
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf) * (1d0 - &
                     Btor2 / B2avg(kf))
             end if
          end select

          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
             j0phi(ktri, tri%ei) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(tri%li(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%ei) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%ei) ** 2
             j0phi(ktri, tri%ei) = clight * r * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
             j0phi(ktri, tri%eo) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(tri%lo(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%eo) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%eo) ** 2
             j0phi(ktri, tri%eo) = clight * r * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
             plot_j0phi = clight * (fs_half%dp_dpsi(kf) * tri%r_Omega + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / tri%r_Omega)
          case (curr_prof_rot)
             plot_j0phi = j0phi_ampere(tri%r_Omega, tri%z_Omega)
          case (curr_prof_ps)
             Btor2 = B0phi_Omega(ktri) ** 2
             plot_j0phi = clight * tri%r_Omega * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          write (fid, '(4(1x, es24.16e3))') &
               j0phi(ktri, 1), j0phi(ktri, 2), j0phi(ktri, 3), plot_j0phi
       end do
    end do
    close(fid)

    ! TODO: replace by real array with kedge index
    ! call check_redundant_edges(cmplx(j0phi, 0d0, dp), .true., 'j0phi')

  contains
    function j0phi_ampere(r, z) result (rotB_phi)
      real(dp), intent(in) :: r, z
      real(dp) :: rotB_phi
      real(dp) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
           dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
      call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
           dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
      rotB_phi = 0.25d0 / pi * clight * (dBrdZ - dBzdR)
    end function j0phi_ampere

  end subroutine compute_j0phi

  subroutine check_curr0
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use magdif_conf, only: longlines
    use magdif_util, only: clight
    integer :: kf, kt, fid_amp, fid_gs, fid_prof
    real(dp) :: cmp_gradp, cmp_amp, cmp_gs, theta, R, Z, dum, B0_R, B0_phi, B0_Z, &
         dB0R_dZ, dB0phi_dR, dB0phi_dZ, dB0Z_dR, J0_R, J0_phi, J0_Z, grad_psi(3)

    open(newunit = fid_prof, file = 'cmp_prof.dat', recl = longlines, status = 'replace')
    open(newunit = fid_amp, file = 'j0_amp.dat', recl = longlines, status = 'replace')
    open(newunit = fid_gs, file = 'j0_gs.dat', recl = longlines, status = 'replace')
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          theta = (dble(kt) - 0.5d0) / dble(mesh%kt_max(kf)) * 2d0 * pi
          ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
          call magdata_in_symfluxcoord_ext(2, dum, fs_half%psi(kf) - fs%psi(0), &
               theta, dum, dum, dum, dum, dum, R, dum, dum, Z, dum, dum)
          call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dum, dum, dB0R_dZ, dB0phi_dR, &
               dum, dB0phi_dZ, dB0Z_dR, dum, dum)
          ! left-hand side of iMHD force balance
          grad_psi = [R * B0_Z, 0d0, -R * B0_R]
          cmp_gradp = fs_half%dp_dpsi(kf) * dot_product(grad_psi, grad_psi) / &
               norm2(grad_psi)
          ! current density via Grad-Shafranov equation
          J0_R = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * B0_R
          J0_Z = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * B0_Z
          J0_phi = clight * (fs_half%dp_dpsi(kf) * R + &
               0.25d0 / pi * fs_half%FdF_dpsi(kf) / R)
          write (fid_gs, '(3(1x, es24.16e3))') J0_R, J0_phi, J0_Z
          cmp_gs = dot_product([J0_phi * B0_Z - J0_Z * B0_phi, J0_Z * B0_R - J0_R * B0_Z, &
               J0_R * B0_phi - J0_phi * B0_R], grad_psi) / norm2(grad_psi) / clight
          ! current density via Ampere's equation
          J0_R = 0.25d0 / pi * clight * (-dB0phi_dZ)
          J0_phi = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          J0_Z = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
          write (fid_amp, '(3(1x, es24.16e3))') J0_R, J0_phi, J0_Z
          cmp_amp = dot_product([J0_phi * B0_Z - J0_Z * B0_phi, J0_Z * B0_R - J0_R * B0_Z, &
               J0_R * B0_phi - J0_phi * B0_R], grad_psi) / norm2(grad_psi) / clight
          write (fid_prof, '(3(1x, es24.16e3))') cmp_gradp, cmp_amp, cmp_gs
       end do
    end do
    close(fid_prof)
    close(fid_amp)
    close(fid_gs)
  end subroutine check_curr0

  subroutine flux_func_cache_check
    use magdif_conf, only: log
    log%msg = 'checking flux_func_cache...'
    if (log%debug) call log%write
    write (log%msg, '("array bounds: fs%psi(", i0, ":", i0, "), ' // &
         ' fs%rad(", i0, ":", i0, "), fs_half%psi(", i0, ":", i0, "), ' // &
         'fs_half%rad(", i0, ":", i0, ")")') lbound(fs%psi, 1), ubound(fs%psi, 1), &
         lbound(fs%rad, 1), ubound(fs%rad, 1), lbound(fs_half%psi, 1), &
         ubound(fs_half%psi, 1), lbound(fs_half%rad, 1), ubound(fs_half%rad, 1)
    if (log%debug) call log%write
    write (log%msg, '("expected sign of psi''(r): ", sp, i0, ss)') equil%cocos%sgn_dpsi
    if (log%debug) call log%write
    write (log%msg, '(i0, " ordering violations for psi")') &
         count((fs%psi(1:) - fs_half%psi) * equil%cocos%sgn_dpsi <= 0d0) + &
         count([(fs_half%psi(1) - fs%psi(0)) * equil%cocos%sgn_dpsi] <= 0d0)
    if (log%debug) call log%write
    write (log%msg, '(i0, " ordering violations for rad")') &
         count(fs%rad(1:) <= fs_half%rad) + count([fs_half%rad(1)] <= [fs%rad(0)])
    if (log%debug) call log%write
  end subroutine flux_func_cache_check

  subroutine write_fluxvar
    use magdif_conf, only: conf, longlines
    integer :: kf, fid

    open(newunit = fid, file = conf%fluxvar_file, recl = longlines, status = 'replace')
    do kf = 0, mesh%nflux
       write (fid, '(5(1x, es24.16e3))') &
            fs%rad(kf), fs%psi(kf), fs%q(kf), fs%p(kf), fs%dp_dpsi(kf)
    end do
    close(fid)
  end subroutine write_fluxvar

end module magdif_mesh
