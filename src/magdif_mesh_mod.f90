module magdif_mesh_mod

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: generate_mesh, write_kilca_convexfile, kilca_vacuum

contains

  subroutine generate_mesh(unprocessed_geqdsk)
    use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp
    use magdif_conf, only: conf, log
    use magdif_util, only: get_equil_filenames, initialize_globals
    use magdif, only: equil, cache_mesh_data, Bnflux, Bnphi

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

    call create_mesh_points(convexfile)
    call connect_mesh_points
    call write_mesh_data
    call cache_mesh_data
    if (conf%kilca_scale_factor /= 0) then
       call compute_kilca_vac_coeff
       call compute_kilca_vacuum
       call check_kilca_vacuum
       call check_RT0
    end if

    if (allocated(Bnflux)) deallocate(Bnflux)
    if (allocated(Bnphi)) deallocate(Bnphi)
    if (allocated(mesh_element_rmp)) deallocate(mesh_element)
    if (allocated(mesh_element)) deallocate(mesh_element)
    if (allocated(mesh_point)) deallocate(mesh_point)
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
       conf%nflux = conf%nflux_unref
       allocate(refined(0:conf%nflux))
       refined = linspace(0d0, 1d0, conf%nflux + 1, 0, 0)
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
    conf%nflux = conf%nflux_unref + 2 * sum(additions - deletions)
    allocate(refined(0:conf%nflux))
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
    refined(fine_hi(nref)+1:conf%nflux) = linspace(refined(fine_hi(nref)), 1d0, &
         conf%nflux - fine_hi(nref), 1, 0)
  end subroutine refine_eqd_partition

  subroutine refine_resonant_surfaces(psi_sample, q_sample, psi2rho_norm, rho_norm_ref)
    use magdif_conf, only: conf, conf_arr, log
    use magdif_util, only: flux_func
    use magdif, only: mesh
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
    integer :: m, kref
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
    mesh%m_res_min = ceiling(minval(abs(q_sample)) * dble(conf%n))
    mesh%m_res_max = floor(maxval(abs(q_sample)) * dble(conf%n))
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
    ! TODO: make separate array variables in mesh type for actual m_[res_]min/max
    allocate(mask(mesh%m_res_min:mesh%m_res_max))
    mask = 0d0 < conf_arr%refinement .and. conf_arr%refinement < 1d0
    allocate(ref_ind(count(mask)))
    call refine_eqd_partition(count(mask), pack(conf_arr%deletions, mask), &
         pack(conf_arr%additions, mask), pack(conf_arr%refinement, mask), &
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
      q_interp_resonant = psi_eqd%interp(abs(q_sample), psi) - dble(m) / dble(conf%n)
    end function q_interp_resonant
  end subroutine refine_resonant_surfaces

  subroutine write_kilca_convexfile(rho_max, convexfile)
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    use magdif, only: equil

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
    use magdif, only: equil, mesh, init_indices, fs, fs_half, flux_func_cache_check
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
         raxis, zaxis
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
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    call psi_interpolator%init(4, psisurf(1:) * psipol_max + psi_axis)
    ! interpolate between psi and rho
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd = rbeg / hypot(theta_axis(1), theta_axis(2))

    call refine_resonant_surfaces(psisurf(1:) * psipol_max + psi_axis, qsaf, &
         psi2rho_norm, rho_norm_ref)
    call fs%init(conf%nflux, .false.)
    call fs_half%init(conf%nflux, .true.)
    fs%rad = rho_norm_ref
    fs%psi = [(interp_psi_pol(raxis + fs%rad(kf) * theta_axis(1), &
         zaxis + fs%rad(kf) * theta_axis(2)), kf = 0, conf%nflux)]
    fs_half%rad = 0.5d0 * (fs%rad(0:conf%nflux-1) + fs%rad(1:conf%nflux))
    fs_half%psi = [(interp_psi_pol(raxis + fs_half%rad(kf) * theta_axis(1), &
         zaxis + fs_half%rad(kf) * theta_axis(2)), kf = 1, conf%nflux)]
    fs%rad = fs%rad * hypot(theta_axis(1), theta_axis(2))
    fs_half%rad = fs_half%rad * hypot(theta_axis(1), theta_axis(2))
    call flux_func_cache_check
    ! dump presumedly optimal values for poloidal resolution
    open(newunit = fid, file = 'optpolres.dat', status = 'replace')
    do kf = 1, conf%nflux - 1
       write (fid, '(2(1x, es24.16e3))') fs%rad(kf), 2d0 * pi * fs%rad(kf) / &
            (fs_half%rad(kf + 1) - fs_half%rad(kf))
    end do
    write (fid, '(2(1x, es24.16e3))') fs%rad(conf%nflux), 2d0 * pi * fs%rad(conf%nflux) / &
         (fs%rad(conf%nflux) - fs%rad(conf%nflux - 1))
    close(fid)

    call init_indices
    ntri = mesh%kt_low(conf%nflux + 1)
    npoint = mesh%kp_low(conf%nflux + 1)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))

    allocate(n_theta(conf%nflux))
    allocate(points(3, npoint))
    allocate(points_s_theta_phi(3, npoint))
    n_theta = conf%nkpol
    s_min = 1d-16
    ! inp_label 2 to use poloidal psi with magdata_in_symfluxcoord_ext
    ! psi is not normalized by psipol_max, but shifted by -psi_axis
    call create_points_2d(2, n_theta, points, points_s_theta_phi, r_scaling_func = psi_ref)
    mesh_point(1)%rcoord = raxis
    mesh_point(1)%zcoord = zaxis
    do kpoi = 2, npoint
       mesh_point(kpoi)%rcoord = points(1, kpoi)
       mesh_point(kpoi)%zcoord = points(3, kpoi)
    end do

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
      if (conf%nflux /= size(psi_eqd)) then
         call log%msg_arg_size('psi_ref', 'nflux', 'size(psi_eqd)', &
              conf%nflux, size(psi_eqd))
         if (log%err) call log%write
         error stop
      end if
      psi_ref = [(interp_psi_pol(raxis + rho_norm_ref(kf) * theta_axis(1), &
           zaxis + rho_norm_ref(kf) * theta_axis(2)) - psi_axis, kf = 1, conf%nflux)]
    end function psi_ref
  end subroutine create_mesh_points

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
    use magdif_conf, only: conf
    use magdif, only: mesh
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
    do kf = 2, conf%nflux
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
    do kf = 2, conf%nflux
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
    do kf = 1, conf%nflux - 1
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
    kf = conf%nflux
    do kt = 1, mesh%kt_max(kf)
       if (mod(kt, 2) == 1) then
          ktri = mesh%kt_low(kf) + kt
          ktri_adj = ktri + mesh%kt_max(kf) + 1
          mesh_element(ktri)%neighbour(1) = ktri_adj
          mesh_element(ktri)%neighbour_edge(1) = 2
       end if
    end do
  end subroutine connect_mesh_points

  subroutine write_mesh_data
    use mesh_mod, only: npoint, ntri, knot, triangle, mesh_point, mesh_element
    use magdif_conf, only: conf, longlines
    use magdif, only: mesh, fs, fs_half, flux_func_cache_check

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
    open(newunit = fid, file = conf%meshdata_file, form = 'unformatted', status = 'replace')
    write (fid) conf%nflux, npoint, ntri, mesh%m_res_min, mesh%m_res_max
    write (fid) fs%psi, fs%rad
    write (fid) fs_half%psi, fs_half%rad
    write (fid) mesh_point
    write (fid) mesh_element
    close(fid)

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace')
    write (fid, '(3(1x, i0))') npoint, ntri, mesh%kp_max(conf%nflux) - 1
    do kpoi = 1, mesh%kp_low(conf%nflux + 1)
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 0
    end do
    do kpoi = mesh%kp_low(conf%nflux + 1) + 1, npoint
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 1
    end do
    do ktri = 1, ntri
       write (fid, '(4(1x, i0))') mesh_element(ktri)%i_knot(:), 0
    end do
    do kp = 1, mesh%kp_max(conf%nflux) - 1
       write (fid, '(4(1x, i0))') mesh%kp_low(conf%nflux) + kp, mesh%kp_low(conf%nflux) + kp + 1, 1
    end do
    close(fid)
  end subroutine write_mesh_data

  ! calculate resonant vacuum perturbation
  subroutine compute_kilca_vacuum
    use mesh_mod, only: ntri, mesh_element_rmp, mesh_point, mesh_element
    use magdif_conf, only: conf
    use magdif_util, only: imun, gauss_legendre_unit_interval
    use magdif, only: equil, Bnflux, Bnphi, check_redundant_edges, check_div_free, &
         write_vector_dof

    integer, parameter :: order = 2
    integer :: ktri, k, ke, pol_modes(2)
    real(dp) :: R, Z, rho, theta, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_phi, B_Z

    call gauss_legendre_unit_interval(order, points, weights)
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    Bnflux = (0d0, 0d0)
    Bnphi = (0d0, 0d0)
    pol_modes = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
    do ktri = 1, ntri
       associate(tri => mesh_element_rmp(ktri), &
            knots => mesh_point(mesh_element(ktri)%i_knot))
         node_R = [knots(:)%rcoord, knots(1)%rcoord]
         node_Z = [knots(:)%zcoord, knots(1)%zcoord]
         do ke = 1, 3
            edge_R = node_R(ke + 1) - node_R(ke)
            edge_Z = node_Z(ke + 1) - node_Z(ke)
            do k = 1, order
               R = node_R(ke) * points(k) + node_R(ke + 1) * points(order - k + 1)
               Z = node_Z(ke) * points(k) + node_Z(ke + 1) * points(order - k + 1)
               rho = hypot(R - equil%rmaxis, Z - equil%zmaxis)
               theta = atan2(Z - equil%zmaxis, R - equil%rmaxis)
               call kilca_vacuum(conf%n, pol_modes, equil%rcentr, rho, theta, B_R, B_phi, B_Z)
               Bnflux(ktri, ke) = Bnflux(ktri, ke) + &
                    weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
            end do
         end do
         ! toroidal flux via zero divergence
         Bnphi(ktri) = imun / conf%n * sum(Bnflux(ktri, :)) / tri%area
       end associate
    end do
    call check_redundant_edges(Bnflux, .false., 'vacuum B_n')
    call check_div_free(Bnflux, Bnphi, conf%n, 1d-9, 'vacuum B_n')
    call write_vector_dof(Bnflux, Bnphi, conf%Bn_vac_file)
  end subroutine compute_kilca_vacuum

  !> Calculate the vacuum perturbation field in cylindrical coordinates from the Fourier
  !> series of all given modes.
  !>
  !> @param tor_mode toroidal mode number, scaled by magdif_config::kilca_scale_factor
  !> (usually magdif_config::n)
  !> @param pol_modes array of poloidal mode numbers
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> magdif_util::g_eqdsk::rcentr)
  !> @param r radial distance \f$ r \f$ from magnetic axis
  !> @param theta geometrical poloidal angle \f$ theta \f$ (coinciding with symmetry flux
  !> coordinates' poloidal angle in this geometry)
  !> @param B_R physical component \f$ B_{R} (r, theta, n) \f$ of the vacuum perturbation
  !> field
  !> @param B_phi physical component \f$ B_{(\varphi)} (r, theta, n) \f$ of the vacuum
  !> perturbation field
  !> @param B_Z physical component \f$ B_{Z} (r, theta, n) \f$ of the vacuum perturbation
  !> field
  subroutine kilca_vacuum(tor_mode, pol_modes, R_0, r, theta, B_R, B_phi, B_Z)
    use magdif_conf, only: conf_arr
    use magdif_util, only: imun, straight_cyl2bent_cyl
    integer, intent(in) :: tor_mode, pol_modes(1:)
    real(dp), intent(in) :: R_0, r, theta
    complex(dp), intent(out) :: B_R, B_phi, B_Z
    complex(dp) :: B_rad, B_pol, B_tor, temp_B_rad, temp_B_pol, temp_B_tor
    complex(dp), dimension(size(pol_modes)) :: fourier_basis
    integer :: k

    B_rad = (0d0, 0d0)
    B_pol = (0d0, 0d0)
    B_tor = (0d0, 0d0)
    fourier_basis = exp(imun * pol_modes * theta)
    do k = 1, ubound(pol_modes, 1)
       call kilca_vacuum_fourier(tor_mode, pol_modes(k), R_0, r, &
            conf_arr%kilca_vac_coeff(abs(pol_modes(k))), temp_B_rad, temp_B_pol, temp_B_tor)
       B_rad = B_rad + temp_B_rad * fourier_basis(k)
       B_pol = B_pol + temp_B_pol * fourier_basis(k)
       B_tor = B_tor + temp_B_tor * fourier_basis(k)
    end do
    call straight_cyl2bent_cyl(B_rad, B_pol, B_tor, theta, B_R, B_phi, B_Z)
  end subroutine kilca_vacuum

  subroutine compute_kilca_vac_coeff
    use magdif_conf, only: conf, conf_arr, log, cmplx_fmt
    use magdif, only: mesh, equil
    integer :: m
    complex(dp) :: B_rad, B_pol, B_tor

    do m = mesh%m_res_min, mesh%m_res_max
       if (abs(conf_arr%kilca_vac_r(m)) <= 0d0) then
          write (log%msg, '("ignoring kilca_vac_r(", i0, "), ' // &
               'resorting to kilca_vac_coeff(", i0, ")")') m, m
          if (log%info) call log%write
          cycle
       end if
       call kilca_vacuum_fourier(conf%n, m, equil%rcentr, conf_arr%kilca_vac_r(m), (1d0, 0d0), &
            B_rad, B_pol, B_tor)
       conf_arr%kilca_vac_coeff(m) = conf_arr%kilca_vac_Bz(m) / B_tor
    end do
    log%msg = 'effective vacuum perturbation field coefficients:'
    if (log%info) call log%write
    do m = mesh%m_res_min, mesh%m_res_max
       write (log%msg, '("kilca_vac_coeff(", i0, ") = ", ' // cmplx_fmt // ')') m, &
            conf_arr%kilca_vac_coeff(m)
       if (log%info) call log%write
    end do
  end subroutine compute_kilca_vac_coeff

  !> Calculate the Fourier coefficient of the vacuum perturbation field for a given
  !> toroidal-poloidal mode.
  !>
  !> @param tor_mode toroidal mode number, scaled by magdif_config::kilca_scale_factor
  !> (usually magdif_config::n)
  !> @param pol_mode poloidal mode number
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> magdif_util::g_eqdsk::rcentr)
  !> @param r radial distance \f$ r \f$ from magnetic axis
  !> @param vac_coeff coefficient of modified Bessel functions (integration constant)
  !> @param B_rad physical component \f$ B_{r} (r, m, n) \f$ of the vacuum perturbation
  !> field
  !> @param B_pol physical component \f$ B_{(\theta)} (r, m, n) \f$ of the vacuum
  !> perturbation field
  !> @param B_tor physical component \f$ B_{z} (r, m, n) \f$ of the vacuum perturbation
  !> field
  subroutine kilca_vacuum_fourier(tor_mode, pol_mode, R_0, r, vac_coeff, &
       B_rad, B_pol, B_tor)
    use fgsl, only: fgsl_double, fgsl_int, fgsl_success, fgsl_sf_bessel_icn_array
    use magdif_conf, only: log
    use magdif_util, only: imun
    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r
    complex(dp), intent(in) :: vac_coeff
    complex(dp), intent(out) :: B_rad, B_pol, B_tor
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(abs(pol_mode)-1, abs(pol_mode)+1, k_z_r, I_m)
    if (status /= fgsl_success .and. log%err) then
       write (log%msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
       call log%write
    end if
    B_rad = 0.5d0 * (I_m(-1) + I_m(1)) * vac_coeff
    B_pol = imun * pol_mode / k_z_r * I_m(0) * vac_coeff
    B_tor = imun * I_m(0) * vac_coeff
  end subroutine kilca_vacuum_fourier

  subroutine check_kilca_vacuum
    use magdif_conf, only: conf, conf_arr, longlines
    use magdif, only: equil, fs_half
    complex(dp) :: B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    real(dp) :: rad
    integer :: kf, fid, abs_pol_mode

    abs_pol_mode = abs(conf%kilca_pol_mode)
    open(newunit = fid, file = 'cmp_vac.dat', recl = 3 * longlines)
    do kf = 1, conf%nflux
       rad = fs_half%rad(kf)
       call kilca_vacuum_fourier(conf%n, -abs_pol_mode, equil%rcentr, rad, &
            conf_arr%kilca_vac_coeff(abs_pol_mode), B_rad_neg, B_pol_neg, B_tor_neg)
       call kilca_vacuum_fourier(conf%n, abs_pol_mode, equil%rcentr, rad, &
            conf_arr%kilca_vac_coeff(abs_pol_mode), B_rad_pos, B_pol_pos, B_tor_pos)
       write (fid, '(13(1x, es24.16e3))') rad, &
           B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    end do
    close(fid)
  end subroutine check_kilca_vacuum

  subroutine check_RT0
    use magdif_conf, only: conf, longlines
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    use magdif, only: equil, fs_half, Bnflux, Bnphi, point_location, interp_RT0
    integer :: fid, kf, kpol, ktri
    real(dp) :: rad, theta, R, Z
    complex(dp) :: B_R, B_Z, B_phi, B_R_interp, B_Z_interp, B_phi_interp
    integer :: pol_modes(2)

    pol_modes = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
    open(newunit = fid, file = 'cmp_RT0.dat', recl = longlines)
    do kf = conf%nflux / 3, conf%nflux / 3 ! 1, nflux
       rad = fs_half%rad(kf)
       do kpol = 1, 2 * conf%nkpol
          theta = (kpol - 0.5d0) / dble(2 * conf%nkpol) * 2d0 * pi
          call kilca_vacuum(conf%n, pol_modes, equil%rcentr, rad, theta, B_R, B_phi, B_Z)
          R = equil%rmaxis + rad * cos(theta)
          Z = equil%zmaxis + rad * sin(theta)
          ktri = point_location(R, Z)
          call interp_RT0(ktri, Bnflux, R, Z, B_R_interp, B_Z_interp)
          B_phi_interp = Bnphi(ktri)
          write (fid, '(14(1x, es23.15e3))') rad, theta, B_R, B_phi, B_Z, &
               B_R_interp, B_phi_interp, B_Z_interp
       end do
    end do
    close(fid)
  end subroutine check_RT0
end module magdif_mesh_mod
