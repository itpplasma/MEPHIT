module magdif_mesh_mod

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: generate_mesh, write_kilca_convexfile, kilca_vacuum

contains

  subroutine generate_mesh(unprocessed_geqdsk)
    use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp
    use magdif_config, only: log_msg, log_info, log_write, kilca_scale_factor
    use magdif_util, only: get_equil_filenames, initialize_globals
    use magdif, only: equil, cache_mesh_data, Bnflux, Bnphi

    character(len = *), intent(in) :: unprocessed_geqdsk
    character(len = 1024) :: gfile, convexfile

    call get_equil_filenames(gfile, convexfile)
    log_msg = 'attempting to read unprocessed G EQDSK file ' // trim(unprocessed_geqdsk)
    if (log_info) call log_write
    call equil%read(trim(unprocessed_geqdsk))
    call equil%classify
    call equil%standardise
    if (kilca_scale_factor /= 0) then
       call equil%scale(kilca_scale_factor)
    end if
    log_msg = 'attempting to write processed G EQDSK file ' // trim(gfile)
    if (log_info) call log_write
    call equil%write(trim(gfile))
    call initialize_globals(equil%rmaxis, equil%zmaxis)

    call create_mesh_points(convexfile)
    call connect_mesh_points
    call write_mesh_data
    call cache_mesh_data
    if (kilca_scale_factor /= 0) then
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
    use magdif_config, only: nflux_unref, nflux, log_msg_arg_size, log_err, log_write
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
       nflux = nflux_unref
       allocate(refined(0:nflux))
       refined = linspace(0d0, 1d0, nflux + 1, 0, 0)
       return
    end if
    if (nref /= size(deletions)) then
       call log_msg_arg_size('refine_eqd_partition', 'nref', 'size(deletions)', nref, &
            size(deletions))
       if (log_err) call log_write
       error stop
    end if
    if (nref /= size(additions)) then
       call log_msg_arg_size('refine_eqd_partition', 'nref', 'size(additions)', nref, &
            size(additions))
       if (log_err) call log_write
       error stop
    end if
    if (nref /= size(refinement)) then
       call log_msg_arg_size('refine_eqd_partition', 'nref', 'size(refinement)', nref, &
            size(refinement))
       if (log_err) call log_write
       error stop
    end if
    if (nref /= size(res)) then
       call log_msg_arg_size('refine_eqd_partition', 'nref', 'size(res)', nref, size(res))
       if (log_err) call log_write
       error stop
    end if
    if (nref /= size(ref_ind)) then
       call log_msg_arg_size('refine_eqd_partition', 'nref', 'size(ref_ind)', nref, &
            size(ref_ind))
       if (log_err) call log_write
       error stop
    end if
    nflux = nflux_unref + 2 * sum(additions - deletions)
    allocate(refined(0:nflux))
    refined = 0d0
    coarse_lo = floor(res * (nflux_unref + 1)) - deletions
    coarse_hi = ceiling(res * (nflux_unref + 1)) + deletions
    ! compute upper and lower array indices of refined regions
    fine_lo = coarse_lo + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 0, nref - 1)]
    fine_hi = coarse_hi + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 1, nref)]
    ref_ind = fine_hi - additions
    ! compute separations between flux surfaces using geometric series
    coarse_sep = 1d0 / dble(nflux_unref)
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
    refined(fine_hi(nref)+1:nflux) = linspace(refined(fine_hi(nref)), 1d0, &
         nflux - fine_hi(nref), 1, 0)
  end subroutine refine_eqd_partition

  subroutine refine_resonant_surfaces(psi_sample, q_sample, psi2rho_norm, rho_norm_ref)
    use magdif_config, only: n, refinement, additions, deletions, sheet_current_factor, &
         read_delayed_config, log_msg, log_debug, log_err, log_write, log_msg_arg_size
    use magdif_util, only: flux_func
    use magdif, only: m_res_min, m_res_max
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
       call log_msg_arg_size('refine_resonant_surfaces', 'size(psi_sample)', &
            'size(q_sample)', size(psi_sample), size(q_sample))
       if (log_err) call log_write
       error stop
    end if
    psi_min = minval(psi_sample)
    psi_max = maxval(psi_sample)
    call psi_eqd%init(4, psi_sample)
    m_res_min = ceiling(minval(abs(q_sample)) * dble(n))
    m_res_max = floor(maxval(abs(q_sample)) * dble(n))
    allocate(psi_res(m_res_min:m_res_max))
    allocate(rho_res(m_res_min:m_res_max))
    log_msg = 'resonance positions:'
    if (log_debug) call log_write
    do m = m_res_min, m_res_max
       psi_res(m) = zeroin(psi_min, psi_max, q_interp_resonant, 1d-9)
       rho_res(m) = psi2rho_norm(psi_res(m))
       write (log_msg, '("m = ", i2, ", psi_m = ", es24.16e3, ", rho_m = ", f19.16)') &
            m, psi_res(m), rho_res(m)
       if (log_debug) call log_write
    end do
    call read_delayed_config(m_res_min, m_res_max)
    allocate(mask(m_res_min:m_res_max))
    mask = 0d0 < refinement .and. refinement < 1d0
    allocate(ref_ind(count(mask)))
    call refine_eqd_partition(count(mask), pack(deletions, mask), pack(additions, mask), &
         pack(refinement, mask), pack(rho_res, mask), rho_norm_ref, ref_ind)

    log_msg = 'refinement positions:'
    if (log_debug) call log_write
    kref = 0
    do m = m_res_min, m_res_max
       if (.not. mask(m)) cycle
       kref = kref + 1
       write (log_msg, '("m = ", i0, ", kf = ", i0, ' // &
            '", rho: ", f19.16, 2(" < ", f19.16))') m, ref_ind(kref), &
            rho_norm_ref(ref_ind(kref) - 1), rho_res(m), rho_norm_ref(ref_ind(kref))
       if (log_debug) call log_write
    end do

    if (allocated(ref_ind)) deallocate(ref_ind)
    if (allocated(mask)) deallocate(mask)
    if (allocated(rho_res)) deallocate(rho_res)
    if (allocated(psi_res)) deallocate(psi_res)
    if (allocated(refinement)) deallocate(refinement)
    if (allocated(deletions)) deallocate(deletions)
    if (allocated(additions)) deallocate(additions)
    if (allocated(sheet_current_factor)) deallocate(sheet_current_factor)

  contains
    function q_interp_resonant(psi)
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant
      q_interp_resonant = psi_eqd%interp(abs(q_sample), psi) - dble(m) / dble(n)
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
    use mesh_mod, only: npoint, mesh_point, ntri, mesh_element, mesh_element_rmp
    use magdif_config, only: nflux, nkpol, kilca_scale_factor, log_msg_arg_size, &
         log_err, log_write
    use magdif_util, only: interp_psi_pol, flux_func
    use magdif, only: equil, kp_low, kt_low, init_indices, fs, fs_half, &
         flux_func_cache_check
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
         raxis, zaxis
    use field_line_integration_mod, only: circ_mesh_scale, o_point, x_point, &
         theta0_at_xpoint, theta_axis
    use points_2d, only: s_min, create_points_2d

    character(len = *), intent(in) :: convexfile
    integer :: kf, kpoi
    integer, dimension(:), allocatable :: n_theta
    real(dp), dimension(:), allocatable :: rho_norm_eqd, rho_norm_ref
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi
    type(flux_func) :: psi_interpolator
    real(dp) :: psi_axis, rho_max

    theta0_at_xpoint = .true.
    circ_mesh_scale = kilca_scale_factor
    if (kilca_scale_factor /= 0) then
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

    call refine_resonant_surfaces(psisurf(1:) * psipol_max + psi_axis, qsaf, psi2rho_norm, &
         rho_norm_ref)
    call fs%init(nflux, .false.)
    call fs_half%init(nflux, .true.)
    fs%rad = rho_norm_ref
    fs%psi = [(interp_psi_pol(raxis + fs%rad(kf) * theta_axis(1), &
         zaxis + fs%rad(kf) * theta_axis(2)), kf = 0, nflux)]
    fs_half%rad = 0.5d0 * (fs%rad(0:nflux-1) + fs%rad(1:nflux))
    fs_half%psi = [(interp_psi_pol(raxis + fs_half%rad(kf) * theta_axis(1), &
         zaxis + fs_half%rad(kf) * theta_axis(2)), kf = 1, nflux)]
    fs%rad = fs%rad * hypot(theta_axis(1), theta_axis(2))
    fs_half%rad = fs_half%rad * hypot(theta_axis(1), theta_axis(2))
    call flux_func_cache_check

    call init_indices
    ntri = kt_low(nflux+1)
    npoint = kp_low(nflux+1)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))

    allocate(n_theta(nflux))
    allocate(points(3, npoint))
    allocate(points_s_theta_phi(3, npoint))
    n_theta = nkpol
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
      if (nflux /= size(psi_eqd)) then
         call log_msg_arg_size('psi_ref', 'nflux', 'size(psi_eqd)', nflux, size(psi_eqd))
         if (log_err) call log_write
         error stop
      end if
      psi_ref = [(interp_psi_pol(raxis + rho_norm_ref(kf) * theta_axis(1), &
           zaxis + rho_norm_ref(kf) * theta_axis(2)) - psi_axis, kf = 1, nflux)]
    end function psi_ref
  end subroutine create_mesh_points

  subroutine connect_mesh_points
    use mesh_mod, only: mesh_element
    use magdif_config, only: nflux
    use magdif_util, only: interleave, calculate_det_3
    use magdif, only: kt_max, kt_low, kp_low

    integer :: kf, k
    integer, dimension(maxval(kt_max)) :: kq, kt, tri_f, tri_i, tri_o
    integer, dimension(4, maxval(kt_max)) :: quad

    kt = [(k, k = 1, maxval(kt_max))]
    kf = 1
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%i_knot(1) = &
         kt_low(kf) + kt(1:kt_max(kf)) + 1
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%i_knot(2) = &
         kt_low(kf) + mod(kt(1:kt_max(kf)), kt_max(kf)) + 2
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%i_knot(3) = 1
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%knot_h = 3
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour(1) = &
         kt_low(kf+1) + 2 * kt(1:kt_max(kf))
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour(2) = &
         mod(kt(1:kt_max(kf)), kt_max(kf)) + 1
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour(3) = &
         mod(kt(1:kt_max(kf)) + kt_max(kf) - 2, kt_max(kf)) + 1
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour_edge(1) = 2
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour_edge(2) = 3
    mesh_element(kt_low(kf)+1:kt_low(kf+1))%neighbour_edge(3) = 2
    kq = (kt + 1) / 2
    do kf = 2, nflux
       ! assign nodes of quadrilaterals in ascending global order
       quad(1,:) = kp_low(kf-1) + kq
       quad(2,:) = kp_low(kf-1) + mod(kq, kt_max(kf) / 2) + 1
       quad(3,:) = kp_low(kf) + kq
       quad(4,:) = kp_low(kf) + mod(kq, kt_max(kf) / 2) + 1
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%i_knot(1) = &
            interleave(quad(3,:), quad(4,:), kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%i_knot(2) = &
            interleave(quad(4,:), quad(2,:), kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%i_knot(3) = quad(1, 1:kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%knot_h = &
            interleave(3, 1, kt_max(kf))
       tri_o = kt_low(kf) + mod(kt, kt_max(kf)) + 1
       tri_i = kt_low(kf) + mod(kt + kt_max(kf) - 2, kt_max(kf)) + 1
       if (kf == 2) then
          tri_f = interleave(kt_low(kf+1) + kt + 1, kt_low(kf-1) + kq, kt_max(kf))
!!$     elseif (kf == nflux) then
!!$        tri_f = interleave(0, kt_low(kf-1) + kt - 1, kt_max(kf))
       else
          tri_f = interleave(kt_low(kf+1) + kt + 1, kt_low(kf-1) + kt - 1, kt_max(kf))
       end if
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour(1) = &
            interleave(tri_f, tri_o, kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour(2) = &
            interleave(tri_o, tri_f, kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour(3) = tri_i(1:kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour_edge(1) = &
            interleave(2, 3, kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour_edge(2) = &
            interleave(3, 1, kt_max(kf))
       mesh_element(kt_low(kf)+1:kt_low(kf)+kt_max(kf))%neighbour_edge(3) = &
            interleave(1, 2, kt_max(kf))
    end do
    call calculate_det_3(mesh_element(1:kt_low(nflux+1)))
  end subroutine connect_mesh_points

  subroutine write_mesh_data
    use mesh_mod, only: npoint, ntri, knot, triangle, mesh_point, mesh_element
    use magdif_config, only: nflux, longlines, meshdata_file
    use magdif, only: kp_max, kp_low, fs, fs_half, m_res_min, m_res_max, &
         flux_func_cache_check

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
    open(newunit = fid, file = meshdata_file, form = 'unformatted', status = 'replace')
    write (fid) nflux, npoint, ntri, m_res_min, m_res_max
    write (fid) fs%psi, fs%rad
    write (fid) fs_half%psi, fs_half%rad
    write (fid) mesh_point
    write (fid) mesh_element
    close(fid)

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace')
    write (fid, '(3(1x, i0))') npoint, ntri, kp_max(nflux+1) - 1
    do kpoi = 1, kp_low(nflux+1)
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 0
    end do
    do kpoi = kp_low(nflux+1) + 1, npoint
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 1
    end do
    do ktri = 1, ntri
       write (fid, '(4(1x, i0))') mesh_element(ktri)%i_knot(:), 0
    end do
    do kp = 1, kp_max(nflux+1) - 1
       write (fid, '(4(1x, i0))') kp_low(nflux+1) + kp, kp_low(nflux+1) + kp + 1, 1
    end do
    close(fid)
  end subroutine write_mesh_data

  ! calculate resonant vacuum perturbation
  subroutine compute_kilca_vacuum
    use mesh_mod, only: ntri, triangle_rmp, mesh_element_rmp, mesh_point
    use magdif_config, only: n, kilca_pol_mode, Bn_vac_file
    use magdif_util, only: imun, gauss_legendre_unit_interval
    use magdif, only: equil, Bnflux, Bnphi, check_redundant_edges, check_div_free, &
         write_vector_dof

    integer, parameter :: order = 2
    integer :: ktri, k
    real(dp) :: R, Z, rho, theta, edge_R, edge_Z, base_R, base_Z, tip_R, tip_Z
    real(dp), dimension(order) :: points, weights
    integer :: pol_modes(2)
    complex(dp) :: B_R, B_phi, B_Z
    type(triangle_rmp) :: tri

    call gauss_legendre_unit_interval(order, points, weights)
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    Bnflux = (0d0, 0d0)
    Bnphi = (0d0, 0d0)
    pol_modes = [kilca_pol_mode, -kilca_pol_mode]
    do ktri = 1, ntri
       tri = mesh_element_rmp(ktri)
       ! flux through edge f
       base_R = mesh_point(tri%lf(1))%rcoord
       base_Z = mesh_point(tri%lf(1))%zcoord
       tip_R = mesh_point(tri%lf(2))%rcoord
       tip_Z = mesh_point(tri%lf(2))%zcoord
       edge_R = tip_R - base_R
       edge_Z = tip_Z - base_Z
       do k = 1, order
          R = base_R * points(k) + tip_R * points(order - k + 1)
          Z = base_Z * points(k) + tip_Z * points(order - k + 1)
          rho = hypot(R - equil%rmaxis, Z - equil%zmaxis)
          theta = atan2(Z - equil%zmaxis, R - equil%rmaxis)
          call kilca_vacuum(n, pol_modes, equil%rcentr, rho, theta, B_R, B_phi, B_Z)
          Bnflux(ktri, tri%ef) = Bnflux(ktri, tri%ef) + &
               weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
       end do
       ! flux through edge i
       base_R = mesh_point(tri%li(1))%rcoord
       base_Z = mesh_point(tri%li(1))%zcoord
       tip_R = mesh_point(tri%li(2))%rcoord
       tip_Z = mesh_point(tri%li(2))%zcoord
       edge_R = tip_R - base_R
       edge_Z = tip_Z - base_Z
       do k = 1, order
          R = base_R * points(k) + tip_R * points(order - k + 1)
          Z = base_Z * points(k) + tip_Z * points(order - k + 1)
          rho = hypot(R - equil%rmaxis, Z - equil%zmaxis)
          theta = atan2(Z - equil%zmaxis, R - equil%rmaxis)
          call kilca_vacuum(n, pol_modes, equil%rcentr, rho, theta, B_R, B_phi, B_Z)
          Bnflux(ktri, tri%ei) = Bnflux(ktri, tri%ei) + &
               weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
       end do
       ! flux through edge o
       base_R = mesh_point(tri%lo(1))%rcoord
       base_Z = mesh_point(tri%lo(1))%zcoord
       tip_R = mesh_point(tri%lo(2))%rcoord
       tip_Z = mesh_point(tri%lo(2))%zcoord
       edge_R = tip_R - base_R
       edge_Z = tip_Z - base_Z
       do k = 1, order
          R = base_R * points(k) + tip_R * points(order - k + 1)
          Z = base_Z * points(k) + tip_Z * points(order - k + 1)
          rho = hypot(R - equil%rmaxis, Z - equil%zmaxis)
          theta = atan2(Z - equil%zmaxis, R - equil%rmaxis)
          call kilca_vacuum(n, pol_modes, equil%rcentr, rho, theta, B_R, B_phi, B_Z)
          Bnflux(ktri, tri%eo) = Bnflux(ktri, tri%eo) + &
               weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
       end do
       ! toroidal flux
       Bnphi(ktri) = imun / n * sum(Bnflux(ktri, :)) / tri%area
    end do
    call check_redundant_edges(Bnflux, .false., 'vacuum B_n')
    call check_div_free(Bnflux, Bnphi, n, 1d-9, 'vacuum B_n')
    call write_vector_dof(Bnflux, Bnphi, Bn_vac_file)
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
    use magdif_config, only: kilca_vac_coeff
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
            kilca_vac_coeff(abs(pol_modes(k))), temp_B_rad, temp_B_pol, temp_B_tor)
       B_rad = B_rad + temp_B_rad * fourier_basis(k)
       B_pol = B_pol + temp_B_pol * fourier_basis(k)
       B_tor = B_tor + temp_B_tor * fourier_basis(k)
    end do
    call straight_cyl2bent_cyl(B_rad, B_pol, B_tor, theta, B_R, B_phi, B_Z)
  end subroutine kilca_vacuum

  subroutine compute_kilca_vac_coeff
    use magdif_config, only: kilca_vac_r, kilca_vac_Bz, kilca_vac_coeff, n, cmplx_fmt, &
         log_msg, log_info, log_write
    use magdif, only: m_res_min, m_res_max, equil
    integer :: m
    complex(dp) :: B_rad, B_pol, B_tor

    do m = m_res_min, m_res_max
       if (abs(kilca_vac_r(m)) <= 0d0) then
          write (log_msg, '("ignoring kilca_vac_r(", i0, "), ' // &
               'resorting to kilca_vac_coeff(", i0, ")")') m, m
          if (log_info) call log_write
          cycle
       end if
       call kilca_vacuum_fourier(n, m, equil%rcentr, kilca_vac_r(m), (1d0, 0d0), &
            B_rad, B_pol, B_tor)
       kilca_vac_coeff(m) = kilca_vac_Bz(m) / B_tor
    end do
    log_msg = 'effective vacuum perturbation field coefficients:'
    if (log_info) call log_write
    do m = m_res_min, m_res_max
       write (log_msg, '("kilca_vac_coeff(", i0, ") = ", ' // cmplx_fmt // ')') m, &
            kilca_vac_coeff(m)
       if (log_info) call log_write
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
    use magdif_config, only: log_msg, log_err, log_write
    use magdif_util, only: imun
    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r
    complex(dp), intent(in) :: vac_coeff
    complex(dp), intent(out) :: B_rad, B_pol, B_tor
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(abs(pol_mode)-1, abs(pol_mode)+1, k_z_r, I_m)
    if (status /= fgsl_success .and. log_err) then
       write (log_msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
       call log_write
    end if
    B_rad = 0.5d0 * (I_m(-1) + I_m(1)) * vac_coeff
    B_pol = imun * pol_mode / k_z_r * I_m(0) * vac_coeff
    B_tor = imun * I_m(0) * vac_coeff
  end subroutine kilca_vacuum_fourier

  subroutine check_kilca_vacuum
    use magdif_config, only: n, nflux, kilca_pol_mode, kilca_vac_coeff, longlines
    use magdif, only: equil, fs_half
    complex(dp) :: B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    real(dp) :: rad
    integer :: kf, fid, abs_pol_mode

    abs_pol_mode = abs(kilca_pol_mode)
    open(newunit = fid, file = 'cmp_vac.dat', recl = 3 * longlines)
    do kf = 1, nflux
       rad = fs_half%rad(kf)
       call kilca_vacuum_fourier(n, -abs_pol_mode, equil%rcentr, rad, &
            kilca_vac_coeff(abs_pol_mode), B_rad_neg, B_pol_neg, B_tor_neg)
       call kilca_vacuum_fourier(n, abs_pol_mode, equil%rcentr, rad, &
            kilca_vac_coeff(abs_pol_mode), B_rad_pos, B_pol_pos, B_tor_pos)
       write (fid, '(13(1x, es24.16e3))') rad, &
           B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    end do
    close(fid)
  end subroutine check_kilca_vacuum

  subroutine check_RT0
    use magdif_config, only: n, kilca_pol_mode, nflux, nkpol, longlines
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    use magdif, only: equil, fs_half, Bnflux, Bnphi, point_location, interp_RT0
    integer :: fid, kf, kpol, ktri
    real(dp) :: rad, theta, R, Z
    complex(dp) :: B_R, B_Z, B_phi, B_R_interp, B_Z_interp, B_phi_interp
    integer :: pol_modes(2)

    pol_modes = [kilca_pol_mode, -kilca_pol_mode]
    open(newunit = fid, file = 'cmp_RT0.dat', recl = longlines)
    do kf = nflux / 3, nflux / 3 ! 1, nflux
       rad = fs_half%rad(kf)
       do kpol = 1, 2 * nkpol
          theta = (kpol - 0.5d0) / dble(2 * nkpol) * 2d0 * pi
          call kilca_vacuum(n, pol_modes, equil%rcentr, rad, theta, B_R, B_phi, B_Z)
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
