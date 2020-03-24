module magdif_mesh_mod

  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90

  implicit none

  private

  public :: generate_mesh, write_kilca_convexfile, kilca_vacuum

  interface
     function mapping(val)
       import :: dp
       real(dp), intent(in) :: val
       real(dp) :: mapping
     end function mapping
  end interface

contains

  subroutine generate_mesh(unprocessed_geqdsk)
    use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp
    use magdif_config, only: log_msg, log_info, log_write, kilca_scale_factor
    use magdif_util, only: get_equil_filenames, initialize_globals
    use magdif, only: equil, cache_mesh_data

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

    if (kilca_scale_factor /= 0) then
       call create_kilca_mesh_points(convexfile)
    else
       call create_mesh_points
    end if
    call connect_mesh_points
    call write_mesh_data
    call cache_mesh_data
    if (kilca_scale_factor /= 0) then
       call compute_kilca_vacuum
    end if

    if (allocated(mesh_element_rmp)) deallocate(mesh_element)
    if (allocated(mesh_element)) deallocate(mesh_element)
    if (allocated(mesh_point)) deallocate(mesh_point)
  end subroutine generate_mesh

  subroutine compute_resonant_surfaces(m_res_min, m_res_max, psi_res)
    use magdif_config, only: n
    use magdif_util, only: linspace, flux_func
    use magdif, only: equil
    use netlib_mod, only: zeroin
    integer, intent(in) :: m_res_min, m_res_max
    real(dp), intent(out) :: psi_res(m_res_min:m_res_max)
    integer :: m
    type(flux_func) :: psi_eqd

    call psi_eqd%init(4, linspace(equil%simag, equil%sibry, equil%nw, 0, 0))
    do m = m_res_min, m_res_max
       psi_res(m) = zeroin(equil%simag, equil%sibry, q_interp_resonant, 1d-9)
    end do
  contains
    function q_interp_resonant(psi)
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant
      q_interp_resonant = psi_eqd%interp(abs(equil%qpsi) - dble(m) / dble(n), psi)
    end function q_interp_resonant
  end subroutine compute_resonant_surfaces

  subroutine refine_eqd_partition(nref, deletions, additions, refinement, res, refined)
    use magdif_config, only: nflux_unref, nflux
    use magdif_util, only: linspace
    integer, intent(in) :: nref
    integer, dimension(nref), intent(in) :: deletions, additions
    real(dp), dimension(nref), intent(in) :: res, refinement
    real(dp), dimension(:), allocatable, intent(out) :: refined
    integer :: kref, k
    integer, dimension(:), allocatable :: coarse_lo, coarse_hi, fine_lo, fine_hi
    real(dp) :: coarse_sep
    real(dp), dimension(:), allocatable :: fine_sep, factor
    real(dp), dimension(:, :), allocatable :: geom_ser

    if (nref < 1) then
       nflux = nflux_unref
       allocate(refined(0:nflux))
       refined = linspace(0d0, 1d0, nflux + 1, 0, 0)
       return
    end if
    allocate(coarse_lo(nref))
    allocate(coarse_hi(nref))
    allocate(fine_lo(nref))
    allocate(fine_hi(nref))
    allocate(fine_sep(nref))
    allocate(factor(nref))
    allocate(geom_ser(maxval(additions) + 1, nref))
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
    if (allocated(coarse_lo)) deallocate(coarse_lo)
    if (allocated(coarse_hi)) deallocate(coarse_hi)
    if (allocated(fine_lo)) deallocate(fine_lo)
    if (allocated(fine_hi)) deallocate(fine_hi)
    if (allocated(fine_sep)) deallocate(fine_sep)
    if (allocated(factor)) deallocate(factor)
    if (allocated(geom_ser)) deallocate(geom_ser)
  end subroutine refine_eqd_partition

  subroutine radial_refinement(rho_norm_ref, psi2rho_norm)
    use magdif_config, only: n, refinement, additions, deletions, sheet_current_factor, &
         read_delayed_config
    use magdif, only: equil
    real(dp), dimension(:), allocatable, intent(out) :: rho_norm_ref
    procedure(mapping), optional :: psi2rho_norm
    integer :: m_res_min, m_res_max, m
    real(dp), dimension(:), allocatable :: psi_res, rho_res
    logical, dimension(:), allocatable :: mask

    m_res_min = ceiling(minval(abs(equil%qpsi)) * dble(n))
    m_res_max = floor(maxval(abs(equil%qpsi)) * dble(n))
    allocate(psi_res(m_res_min:m_res_max))
    call compute_resonant_surfaces(m_res_min, m_res_max, psi_res)
    allocate(rho_res(size(psi_res)))
    if (present(psi2rho_norm)) then
       rho_res = [(psi2rho_norm(psi_res(m)), m = m_res_min, m_res_max)]
    else
       rho_res = sqrt((psi_res - equil%simag) / (equil%sibry - equil%simag))  ! stub
    end if
    call read_delayed_config(m_res_min, m_res_max)
    mask = 0d0 < refinement .and. refinement < 1d0
    call refine_eqd_partition(count(mask), pack(deletions, mask), pack(additions, mask), &
         pack(refinement, mask), pack(rho_res, mask), rho_norm_ref)

    if (allocated(rho_res)) deallocate(rho_res)
    if (allocated(psi_res)) deallocate(psi_res)
    if (allocated(refinement)) deallocate(refinement)
    if (allocated(deletions)) deallocate(deletions)
    if (allocated(additions)) deallocate(additions)
    if (allocated(sheet_current_factor)) deallocate(sheet_current_factor)
  end subroutine radial_refinement

  subroutine create_kilca_mesh_points(convexfile)
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    use mesh_mod, only: npoint, mesh_point, ntri, mesh_element, mesh_element_rmp
    use magdif_config, only: nflux
    use magdif_util, only: interp_psi_pol, linspace, flux_func
    use magdif, only: equil, kp_low, kp_max, kt_low, init_indices, fs, fs_half

    character(len = *), intent(in) :: convexfile
    integer :: fid, k, nlabel, kf, kp, kpoi
    real(dp) :: rho_max, theta
    real(dp), dimension(:), allocatable :: r_eqd, rho_norm_eqd, sample_psi, rho_ref, &
         rho_half
    type(flux_func) :: psi_interpolator

    ! calculate maximal extent from magnetic axis
    rho_max = min(equil%rmaxis - equil%rleft, &
         equil%rleft + equil%rdim - equil%rmaxis, &
         equil%zdim * 0.5d0 + equil%zmid - equil%zmaxis, &
         equil%zdim * 0.5d0 - equil%zmid + equil%zmaxis)

    ! interpolate between rho and psi
    open(newunit = fid, file = 'preload_for_SYNCH.inp', status = 'old')
    read (fid, *)
    read (fid, *) nlabel
    close(fid)
    allocate(r_eqd(nlabel))
    allocate(sample_psi(nlabel))
    allocate(rho_norm_eqd(nlabel))
    r_eqd = linspace(equil%rmaxis, equil%rmaxis + rho_max, nlabel, 0, 0)
    sample_psi = [(interp_psi_pol(r_eqd(k), equil%zmaxis), k = 1, nlabel)]
    rho_norm_eqd = linspace(0d0, 1d0, nlabel, 0, 0)
    call psi_interpolator%init(4, sample_psi)

    call radial_refinement(rho_ref, psi2rho_norm)
    rho_ref = rho_ref * rho_max

    call fs%init(nflux, .false.)
    call fs_half%init(nflux, .true.)
    fs%psi = [(interp_psi_pol(equil%rmaxis + rho_ref(kf), equil%zmaxis), kf = 0, nflux)]
    allocate(rho_half(nflux))
    rho_half = 0.5d0 * (rho_ref(0:nflux-1) + rho_ref(1:nflux))
    fs_half%psi = [(interp_psi_pol(equil%rmaxis + rho_half(kf), equil%zmaxis), &
         kf = 1, nflux)]
    if (allocated(rho_half)) deallocate(rho_half)

    call init_indices
    ntri = kt_low(nflux+1)
    npoint = kp_low(nflux+1)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))

    mesh_point(1)%rcoord = equil%rmaxis
    mesh_point(1)%zcoord = equil%zmaxis
    do kf = 1, nflux
       do kp = 1, kp_max(kf)
          kpoi = kp_low(kf) + kp
          theta = dble(kp-1) / dble(kp_max(kf)) * 2d0 * pi  ! [0, 2\pi)
          mesh_point(kpoi)%rcoord = equil%rmaxis + rho_ref(kf) * cos(theta)
          mesh_point(kpoi)%zcoord = equil%zmaxis + rho_ref(kf) * sin(theta)
       end do
    end do
    call write_kilca_convexfile(rho_max, convexfile)
    if (allocated(rho_ref)) deallocate(rho_ref)
    if (allocated(rho_norm_eqd)) deallocate(rho_norm_eqd)
    if (allocated(sample_psi)) deallocate(sample_psi)
    if (allocated(r_eqd)) deallocate(r_eqd)

  contains
    function psi2rho_norm(psi) result(rho_norm)
      real(dp), intent(in) :: psi
      real(dp) :: rho_norm
      rho_norm = psi_interpolator%interp(rho_norm_eqd, psi)
    end function psi2rho_norm
  end subroutine create_kilca_mesh_points

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
       write (fid, '(2(1x, es23.16))') equil%rmaxis + rho_max * cos(theta), &
            equil%zmaxis + rho_max * sin(theta)
    end do
    close(fid)
  end subroutine write_kilca_convexfile

  subroutine create_mesh_points
    use mesh_mod, only: npoint, mesh_point, ntri, mesh_element, mesh_element_rmp
    use magdif_config, only: nflux, nkpol
    use magdif_util, only: interp_psi_pol, flux_func
    use magdif, only: kp_low, kt_low, init_indices, fs, fs_half
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, raxis, zaxis
    use field_line_integration_mod, only: theta0_at_xpoint, theta_axis
    use points_2d, only: s_min, create_points_2d

    integer :: kf, kpoi
    integer, dimension(:), allocatable :: n_theta
    real(dp), dimension(:), allocatable :: rho_norm_eqd, rho_ref, rho_half
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi
    type(flux_func) :: psi_interpolator
    real(dp) :: psi_axis

    ! calculates points on a fine grid in the core region by integrating along field lines
    call preload_for_SYNCH
    ! loads points that are calculated in preload_for_SYNCH into module variables
    ! and spline interpolates them for use with magdata_in_symfluxcoord_ext
    call load_magdata_in_symfluxcoord
    ! interpolate between rho and psi
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd = rbeg / hypot(theta_axis(1), theta_axis(2))
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    call psi_interpolator%init(4, psisurf * psipol_max + psi_axis)

    call radial_refinement(rho_ref, psi2rho_norm)

    call fs%init(nflux, .false.)
    call fs_half%init(nflux, .true.)
    fs%psi = [(interp_psi_pol(raxis + rho_ref(kf) * theta_axis(1), &
         zaxis + rho_ref(kf) * theta_axis(2)), kf = 0, nflux)]
    allocate(rho_half(nflux))
    rho_half = 0.5d0 * (rho_ref(0:nflux-1) + rho_ref(1:nflux))
    fs_half%psi = [(interp_psi_pol(raxis + rho_half(kf) * theta_axis(1), &
         zaxis + rho_half(kf) * theta_axis(2)), kf = 1, nflux)]
    if (allocated(rho_half)) deallocate(rho_half)

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
    theta0_at_xpoint = .true.
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
    if (allocated(rho_ref)) deallocate(rho_ref)
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
      psi_ref = [(interp_psi_pol(raxis + rho_ref(kf) * theta_axis(1), &
           zaxis + rho_ref(kf) * theta_axis(2)) - psi_axis, kf = 1, nflux)]
    end function psi_ref
  end subroutine create_mesh_points

  subroutine connect_mesh_points
    use mesh_mod, only: mesh_element
    use magdif_config, only: nflux
    use magdif_util, only: interleave, calculate_det_3
    use magdif, only: kt_max, kt_low, kp_low

    integer :: kf, k
    integer, dimension(:), allocatable :: kq, kt, tri_f, tri_i, tri_o
    integer, dimension(:,:), allocatable :: quad

    allocate (kt(maxval(kt_max)))
    kt = [(k, k = 1, maxval(kt_max))]
    allocate (tri_f(maxval(kt_max)))
    allocate (tri_i(maxval(kt_max)))
    allocate (tri_o(maxval(kt_max)))
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
    allocate (kq(maxval(kt_max)))
    kq = (kt + 1) / 2
    allocate (quad(4, maxval(kt_max)))
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
    if (allocated(kt)) deallocate(kt)
    if (allocated(tri_f)) deallocate(tri_f)
    if (allocated(tri_i)) deallocate(tri_i)
    if (allocated(tri_o)) deallocate(tri_o)
    if (allocated(kq)) deallocate(kq)
    if (allocated(quad)) deallocate(quad)
  end subroutine connect_mesh_points

  subroutine write_mesh_data
    use mesh_mod, only: npoint, ntri, knot, triangle, mesh_point, mesh_element
    use magdif_config, only: nflux, longlines, meshdata_file
    use magdif, only: kp_max, kp_low, fs, fs_half

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

    open(newunit = fid, file = meshdata_file, form = 'unformatted', status = 'replace')
    write (fid) nflux, npoint, ntri
    write (fid) fs%psi
    write (fid) fs_half%psi
    write (fid) mesh_point
    write (fid) mesh_element
    close(fid)

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace')
    write (fid, '(3(1x, i0))') npoint, ntri, kp_max(nflux+1) - 1
    do kpoi = 1, kp_low(nflux+1)
       write (fid, '(2(1x, es23.16), 1x, i0)') &
            mesh_point(kpoi)%rcoord, mesh_point(kpoi)%zcoord, 0
    end do
    do kpoi = kp_low(nflux+1) + 1, npoint
       write (fid, '(2(1x, es23.16), 1x, i0)') &
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
    use magdif_config, only: n, R0, kilca_scale_factor, kilca_pol_mode, Bn_vac_file
    use magdif_util, only: imun
    use magdif, only: equil, Bnflux, Bnphi, check_redundant_edges, check_div_free, &
         write_vector_dof

    integer :: ktri
    real(dp) :: r, z, n_r, n_z, rho, theta
    complex(dp) :: Br, Bp, Bz
    type(triangle_rmp) :: tri

    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    do ktri = 1, ntri
       tri = mesh_element_rmp(ktri)
       ! flux through edge f
       n_r = mesh_point(tri%lf(2))%zcoord - mesh_point(tri%lf(1))%zcoord
       n_z = mesh_point(tri%lf(1))%rcoord - mesh_point(tri%lf(2))%rcoord
       r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
       z = sum(mesh_point(tri%lf(:))%zcoord) * 0.5d0
       rho = hypot(r - equil%rmaxis, z - equil%zmaxis)
       theta = atan2(z - equil%zmaxis, r - equil%rmaxis)
       call kilca_vacuum(n, kilca_pol_mode, R0, rho, theta, Br, Bp, Bz)
       Bnflux(ktri, tri%ef) = (Br * n_r + Bz * n_z) * r
       ! flux through edge i
       n_r = mesh_point(tri%li(2))%zcoord - mesh_point(tri%li(1))%zcoord
       n_z = mesh_point(tri%li(1))%rcoord - mesh_point(tri%li(2))%rcoord
       r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
       z = sum(mesh_point(tri%li(:))%zcoord) * 0.5d0
       rho = hypot(r - equil%rmaxis, z - equil%zmaxis)
       theta = atan2(z - equil%zmaxis, r - equil%rmaxis)
       call kilca_vacuum(n, kilca_pol_mode, R0, rho, theta, Br, Bp, Bz)
       Bnflux(ktri, tri%ei) = (Br * n_r + Bz * n_z) * r
       ! flux through edge o
       n_r = mesh_point(tri%lo(2))%zcoord - mesh_point(tri%lo(1))%zcoord
       n_z = mesh_point(tri%lo(1))%rcoord - mesh_point(tri%lo(2))%rcoord
       r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
       z = sum(mesh_point(tri%lo(:))%zcoord) * 0.5d0
       rho = hypot(r - equil%rmaxis, z - equil%zmaxis)
       theta = atan2(z - equil%zmaxis, r - equil%rmaxis)
       call kilca_vacuum(n, kilca_pol_mode, R0, rho, theta, Br, Bp, Bz)
       Bnflux(ktri, tri%eo) = (Br * n_r + Bz * n_z) * r
       ! toroidal flux
       Bnphi(ktri) = imun / n * sum(Bnflux(ktri, :)) / tri%area
    end do
    call check_redundant_edges(Bnflux, .false., 'non-resonant B_n')
    call check_div_free(Bnflux, Bnphi, n, 1d-9, 'non-resonant B_n')
    Bnflux = Bnflux * kilca_scale_factor
    call write_vector_dof(Bnflux, Bnphi, Bn_vac_file)
    if (allocated(Bnflux)) deallocate(Bnflux)
    if (allocated(Bnphi)) deallocate(Bnphi)
  end subroutine compute_kilca_vacuum

  subroutine kilca_vacuum(tor_mode, pol_mode, R_0, r, theta, Br, Bp, Bz)
    use fgsl, only: fgsl_double, fgsl_int, fgsl_success, fgsl_sf_bessel_icn_array
    use magdif_config, only: log_msg, log_err, log_write, kilca_vac_coeff
    use magdif_util, only: imun

    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r, theta
    complex(dp), intent(out) :: Br, Bp, Bz
    complex(dp) :: B_r, B_theta, B_z
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(abs(pol_mode)-1, abs(pol_mode)+1, k_z_r, I_m)
    if (status /= fgsl_success .and. log_err) then
       write (log_msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
       call log_write
    end if
    B_r = 0.5d0 * (I_m(-1) + I_m(1)) * cos(pol_mode * theta) * kilca_vac_coeff
    B_theta = -pol_mode / k_z_r * I_m(0) * sin(pol_mode * theta) * kilca_vac_coeff
    B_z = imun * I_m(0) * cos(pol_mode * theta) * kilca_vac_coeff
    Br = B_r * cos(theta) - B_theta * sin(theta)
    Bp = B_z
    Bz = B_r * sin(theta) + B_theta * cos(theta)
  end subroutine kilca_vacuum

end module magdif_mesh_mod
