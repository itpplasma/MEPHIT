module magdif_mesh_mod

  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90

  implicit none

  private

  public :: generate_mesh, write_kilca_convexfile, kilca_vacuum

contains

  subroutine generate_mesh(unprocessed_geqdsk)
    use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, mesh_element_rmp
    use magdif_config, only: log_msg, log_info, log_write, nflux, kilca_scale_factor
    use magdif_util, only: get_equil_filenames
    use magdif, only: equil, init_indices, kt_low, kp_low, cache_mesh_data

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

    call init_indices
    ntri = kt_low(nflux+1)
    npoint = kp_low(nflux+1)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))

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

    allocate(coarse_lo(nref))
    allocate(coarse_hi(nref))
    allocate(fine_lo(nref))
    allocate(fine_hi(nref))
    allocate(fine_sep(nref))
    allocate(factor(nref))
    allocate(geom_ser(maxval(additions) + 1, nref))
    nflux = nflux_unref + 2 * sum(additions - deletions)
    allocate(refined(nflux))
    refined = 0d0
    coarse_lo = floor(res * nflux_unref) - deletions;
    coarse_hi = ceiling(res * nflux_unref) + deletions;
    ! compute upper and lower array indices of refined regions
    fine_lo = coarse_lo + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 0, nref - 1)]
    fine_hi = coarse_hi + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 1, nref)]
    ! compute separations between flux surfaces using geometric series
    coarse_sep = 1d0 / dble(nflux_unref - 1)
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
    refined(1:fine_lo(1)-1) = linspace(0d0, refined(fine_lo(1)), fine_lo(1) - 1, 0, 1)
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

  subroutine radial_refinement(rho)
    use magdif_config, only: n, refinement, additions, deletions, sheet_current_factor, &
         read_delayed_config
    use magdif, only: equil
    real(dp), dimension(:), allocatable, intent(out) :: rho
    integer :: m_res_min, m_res_max
    real(dp), dimension(:), allocatable :: psi_res, rho_res
    logical, dimension(:), allocatable :: mask

    m_res_min = ceiling(minval(abs(equil%qpsi)) * dble(n))
    m_res_max = floor(maxval(abs(equil%qpsi)) * dble(n))
    allocate(psi_res(m_res_min:m_res_max))
    call compute_resonant_surfaces(m_res_min, m_res_max, psi_res)
    allocate(rho_res(size(psi_res)))
    rho_res = sqrt((psi_res - equil%simag) / (equil%sibry - equil%simag))  ! stub
    call read_delayed_config(m_res_min, m_res_max)
    mask = 0d0 < refinement .and. refinement < 1d0
    call refine_eqd_partition(count(mask), pack(deletions, mask), pack(additions, mask), &
         pack(refinement, mask), pack(rho_res, mask), rho)

    if (allocated(rho_res)) deallocate(rho_res)
    if (allocated(psi_res)) deallocate(psi_res)
    if (allocated(refinement)) deallocate(refinement)
    if (allocated(deletions)) deallocate(deletions)
    if (allocated(additions)) deallocate(additions)
    if (allocated(sheet_current_factor)) deallocate(sheet_current_factor)
  end subroutine radial_refinement

  subroutine create_kilca_mesh_points(convexfile)
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    use mesh_mod, only: knot, mesh_point
    use magdif_config, only: nflux
    use magdif_util, only: interp_psi_pol
    use magdif, only: equil, kp_low, kp_max

    character(len = *), intent(in) :: convexfile
    integer :: kf, kp
    real(dp) :: rho, rho_max, theta
    type(knot) :: node

    ! calculate maximal extent from magnetic axis
    rho_max = min(equil%rmaxis - equil%rleft, &
         equil%rleft + equil%rdim - equil%rmaxis, &
         equil%zdim * 0.5d0 + equil%zmid - equil%zmaxis, &
         equil%zdim * 0.5d0 - equil%zmid + equil%zmaxis)

    mesh_point(1)%rcoord = equil%rmaxis
    mesh_point(1)%zcoord = equil%zmaxis
    mesh_point(1)%psi_pol = interp_psi_pol(equil%rmaxis, equil%zmaxis)
    do kf = 1, nflux
       rho = dble(kf) / dble(nflux+1) * rho_max
       do kp = 1, kp_max(kf)
          theta = dble(kp-1) / dble(kp_max(kf)) * 2d0 * pi  ! [0, 2\pi)
          node%rcoord = equil%rmaxis + rho * cos(theta)
          node%zcoord = equil%zmaxis + rho * sin(theta)
          node%psi_pol = interp_psi_pol(equil%rmaxis + rho * cos(theta), &
               equil%zmaxis + rho * sin(theta))
          node%n_owners = 0
          mesh_point(kp_low(kf) + kp) = node
       end do
    end do
    call write_kilca_convexfile(rho_max, convexfile)
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
    use mesh_mod, only: npoint, mesh_point
    use magdif_config, only: nflux, nkpol
    use magdif_util, only: initialize_globals, interp_psi_pol
    use magdif, only: equil
    use field_line_integration_mod, only: theta0_at_xpoint
    use points_2d, only: s_min, create_points_2d

    integer :: kpoi
    integer, dimension(:), allocatable :: n_theta
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi

    allocate(n_theta(nflux))
    allocate(points(3, npoint))
    allocate(points_s_theta_phi(3, npoint))

    n_theta = nkpol
    s_min = 1d-16
    theta0_at_xpoint = .true.
    call initialize_globals(equil%rmaxis, equil%zmaxis)
    call create_points_2d(n_theta, points, points_s_theta_phi, r_scaling_func = sqr)
    points(:, 1) = [equil%rmaxis, 0d0, equil%zmaxis]
    mesh_point(1)%rcoord = equil%rmaxis
    mesh_point(1)%zcoord = equil%zmaxis
    mesh_point(1)%psi_pol = interp_psi_pol(equil%rmaxis, equil%zmaxis)
    do kpoi = 2, npoint
       mesh_point(kpoi)%rcoord = points(1, kpoi)
       mesh_point(kpoi)%zcoord = points(3, kpoi)
       mesh_point(kpoi)%psi_pol = interp_psi_pol(points(1, kpoi), points(3, kpoi))
    end do

    if (allocated(n_theta)) deallocate(n_theta)
    if (allocated(points)) deallocate(points)
    if (allocated(points_s_theta_phi)) deallocate(points_s_theta_phi)

  contains
    pure function sqr(x) result(x_squared)
      real(dp), dimension(:), intent(in) :: x
      real(dp), dimension(size(x)) :: x_squared
      x_squared = x * x
    end function sqr
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
    use magdif_config, only: nflux, longlines, point_file, tri_file
    use magdif, only: kp_max, kp_low, equil

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

    open(newunit = fid, file = point_file, form = 'unformatted', status = 'replace')
    write (fid) npoint
    write (fid) mesh_point
    close(fid)

    open(newunit = fid, file = tri_file, form = 'unformatted', status = 'replace')
    write (fid) ntri
    write (fid) mesh_element
    write (fid) equil%bcentr
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
