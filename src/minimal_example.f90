program minimal_example

  use magdif_config
  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
  use input_files, only: gfile, convexfile
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       bphicovar, knot, triangle, triangle_rmp, mesh_element_rmp
  use magdif_util, only: imun, interp_psi_pol
  use magdif, only: init_indices, kt_max, kt_low, kp_max, kp_low, Bnflux, Bnphi, &
       get_labeled_edges, check_redundant_edges, check_div_free, write_vector_dof, &
       cache_mesh_data

  implicit none

  interface interleave
     procedure interleave_vv, interleave_vs, interleave_sv, interleave_ss
  end interface interleave

  character(len = 1024) :: unscaled_geqdsk
  real(dp) :: rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr

  integer :: fid, kf, kp, k, ke, ktri, kpoi
  integer, dimension(:), allocatable :: kq, kt, tri_f, tri_i, tri_o
  integer, dimension(:,:), allocatable :: quad
  real(dp) :: rho, rho_max, theta
  integer, parameter :: nrz = 64  ! at most 100 values are read in by field_divB0.f90

  type(knot) :: node
  type(triangle) :: elem
  type(triangle_rmp) :: tri

  real(dp) :: r, z, n_r, n_z
  complex(dp) :: Br, Bp, Bz

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
     if (command_argument_count() >= 2) then
        call get_command_argument(2, unscaled_geqdsk)
     else
        error stop 'expected path to unscaled G EQDSK file as second parameter'
     end if
  else
     error stop 'expected path to magdif config file as first parameter'
  end if
  call read_config

  log_file = '-'
  call log_open

  call initialize_globals
  call scale_geqdsk(kilca_scale_factor, trim(unscaled_geqdsk), trim(gfile), &
       rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr)
  ! convert from SI to Gaussian units
  rdim = rdim * 1d2
  zdim = zdim * 1d2
  rleft = rleft * 1d2
  zmid = zmid * 1d2
  rmaxis = rmaxis * 1d2
  zmaxis = zmaxis * 1d2
  bphicovar = bcentr * 1d5

  ! write configuration to FreeFem++ (for Delaunay triangulation)
  open(newunit = fid, file = 'extmesh.idp', status = 'replace')
  write (fid, '("real rmaxis = ", es23.16, ";")') rmaxis
  write (fid, '("real zmaxis = ", es23.16, ";")') zmaxis
  write (fid, '("real rdim = ", es23.16, ";")') rdim
  write (fid, '("real zdim = ", es23.16, ";")') zdim
  write (fid, '("int nkpol = ", i0, ";")') nkpol
  write (fid, '("int nflux = ", i0, ";")') nflux
  close(fid)

  ! calculate maximal extent from magnetic axis
  rho_max = min(rmaxis - rleft, rleft + rdim - rmaxis, &
       zdim * 0.5d0 + zmid - zmaxis, zdim * 0.5d0 - zmid + zmaxis)

  call init_indices
  ntri = kt_low(nflux+1)
  npoint = kp_low(nflux+1)
  allocate (mesh_element(ntri))
  allocate (mesh_point(npoint))

  mesh_point(1)%rcoord = rmaxis
  mesh_point(1)%zcoord = zmaxis
  mesh_point(1)%psi_pol = interp_psi_pol(rmaxis, zmaxis)
  do kf = 1, nflux
     rho = dble(kf) / dble(nflux+1) * rho_max
     do kp = 1, kp_max(kf)
        theta = dble(kp-1) / dble(kp_max(kf)) * 2d0 * pi  ! [0, 2\pi)
        node%rcoord = rmaxis + rho * cos(theta)
        node%zcoord = zmaxis + rho * sin(theta)
        node%psi_pol = &
             interp_psi_pol(rmaxis + rho * cos(theta), zmaxis + rho * sin(theta))
        node%n_owners = 0
        mesh_point(kp_low(kf) + kp) = node
     end do
  end do

  ! convexfile is read in from field_divB0.inp in initialize_globals
  open(newunit = fid, file = convexfile, status = 'replace')
  rho = rho_max
  do kp = 1, nrz
     theta = dble(kp) / dble(nrz) * 2d0 * pi
     write (fid, '(2(1x, es23.16))') rmaxis + rho * cos(theta), zmaxis + rho * sin(theta)
  end do
  close(fid)

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
  write (fid) bphicovar
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

  allocate(mesh_element_rmp(ntri))
  call cache_mesh_data

  ! calculate resonant vacuum perturbation
  allocate(Bnflux(ntri, 3))
  allocate(Bnphi(ntri))
  do ktri = 1, ntri
        tri = mesh_element_rmp(ktri)
        ! flux through edge f
        n_r = mesh_point(tri%lf(2))%zcoord - mesh_point(tri%lf(1))%zcoord
        n_z = mesh_point(tri%lf(1))%rcoord - mesh_point(tri%lf(2))%rcoord
        r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
        z = sum(mesh_point(tri%lf(:))%zcoord) * 0.5d0
        rho = hypot(r - rmaxis, z - zmaxis)
        theta = atan2(z - zmaxis, r - rmaxis)
        call kilca_vacuum(n, kilca_pol_mode, R0, rho, theta, Br, Bp, Bz)
        Bnflux(ktri, tri%ef) = (Br * n_r + Bz * n_z) * r
        ! flux through edge i
        n_r = mesh_point(tri%li(2))%zcoord - mesh_point(tri%li(1))%zcoord
        n_z = mesh_point(tri%li(1))%rcoord - mesh_point(tri%li(2))%rcoord
        r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
        z = sum(mesh_point(tri%li(:))%zcoord) * 0.5d0
        rho = hypot(r - rmaxis, z - zmaxis)
        theta = atan2(z - zmaxis, r - rmaxis)
        call kilca_vacuum(n, kilca_pol_mode, R0, rho, theta, Br, Bp, Bz)
        Bnflux(ktri, tri%ei) = (Br * n_r + Bz * n_z) * r
        ! flux through edge o
        n_r = mesh_point(tri%lo(2))%zcoord - mesh_point(tri%lo(1))%zcoord
        n_z = mesh_point(tri%lo(1))%rcoord - mesh_point(tri%lo(2))%rcoord
        r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
        z = sum(mesh_point(tri%lo(:))%zcoord) * 0.5d0
        rho = hypot(r - rmaxis, z - zmaxis)
        theta = atan2(z - zmaxis, r - rmaxis)
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
  if (allocated(mesh_element_rmp)) deallocate(mesh_element)
  if (allocated(mesh_element)) deallocate(mesh_element)
  if (allocated(mesh_point)) deallocate(mesh_point)

  call log_close

contains

  ! better future solution: put this in a separate subroutine in field_divB0.f90
  subroutine initialize_globals
    integer :: fid
    open(newunit = fid, file = 'field_divB0.inp', status = 'old')
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *) gfile        ! equilibrium file
    read (fid, *)
    read (fid, *) convexfile   ! convex file for stretchcoords
    close(fid)
  end subroutine initialize_globals

  subroutine scale_geqdsk(gamma, file_in, file_out, &
       rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr)
    integer, intent(in) :: gamma
    character(len = *), intent(in) :: file_in, file_out
    real(dp), intent(out) :: rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr

    character(len = *), parameter :: geqdsk_2000 = '(6a8,3i4)'
    character(len = *), parameter :: geqdsk_2020 = '(5e16.9)'
    character(len = *), parameter :: geqdsk_2022 = '(2i5)'

    character(len = 10) :: text(6)
    integer :: fid_in, fid_out, i, j, idum, nw, nh, nbbbs, limitr
    real(dp) :: xdum, rcentr, simag, sibry, current, unscaled_rcentr, r_shift
    real(dp), allocatable :: fpol(:), pres(:), ffprim(:), pprime(:), psirz(:, :), &
         qpsi(:), rbbbs(:), zbbbs(:), rlim(:), zlim(:)

    write (log_msg, '("attempting to read unscaled G EQDSK file ", a)') file_in
    call log_write
    open(newunit = fid_in, file = file_in, status = 'old')
    read (fid_in, geqdsk_2000) (text(i), i = 1, 6), idum, nw, nh
    allocate(fpol(nw))
    allocate(pres(nw))
    allocate(ffprim(nw))
    allocate(pprime(nw))
    allocate(psirz(nw, nh))
    allocate(qpsi(nw))
    read (fid_in, geqdsk_2020) rdim, zdim, rcentr, rleft, zmid
    read (fid_in, geqdsk_2020) rmaxis, zmaxis, simag, sibry, bcentr
    read (fid_in, geqdsk_2020) current, simag, xdum, rmaxis, xdum
    read (fid_in, geqdsk_2020) zmaxis, xdum, sibry, xdum, xdum
    read (fid_in, geqdsk_2020) (fpol(i), i = 1, nw)
    read (fid_in, geqdsk_2020) (pres(i), i = 1, nw)
    read (fid_in, geqdsk_2020) (ffprim(i), i = 1, nw)
    read (fid_in, geqdsk_2020) (pprime(i), i = 1, nw)
    read (fid_in, geqdsk_2020) ((psirz(i, j), i = 1, nw), j = 1, nh)
    read (fid_in, geqdsk_2020) (qpsi(i), i = 1, nw)
    read (fid_in, geqdsk_2022) nbbbs, limitr
    allocate(rbbbs(nbbbs))
    allocate(zbbbs(nbbbs))
    allocate(rlim(limitr))
    allocate(zlim(limitr))
    read (fid_in, geqdsk_2020) (rbbbs(i), zbbbs(i), i = 1, nbbbs)
    read (fid_in, geqdsk_2020) (rlim(i), zlim(i), i = 1, limitr)
    close(fid_in)

    unscaled_rcentr = rcentr
    r_shift = rcentr * (gamma - 1)
    rcentr = rcentr * gamma
    rleft = rcentr - 0.5d0 * rdim
    rmaxis = rmaxis + r_shift
    simag = simag * gamma
    sibry = sibry * gamma
    fpol = fpol * gamma
    ffprim = ffprim * gamma
    pprime = pprime / gamma
    psirz = psirz * gamma
    qpsi = qpsi / gamma
    rbbbs = rbbbs + r_shift
    rlim = rlim + r_shift

    write (log_msg, '("attempting to write scaled G EQDSK file ", a)') file_out
    call log_write
    open(newunit = fid_out, file = file_out, status = 'replace')
    write (fid_out, geqdsk_2000) (text(i), i = 1, 6), idum, nw, nh
    write (fid_out, geqdsk_2020) rdim, zdim, rcentr, rleft, zmid
    write (fid_out, geqdsk_2020) rmaxis, zmaxis, simag, sibry, bcentr
    write (fid_out, geqdsk_2020) current, simag, xdum, rmaxis, xdum
    write (fid_out, geqdsk_2020) zmaxis, xdum, sibry, xdum, xdum
    write (fid_out, geqdsk_2020) (fpol(i), i = 1, nw)
    write (fid_out, geqdsk_2020) (pres(i), i = 1, nw)
    write (fid_out, geqdsk_2020) (ffprim(i), i = 1, nw)
    write (fid_out, geqdsk_2020) (pprime(i), i = 1, nw)
    write (fid_out, geqdsk_2020) ((psirz(i, j), i = 1, nw), j = 1, nh)
    write (fid_out, geqdsk_2020) (qpsi(i), i = 1, nw)
    write (fid_out, geqdsk_2022) nbbbs, limitr
    write (fid_out, geqdsk_2020) (rbbbs(i), zbbbs(i), i = 1, nbbbs)
    write (fid_out, geqdsk_2020) (rlim(i), zlim(i), i = 1, limitr)
    close(fid_out)

    if (allocated(fpol)) deallocate(fpol)
    if (allocated(pres)) deallocate(pres)
    if (allocated(ffprim)) deallocate(ffprim)
    if (allocated(pprime)) deallocate(pprime)
    if (allocated(psirz)) deallocate(psirz)
    if (allocated(qpsi)) deallocate(qpsi)
    if (allocated(rbbbs)) deallocate(rbbbs)
    if (allocated(zbbbs)) deallocate(zbbbs)
    if (allocated(rlim)) deallocate(rlim)
    if (allocated(zlim)) deallocate(zlim)
  end subroutine scale_geqdsk

  subroutine add_node_owner(kpoint, ktri)
    use mesh_mod, only: n_owners_max
    integer, intent(in) :: kpoint, ktri
    if (mesh_point(kpoint)%n_owners < n_owners_max) then
       mesh_point(kpoint)%i_owner_tri(mesh_point(kpoint)%n_owners + 1) = ktri
       mesh_point(kpoint)%n_owners = mesh_point(kpoint)%n_owners + 1
    else
       write (log_msg, '("Maximal number of owning triangles exceeded at point ", i0)') &
            kpoint
       if (log_warn) call log_write
    end if
  end subroutine add_node_owner

  elemental subroutine calculate_det_3(elem)
    type(triangle), intent(inout) :: elem
    real(dp) :: e1_r, e1_z, e2_r, e2_z
    e1_r = mesh_point(elem%i_knot(1))%rcoord - mesh_point(elem%i_knot(3))%rcoord
    e1_z = mesh_point(elem%i_knot(1))%zcoord - mesh_point(elem%i_knot(3))%zcoord
    e2_r = mesh_point(elem%i_knot(2))%rcoord - mesh_point(elem%i_knot(3))%rcoord
    e2_z = mesh_point(elem%i_knot(2))%zcoord - mesh_point(elem%i_knot(3))%zcoord
    elem%det_3 = abs(e1_r * e2_z - e1_z * e2_r)
  end subroutine calculate_det_3

  pure function interleave_vv(first_v, second_v, num) result(merged)
    integer, intent(in) :: first_v(:), second_v(:), num
    integer :: merged(num)
    integer :: k
    merged = merge(first_v, second_v, [([.true., .false.], k = 1, num / 2)])
  end function interleave_vv

  pure function interleave_vs(first_v, second_s, num) result(merged)
    integer, intent(in) :: first_v(:), second_s, num
    integer :: merged(num)
    integer :: k
    merged = merge(first_v, [(second_s, k = 1, num)], &
         [([.true., .false.], k = 1, num / 2)])
  end function interleave_vs

  pure function interleave_sv(first_s, second_v, num) result(merged)
    integer, intent(in) :: first_s, second_v(:), num
    integer :: merged(num)
    integer :: k
    merged = merge([(first_s, k = 1, num)], second_v, &
         [([.true., .false.], k = 1, num / 2)])
  end function interleave_sv

  pure function interleave_ss(first_s, second_s, num) result(merged)
    integer, intent(in) :: first_s, second_s, num
    integer :: merged(num)
    integer :: k
    merged = [([first_s, second_s], k = 1, num / 2)]
  end function interleave_ss

  subroutine kilca_vacuum(tor_mode, pol_mode, R_0, r, theta, Br, Bp, Bz)
    use fgsl, only: fgsl_double, fgsl_int, fgsl_success, fgsl_sf_bessel_icn_array
    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r, theta
    complex(dp), intent(out) :: Br, Bp, Bz
    complex(dp) :: B_r, B_theta, B_z
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(pol_mode-1, pol_mode+1, k_z_r, I_m)
    if (status /= fgsl_success .and. log_err) then
       write (log_msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
       call log_write
    end if
    B_r = (0.5d0, 0d0) * (I_m(-1) + I_m(1)) * cos(pol_mode * theta) &
         * kilca_vac_coeff
    B_theta = imun * pol_mode / k_z_r * I_m(0) * cos(pol_mode * theta) &
         * kilca_vac_coeff
    B_z = imun * I_m(0) * cos(pol_mode * theta) * kilca_vac_coeff
    Br = B_r * cos(theta) - B_theta * sin(theta)
    Bp = B_z
    Bz = B_r * sin(theta) + B_theta * cos(theta)
  end subroutine kilca_vacuum
end program minimal_example