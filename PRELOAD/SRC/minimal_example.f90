program minimal_example

  use magdif_config
  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
  use input_files, only: gfile, convexfile
  use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, & ! PRELOAD/SRC/mesh_mod.f90
       bphicovar, knot, triangle
  use magdif, only: init_indices, kt_max, kt_low, kp_max, kp_low

  implicit none

  interface interleave
     procedure interleave_vv, interleave_vs, interleave_sv, interleave_ss
  end interface interleave

  character(len = 1024) :: config_file
  real(dp) :: rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr

  integer :: fid, kf, kp, k
  integer, dimension(:), allocatable :: kq, kt, tri_f, tri_i, tri_o
  integer, dimension(:,:), allocatable :: quad
  real(dp) :: rho, rho_max, theta
  integer, parameter :: nrz = 64  ! at most 100 values are read in by field_divB0.f90

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     stop 'Error: expected path to magdif config file as first parameter'
  endif
  call read_config(config_file)

  call initialize_globals
  call read_geqdsk_header(gfile, rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr)
  ! convert from SI to Gaussian units
  rdim = rdim * 1d2
  zdim = zdim * 1d2
  rleft = rleft * 1d2
  zmid = zmid * 1d2
  rmaxis = rmaxis * 1d2
  zmaxis = zmaxis * 1d2
  bphicovar = bcentr * 1d5

  ! calculate maximal extent from magnetic axis
  rho_max = min(rmaxis - rleft, rleft + rdim - rmaxis, &
       zdim * 0.5d0 + zmid - zmaxis, zdim * 0.5d0 - zmid + zmaxis)

  call init_indices
  ntri = kt_low(nflux+1) + kt_max(nflux+1)
  npoint = kp_low(nflux+1) + kp_max(nflux+1)
  allocate (mesh_element(ntri))
  allocate (mesh_point(npoint))

  mesh_point(1)%rcoord = rmaxis
  mesh_point(1)%zcoord = zmaxis
  mesh_point(1)%psi_pol = interp_psi_pol(rmaxis, zmaxis)
  do kf = 1, nflux+1
     rho = dble(kf) / dble(nflux+2) * rho_max
     do kp = 1, kp_max(kf)
        theta = dble(kp-1) / dble(kp_max(kf)) * 2d0 * pi  ! [0, 2\pi)
        mesh_point(kp_low(kf) + kp)%rcoord = rmaxis + rho * cos(theta)
        mesh_point(kp_low(kf) + kp)%zcoord = zmaxis + rho * sin(theta)
        mesh_point(kp_low(kf) + kp)%psi_pol = &
             interp_psi_pol(rmaxis + rho * cos(theta), zmaxis + rho * sin(theta))
        mesh_point(kp_low(kf) + kp)%n_owners = 0
     end do
  end do

  ! convexfile is read in from field_divB0.inp by first call to subroutine field
  ! this is done in initialize_globals
  open(newunit = fid, file = convexfile)
  rho = rho_max
  do kp = 1, nrz
     theta = dble(kp) / dble(nrz) * 2d0 * pi
     write (fid, *) rmaxis + rho * cos(theta), zmaxis + rho * sin(theta)
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
  do kf = 2, nflux+1
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
     elseif (kf == nflux+1) then
        tri_f = interleave(0, kt_low(kf-1) + kt - 1, kt_max(kf))
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

!!$  open(newunit = fid, file = 'mesh_new.asc', recl = longlines)
!!$  do k = 1, ntri
!!$     elem = mesh_element(k)
!!$     write (fid, '(i5,1x,i5,1x,i5,1x,i1,1x,i5,1x,i5,1x,i5,1x,i3,1x,i3,1x,i3)') &
!!$          (elem%i_knot(ke), ke = 1, 3), elem%knot_h, &
!!$          (elem%neighbour(ke), ke = 1, 3), (elem%neighbour_edge(ke), ke = 1, 3)
!!$  end do
!!$  close(fid)

  open(newunit = fid, file = point_file, form = 'unformatted')
  write (fid) npoint
  write (fid) mesh_point
  close(fid)

  open(newunit = fid, file = tri_file, form = 'unformatted')
  write (fid) ntri
  write (fid) mesh_element
  write (fid) bphicovar
  close(fid)

  open(newunit = fid, file = 'inputformaxwell.msh')
  write (fid, *) npoint, ntri, kp_max(nflux+1) - 1
  do k = 1, kp_low(nflux)
     write (fid, *) mesh_point(k)%rcoord, mesh_point(k)%zcoord, 0
  end do
  do k = kp_low(nflux) + 1, npoint
     write (fid, *) mesh_point(k)%rcoord, mesh_point(k)%zcoord, 1
  end do
  do k = 1, ntri
     write (fid, *) mesh_element(k)%i_knot(:), 0
  end do
  do kp = 1, kp_max(nflux+1) - 1
     write (fid, *) kp_low(nflux+1) + kp, kp_low(nflux+1) + kp + 1, 1
  end do
  close(fid)

  if (allocated(mesh_element)) deallocate(mesh_element)
  if (allocated(mesh_point)) deallocate(mesh_point)

contains

  subroutine initialize_globals
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    r = rmaxis
    z = zmaxis
    p = 0d0
    call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  end subroutine initialize_globals

  subroutine read_geqdsk_header(filename, rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr)
    character(len = *), intent(in) :: filename
    real(dp), intent(out) :: rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr

    character(len = *), parameter :: geqdsk_2000 = '(6a8,3i4)'
    character(len = *), parameter :: geqdsk_2020 = '(5e16.9)'

    character(len = 10) :: text(6)
    integer :: fid, k, idum, nw, nh
    real(dp) :: rcentr, simag, sibry

    open(newunit = fid, file = filename)
    read(fid, geqdsk_2000) (text(k), k = 1, 6), idum, nw, nh
    write(*, *) 'Header from G EQDSK file: ', (text(k), k = 1, 6)
    read(fid, geqdsk_2020) rdim, zdim, rcentr, rleft, zmid
    read(fid, geqdsk_2020) rmaxis, zmaxis, simag, sibry, bcentr
    close(fid)
  end subroutine read_geqdsk_header

  function interp_psi_pol(r, z) result(psi_pol)
    use field_eq_mod, only: psif

    real(dp), intent(in) :: r, z
    real(dp) :: psi_pol

    real(dp) :: p, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    p = 0d0
    call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
    psi_pol = psif
  end function interp_psi_pol

  subroutine add_node_owner(kpoint, ktri)
    use mesh_mod, only: n_owners_max
    integer, intent(in) :: kpoint, ktri
    if (mesh_point(kpoint)%n_owners < n_owners_max) then
       mesh_point(kpoint)%i_owner_tri(mesh_point(kpoint)%n_owners + 1) = ktri
       mesh_point(kpoint)%n_owners = mesh_point(kpoint)%n_owners + 1
    else
       write (logfile, *) 'Maximal number of owning triangles exceeded at point ', kpoint
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

end program minimal_example
