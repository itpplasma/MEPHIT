program geomint_mesh

  use magdif_config
  use points_2d
  use from_nrtype, only: dp
  use input_files, only: gfile
  use mesh_mod, only: npoint, bphicovar

  implicit none

  real(dp) :: rdim, zdim, rleft, zmid, rmaxis, zmaxis, bcentr

  integer, dimension(:), allocatable :: n_theta
  real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi

  integer :: k, fid

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     stop 'Error: expected path to magdif config file as first parameter'
  endif
  call read_config

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

  open(newunit = fid, file = 'preload_for_SYNCH.inp')
  write (fid, '(i4)') 3600
  write (fid, '(i4)') 3000
  write (fid, '(i3)') 999
  write (fid, '(i5)') 10000
  close(fid)

  npoint = nflux * nkpol + 1
  allocate(n_theta(nflux))
  allocate(points(3, npoint))
  allocate(points_s_theta_phi(3, npoint))

  n_theta = nkpol
  s_min = 1d-16
  call create_points_2d(n_theta, points, points_s_theta_phi, r_scaling_func = sqr)
  points(:, 1) = [rmaxis, 0d0, zmaxis]
  open(newunit = fid, file = 'points.fmt')
  do k = 1, npoint
     write (fid, '(2(1x, es23.16))') points(1, k), points(3, k)
  end do
  close(fid)

  if (allocated(n_theta)) deallocate(n_theta)
  if (allocated(points)) deallocate(points)
  if (allocated(points_s_theta_phi)) deallocate(points_s_theta_phi)

contains
    subroutine initialize_globals
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    if (kilca_scale_factor /= 0) then
       r = R0 * kilca_scale_factor
    else
       r = R0
    end if
    z = 0d0
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
    write (*, '("Header from G EQDSK file: ", 6a)') (text(k), k = 1, 6)
    read(fid, geqdsk_2020) rdim, zdim, rcentr, rleft, zmid
    read(fid, geqdsk_2020) rmaxis, zmaxis, simag, sibry, bcentr
    close(fid)
  end subroutine read_geqdsk_header

  pure function sqr(x) result(x_squared)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: x_squared
    x_squared = x * x
  end function sqr
end program geomint_mesh
