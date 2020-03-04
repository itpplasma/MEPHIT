program geomint_mesh

  use magdif_config
  use magdif_util, only: get_equil_filenames, initialize_globals, g_eqdsk
  use points_2d, only: s_min, create_points_2d
  use from_nrtype, only: dp
  use mesh_mod, only: npoint
  use field_line_integration_mod, only: theta0_at_xpoint

  implicit none

  character(len = 1024) :: unprocessed_geqdsk, gfile, convexfile

  integer, dimension(:), allocatable :: n_theta
  real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi

  type(g_eqdsk) :: equil

  integer :: k, fid

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
     if (command_argument_count() >= 2) then
        call get_command_argument(2, unprocessed_geqdsk)
     else
        error stop 'expected path to unprocessed G EQDSK file as second parameter'
     end if
  else
     error stop 'expected path to magdif config file as first parameter'
  endif
  call read_config

  log_file = '-'
  call log_open

  call get_equil_filenames(gfile, convexfile)
  log_msg = 'attempting to read unprocessed G EQDSK file ' // trim(unprocessed_geqdsk)
  if (log_info) call log_write
  call equil%read(trim(unprocessed_geqdsk))
  call equil%classify
  call equil%standardise
  log_msg = 'attempting to write scaled G EQDSK file ' // trim(gfile)
  if (log_info) call log_write
  call equil%write(trim(gfile))

  npoint = nflux * nkpol + 1
  allocate(n_theta(nflux))
  allocate(points(3, npoint))
  allocate(points_s_theta_phi(3, npoint))

  n_theta = nkpol
  s_min = 1d-16
  theta0_at_xpoint = .true.
  call initialize_globals(equil%rmaxis, equil%zmaxis)
  call create_points_2d(n_theta, points, points_s_theta_phi, r_scaling_func = sqr)
  points(:, 1) = [equil%rmaxis, 0d0, equil%zmaxis]
  open(newunit = fid, file = 'points.fmt')
  do k = 1, npoint
     write (fid, '(2(1x, es23.16))') points(1, k), points(3, k)
  end do
  close(fid)

  if (allocated(n_theta)) deallocate(n_theta)
  if (allocated(points)) deallocate(points)
  if (allocated(points_s_theta_phi)) deallocate(points_s_theta_phi)

contains

  pure function sqr(x) result(x_squared)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: x_squared
    x_squared = x * x
  end function sqr
end program geomint_mesh
