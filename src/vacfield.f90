! vacfield.f90
! writes out quantities from vacuum perturbation field
! into ascii files
program vacfield

  use iso_fortran_env, only: dp => real64
  use mesh_mod, only: npoint, ntri, knot, mesh_point, mesh_element, mesh_element_rmp
  use magdif_conf, only: conf, magdif_config_read
  use magdif_util, only: initialize_globals, gauss_legendre_unit_interval, imun
  use magdif, only: fs, fs_half

  implicit none

  character(len = 1024) :: config_file
  integer, parameter :: order = 2
  integer :: fid, ktri, ke, k
  real(dp) :: R, Z, edge_R, edge_Z, node_R(4), node_Z(4)
  real(dp), dimension(order) :: points, weights
  complex(dp) :: B_Rn, B_Zn, Bnflux(3), Bnphiflux

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  call magdif_config_read(conf, config_file)

  call gauss_legendre_unit_interval(order, points, weights)

  open(newunit = fid, file = conf%meshdata_file, form = 'unformatted', status = 'old')
  read (fid) conf%nflux, npoint, ntri
  call fs%init(conf%nflux, .false.)
  call fs_half%init(conf%nflux, .true.)
  read (fid) fs%psi, fs%rad
  read (fid) fs_half%psi, fs_half%rad
  allocate(mesh_point(npoint))
  read (fid) mesh_point
  allocate(mesh_element(ntri))
  read (fid) mesh_element
  close(fid)
  allocate(mesh_element_rmp(ntri))

  R = mesh_point(1)%rcoord
  Z = mesh_point(1)%zcoord
  call initialize_globals(R, Z)

  open(newunit = fid, file = trim(conf%Bn_vac_file), status = 'replace')
  do ktri = 1, ntri
     Bnflux = (0d0, 0d0)
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
             call spline_bpol_n(conf%n, R, Z, B_Rn, B_Zn)
             ! project vector onto outward edge normal
             Bnflux(ke) = Bnflux(ke) + weights(k) * (B_Rn * edge_Z - B_Zn * edge_R) * R
          end do
       end do
     end associate
     ! toroidal flux via zero divergence
     Bnphiflux = imun / conf%n * sum(Bnflux)
     write (fid, '(8(1x, es24.16e3))') Bnflux, Bnphiflux
  end do
  close(fid)

  if (allocated(mesh_point)) deallocate(mesh_point)
  if (allocated(mesh_element)) deallocate(mesh_element)
  if (allocated(mesh_element_rmp)) deallocate(mesh_element_rmp)

end program vacfield
