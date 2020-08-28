! vacfield.f90
! writes out quantities from vacuum perturbation field
! into ascii files
program vacfield

  use iso_fortran_env, only: dp => real64
  use iso_c_binding, only: c_long
  use mesh_mod, only: npoint, ntri, knot, mesh_point, mesh_element
  use hdf5_tools, only: HID_T, h5_init, h5_open, h5_get, h5_close, h5_deinit
  use magdif_conf, only: conf, magdif_config_read
  use magdif_util, only: initialize_globals, gauss_legendre_unit_interval, imun
  use magdif_mesh, only: mesh

  implicit none

  character(len = 1024) :: config_file
  integer, parameter :: order = 2
  integer :: fid, ktri, ke, k
  integer(HID_T) :: h5id_magdif
  integer(c_long) :: length
  integer, dimension(:, :), allocatable :: tri_node
  real(dp) :: R, Z, edge_R, edge_Z, node_R(4), node_Z(4)
  real(dp), dimension(order) :: points, weights
  complex(dp) :: B_Rn, B_Zn, Bnflux(3), Bnphiflux

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  call h5_init
  call magdif_config_read(conf, config_file)

  call gauss_legendre_unit_interval(order, points, weights)

  call h5_open('magdif.h5', h5id_magdif)
  call h5_get(h5id_magdif, 'mesh/R_O', mesh%R_O)
  call h5_get(h5id_magdif, 'mesh/Z_O', mesh%Z_O)
  call h5_get(h5id_magdif, 'mesh/npoint', npoint)
  call h5_get(h5id_magdif, 'mesh/ntri', ntri)
  allocate(mesh_point(npoint))
  allocate(mesh_element(ntri))
  allocate(tri_node(3, ntri))
  call h5_get(h5id_magdif, 'mesh/node_R', mesh_point%rcoord)
  call h5_get(h5id_magdif, 'mesh/node_Z', mesh_point%zcoord)
  call h5_get(h5id_magdif, 'mesh/tri_node', tri_node)
  call h5_close(h5id_magdif)
  ! intermediate until mesh_mod is refactored into magdif_mesh
  do ktri = 1, ntri
     mesh_element(ktri)%i_knot = tri_node(:, ktri)
  end do
  deallocate(tri_node)

  R = mesh%R_O
  Z = mesh%Z_O
  call initialize_globals(R, Z)

  open(newunit = fid, file = conf%Bn_vac_file, access = 'stream', status = 'replace')
  length = 8 * ntri
  write (fid) length
  do ktri = 1, ntri
     Bnflux = (0d0, 0d0)
     associate(knots => mesh_point(mesh_element(ktri)%i_knot))
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
     write (fid) Bnflux, Bnphiflux
  end do
  close(fid)

  if (allocated(mesh_point)) deallocate(mesh_point)
  if (allocated(mesh_element)) deallocate(mesh_element)
  call h5_deinit
end program vacfield
