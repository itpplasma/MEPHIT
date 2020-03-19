! vacfield.f90
! writes out quantities from vacuum perturbation field
! into ascii files
program vacfield

  use from_nrtype, only: dp
  use mesh_mod, only: npoint, ntri, knot, mesh_point, mesh_element, mesh_element_rmp
  use magdif_config, only: n_tor_out => n, nflux, meshdata_file, Bn_vac_file, &
       config_file, read_config
  use magdif_util, only: initialize_globals, imun
  use magdif, only: fs, fs_half

  implicit none

  integer :: fid, ktri, ke
  real(dp) :: R, Z
  type(knot), dimension(3) :: knots
  real(dp), dimension(3) :: rr, zz, lr, lz
  complex(dp), dimension(3) :: Bnflux
  complex(dp) :: B_Rn, B_Zn, Bnphiflux

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  call read_config

  open(newunit = fid, file = meshdata_file, form = 'unformatted', status = 'old')
  read (fid) nflux, npoint, ntri
  call fs%init(nflux, .false.)
  call fs_half%init(nflux, .true.)
  read (fid) fs%psi
  read (fid) fs_half%psi
  allocate(mesh_point(npoint))
  read (fid) mesh_point
  allocate(mesh_element(ntri))
  read (fid) mesh_element
  close(fid)
  allocate(mesh_element_rmp(ntri))

  R = mesh_point(1)%rcoord
  Z = mesh_point(1)%zcoord
  call initialize_globals(R, Z)

  open(newunit = fid, file = trim(Bn_vac_file), status = 'replace')
  do ktri = 1, ntri
     knots = mesh_point(mesh_element(ktri)%i_knot)
     ! edge midpoints
     rr(1) = (knots(1)%rcoord + knots(2)%rcoord) * 0.5d0
     rr(2) = (knots(2)%rcoord + knots(3)%rcoord) * 0.5d0
     rr(3) = (knots(3)%rcoord + knots(1)%rcoord) * 0.5d0
     zz(1) = (knots(1)%zcoord + knots(2)%zcoord) * 0.5d0
     zz(2) = (knots(2)%zcoord + knots(3)%zcoord) * 0.5d0
     zz(3) = (knots(3)%zcoord + knots(1)%zcoord) * 0.5d0
     ! edge vector components (counterclockwise)
     lr(1) = knots(2)%rcoord - knots(1)%rcoord
     lr(2) = knots(3)%rcoord - knots(2)%rcoord
     lr(3) = knots(1)%rcoord - knots(3)%rcoord
     lz(1) = knots(2)%zcoord - knots(1)%zcoord
     lz(2) = knots(3)%zcoord - knots(2)%zcoord
     lz(3) = knots(1)%zcoord - knots(3)%zcoord
     do ke = 1, 3
        call spline_bpol_n(n_tor_out, rr(ke), zz(ke), B_Rn, B_Zn)
        ! project vector onto outward edge normal
        Bnflux(ke) = (lz(ke) * B_Rn - lr(ke) * B_Zn) * rr(ke)
     end do
     Bnphiflux = imun / n_tor_out * sum(Bnflux)
     write (fid, '(8(1x, es23.16))') Bnflux, Bnphiflux
  end do
  close(fid)

  if (allocated(mesh_point)) deallocate(mesh_point)
  if (allocated(mesh_element)) deallocate(mesh_element)
  if (allocated(mesh_element_rmp)) deallocate(mesh_element_rmp)

end program vacfield
