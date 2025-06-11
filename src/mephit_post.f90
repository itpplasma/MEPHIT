program mephit_debug

  use iso_fortran_env, only: dp => real64
  use hdf5_tools, only: HID_T, h5_init, h5_open, h5_get, h5_close, h5_deinit
  use mephit_util, only: ev2erg, resample1d

  implicit none

  integer(HID_T) :: h5id_root
  character(len = 1024) :: suffix, str_m, datafile
  integer :: fid, k, m, m_res_min, m_res_max, m_max, m_res, n, nflux, sgn_dpsi
  real(dp) :: m_i, Z_i, R_O
  real(dp), dimension(:), allocatable :: psi, half_psi, q, half_q, &
    F, Phi0, avg_R2gradpsi2, dens_e, temp_e, temp_i, nu_e, nu_i
  complex(dp), dimension(:), allocatable :: Bmnpsi_over_B0phi
  complex(dp), dimension(:, :), allocatable :: Phimn, Phimn_aligned, Bmn_psi

  if (command_argument_count() < 2) then
    error stop 'expected filename and resonant mode number as parameters'
  end if
  call get_command_argument(1, datafile)
  call get_command_argument(2, str_m)
  read (str_m, *) m
  print '("datafile ''", a, "''")', trim(datafile)
  print '("resonant mode number m = ", i0)', m

  call h5_init
  call h5_open(trim(datafile), h5id_root)
  call h5_get(h5id_root, 'config/n', n)
  call h5_get(h5id_root, 'config/m_i', m_i)
  call h5_get(h5id_root, 'config/Z_i', Z_i)
  call h5_get(h5id_root, 'mesh/R_O', R_O)
  call h5_get(h5id_root, 'mesh/m_res_min', m_res_min)
  call h5_get(h5id_root, 'mesh/m_res_max', m_res_max)
  call h5_get(h5id_root, 'mesh/nflux', nflux)
  allocate(psi(0:nflux), q(0:nflux), F(0:nflux), Phi0(0:nflux), avg_R2gradpsi2(0:nflux), &
    dens_e(0:nflux), temp_e(0:nflux), temp_i(0:nflux), nu_e(0:nflux), nu_i(0:nflux))
  call h5_get(h5id_root, 'cache/fs/psi', psi)
  call h5_get(h5id_root, 'cache/fs/q', q)
  m_res = -abs(m) * int(sign(1d0, q(nflux)))
  allocate(half_psi(1:nflux), half_q(1:nflux))
  call h5_get(h5id_root, 'cache/fs_half/psi', half_psi)
  call h5_get(h5id_root, 'cache/fs_half/q', half_q)
  call h5_get(h5id_root, 'cache/fs/F', F)
  call h5_get(h5id_root, 'equil/cocos/sgn_dpsi', sgn_dpsi)
  call h5_get(h5id_root, 'equil/profiles/Phi0/y', Phi0)
  call h5_get(h5id_root, 'equil/profiles/dens_e/y', dens_e)
  call h5_get(h5id_root, 'equil/profiles/temp_e/y', temp_e)
  call h5_get(h5id_root, 'equil/profiles/temp_i/y', temp_i)
  call h5_get(h5id_root, 'equil/profiles/nu_e/y', nu_e)
  call h5_get(h5id_root, 'equil/profiles/nu_i/y', nu_i)
  call h5_get(h5id_root, 'mesh/avg_R2gradpsi2', avg_R2gradpsi2)
  allocate(Phimn(0:nflux, m_res_min:m_res_max), Phimn_aligned(0:nflux, m_res_min:m_res_max))
  call h5_get(h5id_root, 'iter/Phi_mn', Phimn)
  call h5_get(h5id_root, 'iter/Phi_aligned_mn', Phimn_aligned)
  call h5_get(h5id_root, 'iter/Bmn/m_max', m_max)
  allocate(Bmn_psi(-m_max:m_max, 1:nflux), Bmnpsi_over_B0phi(0:nflux))
  call h5_get(h5id_root, 'iter/Bmn/coeff_rad', Bmn_psi)
  call resample1d(half_psi, Bmn_psi(m_res, :)%Re / half_q * sgn_dpsi, psi, Bmnpsi_over_B0phi%Re, 3)
  call resample1d(half_psi, Bmn_psi(m_res, :)%Im / half_q * sgn_dpsi, psi, Bmnpsi_over_B0phi%Im, 3)
  call h5_close(h5id_root)
  call h5_deinit

  open(newunit = fid, file = 'response_current.dat')
  write (fid, *) 1, -m_res, n, nflux
  write (fid, *) m_i, Z_i, R_O
  write (fid, *) psi(1:)
  write (fid, *) -q(1:)
  write (fid, *) F(1:)
  write (fid, *) Phi0(1:)
  write (fid, *) avg_R2gradpsi2(1:)
  write (fid, *) dens_e(1:)
  write (fid, *) temp_e(1:) * ev2erg
  write (fid, *) temp_i(1:) * ev2erg
  write (fid, *) nu_e(1:)
  write (fid, *) nu_i(1:)
  close(fid)

  open(newunit = fid, file = 'bpsi_over_bphi.dat')
  do k = 1, nflux
    write (fid, *) psi(k), Bmnpsi_over_B0phi(k)%Re, Bmnpsi_over_B0phi(k)%Im
  end do
  close(fid)

  open(newunit = fid, file = 'Phi_m.dat')
  do k = 1, nflux
    write (fid, *) psi(k), Phimn(k, abs(m))%Re, Phimn(k, abs(m))%Im
  end do
  close(fid)

  open(newunit = fid, file = 'Phi_m_aligned.dat')
  do k = 1, nflux
    write (fid, *) psi(k), Phimn_aligned(k, abs(m))%Re, Phimn_aligned(k, abs(m))%Im
  end do
  close(fid)

end program mephit_debug
