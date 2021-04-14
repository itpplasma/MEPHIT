program vacfield

  use iso_fortran_env, only: dp => real64
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use magdif_conf, only: log, magdif_log
  use magdif_pert, only: AUG_coils_read, AUG_coils_write_Nemov, AUG_coils_read_Nemov, &
       AUG_coils_write_GPEC, AUG_coils_read_GPEC, AUG_coils_write_Fourier

  implicit none

  integer :: argc, nseg, ncoil, nwind
  character(len = 1024) :: in_type, out_type, in_dir, out_dir
  real(dp), allocatable :: XYZ(:, :, :)

  argc = command_argument_count()
  if (argc >= 4) then
     call get_command_argument(1, in_type)
     call get_command_argument(2, out_type)
     call get_command_argument(3, in_dir)
     call get_command_argument(4, out_dir)
  else
     error stop 'Expected 4 command line arguments, see documentation'
  endif
  log = magdif_log('-', 4, .false.)

  call h5_init
  h5overwrite = .true.
  if (in_type == 'AUG') then
     call AUG_coils_read(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else if (in_type == 'GPEC') then
     call AUG_coils_read_GPEC(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else if (in_type == 'Nemov') then
     call AUG_coils_read_Nemov(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else
     write (log%msg, '("unknown input type ", a)') trim(in_type)
     if (log%err) call log%write
     error stop
  end if
  if (out_type == 'Fourier') then
     ! use default values for now
     call AUG_coils_write_Fourier(trim(out_dir), ncoil, nseg, nwind, XYZ, &
          64, 75.0d0, 267.0d0, -154.0d0, 150.4d0, 150, 512, 300)
  else if (out_type == 'GPEC') then
     call AUG_coils_write_GPEC(trim(out_dir), ncoil, nseg, nwind, XYZ)
  else if (out_type == 'Nemov') then
     call AUG_coils_write_Nemov(trim(out_dir), ncoil, nseg, XYZ)
  else
     write (log%msg, '("unknown output type ", a)') trim(out_type)
     if (log%err) call log%write
     error stop
  end if
  if (allocated(XYZ)) deallocate(XYZ)
  call h5_deinit

end program vacfield
