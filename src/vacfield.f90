program vacfield

  use iso_fortran_env, only: dp => real64
  use magdif_conf, only: log, magdif_log
  use magdif_pert, only: AUG_coils_read, AUG_coils_write_Nemov, AUG_coils_read_Nemov, &
       AUG_coils_write_GPEC, AUG_coils_read_GPEC

  implicit none

  integer :: argc, nseg, ncoil, nwind, ind
  character(len = 1024) :: executable, operation, input, output
  real(dp), allocatable :: XYZ(:, :, :)

  argc = command_argument_count()
  if (argc >= 2) then
     call get_command_argument(1, operation)
     call get_command_argument(2, input)
     if (argc >= 3) then
        call get_command_argument(3, output)
     else
        output = input
     end if
  else
     call get_command_argument(0, executable)
     print '("Convert RMP coil geometry file formats.")'
     print '("Syntax: ", a, " operation source-directory [target-directory]")', trim(executable)
     print '("where operation is one of: AUG2GPEC, AUG2Nemov, GPEC2Nemov, Nemov2GPEC")'
     print '("When target-directory is omitted, it is assumed to be the same as source-directory.")'
     error stop 'Not enough command line arguments'
  endif
  log = magdif_log('-', 4, .false.)

  input = adjustl(input)
  ind = len_trim(input)
  if (input(ind:ind) == '/') input(ind:ind) = ' '
  ind = len_trim(output)
  output = adjustl(output)
  if (output(ind:ind) == '/') output(ind:ind) = ' '
  if (operation == 'AUG2GPEC') then
     call AUG_coils_read(trim(input), ncoil, nseg, nwind, XYZ)
     call AUG_coils_write_GPEC(trim(output), ncoil, nseg, nwind, XYZ)
  else if (operation == 'AUG2Nemov') then
     call AUG_coils_read(trim(input), ncoil, nseg, nwind, XYZ)
     call AUG_coils_write_Nemov(trim(output), ncoil, nseg, XYZ)
  else if (operation == 'GPEC2Nemov') then
     call AUG_coils_read_GPEC(trim(input), ncoil, nseg, nwind, XYZ)
     call AUG_coils_write_Nemov(trim(output), ncoil, nseg, XYZ)
  else if (operation == 'Nemov2GPEC') then
     call AUG_coils_read_Nemov(trim(input), ncoil, nseg, nwind, XYZ)
     call AUG_coils_write_GPEC(trim(output), ncoil, nseg, nwind, XYZ)
  else
     write (log%msg, '("unknown operation ", a)') trim(operation)
     if (log%err) call log%write
     error stop
  end if
  if (allocated(XYZ)) deallocate(XYZ)

end program vacfield
