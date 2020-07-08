program standardise_equilibrium
  use magdif_config, only: log_file, log_err, log_warn, log_info, log_debug, log_open, log_close
  use magdif_util, only: g_eqdsk

  implicit none

  character(len = 1024) :: unprocessed_geqdsk, standardised_geqdsk
  type(g_eqdsk) :: equil

  if (command_argument_count() >= 1) then
     call get_command_argument(1, unprocessed_geqdsk)
     if (command_argument_count() >= 2) then
        call get_command_argument(2, standardised_geqdsk)
     else
        error stop 'expected path to standardised G EQDSK file as second parameter'
     end if
  else
     error stop 'expected path to unprocessed G EQDSK file as first parameter'
  end if

  log_file = '-'
  log_err = .true.
  log_warn = .true.
  log_info = .true.
  log_debug = .true.
  call log_open
  call equil%read(trim(unprocessed_geqdsk))
  call equil%classify
  call equil%standardise
  call equil%write(trim(standardised_geqdsk))
  call log_close
end program standardise_equilibrium
