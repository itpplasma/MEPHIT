program standardise_equilibrium
  use mephit_conf, only: logger
  use mephit_util, only: g_eqdsk

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

  call logger%init('-', 4, .false.) ! log_level = debug
  call equil%read(trim(unprocessed_geqdsk), '')
  call equil%classify
  call equil%standardise
  call equil%write(trim(standardised_geqdsk))
  call equil%deinit
  call logger%deinit
end program standardise_equilibrium
