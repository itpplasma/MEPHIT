program magdif_test
  use magdif_config, only: config_file, read_config, runmode, runmode_single, &
       runmode_direct, runmode_precon
  use magdif, only: magdif_init, magdif_cleanup, magdif_single, magdif_direct

  implicit none

  if (command_argument_count() == 1) then
     call get_command_argument(1, config_file)
  else
     stop 'Error: expected path to config file as first parameter'
  endif

  call read_config(config_file)
  call magdif_init
  select case (runmode)
  case (runmode_single)
     call magdif_single
  case (runmode_direct)
     call magdif_direct
  case (runmode_precon)
     call magdif_direct
  case default
     stop 'Error: unknown runmode'
  end select
  call magdif_cleanup

end program magdif_test
