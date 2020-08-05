program magdif_test
  use magdif_conf, only: conf, magdif_config_read, conf_arr, runmode_single, &
       runmode_direct, runmode_precon
  use magdif, only: magdif_init, magdif_cleanup, magdif_single, magdif_iterated

  implicit none

  character(len = 1024) :: config_file, bin_dir

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  if (command_argument_count() >= 2) then
     call get_command_argument(2, bin_dir)
  end if

  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  call magdif_init(bin_dir)
  select case (conf%runmode)
  case (runmode_single)
     call magdif_single
  case (runmode_direct)
     call magdif_iterated
  case (runmode_precon)
     call magdif_iterated
  case default
     error stop 'unknown runmode'
  end select
  call magdif_cleanup

end program magdif_test
