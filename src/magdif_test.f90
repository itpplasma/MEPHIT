program magdif_test
  use magdif_conf, only: conf, magdif_config_read, conf_arr
  use magdif, only: magdif_init, magdif_cleanup, magdif_iterate, magdif_postprocess
  use hdf5_tools, only: h5_init, h5_deinit

  implicit none

  character(len = 1024) :: config_file

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif

  call h5_init
  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  call magdif_init
  call magdif_iterate
  call magdif_postprocess
  call magdif_cleanup
  call h5_deinit

end program magdif_test
