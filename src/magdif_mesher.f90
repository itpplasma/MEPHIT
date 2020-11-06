program magdif_preprocess

  use hdf5_tools, only: h5_init, h5_deinit
  use magdif_conf, only: conf, magdif_config_read, conf_arr, log, magdif_log
  use magdif_mesh, only: generate_mesh, write_mesh_cache

  implicit none

  character(len = 1024) :: config_file, unprocessed_geqdsk

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
     if (command_argument_count() >= 2) then
        call get_command_argument(2, unprocessed_geqdsk)
     else
        error stop 'expected path to unprocessed G EQDSK file as second parameter'
     end if
  else
     error stop 'expected path to magdif config file as first parameter'
  endif

  call h5_init
  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  log = magdif_log('-', conf%log_level, conf%quiet)
  call generate_mesh(unprocessed_geqdsk)
  call write_mesh_cache
  call h5_deinit
end program magdif_preprocess
