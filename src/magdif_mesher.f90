program magdif_preprocess

  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use magdif_conf, only: conf, magdif_config_read, log, magdif_log
  use magdif_mesh, only: generate_mesh, write_mesh_cache
  use magdif_pert, only: generate_vacfield

  implicit none

  character(len = 1024) :: config_file

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to magdif config file as first parameter'
  endif

  call h5_init
  h5overwrite = .true.
  call magdif_config_read(conf, config_file)
  log = magdif_log('-', conf%log_level, conf%quiet)
  call generate_mesh
  call write_mesh_cache
  call generate_vacfield
  call h5_deinit
end program magdif_preprocess
