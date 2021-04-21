program magdif_test
  use magdif_conf, only: conf, magdif_config_read
  use magdif, only: magdif_init, magdif_cleanup, magdif_iterate, magdif_postprocess, &
       freefem_pipe
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite

  implicit none

  character(len = 1024) :: config_file

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  if (command_argument_count() >= 2) then
     call get_command_argument(2, freefem_pipe)
  else
     error stop 'expected path to named pipe as second parameter'
  endif

  call h5_init
  h5overwrite = .true.
  call magdif_config_read(conf, config_file)
  call magdif_init
  call magdif_iterate
  call magdif_postprocess
  call magdif_cleanup
  call h5_deinit

end program magdif_test
