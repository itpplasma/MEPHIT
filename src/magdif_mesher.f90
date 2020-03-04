program magdif_mesher

  use magdif_config, only: config_file, read_config, log_file, log_open, log_close
  use magdif_mesh_mod, only: generate_mesh

  implicit none

  character(len = 1024) :: unprocessed_geqdsk

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

  call read_config
  log_file = '-'
  call log_open

  call generate_mesh(unprocessed_geqdsk)

  call log_close

end program magdif_mesher
