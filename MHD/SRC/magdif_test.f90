program magdif_test

  use magdif, only: config_file, magdif_init, magdif_cleanup, &
       compute_presn, compute_j0phi, compute_currn

  implicit none

  if (command_argument_count() == 1) then
     call get_command_argument(1, config_file)
  else
     stop 'Error: expected path to config file as first parameter'
  endif

  call magdif_init
  call compute_presn
  call compute_j0phi
  call compute_currn
  call magdif_cleanup
   
end program magdif_test
