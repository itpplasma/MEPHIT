program magdif_test

  use magdif, only: magdif_init, magdif_cleanup, compute_presn

  implicit none

  call magdif_init
  call compute_presn
  call magdif_cleanup
   
end program magdif_test
