program magdif_test

  use magdif, only : init, test_pressure_profile, test_current_new, test_Bnp
  use mesh_mod, only: mesh_element, mesh_point

  implicit none

  integer :: k,l

  call init
  
!  call test_pressure_profile
!  call test_current_new

!  call test_fdm_simple
!  call test_triangle_strip
!  call test_current
!  call test_dpsidr
  call test_Bnp
  
end program magdif_test
