module fixed_coord
! to change the coordinate for minimum finding
  real*8 :: fixed 
  integer :: num_fixed
end module fixed_coord
! ------------------------------------------
module psi4root
  real*8 :: psi_root, R_o, Z_o, theta, R_x, Z_x, R0, Z0
end module psi4root
! -------------------------------------------
module const
  real*8, parameter :: pi = 3.14159265358979d0, twopi = 2.d0*pi
end module const
! -------------------------------------------
