subroutine twostreights(R11, Z11, R12, Z12, R21, Z21, R22, Z22, R, Z)
! calculates the point R,Z of intersection of two streight lines
! defined by two point each, 1st index - line number, 2nd - point number
  implicit none 
  real(kind=8), intent(IN) :: R11, Z11, R12, Z12, R21, Z21, R22, Z22  
  real(kind=8), intent(OUT) :: R, Z
  real(kind=8) :: dr1, dz1, rz1, zr1, dr2, dz2, alph
  dr1 = R12 - R11
  dz1 = Z12 - Z11
  rz1 = dr1/dz1
  zr1 = dz1/dr1

  dr2 = R22 - R21
  dz2 = Z22 - Z21
  if (abs(dr2) .ge. abs(dz2)) then
     alph = dz2/dr2
     R = (Z21 - Z11 + R11*zr1 - R21*alph)/(zr1 - alph)
     Z = Z11 + (R - R11)*zr1
  else
     alph = dr2/dz2
     Z = (R21 - R11 + Z11*rz1 - Z21*alph)/(rz1 - alph)
     R = R11 + (Z - Z11)*rz1           
  endif
  return
end subroutine twostreights
