FUNCTION rtbis(func,x1,x2,xacc)
  USE from_nrtype ! ; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x1,x2,xacc
  REAL(DP) :: rtbis, rhohuge=1.d50
  INTERFACE
     FUNCTION func(x)
       USE from_nrtype
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=40
  INTEGER(I4B) :: j
  REAL(DP) :: dx,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0) then
!     print *,x1, x2, func(x1), func(x2)
!     print *, 'rtbis: root must be bracketed' ! call nrerror('rtbis: root must be bracketed')
     rtbis=rhohuge  ! attention, just very large number
     return
  end if
  if (f < 0.0) then
     rtbis=x1
     dx=x2-x1
  else
     rtbis=x2
     dx=x1-x2
  end if
  do j=1,MAXIT
     dx=dx*0.5_sp
     xmid=rtbis+dx
     fmid=func(xmid)
     if (fmid <= 0.0) rtbis=xmid
     if (abs(dx) < xacc .or. fmid == 0.0) RETURN
  end do
  stop 'rtbis: too many bisections'
!  call nrerror('rtbis: too many bisections')
END FUNCTION rtbis
