module lsq
  integer, parameter :: npoly=7
  real coef(0:npoly)
end module lsq
!--------------------------------------------------------
subroutine leastsq(npoint,xi,f)
! least square polynom (order NPOLY, see above)
  use lsq
  implicit none
  integer :: npoint
  real, dimension(npoint) :: xi, f
  integer, parameter :: nx=2*npoly, npp1=npoly+1
  integer, dimension(npp1) :: indx
  real, dimension(0:nx) :: sx
  real, dimension(0:npoly) :: sf
  real, dimension(npp1,npp1) :: a
  real, dimension(npp1) :: b
  real :: d,val_lsq
  integer :: i,k,l

  do k = 1, nx
     sx(k) = 0.
     do i=1,npoint
        sx(k) = sx(k)+ xi(i)**k
     enddo
  enddo
  sx(0) = npoint

  do l = 0, npoly
     sf(l) = 0.
     do i=1,npoint
        sf(l) = sf(l) + f(i)*xi(i)**l
     enddo
  enddo

  a = 0.
  b = 0.
  d = 0.
  do i=1,npp1
     do k=1,npp1
        a(i,k) = sx(i+k-2)
     enddo
     b(i) = sf(i-1)
  enddo

  call ludcmp(a,npp1,npp1,indx,d)
  call lubksb(a,npp1,npp1,indx,b)

  do i=1,npp1
     coef(i-1) = b(i)
  enddo
!
  do i=1,npoint
    f(i)=val_lsq(xi(i))
  enddo
!
  return
end subroutine leastsq
! ----------------------------------------------------------------
function val_lsq(x)
  use lsq
  implicit none
  integer :: i
  real :: x, val_lsq

  val_lsq = 0.
  do i=npoly,1,-1
     val_lsq = val_lsq + coef(i)
     val_lsq = val_lsq*x
  enddo
  val_lsq = val_lsq + coef(0)

  return
end function val_lsq
! ----------------------------------------------------------------
function der_lsq(x)
  use lsq
  implicit none
  integer :: i
  real :: x, der_lsq

  der_lsq = 0.
  do i=npoly,2,-1
     der_lsq = der_lsq + coef(i)*i
     der_lsq = der_lsq*x
  enddo
  der_lsq = der_lsq + coef(1)

  return
end function der_lsq
