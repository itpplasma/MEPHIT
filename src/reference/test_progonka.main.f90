!
  implicit none
!
  integer :: npoi,i
!
  double precision :: xbeg,xend,hx,x
!
  double precision, dimension(:),   allocatable :: xarr
  double complex,   dimension(:),   allocatable :: farr,garr,qarr
  double complex,   dimension(:),   allocatable :: a2_in,a2_out,a0
  double complex,   dimension(:),   allocatable :: b2_in,b2_out,b0
  double complex,   dimension(:),   allocatable :: c2_in,c2_out,c0
  double complex,   dimension(:),   allocatable :: d2_in,d2_out,d0
!
  xbeg=0.d0
  xend=3.1415926
  xbeg=-xend !/2.d0*3.d0
  xend=-xbeg
!
  npoi=300
  allocate(xarr(npoi),farr(npoi),garr(npoi),qarr(npoi))
  allocate(a2_in(npoi),a2_out(npoi),a0(npoi))
  allocate(b2_in(npoi),b2_out(npoi),b0(npoi))
  allocate(c2_in(npoi),c2_out(npoi),c0(npoi))
  allocate(d2_in(npoi),d2_out(npoi),d0(npoi))
!
  hx=(xend-xbeg)/dble(npoi-1)
!
  do i=1,npoi
    xarr(i)=(xbeg+hx*dble(i-1))**3
    qarr(i)=exp((0.d0,1.d0)*xarr(i))
    a2_in(i)=1.d0/(1.d0+xarr(i)**2)
    a2_out(i)=a2_in(i)
    a0(i)=-1.d0 !+xarr(i)
    b2_in(i)=-1.d0
    b2_out(i)=b2_in(i)
    b0(i)=1.d0+xarr(i)
  enddo
!
  c2_in=-a2_in
  c2_out=-a2_out
  c0=-a0
  d2_in=b2_in
  d2_out=b2_out
  d0=b0
!
  call progonka(npoi,xarr,qarr,  &
                a2_in,a2_out,a0, &
                b2_in,b2_out,b0, &
                c2_in,c2_out,c0, &
                d2_in,d2_out,d0, &
                farr,garr)
!
  do i=1,npoi
    write (12,*) xarr(i),dble(qarr(i)),dble(farr(i)),dble(garr(i))
    write (14,*) xarr(i),aimag(qarr(i)),aimag(farr(i)),aimag(garr(i))
  enddo
!
  end
