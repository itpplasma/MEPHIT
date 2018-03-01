!"Light" fortran version of Pierre L'Ecyuer's RNGSTREAM
! Operations Research © 2002 INFORMS
! Vol. 50, No. 6, November–December 2002, pp. 1073–1075

module rngstream_short
  real*8, parameter :: m1 = 4294967087.d0
  real*8, parameter :: m2 = 4294944443.d0
end module rngstream_short
!-------------------------------------------------------------------------
! Return (a*s + c) MOD m; a, s, c and m must be < 2^35
!
function MultModM(a_in, s, c, m)
  implicit none
  real*8, parameter ::  two17 = 131072.d0
  real*8, parameter ::  two53 = 9007199254740992.d0
  real*8 a, s, c, m, MultModM, a_in
  real*8 v
  integer*4 a1

  a = a_in
  v = a*s + c

  if (v .ge. two53 .or. v .le. -two53) then
     a1 = int(a / two17)
     a  = a - a1*two17
     v  = a1*s
     a1 = int(v / m)  
     v  = v - a1*m
     v = v*two17 + a*s + c
  endif

  a1 = int(v / m)
! in case v < 0)
  v = v - a1*m 
  if (v .lt. 0.d0) then 
     MultModM = v + m  
  else
     MultModM = v
  endif
  return
end function MultModM
!-------------------------------------------------------------------------
! Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
! Works also when v = s.
subroutine MatVecModM(A, s, v, m)
  implicit none
  real*8 MultModM, m
  real*8 A(0:2,0:2)
  integer i
  real*8, dimension(0:2) :: s, v, x              ! Necessary if v = s

  do i=0,2
     x(i) = MultModM(A(i,0), s(0), 0.d0, m)
     x(i) = MultModM(A(i,1), s(1), x(i), m)
     x(i) = MultModM(A(i,2), s(2), x(i), m)
  enddo
  v=x
  return
end subroutine   MatVecModM
!-------------------------------------------------------------------------
! constructor
subroutine RngStream(nextseed, Cg)
!    use lpk_include, only: m1, m2, Cg
! Information on a stream. The arrays {Cg, Bg, Ig} contain the current
! state of the stream, the starting state of the current SubStream, and the
! starting state of the stream. This stream generates antithetic variates
! if anti = true. It also generates numbers with extended precision (53
! bits if machine follows IEEE 754 standard) if incPrec = true. nextSeed
! will be the seed of the next declared RngStream. 
  use rngstream_short
  implicit none 
  real*8, parameter :: A1p127(0:2,0:2) = reshape((/ 2427906178.d0, 226153695.d0,  1988835001.d0,   &
                                                    3580155704.d0, 1230515664.d0, 986791581.d0,    &
                                                    949770784.d0,  3580155704.d0, 1230515664.d0 /), (/3,3/))
 
  real*8, parameter :: A2p127(0:2,0:2) = reshape((/ 1464411153.d0, 32183930.d0,   2824425944.d0,   &
                                                    277697599.d0,  1464411153.d0, 32183930.d0,     &
                                                    1610723613.d0, 1022607788.d0, 2093834863.d0 /), (/3,3/))
  real*8, dimension(0:5) :: nextSeed, Cg

  Cg = nextSeed 

  call MatVecModM (A1p127, nextSeed(0:2), nextSeed(0:2), m1)
  call MatVecModM (A2p127, nextSeed(3:5), nextSeed(3:5), m2)
  return
end subroutine RngStream
!-------------------------------------------------------------------------
! Generate the next uniformly distributed (0,1) random number.
function ran_u_01(Cg)
  use rngstream_short

  implicit none  
  real*8, parameter :: norm = 1.d0/(m1 + 1.d0) 
  integer*4  k
  real*8 p1, p2, ran_u_01
  real*8, parameter :: a12 = 1403580.d0, a21 = 527612.d0,       &
                       a13n = 810728.d0, a23n = 1370589.d0
  real*8 Cg(0:5)
! Component 1 
  p1 = a12*Cg(1) - a13n*Cg(0)
  k = int(p1 / m1)
  p1 = p1 - k*m1
  if (p1 .lt. 0.d0) p1 = p1 + m1
  Cg(0) = Cg(1)
  Cg(1) = Cg(2)
  Cg(2) = p1
! Component 2 
  p2 = a21*Cg(5) - a23n*Cg(3)
  k = int(p2 / m2)
  p2 = p2 - k*m2
  if (p2 .lt. 0.d0) p2 = p2 + m2
  Cg(3) = Cg(4)
  Cg(4) = Cg(5)
  Cg(5) = p2
! Combination 
  if (p1 .gt. p2) then
     ran_u_01 = (p1 - p2)*norm
  else 
     ran_u_01 = (p1 - p2 + m1)*norm
  endif
  return
end function ran_u_01
!------------------------------------------------------------------
! Generate uniformly distributed (-sqrt(3),+sqrt(3)) random number
! (unit variance)
function ran_u_sq3(Cg)
  real*8 ran_u_01, ran_u_sq3
  real*8 Cg(0:5)
  ran_u_sq3 = 3.46410161513775458705d0*(ran_u_01(Cg) - 5.d-1)
  return
end function ran_u_sq3
!------------------------------------------------------------------
! Generate \pm 1 random number (unit variance)
function ran_pm1(Cg)
  real*8 ran_u_01, ran_pm1
  real*8 Cg(0:5)
  ran_pm1 = 1.d0
  if(ran_u_01(Cg) .lt.  5.d-1) ran_pm1 = -1.d0
  return
end function ran_pm1
!----------------------------------------------------------------------
! Generate gaussian random number
function gauss_openmp(Cg)
  implicit none
  real*8 ran_u_01
  real*8 Cg(0:5)
  real*8 gauss_openmp
  integer :: iset=0
  real*8 fac, gset, rsq, v1, v2
  save iset, gset
!$OMP THREADPRIVATE (iset, gset)

  if (iset .eq. 0) then
1    v1 = 2.d0*ran_u_01(Cg) - 1.d0
     v2 = 2.d0*ran_u_01(Cg) - 1.d0
     rsq = v1**2 + v2**2
     if(rsq.ge.1.d0 .or. rsq.eq.0.d0) goto 1
     fac = sqrt(-2.d0*log(rsq)/rsq)
     gset = v1*fac
     gauss_openmp = v2*fac
     iset = 1
  else 
     gauss_openmp = gset
     iset = 0
  endif

  return
end function gauss_openmp
!-------------------------------------------------------------------------------------------
! No-OpenMP generator from Numerical Recepies 
function ran2(idum)
  integer :: idum, j, k
  real*8 :: ran2
  integer, parameter :: im1 = 2147483563, im2 = 2147483399, imm1 = im1-1, ia1 = 40014, ia2 = 40692,      &
                        iq1 = 53668, iq2 = 52774, ir1 = 12211, ir2 = 3791, ntab = 32, ndiv = 1+imm1/ntab
  real*8, parameter :: am = 1.d0/im1, eps = 1.2d-7, rnmx = 1.d0-eps
  integer, save :: idum2 = 123456789, iv(1:ntab) = 0, iy = 0

  if (idum .le. 0) then
     idum=max(-idum,1)
     idum2=idum
     do  j=ntab+8,1,-1
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if (idum.lt.0) idum=idum+im1
        if (j.le.ntab) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/iq1
  idum=ia1*(idum-k*iq1)-k*ir1
  if (idum.lt.0) idum=idum+im1
  k=idum2/iq2
  idum2=ia2*(idum2-k*iq2)-k*ir2
  if (idum2.lt.0) idum2=idum2+im2
  j=1+iy/ndiv
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+imm1
  ran2=min(am*iy,rnmx)
  return
end function ran2
