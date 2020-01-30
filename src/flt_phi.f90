program flt_phi
  use period_mod

  implicit none

  integer, parameter :: nequat = 2
  real*8, parameter :: pi = 3.14159265358979d0
  real*8, dimension (:) :: y(nequat)
  real*8 :: rmn, rmx, zmn, zmx, raxis, zaxis, rmnplot, rmxplot, phi0, zet0
  real*8 :: hn1, hn1min, relerr, ah, r, z, phi, phiout, rrr, zzz, ppp
  real*8 :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  real*8 :: epsB, divB
  integer idir, nplot, nturns, neqn, nphi, npoint, i, j, iunit, nbad, nok

  external fff, rkqs

  open(11, file='flt.inp')
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*) rmn
  read(11,*) rmx
  read(11,*) zmn
  read(11,*) zmx
  read(11,*) raxis
  read(11,*) zaxis
  read(11,*) per_phi
  per_phi = per_phi*2.*pi
  read(11,*)
  read(11,*)
  read(11,*) idir
  read(11,*) nplot, nturns
  read(11,*) rmnplot
  read(11,*) rmxplot
  read(11,*) phi0, zet0
  close(11)
  phi0 = phi0*2.*pi
  ah = (rmxplot - rmnplot)/(nplot-1)
  if(nplot .gt. 50) STOP 'too large number of plots'

! parameters for ODEINT
  hn1 = 0.1*idir
  hn1min = 0.
  relerr = 1.e-12
  neqn = nequat
  
  do i=1, nplot
     iunit = 50 + (i-1)*idir
     npoint = 0
     phi = phi0

     y(1) = rmnplot + (i-1)*ah
     y(2) = zet0

     r = y(1)
     z = y(2)

!     write(*,*)'surface #',i, r, z
     do j=1,nturns

        phiout = phi + per_phi*idir

        call odeint(y,neqn,phi,phiout,relerr,hn1,hn1min,nok,nbad,fff,rkqs) 

        npoint = npoint + 1
        phi = phiout
        rrr = y(1)  ! cylindic R
        ppp = phi   ! cylindic $\phi$
        zzz = y(2)  ! cylindic Z
        if (ppp .ge. per_phi) ppp = ppp - (int(ppp/per_phi))*per_phi
        if (ppp .lt. 0.) ppp = ppp + (int(abs(ppp)/per_phi) +1)*per_phi
        call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ                &
             ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

        divB = Br/rrr + dBrdR + dBzdZ + dBpdp/rrr
        epsB = divB*rrr/Bp
        write(iunit,101)rrr,zzz,phi,epsB,Br,Bp,Bz
!        write(11,101)r,z,rrr,zzz,phi,epsB,Br,Bp,Bz
        if( y(1).lt.rmn .or. y(1).gt.rmx .or.            &
                    y(2).lt.zmn .or. y(2).gt.zmx ) exit
     enddo
!     write(*,*) npoint, 'points within the domain'
  enddo
 101  format(1000e15.7)
end program flt_phi
! -----------------------------------------------------------------
subroutine fff(phi,y,yp)
  use period_mod

  implicit none

  integer, parameter :: nequat = 2
  real*8 y(nequat),yp(nequat)
  real*8 :: phi, rrr, zzz, ppp
  real*8 :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  rrr = y(1)
  zzz = y(2)
  ppp = phi
  if (ppp .ge. per_phi) ppp = ppp - (int(ppp/per_phi))*per_phi
  if (ppp .lt. 0.) ppp = ppp + (int(abs(ppp)/per_phi) +1)*per_phi

  call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
       ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  yp(1) = Br*rrr/Bp
  yp(2) = Bz*rrr/Bp

  return
end subroutine fff
