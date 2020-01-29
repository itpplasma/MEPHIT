module collis_mod
  use constants, only : nsorts
  double precision, dimension(nsorts) :: dpp, dhh, fpeff
  double precision, parameter :: dlam_step=5.d-2, dlam_max=5.d-1
  double precision, parameter, dimension(nsorts) :: pmin =(/ 1.5d-1, 1.0d-2 /)
  double precision, dimension(nsorts,nsorts) :: velrat, coulog
  double precision, dimension(nsorts) :: vth, frecol
  logical :: Lorentz=.false., no_cons_laws=.true.
end module collis_mod
! --------------------------------------------------------------------
subroutine collis_time(dens_cl, temp_cl, vel_cl, v3, vt, is_test, dt_coll, col_freq, col_freq_e)
  use constants, only : nsorts
  use collis_mod
  implicit none
  integer ::  is_test, is_field, isort1st
  real*8, dimension(nsorts) :: dens_cl, temp_cl, vel_cl
  real*8 :: v3, vt, dt_coll, v3mrf, vmod, u, col_freq, col_freq_e

  call velrat_vth_lcoul(temp_cl, dens_cl) 
  isort1st = 1
!  if(is_test.eq.1) isort1st = 1
  do is_field = nsorts, isort1st,-1  
     v3mrf = v3 - vel_cl(is_field) !  in moving reference frame
     vmod = sqrt(v3mrf**2 + vt**2)
     u = max(vmod/vth(is_test), pmin(is_test))
     call coleff(is_test, is_field, u)
  enddo
  col_freq = frecol(is_test)
  dt_coll = dlam_step**2/(2.d0*dhh(nsorts))/col_freq
  col_freq_e = frecol(1) ! for mesostep
  return
end subroutine collis_time
! --------------------------------------------------------------------
subroutine collis_step(vel_cl,v3,vt,weight,w_c_V, w_c_T,is_test,dt_fact,cg_rng,deltat,deltav)
  use constants, only : nsorts, sigma_random, one3rd, nsorts, amass, erg2ev, rmass
  use collis_mod
  implicit none
  integer ::  is_test, is_field, isort1st
  real*8, dimension(nsorts) :: vel_cl, deltat, deltav, dlambda, du1st, w_c_V, w_c_T
  real*8 :: v3, vt, dt_fact, v3mrf, vmod, u, dtau, weight
  real*8 :: alambda, v3_n, vt_n, delta_p2
  real*8, parameter :: epsv = 1.d-4, alammx = 9.9d-1
  real*8 :: cg_rng(0:5), ran_u_sq3, ran_pm1

  dtau = dt_fact*frecol(is_test)
  isort1st = 2 ! attention, i-e collisions switched off
  if(is_test.eq.1) isort1st = 1
  do is_field = nsorts, isort1st,-1  
     v3mrf = v3 - vel_cl(is_field) !  in moving reference frame
     vmod = sqrt(v3mrf**2 + vt**2)
     u = vmod/vth(is_test)
     alambda = v3mrf/max(vmod,epsv)

     dlambda(is_field) = sqrt(2.d0*(1.d0- alambda**2)*dhh(is_field)*dtau)*ran_pm1(cg_rng)*sigma_random    &
               - 2.d0*alambda*dhh(is_field)*dtau

     if(abs(dlambda(is_field)) .gt. dlam_max) dlambda(is_field) = dlam_max*sign(1.d0,dlambda(is_field))

     du1st(is_field) = sqrt(2.d0*dpp(is_field)*dtau)*ran_u_sq3(cg_rng) + fpeff(is_field)*dtau

     if(Lorentz) du1st(is_field) = 0.d0 ! no change in |v|

     u = abs(u + du1st(is_field)) ! reflection from 0
     vmod = u*vth(is_test)
     alambda = alambda + dlambda(is_field)
     if(alambda .gt. alammx) alambda = 2.d0*alammx - alambda 
     if(alambda .lt. -alammx) alambda = -2.d0*alammx - alambda
     v3mrf = alambda*vmod 
     vt_n = sqrt(vmod**2 - v3mrf**2) 
     v3_n = v3mrf + vel_cl(is_field) 

     deltav(is_field) = - (v3_n - v3)*weight/w_c_V(is_field)*rmass(is_test,is_field)
     delta_p2 = (v3_n**2 + vt_n**2 - (v3**2 + vt**2))*weight
     deltat(is_field) = - ( delta_p2*amass(is_test)                            &
                          + amass(is_field)*w_c_V(is_field)*deltav(is_field)   &
                          * (deltav(is_field) + 2.d0*vel_cl(is_field))         &
                           )/w_c_T(is_field)*erg2ev*one3rd
     v3 = v3_n
     vt = vt_n
     if(Lorentz .or. no_cons_laws) then ! no evolution of profiles
        deltav(is_field) = 0.d0
        deltat(is_field) = 0.d0
     endif
  enddo

  return
end subroutine collis_step
! --------------------------------------------------------------------
subroutine coleff(itest, ifield, p)
!
!  Computes the local values of dimensionless contravariant components
!  of collisional diffusion tensor and friction force for nonrelativistic
!  plasma.
!
!     Input variables:
!        formal: itest, ifield - sorts of colliding particles (alpha & betha in DK)
!                 p        - dimensionless momentum module p/p_T
!        common: 
!                velrat  - ratio of test species thermal velocity to
!                           background species thermal velocity
!     Output variables:
!        formal: dpp     - dimensionless momentum module diffusion
!                          coefficient
!                dhh     - dimensionless pitch angle diffusion coeff.
!                fpeff   - effective dimensionless drag force (prop. to linear
!                          deviation in Fokker-Planck eq.)
!
  use constants, only: amass, charge, rmass, c_c_2
  use collis_mod
  implicit none
  real*8 :: p, xbeta, a, b, c, a1, f1, f2, pm2, pm
  integer :: itest, ifield
  
  xbeta = p*velrat(itest,ifield)

  call onseff(xbeta, a, b, c, a1)

  pm = 1.d0/p
  pm2 = pm**2
  f1 = c_c_2(ifield,itest)    ! *coulog(itest,ifield)/coulog(itest,itest) attention
  f2 = f1*5.d-1*velrat(ifield,itest)**2
  dpp(ifield) = f2*a*pm*pm2    
  fpeff(ifield) = f2*(a1*xbeta - a)*pm2*pm2 - f1*rmass(itest,ifield)*b*pm2  
  dhh(ifield) = f1*velrat(itest,ifield)*c*pm2

  return
end subroutine coleff
!----------------------------------------------------------------------------------------------------------------
subroutine onseff(w, a, b, c, a1)
!  see Dnestrovskii & Kostomarov book
  implicit none
  integer :: i, j
  real*8 :: w, a, b, c, ex, a1, hw, u, dw1, dw2, hwm1, w2, w3, wmx, wmn
  real*8, parameter :: oneoversqp = 0.56418958354775628695d0 ! \frac{1}{\sqrt{\pi}}
  real*8, parameter :: frac43sp = 0.75225277806367504926d0 ! \frac{4}{3\sqrt{\pi}} 
  real*8, parameter :: frac13sp = 0.18806319451591876232d0 ! \frac{1}{3\sqrt{\pi}}
  real*8, parameter :: wmin=1.d-2, wmax = 6.d0, epsw=1.d-5
  integer, parameter :: npoint=10001
  real*8, dimension(npoint) :: ay, by, cy, a1y   
  logical :: firstcall=.true.
  save firstcall, ay, by, cy, a1y, hw, hwm1, wmx, wmn

  if(firstcall) then
     firstcall = .false.
     hw = (wmax - wmin)/dble(npoint-1)
     hwm1 = 1.d0/hw
     do i=1,npoint
        u = wmin + (i-1)*hw
        ex = exp(-u**2)
        ay(i) = erf(u) - 2.d0*ex*u*oneoversqp
        a1y(i) = 4.d0*ex*u**2*oneoversqp
        by(i) = ay(i)
        cy(i) = oneoversqp*ex + (2.d0/u - 1.d0/u**3)*2.5d-1*ay(i)
     enddo
     wmx = wmax - epsw
     wmn = wmin + epsw
  endif

  w2 = w**2
  w3 = w2*w
  if(w .lt. wmn) then
     a = frac43sp*w3
     a1 = 4.d0*(1.d0-w2)*w2*oneoversqp
     b = a
     c = frac13sp*(2.d0- w2)
  elseif(w .gt. wmx) then
     ex = 0.d0 ! exp(-w2)
     a = 1.d0 ! erf(w) - 2.*ex*w*oneoversqp
     a1 = 0.d0 ! 4*ex*w2*oneoversqp
     b = a
     c = (2.d0/w - 1.d0/w3)*2.5d-1*a ! + oneoversqp*ex
  else
     dw1 = (w - wmin)*hwm1
     j = int(dw1) + 1
     dw2 = dw1 - j + 1
     a = ay(j) + (ay(j+1) - ay(j))*dw2 
     a1 = a1y(j) + (a1y(j+1) - a1y(j))*dw2  
     b = a
     c = cy(j) + (cy(j+1) - cy(j))*dw2 
  endif
  return
end subroutine onseff
! ------------------------------------------------------------------------
subroutine velrat_vth_lcoul(temp_l, dens_l)
!   Input variables:
!        formal: temp_l  - temperatures of all species, eV
!                dens_l  - densities of all species, cm-3     
!    Output variables:
!        common: vth(i)  -     thermal velocity of species (i),
!                              vth(i) = sqrt(2T(i)/m(i))
!                frecol(i) - collision frequency of species (i),
!                         frecol=4*Pi*n*e**4*Lambda/(sqrt(2)*m**2*vth**3)
!                velrat(j,i) - ratio of thermal velocity of species j to
!                              species i thermal velocity

  use constants
  use collis_mod
  implicit none
  real*8 :: dens, ch2_6
  real*8, dimension(nsorts) :: temp_l, dens_l 
  real*8, dimension(nsorts) :: vth02, frecol0
  integer i,j
  logical :: firststep=.true.
  save firststep, ch2_6, vth02, frecol0
!  external clog

  if(firststep) then
     firststep = .false.
     do j=1,nsorts
        do i=1,nsorts
           rmass(i,j) = amass(i)/amass(j)
           c_c_2(i,j) = (charge(i)/charge(j))**2
        enddo
     enddo
     vth02(1:nsorts) = 2.d0*ev2erg/amass(1:nsorts) 
     frecol0(1:nsorts) = 4.d0*pi*(echarge*charge(1:nsorts))**4/(amass(1:nsorts)**2)
     ch2_6 = charge(2)**6*2.d0
  endif
  dens = sum(dens_l(2:nsorts)*charge(2:nsorts))  ! ions used for Spitzer problem only

  coulog(1,1) = 24.d0 - log(dens/temp_l(1)**2)*5.d-1
  coulog(1,2) = coulog(1,1)
  coulog(2,1) = coulog(1,1)
  coulog(2,2) = 23.d0 - log(ch2_6*dens/temp_l(2)**3)*5.d-1
  do i=1,nsorts
!     coulog(i,i) = clog(i,i,temp_l(i),dens) !interpolation
     vth(i) = sqrt(temp_l(i)*vth02(i))
     frecol(i) = frecol0(i)*dens*coulog(i,i)/(vth(i)**3) 
  enddo


  do i=1,nsorts
     do j=1,nsorts
        velrat(j,i) = vth(j)/vth(i) ! sqrt( amass(i)*temp_l(j) / (temp_l(i)*amass(j)) )
     enddo
  enddo
end subroutine velrat_vth_lcoul
! ------------------------------------------------------------------------------------------------
!!$function clog(itest, ifield, temp_f, dens_f)
!!$! interpolation of Coulomb logarithm
!!$  use constants, only : charge, nsorts
!!$  implicit none
!!$
!!$  real*8, parameter :: tmin=1.d0, tmax=2.d4, dmin=1.e11, dmax=1.e15
!!$  integer, parameter :: ndens_clog = 201, ntemp_clog = 501
!!$
!!$  real*8, dimension(nsorts, ndens_clog, ntemp_clog ) :: clog_points
!!$  integer :: itest, ifield, j, k, isort, jt, jn
!!$  real*8 :: temp_f, htemp, htempm1, temp_p, dt1, dt2, clog, ch2_6
!!$  real*8 :: dens_f, hdens, hdensm1, dens_p, dn1, dn2, c1, c2
!!$  logical :: firstcall=.true.
!!$  save firstcall, htempm1, htemp, hdens, hdensm1, clog_points, ch2_6
!!$
!!$  if(firstcall) then
!!$     firstcall = .false.
!!$     ch2_6 = charge(2)**6*2.d0
!!$     htemp = (tmax - tmin)/dble(ntemp_clog-1)
!!$     htempm1 = 1.d0/htemp
!!$     hdens = (dmax - dmin)/dble(ndens_clog-1)
!!$     hdensm1 = 1.d0/hdens
!!$
!!$     do k=1, ntemp_clog
!!$        temp_p = tmin + dble(k-1)*htemp
!!$        do j=1, ndens_clog 
!!$           dens_p = dmin + dble(j-1)*hdens
!!$           clog_points(1,j,k) = 24.d0- log(sqrt(dens_p)/temp_p)
!!$           clog_points(2,j,k) = 23.d0- log(charge(2)**2/temp_p*sqrt(2.d0*dens_p*charge(2)**2/temp_p))
!!$        enddo
!!$     enddo
!!$  end if
!!$
!!$  if(itest.eq.2 .and. ifield.eq.2) then
!!$     isort = 2
!!$  else
!!$     isort = 1
!!$  endif
!!$
!!$!  if(temp_f.le.tmin .or. temp_f.ge.tmax .or. dens_f.le.dmin .or. dens_f.ge.dmax) then
!!$       if(isort .eq. 1) then
!!$          clog = 24.d0- log(dens_f/temp_f**2)*5.d-1
!!$       else if(isort .eq. 2) then
!!$          clog = 23.d0- log(ch2_6*dens_f/temp_f**3)*5.d-1
!!$       else
!!$          print *,'Wrong sort of particles'
!!$          stop
!!$       endif
!  else
!!$     dt1 = (temp_f - tmin)*htempm1
!!$     jt = int(dt1) + 1
!!$     dn1 = (dens_f - dmin)*hdensm1
!!$     jn = int(dn1) + 1
!!$     dt2 = dt1 - jt + 1
!!$     dn2 = dn1 - jn + 1
!!$     c1 = clog_points(isort,jn,jt) + (clog_points(isort,jn,jt+1) - clog_points(isort,jn,jt))*dt2 
!!$     c2 = clog_points(isort,jn+1,jt) + (clog_points(isort,jn+1,jt+1) - clog_points(isort,jn+1,jt))*dt2
!!$     clog = c1 + (c2 - c1)*dn2
!!$!     clog = clog_points(isort,jn,jt)
!!$  endif
!!$  return
!!$
!!$end function clog
