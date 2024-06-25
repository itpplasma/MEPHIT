  module rhs_surf_mod
    double precision :: dr_dphi, dz_dphi
  end module rhs_surf_mod
!
  module field_line_integration_mod
    ! If a circular mesh with large aspect ratio is used, this gives the scale factor.
    ! In this case, o_point and x_point are input instead of output variables.
    ! Set to 0 for regular toroidal geometry.
    integer :: circ_mesh_scale = 0
    ! set this to true to let theta start at the line between O- and X-Point
    logical :: theta0_at_xpoint = .true.
    double precision, dimension(2) :: o_point, x_point, theta_axis
    double precision :: theta0
  end module field_line_integration_mod
!
  subroutine field_line_integration_for_SYNCH(nstep,nsurfmax,nlabel,ntheta,    &
                                              rmn,rmx,zmn,zmx,raxis,zaxis,     &
                                              rbeg,rsmall,qsaf,psisurf,phitor, &
                                              circumf,R_st,Z_st,bmod_st,sqgnorm_st)
!
  use field_eq_mod,  only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2, &
                            icall_eq,nrad,nzet,rad,zet,rtf,btf
  !use theta_rz_mod, only : nsqp,hsqpsi,spllabel
  use rhs_surf_mod, only : dr_dphi, dz_dphi
  use field_line_integration_mod, only: circ_mesh_scale, o_point, x_point, &
       theta0_at_xpoint, theta_axis, theta0
  use magdata_in_symfluxcoor_mod, only: btor, rbig
  use field_sub, only : field_eq
!
  implicit none
!
  integer, parameter :: neq=5
  integer, parameter :: niter_axis=20  !number of iterations for finding axis
  integer, parameter :: niter=50    !number of iterations for Newton method
  integer, parameter :: nstep_min=10   !minimum number of steps
!
  integer, intent(in) :: nstep,nsurfmax,nlabel,ntheta
  integer :: i,j,nsurf,nmap,isurf,iter
!
  double precision, parameter :: pi = 3.14159265358979d0
  double precision, intent(out) :: rmn,rmx,zmn,zmx,raxis,zaxis
  double precision :: relerr,phiout
  double precision :: phi,rrr,ppp,zzz
  double precision :: aiota,hbr
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: psi_axis,h,sig,sig_start,sig_end,phi_sep,sigma,min_d,new_d,r_sep, alpha, beta
!
  double precision, dimension(2) :: prev_ymet, ymet_axis

  double precision, dimension(neq)           :: ymet
  double precision, dimension(nlabel), intent(out) :: rbeg,rsmall,qsaf,psisurf,phitor,circumf
  double precision, dimension(nlabel,ntheta), intent(out) :: R_st,Z_st,bmod_st,sqgnorm_st
!
  external :: rhs_axis, rhs_surf  !, rhs_surf_theta
  double precision, external :: cross_2d_sign
!
!
  relerr = 1.d-9
  nmap=10         !number of maps for finding magnetic axis
!
! Computation box:
  rmn=rad(1)
  rmx=rad(nrad)
  zmn=zet(1)
  zmx=zet(nzet)
!
  if (circ_mesh_scale /= 0) then
     raxis = o_point(1)
     zaxis = o_point(2)
  else
     rrr=0.5d0*(rmn+rmx)
     ppp=0.d0
     zzz=0.5d0*(zmn+zmx)
     !
     ! Search for the magnetic axis
     !
     phi=0.d0
     phiout=2.d0*pi
     ymet=0.d0
     ymet(1)=rrr
     ymet(2)=zzz
     do iter=1,niter_axis
        ymet(3:4)=0.d0
        do i=1,nmap
           call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_axis)
        enddo
        ymet(1:2)=ymet(3:4)/(phiout-phi)/dfloat(nmap)
     enddo
     raxis=ymet(1)
     zaxis=ymet(2)
  end if
  o_point = [raxis, zaxis]
  print *, 'O-point found ', raxis, zaxis
!
  call field_eq(raxis,ppp,zaxis,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  psi_axis=psif
  print *,'toroidal field = ',btf
  btor = btf
  rbig = rtf
!
! End of search for the magnetic axis
!
  hbr=(rmx-raxis)/nsurfmax
!
  rrr=raxis+hbr
  zzz=zaxis
  call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Direction of poloidal field:
  sigma=sign(1.d0,Bz*Bp)
!
!-------------------------------------------------------------------------------
!
! Scan of flux surfaces
!
  h=2.d0*pi/nstep_min
!
  if (circ_mesh_scale /= 0) then
     h = h / dble(circ_mesh_scale)  ! for odeint integration interval
     r_sep = x_point(1)
     nsurf = floor((r_sep - raxis) / hbr) - 1
  else
     surf: do isurf=1,nsurfmax
        phi=0.d0
        phiout=h
        ymet(1)=raxis+hbr*isurf
        ymet(2)=zaxis
        ymet(3)=0.d0
        ymet(4)=0.d0
        ymet(5)=0.d0
        call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
        sig=ymet(2)-zaxis
        do while(sig*(ymet(2)-zaxis).gt.0.d0)
           call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
           if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
                ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
              nsurf=isurf-1
              exit surf
           endif
        enddo
        sig=ymet(2)-zaxis
        do while(sig*(ymet(2)-zaxis).gt.0.d0)
           call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
           if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
                ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
              nsurf=isurf-1
              exit surf
           endif
        enddo
!
        r_sep = raxis+hbr*isurf
     enddo surf
!
     nsurf=nsurf-1  !last point is bad, remove it
  end if
!
  print *,'Separatrix found'
!
!------------------------------------------------------------------------------
!
! Re-define start points step size in R for data storage
  hbr=hbr*dfloat(nsurf)/dfloat(nlabel)
!
! find x-point
  if (circ_mesh_scale /= 0) then
  else
     phi = 0.d0
     phiout = h * sigma
!
     ymet(1) = raxis+hbr*dfloat(nlabel)
     ymet(2) = zaxis
     ymet(3) = 0.d0
     ymet(4) = 0.d0
     ymet(5) = 0.d0
     theta_axis = ymet(1:2) - [raxis, zaxis]

     prev_ymet = ymet(1:2)
     min_d = huge(0.d0)

     sig_start = 1.d0
     sig_end = 1.d0
     do while (sig_start >= sig_end)
        call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
        new_d = sum((prev_ymet - ymet(1:2)) * (prev_ymet - ymet(1:2)))
        if (new_d < min_d) then
           min_d  = new_d
           x_point = ymet(1:2)
        end if
        prev_ymet = ymet(1:2)
        sig_start = sig_end
        sig_end = cross_2d_sign(theta_axis, ymet(1:2) - [raxis, zaxis])
     enddo
  end if
!
  print *,'X-Point found', x_point
  if (theta0_at_xpoint) then
    theta_axis = [x_point(1) - raxis, x_point(2) - zaxis]
    !theta0 = modulo(atan2(theta_axis(2),atan2(theta_axis(1))),2.d0*pi)

  else
    theta_axis = [hbr*dfloat(nlabel), 0.d0]
  end if
  theta0 = atan2(theta_axis(2), theta_axis(1))
!------------------------------------------------------------------------------
!
! Computation of flux functions: effective radius, safety factor, poloidal and toroidal fluxes, circumference
!
  do isurf=1,nlabel
    phi=0.d0
    phiout=h
    phi_sep=phi
    ymet(1:2)=[raxis, zaxis] + theta_axis * (isurf * 1.d0) / nlabel
    ymet(3)=0.d0
    ymet(4)=0.d0
    ymet(5)=0.d0

    sig_start = 1.d0 * sigma
    sig_end = 1.d0 * sigma
    do while (sig_start * sigma >= sig_end * sigma)
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
      sig_start = sig_end
      sig_end = cross_2d_sign(theta_axis, ymet(1:2) - [raxis, zaxis])
    enddo

! Newton method
    do iter=1,niter
      ymet_axis = ymet(1:2) - [raxis, zaxis]
      alpha = atan2(dz_dphi, dr_dphi) - theta0
      beta = atan2(ymet_axis(2), ymet_axis(1)) - theta0
!
      phiout = norm2(ymet_axis) * abs(sin(beta) / sin(alpha)) / norm2([dr_dphi, dz_dphi]) &
           * cross_2d_sign(ymet_axis, theta_axis) * sigma
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
    enddo
!
    aiota=2.d0*pi/phi_sep
    rrr=ymet(1)
    zzz=ymet(2)
!
    call  field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    rbeg(isurf) = hypot(rrr - raxis, zzz - zaxis)
    rsmall(isurf)=sqrt(abs(ymet(3))/pi)
    qsaf(isurf) = sigma / aiota
    psisurf(isurf)=psif-psi_axis
    phitor(isurf)=ymet(4)/(2.d0*pi)
    circumf(isurf)=ymet(5)
  enddo
!
  print *,'Flux functions done'
!
!------------------------------------------------------------------------------
!
! Compute 2D functions:
!
  do isurf=1,nlabel
    phi = 0.d0
    phiout = 2.d0*pi*qsaf(isurf)/ntheta
!
    ymet(1:2)=[raxis, zaxis] + theta_axis * (isurf * 1.d0) / nlabel
    ymet(3) = 0.d0
    ymet(4) = 0.d0
    ymet(5) = 0.d0
!
    do j=1,ntheta
!
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
!
      rrr=ymet(1)
      zzz=ymet(2)
!
      call field_eq(rrr,phi,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
                   ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
      R_st(isurf,j) = rrr
      Z_st(isurf,j) = zzz
      bmod_st(isurf,j) = sqrt(Br**2+Bp**2+Bz**2)
      sqgnorm_st(isurf,j) = rrr/abs(Bp)
    enddo
!
  enddo
!
  print *,'2D functions done'
!
!-------------------------------------------------------------------------------
!
  end subroutine field_line_integration_for_SYNCH
! -----------------------------------------------------------------
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_axis(phi,y,dy)
  use field_sub, only : field_eq
!
  implicit none
!
  integer, parameter :: ndim = 5
!
  double precision, dimension(ndim) :: y,dy
  double precision :: R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                      dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
  R=y(1)
  Z=y(2)
!
  call field_eq(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  dy(1)=Br*R/Bp
  dy(2)=Bz*R/Bp
  dy(3)=y(1)
  dy(4)=y(2)
  dy(5)=0.d0
!
  return
  end subroutine rhs_axis
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!-------------------------------------------------------------------------------
!
  subroutine rhs_surf(phi,y,dy)
!
  use rhs_surf_mod , only: dr_dphi, dz_dphi
  use field_sub, only : field_eq
!
  implicit none
!
  integer, parameter :: ndim = 5
!
  double precision, dimension(ndim) :: y,dy
  double precision :: R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                      dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
  R=y(1)
  Z=y(2)
!
  call field_eq(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  dy(1)=Br*R/Bp
  dy(2)=Bz*R/Bp
  dy(3)=y(1)*dy(2)
  dy(4)=y(1)*y(2)*Br
  dy(5)=hypot(dy(1), dy(2))
!
  dr_dphi=dy(1)
  dz_dphi=dy(2)
!
  return
  end subroutine rhs_surf

pure double precision function cross_2d_sign(a, b)
  ! compute the sign of the 2d cross product a x b
  double precision, intent(in), dimension(2) :: a, b

  cross_2d_sign = sign(1.d0, a(1) * b(2) - a(2) * b(1))
end function cross_2d_sign
