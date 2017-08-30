program axis
  use period_mod
  use field_eq_mod, only : psif, btf, rtf
  use fixed_coord, only : fixed, num_fixed
  use psi4root, only : psi_root, R_o,  Z_o, theta, R_x, Z_x, R0, Z0
  use const, only : pi, twopi
  use check_cw
  use inside_of_wall
  implicit none
  integer, parameter :: nequat = 2, nr=200, nz=200
  real*8, parameter :: eps_separ=1.d-4, rhohuge=1.d40 ! should be less than in rtbis
  real*8 :: y(nequat)
  real*8 :: rmn, rmx, zmn, zmx, raxis, zaxis, rmnplot, rmxplot, phi0, zet0
  real*8 :: hn1, hn1min, relerr, ah, r, z, phi, phiout, rrr, zzz, ppp
  real*8 :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  real*8 :: epsB, divB, dist_max, r_maxdist, z_maxdist 
  integer idir, nplot, nturns, neqn, nphi, npoint, i, j, ifile, nbad, nok, icontinue, k
  real*8 :: hr, hz, psipol(nr), r_psimax, z_psimax, psimax, theta_o_x
  real*8 :: a_min, b_min, c_min, tlr_min, f_min, x_min, z_min, brent, x1_min, x2_min, z1_min
  real*8 :: delta_xo, a_root, b_root, rho, tlr_root, rtbis, f_root
  real*8, dimension(:), allocatable :: R_o_x, Z_o_x, psif_o_x, R_0_pfr, Z_0_pfr, psif_0_pfr
  real*8, dimension(:,:), allocatable :: R_coord, Z_coord
  integer, dimension(:,:), allocatable :: ipoi_core, ipoi_pfr, ipoi_sol ! number of point in the list if inside
  real*8 :: R_dpr_o_x, Z_dpr_o_x, R_dpl_o_x, Z_dpl_o_x ! points of intersection (polar centers for PFR)
  real*8 :: r_tmp, z_tmp, r_wall, z_wall, dummy
  real*8 :: R_dpr_1, Z_dpr_1, R_dpr_2, Z_dpr_2, Z_cw_o_x, R_cw_o_x, theta_dpr, R_cw_dpr, Z_cw_dpr
  real*8 :: R_dpl_1, Z_dpl_1, R_dpl_2, Z_dpl_2, theta_dpl, R_cw_dpl, Z_cw_dpl
! mesh resolution:
!  integer, parameter :: iterat=200, nt_core=100, nr_core=20, nr_pfr=30, nt_pfr=20, nr_sol=20 
!
! case: "finer":
!  integer, parameter :: iterat=200, nt_core=200, nr_core_max=150, nr_pfr=30, nt_pfr=20, nr_sol=60 

  integer, parameter :: iterat=200, nt_core=100, nr_core_max=100, nr_pfr=30, nt_pfr=20, nr_sol=60
  
  integer, parameter :: nt_sol=nt_core+2*nt_pfr+1
  real*8, dimension(1:nt_sol) :: R_sep, Z_sep
  real*8, dimension(:), allocatable :: R_0_sol, Z_0_sol, psif_0_sol
  real*8 :: R1, Z1, R2, Z2, psif1, psif2, R_cw_bis, Z_cw_bis, delta_xcw, epsxcw=1.d-4
  real*8 :: vec_norm, th_mesh_dr=-2.d0, th_mesh_dl=-2.d0, th_mesh_co=-1.5d0, rho_mesh_co=-1.d0
  real*8, dimension(1:nt_core) :: delta_theta
  logical :: out_of_cw, wall_touch
  integer :: ibou_b2_pfr, ibou_b2_sol, number_of_points
  integer, parameter :: n_rational = 5 ! for refinement at rational surfaces
  real*8 :: d_rational(n_rational) ! TODO: don't hardcode this
  real*8 :: tau_rad, tau_rad_old
  integer :: nr_core
  
  external fff, rkqs, f_min, f_root

  ! for case "finer"
!  d_rational(1) = .5675d0  ! q=1.5
!  d_rational(2) = .694d0  ! q=2 
!  d_rational(3) = .7585d0  ! q=2.5
!  d_rational(4) = .808d0  ! q=3
!  d_rational(5) = .850d0  ! q=3.5
!  d_rational(6) = .8905d0 ! q=4
!  d_rational(7) = .9413d0  ! q=4.5
!  d_rational(8) = .971d0  ! q=5
!  d_rational(9) = .985d0  ! q=5.5
  !  d_rational(10) = .997d0 ! q=6
  
  d_rational(1) = .335d0  ! q=1.5
  d_rational(2) = .444d0  ! q=2 
  d_rational(3) = .505d0  ! q=2.5
  d_rational(4) = .558d0  ! q=3
  d_rational(5) = .608d0  ! q=3.5
  ! d_rational(6) = .667d0 ! q=4.0
  
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
  if(nplot .gt. 10) STOP 'too large number of plots'
! parameters for ODEINT
  hn1 = 0.1*idir
  hn1min = 0.
  relerr = 1.e-12
  neqn = nequat
  dummy = 0.d0

!!$rmn = 511.d0
!!$rmx = 512.d0
!!$zmn = -343.d0
!!$zmx = -342.d0
!!$
!!$  hr = (rmx -  rmn)/dfloat(nr - 1)
!!$  hz = (zmx -  zmn)/dfloat(nz - 1)  
!!$  psimax = -1.d30
!!$  do i= 1,nz
!!$     zzz = zmn + (i-1)*hz
!!$     ppp = 0.d0
!!$     do j= 1,nr
!!$        rrr = rmn + (j-1)*hr
!!$        call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ                &
!!$             ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ) 
!!$!        psipol(j) = psif
!!$        write(222,102) rrr, zzz, psif
!!$        if(psif .gt. psimax) then
!!$           psimax = psif
!!$           r_psimax = rrr
!!$           z_psimax = zzz
!!$        endif
!!$     enddo
!!$!     write(222,102)(psipol(k),k=200,300) !1,nr) 
!!$  enddo
!!$  close(222)

!stop

!  raxis = r_psimax
!  zaxis = z_psimax

! O-point:
  do i=1, nplot
     ifile = 50 + (i-1)*idir
     npoint = 0
     phi = phi0

     y(1) = raxis ! rmnplot + (i-1)*ah
     y(2) = zaxis ! zet0

     r = y(1)
     z = y(2)
     
     write(*,*) i, raxis, zaxis 

     dist_max = 0.d0
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
!        write(ifile,101)rrr,zzz,phi,epsB,Br,Bp,Bz, psif
        if( y(1).lt.rmn .or. y(1).gt.rmx .or.            &
             y(2).lt.zmn .or. y(2).gt.zmx ) exit

        raxis = raxis + rrr
        zaxis = zaxis + zzz

     enddo
     raxis = raxis/dfloat(nturns+1)
     zaxis = zaxis/dfloat(nturns+1)
     close(ifile)
  enddo
  R_o = raxis
  Z_o = zaxis
  print *, "Axis: "
  write(*,*) raxis, zaxis, psif

! X-point, coordinate-wise minimization (min and minus max):
  x_min = -341.d0  !60.d0
  tlr_min = 1.d-14
  do i=1, iterat
     do j=1,2
!        num_fixed = modulo(j,2) + 1
        if(j .eq. 1) then
           num_fixed = 2
           a_min = 138.d0 
           b_min = 143.d0 
           c_min = 147.d0 
        else
           num_fixed = 1
           a_min = -93.d0 !50.d0
           b_min = -95.d0 !60.d0
           c_min = -97.d0 !70.d0
        endif
           fixed = x_min
           z_min = brent(a_min, b_min, c_min, f_min, tlr_min, x_min)
        if(j.eq.1) then
           x1_min = x_min
           z1_min = -z_min
        else
           x2_min = x_min
        endif
     enddo
!     print *, x1_min, x2_min, z_min
     write(111,101) x1_min, x2_min, z_min

  enddo
  print *, "X-point: "
  print *, x1_min, x2_min, z_min
  R_x = x1_min 
  Z_x = x2_min

! core region, lines of psi=const:
  R0 = R_o ! for f_root (will be changed for PFR)
  Z0 = Z_o
  delta_xo = 2.d0/dfloat(nr_core_max-1)
  allocate(R_o_x(nr_core_max+1), Z_o_x(nr_core_max+1), psif_o_x(nr_core_max+1), &
       R_coord(nt_core,nr_core_max+1), Z_coord(nt_core,nr_core_max+1), ipoi_core(nt_core+1,nr_core_max+1))
! nt_core+1 to close the array ipoint after period

  do i=1,nr_core_max
     if (i==1) then
        tau_rad = 0.0
        tau_rad_old = 0.0
     else
        tau_rad = tau_rad_old + delta_xo/&
             (1.d0+5.d0*sum(exp(-(tau_rad_old-d_rational)**2/8d-3**2)))
        tau_rad_old = tau_rad
     end if
     if(tau_rad > 1.d0) then
        R_o_x(i) = R_x  
        Z_o_x(i) = Z_x 
    else
        R_o_x(i) = R_x + (1.d0-tau_rad)*(R_o - R_x)
        Z_o_x(i) = Z_x + (1.d0-tau_rad)*(Z_o - Z_x)
        
    !   R_o_x(i) = R_x + sqrt(1.d0-tau_rad)*(R_o - R_x)  
    !   Z_o_x(i) = Z_x + sqrt(1.d0-tau_rad)*(Z_o - Z_x)
     end if
     R_coord(1,i) = R_o_x(i)
     Z_coord(1,i) = Z_o_x(i)
     ppp = 0.d0
     rrr = R_o_x(i)
     zzz = Z_o_x(i)
     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     psif_o_x(i) = psif
     write(112,101)R_o_x(i), Z_o_x(i), psif_o_x(i), sqrt((R_coord(1,i)-R_o)**2 + (Z_coord(1,i)-Z_o)**2)
     if(tau_rad > 1.d0) then
        nr_core = i
        exit
     end if
  enddo

  theta_o_x = atan2(Z_x-Z_o,R_x-R_o)
  R_o_x(nr_core+1) = R_o_x(nr_core)
  Z_o_x(nr_core+1) = Z_o_x(nr_core)
  psif_o_x(nr_core+1) = psif_o_x(nr_core)
  R_coord(1,nr_core+1) = R_coord(1,nr_core)
  Z_coord(1,nr_core+1) = Z_coord(1,nr_core)

  rho = sqrt((R_x-R_o)**2 + (Z_x-Z_o)**2)
!  rho = rho*(1.d0 - eps_separ)
  R_o_x(nr_core) = rho*cos(theta_o_x) + R_o ! R_x + eps_separ*(R_x - R_o)
  Z_o_x(nr_core) = rho*sin(theta_o_x) + Z_o ! Z_x + eps_separ*(Z_x - Z_o)
  R_coord(1,nr_core) = R_o_x(nr_core)
  Z_coord(1,nr_core) = Z_o_x(nr_core)
  ppp = 0.d0
  rrr = R_o_x(nr_core)
  zzz = Z_o_x(nr_core)
  call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  psif_o_x(nr_core) = psif
  write(112,101)R_o_x(nr_core), Z_o_x(nr_core), psif_o_x(nr_core), rho 
  close(112)

print*, 'theta_o_x=', theta_o_x
!  delta_theta(:) = twopi/dfloat(nt_core)
! no data overlap after one period assumed:
  do j=1, nt_core/2
     delta_theta(j) = pi*(exp(-th_mesh_co*dfloat(j)/dfloat(nt_core/2)) - 1.d0) &
                                                /(exp(-th_mesh_co) - 1.d0)
  enddo
  i=0
  do j=nt_core/2+1,nt_core-1
     i = nt_core - j + 1
     delta_theta(j) =  delta_theta(j-1) + (delta_theta(i)-delta_theta(i-1)) 
  enddo
  tlr_root = 1.d-7
  do i=2,nr_core !-1
!     ifile = 40 + (i-1)
!     write(112,*)R_coord(1,i),Z_coord(1,i),'1'
     psi_root = psif_o_x(i)
     a_root = 0.d0
     b_root = sqrt((R_coord(1,i+1)-R_o)**2 + (Z_coord(1,i+1)-Z_o)**2)
     theta = theta_o_x
     do j=2,nt_core
        theta = theta_o_x + delta_theta(j-1) ! *(j-1)
        if (theta .ge. twopi) theta = theta - (int(theta/twopi))*twopi
        if (theta .lt. 0.) theta = theta + (int(abs(theta)/twopi) +1)*twopi
        rho = rtbis(f_root, a_root, b_root, tlr_root)
        R_coord(j,i) = rho*cos(theta) + R_o
        Z_coord(j,i) = rho*sin(theta) + Z_o
!        write(ifile,*)R_coord(j,i),Z_coord(j,i),j
        if(i .eq. nr_core) then
           R_sep(j+nt_pfr) = R_coord(j,i)
           Z_sep(j+nt_pfr) = Z_coord(j,i)
        endif
     enddo
  enddo

  open(212,file='points.fmt')
  number_of_points = 0
  ppp = 0.d0
  rrr = R_coord(1,1)
  zzz = Z_coord(1,1)
  call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  write(212,101)R_coord(1,1), Z_coord(1,1), psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp 
  number_of_points =   number_of_points + 1
  ipoi_core(1:nt_core,1) = 1
  do i=2, nr_core
     do j=1, nt_core
        rrr = R_coord(j,i)
        zzz = Z_coord(j,i)
        call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
             ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
        write(212,101)R_coord(j,i), Z_coord(j,i), psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
        number_of_points =   number_of_points + 1 
        ipoi_core(j,i) = number_of_points
     enddo
  enddo
  ipoi_core(nt_core+1,1:nr_core) = ipoi_core(1,1:nr_core) 
  open(211,file='index_core.pnt')
  write(211,*) number_of_points, nr_core, nt_core+1  
  write(211,999) ipoi_core(1:nt_core+1,1)
  do i=2, nr_core
     write(211,999) ipoi_core(1:nt_core+1,i)
  enddo
  close(211)
! --------------------------------------------------------------
  call  precalc_wall
! ------------------------------------------------------------------------------------------
! PFR, lines psi=const:
! -----------------------------------------------
! intersection of the X-O line with CW, R_cw_o_x, Z_cw_o_x:
  call in_out_cw(R_o, Z_o, R_x, Z_x, out_of_cw, R_cw_o_x, Z_cw_o_x )
  delta_xo = 1.d0/dfloat(nr_pfr-1)
  allocate(R_0_pfr(nr_pfr), Z_0_pfr(nr_pfr), psif_0_pfr(nr_pfr))
  deallocate(R_coord, Z_coord)
  allocate(R_coord(-nt_pfr:nt_pfr,nr_pfr), Z_coord(-nt_pfr:nt_pfr,nr_pfr), ipoi_pfr(-nt_pfr:nt_pfr,nr_pfr))
  ipoi_pfr(:,:) = 0 ! all points are outside
  ibou_b2_pfr = nr_pfr
  do i=1,nr_pfr
     R_0_pfr(i) = R_cw_o_x + delta_xo*(R_x - R_cw_o_x)*dfloat(i-1)  
     Z_0_pfr(i) = Z_cw_o_x + delta_xo*(Z_x - Z_cw_o_x)*dfloat(i-1)
     R_coord(0,i) = R_0_pfr(i)
     Z_coord(0,i) = Z_0_pfr(i)
     ppp = 0.d0
     rrr = R_0_pfr(i)
     zzz = Z_0_pfr(i)

     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     psif_0_pfr(i) = psif
     write(114,101)R_0_pfr(i), Z_0_pfr(i), psif_0_pfr(i), sqrt((rrr-R_cw_o_x)**2 + (zzz-Z_cw_o_x)**2)
     if( PointIsInside(rrr, zzz) .and. i.ne.nr_pfr) then
        write(212,101)R_coord(0,i), Z_coord(0,i), psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
        number_of_points = number_of_points + 1 
        ipoi_pfr(0,i) = number_of_points
     else
        wall_touch=.true.
        ipoi_pfr(0,i) = 0
     endif
  enddo
  R_sep(nt_pfr+1) = R_0_pfr(nr_pfr)
  Z_sep(nt_pfr+1) = Z_0_pfr(nr_pfr)
  R_sep(nt_pfr+nt_core+1) = R_sep(nt_pfr+1)
  Z_sep(nt_pfr+nt_core+1) = Z_sep(nt_pfr+1)
  tlr_root = 1.d-7
!--------------------------------------------------------------------
! polar center for the right half of PFR:
  R_dpr_1 = 1.605d2
  Z_dpr_1 = -1.2d2  
  R_dpr_2 = 1.71146d2 !1.6034d2
  Z_dpr_2 = 4.57d0 ! -1.6d1 
  call twostreights(R_o, Z_o, R_x, Z_x, R_dpr_1, Z_dpr_1, R_dpr_2, Z_dpr_2, R_dpr_o_x, Z_dpr_o_x)
  R0 = R_dpr_o_x !R_cw_o_x ! for f_root
  Z0 = Z_dpr_o_x !Z_cw_o_x
  write(11,*) R_x, Z_x
  write(11,*) R_dpr_o_x, Z_dpr_o_x
  write(11,*) R_dpr_1, Z_dpr_1
  theta_dpr = atan2(Z_dpr_1-Z_dpr_o_x,R_dpr_1-R_dpr_o_x)
  do j=1, nt_pfr
     delta_theta(j) = abs(theta_dpr - theta_o_x)*(exp(th_mesh_dr*dfloat(j)/dfloat(nt_pfr)) - 1.d0) &
                                                /(exp(th_mesh_dr) - 1.d0)
  enddo
  do i=2,nr_pfr
     psi_root = psif_0_pfr(i)
!     b_root = sqrt((R_coord(0,i)-R0)**2 + (Z_coord(0,i)-Z0)**2) !- delta_xo*(R_x - R_cw_o_x)*5.d-1
! last term is (small) shift because of possible convexity of the contour \Psi=const we are looking for
     b_root = sqrt((R_coord(0,nr_pfr)-R0)**2 + (Z_coord(0,nr_pfr)-Z0)**2) ! i nehuj vyebyvat'sja
     wall_touch=.false.
     do j=1, nt_pfr
        theta = theta_o_x + delta_theta(j)
        if (theta .ge. twopi) theta = theta - (int(theta/twopi))*twopi
        if (theta .lt. 0.) theta = theta + (int(abs(theta)/twopi) +1)*twopi
! intersection of new theta=const with CW:
        call in_out_cw(R_dpr_o_x, Z_dpr_o_x, rrr, zzz, out_of_cw, R_cw_dpr, Z_cw_dpr, theta)        
        a_root = sqrt((R_cw_dpr-R0)**2 + (Z_cw_dpr-Z0)**2)
        rho = rtbis(f_root, a_root, b_root, tlr_root)
        if(rho .gt. rhohuge) cycle
!        call in_out_cw(R_dpr_o_x, Z_dpr_o_x, rrr, zzz, out_of_cw, r_wall, z_wall, theta) 
        r_tmp = rho*cos(theta) + R0 
        z_tmp = rho*sin(theta) + Z0 
        if(i .eq. nr_pfr) then
           R_sep(nt_pfr-j+1) = r_tmp
           Z_sep(nt_pfr-j+1) = z_tmp
        endif
!!$        if( j .eq. nt_pfr ) then
!!$           R_coord(j,i) = r_tmp
!!$           Z_coord(j,i) = z_tmp
!!$           call field(r_tmp, ppp, z_tmp, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
!!$                ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
!!$           write(212,101)r_tmp, z_tmp, psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
!!$           number_of_points = number_of_points + 1 
!!$           ipoi_pfr(j,i) = -number_of_points
!!$        else
        if( PointIsInside(r_tmp, z_tmp) .and. j.ne.nt_pfr) then
           call field(r_tmp, ppp, z_tmp, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
                ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
           write(212,101)r_tmp, z_tmp, psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
           number_of_points = number_of_points + 1 
           ipoi_pfr(j,i) = number_of_points
           R_coord(j,i) = r_tmp
           Z_coord(j,i) = z_tmp
        else
           wall_touch=.true.
           ipoi_pfr(j,i) = 0
        endif
     enddo
  enddo
! -----------------------------------------------
! polar center for left half of PFR:
  R_dpl_1 = 1.23d+02 
  Z_dpl_1 =-1.10d+02 
  R_dpl_2 = 1.71146d2 
  Z_dpl_2 = 4.75d0 
  call twostreights(R_o, Z_o, R_x, Z_x, R_dpl_1, Z_dpl_1, R_dpl_2, Z_dpl_2, R_dpl_o_x, Z_dpl_o_x)
  R0 = R_dpl_o_x !R_cw_o_x ! for f_root
  Z0 = Z_dpl_o_x !Z_cw_o_x
  write(11,*) R_dpl_o_x, Z_dpl_o_x
  write(11,*) R_dpl_1, Z_dpl_1
  theta_dpl = atan2(Z_dpl_1-Z_dpl_o_x,R_dpl_1-R_dpl_o_x)
!!$  do j=1,nt_pfr
!!$     delta_theta(j) = abs(theta_dpl - theta_o_x)/dfloat(nt_pfr)
!!$  enddo
  do j=1, nt_pfr
     delta_theta(j) = abs(theta_dpl - theta_o_x)*(exp(th_mesh_dl*dfloat(j)/dfloat(nt_pfr)) - 1.d0) &
                                                /(exp(th_mesh_dl) - 1.d0)
  enddo
  do i=2,nr_pfr
     psi_root = psif_0_pfr(i)
!     b_root = sqrt((R_coord(0,i)-R0)**2 + (Z_coord(0,i)-Z0)**2) - delta_xo*(R_x - R_cw_o_x)*5.d-1
     b_root = sqrt((R_coord(0,nr_pfr)-R0)**2 + (Z_coord(0,nr_pfr)-Z0)**2) ! i nehuj vyebyvat'sja
     wall_touch=.false.
     do j=-1,-nt_pfr,-1
        theta = theta_o_x - delta_theta(-j)
        if (theta .ge. twopi) theta = theta - (int(theta/twopi))*twopi
        if (theta .lt. 0.) theta = theta + (int(abs(theta)/twopi) +1)*twopi
! intersection of new theta=const with CW:
        call in_out_cw(R_dpl_o_x, Z_dpl_o_x, rrr, zzz, out_of_cw, R_cw_dpl, Z_cw_dpl, theta)  

        a_root = sqrt((R_cw_dpl-R0)**2 + (Z_cw_dpl-Z0)**2)

        rho = rtbis(f_root, a_root, b_root, tlr_root)
        if(rho .gt. rhohuge) cycle
!        call in_out_cw(R_dpr_o_x, Z_dpr_o_x, rrr, zzz, out_of_cw, r_wall, z_wall, theta) 
        r_tmp = rho*cos(theta) + R0 
        z_tmp = rho*sin(theta) + Z0 
        if(i .eq. nr_pfr) then
!print *,nt_pfr + nt_core -j +1
           R_sep(nt_pfr+nt_core-j+1) = r_tmp
           Z_sep(nt_pfr+nt_core-j+1) = z_tmp
        endif
!!$        if( j .eq. -nt_pfr) then
!!$           call field(r_tmp, ppp, z_tmp, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
!!$                ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
!!$           write(212,101)r_tmp, z_tmp, psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
!!$           number_of_points = number_of_points + 1 
!!$           ipoi_pfr(j,i) = -number_of_points
!!$           R_coord(j,i) = r_tmp
!!$           Z_coord(j,i) = z_tmp
!!$        else
        if( PointIsInside(r_tmp, z_tmp)  .and. j.ne.-nt_pfr) then
           call field(r_tmp, ppp, z_tmp, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
                ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
           write(212,101)r_tmp, z_tmp, psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
           number_of_points = number_of_points + 1 
           ipoi_pfr(j,i) = number_of_points
           R_coord(j,i) = r_tmp
           Z_coord(j,i) = z_tmp
        else
           wall_touch=.true.
           ipoi_pfr(j,i) = 0
        endif
     enddo
  enddo
!!$  number_of_points = number_of_points - 1 ! remove X-point from the list 
!!$  ipoi_pfr(0,nr_pfr) = 0                  ! 
!!$  do i=1,nr_pfr
!!$     do j=-nt_pfr,nt_pfr
!!$        if(ipoi_pfr(j,i) .ne. 0) then
!!$           rrr = R_coord(j,i)
!!$           zzz = Z_coord(j,i)
!!$           call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
!!$                ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
!!$           write(212,101)R_coord(j,i), Z_coord(j,i), psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
!!$        endif
!!$     enddo
!!$  enddo
  ipoi_pfr(0,nr_pfr) = ipoi_core(nt_core+1,nr_core) ! restore X-point
  open(211,file='index_pfr.pnt')
  write(211,*) number_of_points, nr_pfr, nt_pfr
  do i=1, nr_pfr
     write(211,999) ipoi_pfr(-nt_pfr:nt_pfr,i)
  enddo
  close(211)
! -------------------------------------------------------------------------
! SOL, radial mesh (given on bisectrix from X-point)
  R0 = R_sep(nt_pfr+1) !  R_x 
  Z0 = Z_sep(nt_pfr+1) ! Z_x
  R1 = R_sep(nt_pfr)
  Z1 = Z_sep(nt_pfr)  
  R2 = R_sep(nt_pfr+2)
  Z2 = Z_sep(nt_pfr+2) 
  theta = (atan2(Z1-Z0,R1-R0) + atan2(Z2-Z0,R2-R0))*5.d-1
  if (theta .ge. twopi) theta = theta - (int(theta/twopi))*twopi
  if (theta .lt. 0.) theta = theta + (int(abs(theta)/twopi) +1)*twopi
! intersection of bisectrix passed from X-point with CW:
  call in_out_cw(R_x, Z_x, rrr, zzz, out_of_cw, R_cw_bis, Z_cw_bis, theta)  
  delta_xcw = 1.d0/dfloat(nr_sol-1)
  allocate(R_0_sol(nr_sol), Z_0_sol(nr_sol), psif_0_sol(nr_sol))
  deallocate(R_coord, Z_coord)
  allocate(R_coord(1:nt_sol,nr_sol), Z_coord(1:nt_sol,nr_sol), ipoi_sol(1:nt_sol,nr_sol))
  ipoi_sol(:,:) = 0 ! all points are outside
  do i=1,nr_sol
     R_0_sol(i) = R_x - sqrt(delta_xcw*dfloat(i-1))*(R_x - R_cw_bis)  
     Z_0_sol(i) = Z_x - sqrt(delta_xcw*dfloat(i-1))*(Z_x - Z_cw_bis)
     ppp = 0.d0
     rrr = R_0_sol(i)
     zzz = Z_0_sol(i)
!     if( PointIsInside(rrr, zzz) ) ibou_b2_sol = i

     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     psif_0_sol(i) = psif
     write(113,101)R_0_sol(i), Z_0_sol(i), psif_0_sol(i)
  enddo
! push separatrix contour away from X-point along bisectrix
  R_sep(nt_pfr+1) = R_x - epsxcw*(R_x - R_cw_bis)
  Z_sep(nt_pfr+1) = Z_x - epsxcw*(Z_x - Z_cw_bis)
  R_sep(nt_pfr+nt_core+1) = R_x + epsxcw*(R_x - R_cw_bis)
  Z_sep(nt_pfr+nt_core+1) = Z_x + epsxcw*(Z_x - Z_cw_bis)

  ppp = 0.d0
  tht: do j=2,nt_sol-1
     R_coord(j,1) = R_sep(j) ! correct afterwords at X-point
     Z_coord(j,1) = Z_sep(j) ! correct afterwords at X-point
     ipoi_sol(j,1) = 0 ! to prevent writing of known points
     R1 = R_coord(j,1)
     Z1 = Z_coord(j,1)
     call field(R1, ppp, Z1, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ                &
          , dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

     vec_norm = 1.d-1/sqrt(Br**2+Bz**2)
     psif1 = psif
     R2 = R1 - Bz*vec_norm
     Z2 = Z1 + Br*vec_norm
     call field(R2, ppp, Z2, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ                &
          , dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     psif2 = psif
     rh: do i=2, nr_sol
        do k=1,10000

           if((psif2-psif_0_sol(i))*(psif1-psif_0_sol(i)) .lt. 0.d0) go to 11
           if( .not.PointIsInside(R1, Z1) .or. .not.PointIsInside(R2, Z2) ) exit rh
           R1 = R2
           Z1 = Z2
           psif1 = psif2

           call field(R1, ppp, Z1, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ                &
                , dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

           vec_norm = 1.d-1/sqrt(Br**2+Bz**2)
           psif2 = psif
           R2 = R1 - Bz*vec_norm
           Z2 = Z1 + Br*vec_norm
        enddo
        print *,'No point found in SOL'
        stop
11      continue
        R_coord(j,i) = (R1 + R2)*5.d-1     
        Z_coord(j,i) = (Z1 + Z2)*5.d-1   

        R1 = R_coord(j,i)
        Z1 = Z_coord(j,i)
        call field(R1, ppp, Z1, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ                &
             , dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
        if(PointIsInside(R1, Z1)) then
           number_of_points = number_of_points +1
           ipoi_sol(j,i) = number_of_points
           write(212,101)R_coord(j,i), Z_coord(j,i), psif, sqrt(Br**2 + Bp**2 + Bz**2), Bp
        endif
        vec_norm = 1.d-1/sqrt(Br**2+Bz**2)
        psif1 = psif
        R2 = R1 - Bz*vec_norm
        Z2 = Z1 + Br*vec_norm
        
        call field(R2, ppp, Z2, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ                &
             , dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
        psif2 = psif  
     enddo rh
  enddo tht
  write(212,*) btf, rtf
  close(212)

! separatrix points:
  ipoi_sol(1:nt_pfr+1,1) = ipoi_pfr(nt_pfr:0:-1,nr_pfr)
  ipoi_sol(nt_pfr+2:nt_pfr+1+nt_core,1) = ipoi_core(2:nt_core+1,nr_core) 
  ipoi_sol(nt_pfr+2+nt_core:nt_sol,1) = ipoi_pfr(-1:-nt_pfr:-1,nr_pfr)
  open(211,file='index_sol.pnt')
  write(211,*) number_of_points, nr_sol, nt_sol
  do i=1, nr_sol
     write(211,999) ipoi_sol(1:nt_sol,i)
  enddo
  close(211)

101 format(1000(e21.14,x))
102 format(1000(e15.8,x))
999 format(1000(i6,x))

end program axis
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
! -----------------------------------------------------------------
function f_min(s)
  use fixed_coord, only : fixed, num_fixed
  use field_eq_mod, only : psif
  implicit none
  real*8 :: f_min, s
  real*8 :: rrr, zzz, ppp
  real*8 :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

  ppp = 0.d0
  if(num_fixed .eq. 1) then
     rrr = fixed
     zzz = s
     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     f_min = psif
  else
     rrr = s
     zzz = fixed   
     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     f_min = -psif
  endif
!print *,num_fixed, rrr, zzz, psif

  return
end function f_min
! --------------------------------------------
function f_root(rho)
  use field_eq_mod, only : psif
  use psi4root, only : psi_root, R0,  Z0, theta
  implicit none
  real*8 :: f_root, rho
  real*8 :: rrr, zzz, ppp
  real*8 :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

  ppp = 0.d0
  rrr = rho*cos(theta) + R0
  zzz = rho*sin(theta) + Z0
  call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
       ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  f_root = psif - psi_root
!print *, rho, theta, rrr, zzz, psif, psi_root
  return
end function f_root
