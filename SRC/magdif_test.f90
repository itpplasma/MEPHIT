program magdif_test

  use from_nrtype
  use constants,     only : pi, clight
  use mesh_mod,      only : ntri, npoint, mesh_point, mesh_element,&
                            mesh_element_rmp, knot, triangle, bphicovar
  use magdif,        only : n, presn, curr, assemble_sparse, assemble_system,&
       common_triangles, assign_curr3, assign_currents, dpsidr,&
       get_edge_coord, get_source, pres0, psi, B2avg, nt_core, nr_max, gen_B2avg,&
       Bnflux, Bnp, set_Bnflux_notheta, q, init_safety_factor, assemble_system2
  use sparse_mod,    only : sparse_matmul, sparse_solve
  use for_macrostep, only : t_min, d_min

  implicit none


  integer :: k ! counter
  real(dp) :: dummy_re, dummy_im ! buffer variables
  real(dp) :: dummy8(8)
  real(dp) :: psimin, psimax, dens, temp
  real(dp), parameter ::  ti0=3.d3, di0 = 5.d13 ! from readcarre_m.f90 TODO read

  !real(dp), parameter ::  ti0=3.d2, di0 = 5.d12 ! from readcarre_m.f90 TODO read
  print *, "MAGDIF TEST START"

  n = 2 ! TODO: read this from config file


  open(1,file='points.dat',form='unformatted')
  read(1) npoint
  print *, npoint
  allocate(mesh_point(npoint))
  allocate(presn(npoint))
  read(1) mesh_point
  close(1)

  open(1,file='triangles.dat',form='unformatted')
  read(1) ntri
  print *, ntri
  allocate(mesh_element(ntri))
  allocate(mesh_element_rmp(ntri))
  allocate(Bnflux(ntri,3))
  allocate(Bnp(ntri))
  allocate(curr(ntri,3))
  read(1) mesh_element
  read(1) bphicovar!, fulltime, time_l
  close(1)
  open(1,file='VACFIELD/hpsi_vac.dat')
!open(1,file='hpsi.dat')
  do k=1,ntri
     read (1,*) dummy_re,dummy_im
!     mesh_element_rmp(k)%bnorm_vac=cmplx(dummy_re,dummy_im,dp)
mesh_element_rmp(k)%bnorm_vac = 1d0 ! check with constant source
!mesh_element_rmp(k)%bnorm_vac = 3.d0*172.74467899999999*abs(bphicovar) &
!                  /sum(mesh_point(mesh_element(k)%i_knot(:))%rcoord**2 &
!                  *mesh_point(mesh_element(k)%i_knot(:))%b_mod) ! check 2
  end do
  close(1)

  open(1,file='VACFIELD/Bn_flux.dat')
  do k=1,ntri
     read(1,*) dummy8
!     Bnflux(k,1) = dummy8(1) + (0d0,1d0)*dummy8(2)
!     Bnflux(k,2) = dummy8(3) + (0d0,1d0)*dummy8(4)
!     Bnflux(k,3) = dummy8(5) + (0d0,1d0)*dummy8(6)
!     Bnp(k) = dummy8(7) + (0d0,1d0)*dummy8(8)

call set_Bnflux_notheta(k)
!print *, Bnflux(k,1)

     if(abs((sum(Bnflux(k,:))+&
          (0d0,1d0)*n*Bnp(k)*mesh_element(k)%det_3/2d0)/sum(Bnflux(k,:))) > 1d-10) then
        stop 'Bnvac not divergence free'
     end if
  end do
  close(1)


open(1,file='Bn_flux_test.dat')
do k=1,ntri
write(1,*)  real(Bnflux(k,1)), aimag(Bnflux(k,1)),&
         real(Bnflux(k,2)), aimag(Bnflux(k,2)),&
         real(Bnflux(k,3)), aimag(Bnflux(k,3)),&
         real(Bnp(k)), aimag(Bnp(k))
end do
close(1)

curr = 0d0

  psimin = minval(mesh_point%psi_pol)
  psimax = maxval(mesh_point%psi_pol)

  print *, 'psimin = ', psimin, '  psimax = ', psimax

  allocate(pres0(nr_max+1))
  allocate(psi(nr_max+1))

  psi(1) = mesh_point(1)%psi_pol
  dens = (psi(1) - psimin)/psimax*di0 + d_min !dens
  temp = (psi(1) - psimin)/psimax*ti0 + t_min !temp
  pres0(1) = dens*temp*1.602176d-12
  do k=1,nr_max
     psi(k+1) = mesh_point(1 + (k-1)*nt_core + 1)%psi_pol
     dens = (psi(k+1) - psimin)/psimax*di0 + d_min !dens
     temp = (psi(k+1) - psimin)/psimax*ti0 + t_min !temp
     pres0(k+1) = dens*temp*1.602176d-12
  end do

  allocate(B2avg(nr_max))
  allocate(q(nr_max))

  call init_safety_factor

  open(1,file='q.out')
  do k=1,nr_max
     write(1,*) (psi(k+1)+psi(k))/2d0, q(k)
  end do
  close(1)

!  call test_fdm_simple
!  call test_triangle_strip
  call test_pressure_profile
!  call test_current
  call test_current_new
!  call test_dpsidr

contains

  subroutine test_fdm_simple
    integer, parameter :: ns = 10
    real(dp) :: ds
    real(dp), dimension(ns) :: s

    integer :: m
    complex(dp) :: qmn, pmn
    complex(dp), dimension(ns) :: qn, pn, a, b, c, d, du, x

    integer :: k, nz
    integer, dimension(2*ns) :: irow, icol
    complex(dp), dimension(2*ns) :: aval

    ds = 2d0*pi/ns

    s = (/(k*ds, k=0,ns-1)/)
    m = 1

    qmn = 1d0
    pmn = 1d0/((0d0,1d0)*(m+n))
    qn = qmn*exp((0d0,1d0)*m*s)
    pn = pmn*exp((0d0,1d0)*m*s)

    a = 1d0/ds
    b = (0d0,1d0)*n
    c = 1d0

    call assemble_system(ns, a, b, c, qn, d, du, x)
    call assemble_sparse(ns, d, du, nz, irow, icol, aval)

    call sparse_solve(ns,ns,nz,irow,icol,aval,x)

    open(1,file='test_fdm_simple_x.out', recl=1024)
    do k=1,ns
       write(1,*) real(x(k)), aimag(x(k))
    end do
    close(1)
    open(1,file='test_fdm_simple_pn.out', recl=1024)
    do k=1,ns
       write(1,*) real(pn(k)), aimag(pn(k))
    end do
    close(1)
  end subroutine test_fdm_simple

  subroutine test_triangle_strip
    integer :: k, npmin, npmax, nrow
    complex(dp), allocatable, dimension(:) :: a, b, c, q, d, du, rhs, x, y

    integer, allocatable, dimension(:) :: irow, icol
    complex(dp), allocatable, dimension(:) :: aval
    integer :: nz

    npmin = 52
    npmax = 101
    nrow = npmax-npmin+1

    open(1,file='test_triangle_strip.out')
    do k=npmin,npmax
       write(1,*) mesh_point(k)%rcoord, mesh_point(k)%zcoord, mesh_point(k)%psi_pol
    end do
    close(1)

    allocate(a(nrow), b(nrow), c(nrow), q(nrow))
    allocate(d(nrow), du(nrow), rhs(nrow), x(nrow), y(nrow))
    allocate(irow(2*nrow), icol(2*nrow), aval(2*nrow))

    do k=npmin, npmax
       a(k-npmin+1) = 1d0/sqrt(mesh_point(k)%rcoord**2 + mesh_point(k)%zcoord**2)
    end do

    b = (0d0,1d0)*n
    c = 1d0
    q = 1d0

    call assemble_system(nrow, a, b, c, q, d, du, x)
    call assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    rhs = x

    call sparse_solve(nrow,nrow,nz,irow,icol,aval,x)

    ! Test numerical error
    call sparse_matmul(nrow,nrow,irow,icol,aval,x,y)
    print *, sqrt(sum((y-rhs)**2))

    deallocate(a,b,c,q,d,du,rhs,x,y)
    deallocate(irow,icol,aval)
  end subroutine test_triangle_strip

  subroutine test_pressure_profile
    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: rold, zold, Brold, Bpold, Bzold ! previous values in loop
    complex(dp), dimension(nt_core) :: qn, a, b, c, d, du, x, Mq
    integer :: kp, kl
    integer :: nz
    integer, dimension(2*nt_core) :: irow, icol
    complex(dp), dimension(2*nt_core) :: aval
    real(dp) :: psi_loc, dens_prof, ddens_dpsi, temp_prof, dtemp_dpsi
    type(knot) :: oldknot, curknot
    integer :: common_tri(2)

    r = 0d0; p = 0d0; z = 0d0

    n = 2 ! TODO: read from input file

    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    open(1,file='test_pressure_x.out', recl=1024)
    do kl = 1, nr_max
      do kp=1, nt_core
        rold = r; zold = z
        Brold = Br; Bpold = Bp; Bzold = Bz

        oldknot = curknot
        curknot = mesh_point(1 + (kl-1)*nt_core + kp)
        r = curknot%rcoord; p = 0d0; z = curknot%zcoord

        call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
             ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

        psi_loc = curknot%psi_pol
        dens_prof=(psi_loc - psimin)/psimax*di0 + d_min
        ddens_dpsi=1.d0/psimax*di0
        temp_prof=(psi_loc - psimin)/psimax*ti0 + t_min
        dtemp_dpsi=1.d0/psimax*ti0

        c(kp) = -(ddens_dpsi*temp_prof + dens_prof*dtemp_dpsi)*1.602176d-12 ! -p0'(psi)

        if(kp==1) cycle ! first point

        call common_triangles(oldknot, curknot, common_tri)
        x(kp-1) = c(kp-1)*(oldknot%b_mod+curknot%b_mod)/2d0* &
             sum(mesh_element_rmp(common_tri(:))%bnorm_vac)/2d0

!qn(kp-1) = (oldknot%b_mod*rold+curknot%b_mod*r)/2d0* &
!             sum(mesh_element_rmp(common_tri(:))%bnorm_vac)/2d0

        a(kp-1)=((Br+Brold)/2d0*(r-rold) &
             +(Bz+Bzold)/2d0*(z-zold))/((r-rold)**2+(z-zold)**2)

!a(kp-1)=(r+rold)/2d0*1d3/sqrt((r-rold)**2+(z-zold)**2)
!a(kp-1)=(Br*(r-rold)+Bz*(z-zold))/((r-rold)**2+(z-zold)**2)
!a(kp-1) = 1e3*z + (0d0,1d0)*1d3*r
        b(kp-1)=(0d0,1d0)*n*(Bp/r+Bpold/rold)/2d0
!b(kp-1)=(0d0,1d0)*n*(Bp+Bpold)/2d0
!b(kp-1) = (0d0,1d0)*n*1d0/2d0
      end do ! kp

       ! once more for the last point
       kp = 1
       rold = r; zold = z
       Brold = Br; Bpold = Bp; Bzold = Bz

       oldknot = curknot
       curknot = mesh_point(1 + (kl-1)*nt_core + kp)

       r = curknot%rcoord; p = 0d0; z = curknot%zcoord

       call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

       call common_triangles(oldknot, curknot, common_tri)
       x(nt_core) = c(nt_core)*(oldknot%b_mod+curknot%b_mod)/2d0* &
            sum(mesh_element_rmp(common_tri(:))%bnorm_vac)/2d0

!qn(nt_core) = (oldknot%b_mod*rold+curknot%b_mod*r)/2d0* &
!             sum(mesh_element_rmp(common_tri(:))%bnorm_vac)/2d0

       a(nt_core)=((Br+Brold)/2d0*(r-rold) &
            +(Bz+Bzold)/2d0*(z-zold))/((r-rold)**2+(z-zold)**2)
!a(nt_core)=(r+rold)/2d0*1d3/sqrt((r-rold)**2+(z-zold)**2)
!a(nt_core)=(Br*(r-rold)+Bz*(z-zold))/((r-rold)**2+(z-zold)**2)
!a(nt_core) = 1e3*z + (0d0,1d0)*1d3*r
       b(nt_core)=(0d0,1d0)*n*(Bp/r+Bpold/rold)/2d0
!b(nt_core)=(0d0,1d0)*n*(Bp+Bpold)/2d0
!b(nt_core) = (0d0,1d0)*n*1d0/2d0

       call assemble_system(nt_core, a, b, c, qn, d, du, Mq)
       call assemble_sparse(nt_core, d, du, nz, irow, icol, aval)
       call sparse_solve(nt_core,nt_core,nz,irow,icol,aval,x)

       if(kl==1) then
       write(1,*) mesh_point(1)%psi_pol, dens_prof, temp_prof,&
            temp_prof*dens_prof*1.602176d-12, real(a(nt_core)), Bp,&
            real(sum(x)/size(x)), aimag(sum(x)/size(x))!, 0d0, 0d0
       end if
       do kp=1,nt_core
          presn((kl-1)*nt_core+1 + kp) = x(kp)
          write(1,*) mesh_point((kl-1)*nt_core+1 + kp)%psi_pol, dens_prof,&
            temp_prof,temp_prof*dens_prof*1.602176d-12, real(a(nt_core)),&
            Bp, real(x(kp)), aimag(x(kp))!, real(a(kp)), aimag(a(kp))
       end do
    end do ! kl
    do kp=1,(npoint - nt_core*nr_max - 1)
       write(1,*) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine test_pressure_profile

  subroutine test_current
    real(dp) :: psimin, psimax

    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    integer :: kl, kt
    integer :: nz
    integer, dimension(4*nt_core) :: irow, icol
    complex(dp), dimension(4*nt_core) :: aval

    complex(dp), dimension(2*nt_core) :: qn, a, b, c, d, du, x, Mq
    real(dp) :: psi_inner, dpsi

    type(triangle) :: elem
    type(knot), dimension(3) :: knots

    real(dp), dimension(3) :: rr, zz, lr, lz
    complex(dp) :: jr, jz, jpar, jperp

    r = 0d0
    p = 0d0
    z = 0d0

    psimin = minval(mesh_point%psi_pol)
    psimax = maxval(mesh_point%psi_pol)
    print *,'psimin = ',psimin,'  psimax = ',psimax

    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    open(1,file='test_current_x.out', recl=1024)

    a = -1d0
    c = 1d0
    qn = 1d0

    ! first triangle loop
    psi_inner = mesh_point(1)%psi_pol
    dpsi = mesh_point(nt_core+2)%psi_pol - psi_inner

    call gen_B2avg(1)

    kl = 1
    do kt=1,nt_core
       call assign_curr3(kt, (pres0(kl+1)-pres0(kl))/dpsi, B2avg(kl))
       elem = mesh_element(kt)
       knots = mesh_point(elem%i_knot)
       !r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
       r = sum(knots%rcoord)/3d0
       p = 0d0
       !z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
       z = sum(knots%zcoord)/3d0
       call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
       b(kt) = -(0d0,1d0)*n*mesh_element(kt)%det_3/2d0*Bp/dpsi
       x(kt) = -curr(kt, mod(elem%knot_h,3)+1) + get_source(kt,kl)
       !x(kt) = 1d0
       !print *, x(kt)
!write(1,*) kt, real(b(kt)), aimag(b(kt)), real(x(kt)), aimag(x(kt))
    end do

    call assemble_system(nt_core, a(1:nt_core), b(1:nt_core), c(1:nt_core),&
         qn(1:nt_core), d(1:nt_core), du(1:nt_core), Mq(1:nt_core))
    call assemble_sparse(nt_core, d(1:nt_core), du(1:nt_core), nz,&
         irow(1:2*nt_core), icol(1:2*nt_core), aval(1:2*nt_core))
    call sparse_solve(nt_core,nt_core,nz,irow(1:2*nt_core),icol(1:2*nt_core),&
         aval(1:2*nt_core),x(1:nt_core))

    do kt=1,nt_core
       call assign_currents(kt, x(kt), x(mod(kt,nt_core)+1))
       !write(1,*) kt, real(b(kt)), aimag(b(kt)), real(x(kt)), aimag(x(kt))
       !write(1,*) psi_inner+dpsi, dpsi, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do

    ! TODO: usual triangle loops
    psi_inner = mesh_point(nt_core+2)%psi_pol
    do kl = 2, nr_max
       call gen_B2avg(kl)

       dpsi = mesh_point(kl*nt_core+2)%psi_pol - psi_inner

       do kt=1,2*nt_core
          call assign_curr3(nt_core+(kl-2)*2*nt_core+kt,&
               (pres0(kl+1)-pres0(kl))/dpsi, B2avg(kl))
          ! TODO
          elem = mesh_element(nt_core+(kl-2)*2*nt_core+kt)
          knots = mesh_point(elem%i_knot)
          !r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
          r = sum(knots%rcoord)/3d0
          p = 0d0
          !z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
          z = sum(knots%zcoord)/3d0
          call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
          b(kt) = -(0d0,1d0)*n*elem%det_3/2d0*Bp/dpsi
          x(kt) = -curr(kt, mod(elem%knot_h,3)+1) + get_source(kt,kl)
       end do !kt

       call assemble_system(2*nt_core, a, b, c, qn, d, du, Mq)
       call assemble_sparse(2*nt_core, d, du, nz, irow, icol, aval)
       call sparse_solve(2*nt_core,2*nt_core,nz,irow,icol,aval,x)

       do kt=1,2*nt_core          !print *, kt
         call assign_currents(nt_core+(kl-2)*2*nt_core+kt, x(kt), x(mod(kt,2*nt_core)+1))
       end do ! kt over triangles in strip

       psi_inner = mesh_point(kl*nt_core+2)%psi_pol
    end do ! kl over radial positions

    ! triangle loops outside region of interest
    do kt=1,(ntri - nt_core*(1+2*(nr_max-1)))
    !   write(1,*) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do

    close(1)

    ! write out currents
    open(1,file='currents.dat', recl=1024)
    do kt=1,ntri
       write(1,*) real(curr(kt,1)), real(curr(kt,2)), real(curr(kt,3)),&
            aimag(curr(kt,1)), aimag(curr(kt,2)), aimag(curr(kt,3))
    end do
    close(1)
    open(1,file='currents_toplot.dat', recl=1024)
    do kt=1,ntri
       knots = mesh_point(mesh_element(kt)%i_knot)
       call get_edge_coord(knots, rr, zz, lr, lz)
       r = sum(knots%rcoord)/3d0
       p = 0d0
       z = sum(knots%zcoord)/3d0
       call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
       jr = sum(curr(kt,:)*lr(:))/(3d0*r)
       jz = sum(curr(kt,:)*lz(:))/(3d0*r)
       jpar = (Br*jr + Bz*jz)/sqrt(Br**2+Bz**2)
       jperp = (Bz*jr - Br*jz)/sqrt(Br**2+Bz**2)
       write(1,*) real(jpar), aimag(jpar), real(jperp), aimag(jperp)
    end do
    close(1)
  end subroutine test_current

  subroutine test_dpsidr
    integer :: kt

    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

    type(triangle) :: elem
    type(knot) :: knots(3)

    open(1,file='test_dpsidr.out', recl=1024)
    do kt=1,ntri
       elem = mesh_element(kt)
       knots = mesh_point(elem%i_knot)
       !r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
       r = sum(knots%rcoord)/3d0
       p = 0d0
       !z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
       z = sum(knots%zcoord)/3d0
       call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

       write(1,*) abs(dpsidr(kt)), sqrt((r*Br)**2+(r*Bz)**2)
    end do
    close(1)
  end subroutine test_dpsidr

  subroutine test_current_new
    integer :: kl, kt, k, nt_loop
    type(triangle) :: elem
    integer, dimension(3) :: iknot_s
    type(knot), dimension(3) :: knots, knots_s
    real(8), dimension(3) :: r, z, lr, lz
    real(dp), dimension(3) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,Bmod

    complex(dp), dimension(2*nt_core) :: a1, a2, b1, b2, d, du, x

    integer :: nz
    integer, dimension(4*nt_core) :: irow, icol
    complex(dp), dimension(4*nt_core) :: aval

    complex(dp) :: source2, Bnpar

    real(dp) :: c1, c2, det, h0n1, h0n2, h0p1, h0p2, gradpsi1, gradpsi2,&
         gradpsisq

    integer :: it

    !a1 = 10d0
    !a2 = -10d0
    !b1 = (0d0,0.5d0)
    !b2 = (0d0,0.5d0)

    nt_loop = nt_core

    open(1,file='currents_toplot.dat', recl=1024)
    do kl = 1,nr_max
       do kt = 1,nt_loop
          if(kl==1) then
             it = kt
          else
             it = nt_core+(kl-2)*2*nt_core+kt
          end if

          elem = mesh_element(it)
          knots = mesh_point(elem%i_knot)
          iknot_s(1) = elem%i_knot(mod(elem%knot_h+1,3)+1)
          iknot_s(2) = elem%i_knot(elem%knot_h)
          iknot_s(3) = elem%i_knot(mod(elem%knot_h,3)+1)
          knots_s = mesh_point(iknot_s)
          call get_edge_coord(knots_s, r, z, lr, lz)
          do k=1,3
             call field(r(k),0d0,z(k),Br(k),Bp(k),Bz(k),dBrdR(k),&
                  dBrdp(k),dBrdZ(k),dBpdR(k),dBpdp(k),dBpdZ(k),&
                  dBzdR(k),dBzdp(k),dBzdZ(k))
             Bmod(k) = sqrt(Br(k)**2+Bp(k)**2+Bz(k)**2)
          end do

          h0n1 = (-Br(1)*lz(1)+Bz(1)*lr(1))/Bmod(1)
          h0n2 = (-Br(2)*lz(2)+Bz(2)*lr(2))/Bmod(2)

          h0p1 = Bp(1)/Bmod(1)
          h0p2 = Bp(2)/Bmod(2)

          a1(kt) = r(1)*h0n1
          a2(kt) = r(2)*h0n2
          b1(kt) = (0d0,1d0)*n*elem%det_3/4d0*h0p1
          b2(kt) = (0d0,1d0)*n*elem%det_3/4d0*h0p2

          c1 = (lr(2)**2+lz(2)**2+lr(1)*lr(1)+lz(1)*lz(2))
          c2 = (lr(1)**2+lz(1)**2+lr(1)*lr(1)+lz(1)*lz(2))

          det = (lr(1)**2+lz(1)**2)*(lr(2)**2+lz(2)**2)-(lr(1)*lr(1)+lz(1)*lz(2))**2
          gradpsi1 = c1*(knots_s(2)%psi_pol-knots_s(1)%psi_pol)/det
          gradpsi2 = c2*(knots_s(3)%psi_pol-knots_s(2)%psi_pol)/det
          gradpsisq = gradpsi1**2*(lr(1)**2+lz(1)**2)+&
               2*gradpsi1*gradpsi2*(lr(1)*lr(2)+lz(1)*lz(2))+&
               gradpsi2**2*(lr(2)**2+lz(2)**2)

          x(kt) = -(0d0,1d0)*elem%det_3/2d0*clight*(&
               (presn(iknot_s(2))-presn(iknot_s(1)))*&
               gradpsi1/(r(1)*Bmod(1)**2)&
               +(presn(iknot_s(3))-presn(iknot_s(2)))*&
               gradpsi2/(r(2)*Bmod(2)**2))

          Bnpar = c1/det*Bnflux(it,1)/r(1)*h0n1+&
               c2/det*Bnflux(it,2)/r(2)*h0n2&
               + (h0p1+h0p2)*Bnp(it)/(r(1)+r(2))

          !print *, Bnpar

          source2 = -(0d0,1d0)*elem%det_3/2d0*Bnpar*clight*&
               (pres0(kl+1)-pres0(kl))/(psi(kl+1)-psi(kl))*gradpsisq*&
               3d0/sum(r*Bmod**3)

          x(kt) = 1d0
          !x(kt) = x(kt) + source2

          !print *, a1(kt), a2(kt)
       end do ! kt = 1,nt_loop

       call assemble_system2(nt_loop, a1(1:nt_loop), a2(1:nt_loop),&
            b1(1:nt_loop), b2(1:nt_loop), d(1:nt_loop), du(1:nt_loop))
       call assemble_sparse(nt_loop, d(1:nt_loop), du(1:nt_loop), nz,&
            irow(1:2*nt_loop), icol(1:2*nt_loop), aval(1:2*nt_loop))
       call sparse_solve(nt_loop,nt_loop,nz,irow(1:2*nt_loop),icol(1:2*nt_loop),&
            aval(1:2*nt_loop),x(1:nt_loop))

       do kt = 1, nt_loop-1
          write(1,*) real((x(kt)+x(kt+1))/2d0), aimag((x(kt)+x(kt+1))/2d0)
       end do
       write(1,*) real((x(nt_loop)+x(1))/2d0), aimag((x(nt_loop)+x(1))/2d0)

       if(kl==1) nt_loop = 2*nt_core
       !exit
    end do ! kl = 1,nr_max
    do kt = 1, ntri - (nt_core + 2*nt_core*(nr_max-1))
       write(1,*) 0d0, 0d0
    end do
    close(1)
  end subroutine test_current_new
end program magdif_test
