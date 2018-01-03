module magdif
  use from_nrtype
  use constants, only: clight, pi
  use sparse_mod, only: remap_rc, sparse_solve, sparse_matmul
  use mesh_mod, only: knot, triangle, mesh_point, mesh_element, mesh_element_rmp, ntri,&
       bphicovar, npoint
  use gettormode_mod, only: B_Rn,B_Zn,n_tor_out
  use for_macrostep, only : t_min, d_min
  
  implicit none

  integer, parameter :: nt_core=200
  !integer, parameter :: nr_core=75
  !integer, parameter :: nr_max=74
  integer, parameter :: nr_max=74

  integer :: n ! harmonic index of perturbation
  real(dp), allocatable :: pres0(:)
  real(dp), allocatable :: psi(:)
  real(dp), allocatable :: q(:), dqdpsi(:)
  real(dp), allocatable :: B2avg(:)
  complex(dp), allocatable :: presn(:)
  complex(dp), allocatable :: curr(:,:)
  complex(dp), allocatable :: Bnflux(:,:), Bnp(:) ! fluxes R*Bn*n and physical toroidal comp.

  complex(dp), parameter :: imun = (0d0,1d0)
  real(dp), parameter :: R0 = 172.74467899999999
  
  real(dp), parameter ::  ti0=3.d3, di0 = 5.d13 ! from readcarre_m.f90 TODO read
  real(dp) :: psimin, psimax, dens, temp
  
contains

  subroutine init
    
    integer :: k ! counter
    real(dp) :: dummy8(8)
  
    n = 2 ! TODO: read this from config file
    call read_mesh

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
    
    call read_hpsi

    allocate(B2avg(nr_max))
    allocate(q(nr_max))
    allocate(dqdpsi(nr_max))

    call init_safety_factor

    open(1,file='q.out')
    do k=2,nr_max-1
       write(1,*) (psi(k+1)+psi(k))/2d0, q(k), dqdpsi(k),(q(k+1)-q(k))/((psi(k+1)-psi(k-1))/2d0)
    end do
    close(1)
    
    call set_Bnfluxes_test_notheta
    
    open(1,file='VACFIELD/Bn_flux.dat')
    do k=1,ntri
       read(1,*) dummy8
       !     Bnflux(k,1) = dummy8(1) + (0d0,1d0)*dummy8(2)
       !     Bnflux(k,2) = dummy8(3) + (0d0,1d0)*dummy8(4)
       !     Bnflux(k,3) = dummy8(5) + (0d0,1d0)*dummy8(6)
       !     Bnp(k) = dummy8(7) + (0d0,1d0)*dummy8(8)

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
  end subroutine init

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
         gradpsisq, absgrpsi(3)

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
             absgrpsi(k) = sqrt(r(k)**2*(Br(k)**2+Bz(k)**2))
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
          !gradpsisq = ((absgrpsi(1)+absgrpsi(2))/2d0)**2
          !print *, gradpsisq, ((absgrpsi(1)+absgrpsi(2))/2d0)**2
          stop 'TODO: grad psi'
          
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

          !x(kt) = 1d0
          x(kt) = x(kt) + source2

          !print *, a1(kt), a2(kt)
          !print *, b1(kt), b2(kt)
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
  
  subroutine assemble_system(nrow, a, b, c, q, d, du, Mq)
    integer, intent(in) :: nrow                ! number of system rows
    complex(dp), intent(in) :: a(:), b(:), c(:) ! system coefficients
    complex(dp), intent(in) :: q(:)             ! right-hand side without mass matrix applied
    complex(dp), intent(out) :: d(nrow)         ! diagonal of stiffness matrix
    complex(dp), intent(out) :: du(nrow)        ! superdiagonal of stiffness matrix + A(n,1)
    complex(dp), intent(out) :: Mq(:)           ! right-hand side with mass matrix applied

    integer :: k

    do k=1,nrow-1
       d(k) = -a(k) + b(k)/2d0
       du(k) = a(k) + b(k)/2d0
       Mq(k) = (c(k+1)*q(k+1)+c(k)*q(k))/2d0
    end do
    
    d(nrow) = -a(nrow) + b(nrow)/2d0
    du(nrow) = a(1) + b(1)/2d0
    Mq(nrow) = (c(1)*q(1)+c(nrow)*q(nrow))/2d0
  end subroutine assemble_system
  
subroutine assemble_system2(nrow, a1, a2, b1, b2, d, du)
    integer, intent(in) :: nrow                ! number of system rows
    complex(dp), intent(in) :: a1(:), a2(:), b1(:), b2(:) ! system coefficients
    complex(dp), intent(out) :: d(nrow)         ! diagonal of stiffness matrix
    complex(dp), intent(out) :: du(nrow)        ! superdiagonal of stiffness matrix + A(n,1)
    
    d = a1 + b1/2d0
    du(1:(nrow-1)) = a2(1:(nrow-1)) + b2(1:(nrow-1))/2d0
    du(nrow) = a2(1) + b2(1)/2d0
  end subroutine assemble_system2
  
  subroutine assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    integer, intent(in)  :: nrow                          ! number of system rows
    complex(dp), intent(in)  :: d(nrow)                   ! diagnonal
    complex(dp), intent(in)  :: du(nrow)                  ! upper diagonal
    integer, intent(out) :: nz                            ! number of system rows
    integer, intent(out) :: irow(2*nrow), icol(2*nrow)    ! matrix index representation
    complex(dp), intent(out) :: aval(2*nrow)              ! matrix values
    
    integer :: k
    
    irow(1) = 1
    icol(1) = 1
    aval(1) = d(1)
 
    irow(2) = nrow
    icol(2) = 1
    aval(2) = du(nrow)

    do k = 2,nrow
       ! off-diagonal 
       irow(2*k-1) = k-1
       icol(2*k-1) = k
       aval(2*k-1) = du(k-1)

       ! diagonal 
       irow(2*k) = k
       icol(2*k) = k
       aval(2*k) = d(k)
    end do

    nz = 2*nrow
  end subroutine assemble_sparse

!!$  subroutine solve_full(n, d, du, alpha, Mq)
!!$    integer, intent(in) :: n                   ! number of system rows
!!$    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
!!$    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
!!$    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
!!$    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
!!$    ! TODO: implement
!!$  end subroutine solve_full
!!$  
!!$  subroutine solve_cycl_tridiag(n, d, du, alpha, Mq)
!!$    integer, intent(in) :: n                   ! number of system rows
!!$    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
!!$    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
!!$    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
!!$    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
!!$    ! TODO: implement
!!$  end subroutine solve_cycl_tridiag
!!$
!!$  subroutine solve_sparse(n, d, du, alpha, Mq)
!!$    integer, intent(in) :: n                   ! number of system rows
!!$    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
!!$    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
!!$    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
!!$    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
!!$    ! TODO: implement
!!$  end subroutine solve_sparse

  subroutine common_triangles(knot1, knot2, common_tri)
    type(knot), intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)

    integer :: k, l, kcom
    kcom = 0
    do k=1,knot1%n_owners
       do l=1,knot2%n_owners
          if (knot1%i_owner_tri(k)==knot2%i_owner_tri(l)) then
             kcom = kcom+1
             if (kcom>2) stop "Error: more than two common triangles for knots"
             common_tri(kcom) = knot1%i_owner_tri(k)
          end if
       end do
    end do
  end subroutine common_triangles

  subroutine get_edge_coord(knots, r, z, lr, lz)

    type(knot), intent(in) :: knots(3)
    real(dp), dimension(3) :: r, z, lr, lz

    r(1) = (knots(1)%rcoord + knots(2)%rcoord)/2d0
    r(2) = (knots(2)%rcoord + knots(3)%rcoord)/2d0
    r(3) = (knots(3)%rcoord + knots(1)%rcoord)/2d0
    
    z(1) = (knots(1)%zcoord + knots(2)%zcoord)/2d0
    z(2) = (knots(2)%zcoord + knots(3)%zcoord)/2d0
    z(3) = (knots(3)%zcoord + knots(1)%zcoord)/2d0
    
    lr(1) = knots(2)%rcoord - knots(1)%rcoord
    lr(2) = knots(3)%rcoord - knots(2)%rcoord
    lr(3) = knots(1)%rcoord - knots(3)%rcoord
    
    lz(1) = knots(2)%zcoord - knots(1)%zcoord
    lz(2) = knots(3)%zcoord - knots(2)%zcoord
    lz(3) = knots(1)%zcoord - knots(3)%zcoord
  end subroutine get_edge_coord
  
  subroutine assign_curr3(ktri, p0pr, B2a)
    
    integer, intent(in) :: ktri
    real(dp), intent(in) :: p0pr, B2a
    
    type(triangle) :: elem
    type(knot) :: knots(3)
    
    real(dp), dimension(3) :: r,z,lr,lz
    real(dp) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: B2, j0parB
    
    integer :: indl3

    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)
    
    call get_edge_coord(knots, r, z, lr, lz)
    indl3 = mod(elem%knot_h,3)+1

    call field(r(indl3),0d0,z(indl3),Br,Bp,Bz,dBrdR,&
         dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,&
         dBzdR,dBzdp,dBzdZ)
    
    curr(ktri, indl3) = clight/(Br**2+Bp**2+Bz**2)*(&
         Bp*r(indl3)*(presn(elem%i_knot(mod(indl3,3)+1))-presn(elem%i_knot(indl3)))&
         -imun*n*(presn(elem%i_knot(mod(indl3,3)+1))+presn(elem%i_knot(indl3)))/2d0*&
         (Br*lr(indl3)+Bz*lz(indl3)))
    
    B2 = Br**2+Bp**2+Bz**2
    j0parB = -clight*p0pr*r(indl3)*(1d0/B2 - 1d0/B2a)
    curr(ktri, indl3) = curr(ktri, indl3) - j0parB*Bnflux(ktri,indl3)/B2

  end subroutine assign_curr3
    
    
  subroutine assign_currents(ktri, I1, I2)
    ! assigns currents in the correct direction
    
    integer, intent(in) :: ktri
    complex(dp), intent(in) :: I1, I2

    type(triangle) :: elem
    type(knot) :: knots(3)

    integer :: indl3

    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)    
    indl3 = mod(elem%knot_h,3)+1

    if (minval(elem%i_knot) == elem%i_knot(elem%knot_h)) then
       ! type 1 triangle
       curr(ktri, elem%knot_h) = I2
       curr(ktri, mod(indl3,3)+1) = -I1
    else       
       ! type 2 triangle
       curr(ktri, elem%knot_h) = -I1
       curr(ktri, mod(indl3,3)+1) = I2
    end if
  end subroutine assign_currents

  function absgradpsi(r,z)
    real(dp) , intent(in) :: r,z
    real(dp) absgradpsi
    
    real(dp) :: p,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

    p = 0d0
    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    absgradpsi = sqrt((Br**2+Bz**2)*r**2)
  end function absgradpsi

  function get_source(ktri,kl)
    integer, intent(in) :: ktri,kl
    real(dp) :: p0pr
    complex(dp) get_source

    type(triangle) :: elem
    type(knot), dimension(3) :: knots, knots_s

    integer, dimension(3) :: iknot_s
    complex(dp) pres_source, cur_source
    real(dp), dimension(3) :: r, z, lr, lz
    real(dp), dimension(3) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    integer :: k

    p0pr = (pres0(kl+1)-pres0(kl))/(psi(kl+1)-psi(kl))

    elem = mesh_element(ktri)
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
    end do
    
    pres_source = r(1)/2d0*(presn(iknot_s(1))-presn(iknot_s(2)))/&
         (knots_s(1)%psi_pol - knots_s(2)%psi_pol)&
         + r(2)/2d0*(presn(iknot_s(3))-presn(iknot_s(2)))/&
         (knots_s(3)%psi_pol - knots_s(2)%psi_pol)

    pres_source = pres_source - Bnp(ktri)*2d0/(Bp(1)+Bp(2))*p0pr
    pres_source = clight*pres_source

    cur_source = j0phi(r(1),0d0,z(1),p0pr,B2avg(kl))*(Bnp(ktri)/Bp(1)&
         + Bnflux(ktri,mod(elem%knot_h+1,3)+1)/(knots_s(2)%psi_pol-knots_s(1)%psi_pol))&
         + j0phi(r(2),0d0,z(2),p0pr,B2avg(kl))*(Bnp(ktri)/Bp(2)&
         + Bnflux(ktri,elem%knot_h)/(knots_s(3)%psi_pol-knots_s(2)%psi_pol))
!    print *, Bnp(ktri)/Bp(1)
!pres_source = 0d0    
!cur_source = 0d0
    
    get_source = pres_source + cur_source
    get_source = -(0d0,1d0)*n*get_source*elem%det_3/2d0 ! weigh with surface area
!    get_source = 0d0
  end function get_source
  
  function j0phi(r,p,z,p0pr,B2a)
    real(dp), intent(in) :: r,p,z,p0pr,B2a
    real(dp) :: j0phi, j0phipar, j0phiperp
    real(dp) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: B2, Bpol2, Btor2

    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    B2 = (Br**2+Bp**2+Bz**2)
    Bpol2 = Br**2+Bz**2
    Btor2 = Bp**2
    
    j0phipar = clight*p0pr*Btor2*(1d0/B2 - 1d0/B2a)
    j0phiperp = clight*p0pr*Bpol2/B2
!print *, j0phipar, j0phiperp, Btor2, Bpol2
    j0phi = r*(j0phipar + j0phiperp)
  end function j0phi

  subroutine gen_B2avg(kl)
    integer, intent(in) :: kl

    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    integer :: kt
    type(triangle) :: elem
    type(knot), dimension(3) :: knots

    B2avg(kl) = 0d0
    
    if(kl==1) then
       do kt=1,nt_core
          elem = mesh_element(kt)
          knots = mesh_point(elem%i_knot)
          
          r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
          z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
          call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
          B2avg = B2avg + Br**2 + Bp**2 + Bz**2
       end do !kt
       B2avg = B2avg/nt_core
    else
       do kt=1,2*nt_core
          elem = mesh_element(nt_core+(kl-2)*2*nt_core+kt)
          knots = mesh_point(elem%i_knot)
          
          r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
          z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
          call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
          B2avg = B2avg + Br**2 + Bp**2 + Bz**2
       end do !kt
       B2avg = B2avg/(2*nt_core)
    end if
  end subroutine gen_B2avg

  function get_layer(ktri)
    integer, intent(in) :: ktri
    integer :: get_layer

    if(ktri<=nt_core) then
       get_layer = 1
    else
       get_layer = (ktri-1-nt_core)/(2*nt_core)+2
    end if
  end function get_layer
  
  subroutine set_Bnflux3_Bp_notheta(ktri)
    ! generates Bnflux and Bnp without theta component
    ! based on mesh_element_rmp(k)%bnorm_vac
    integer :: ktri, kl
    type(triangle) :: elem
    integer, dimension(3) :: iknot_s
    type(knot), dimension(3) :: knots, knots_s
    complex(dp) :: Bnflux3
    real(dp), dimension(3) :: r, z, lr, lz
    integer :: common_tri(2)
    
    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)

    iknot_s(1) = elem%i_knot(mod(elem%knot_h+1,3)+1)
    iknot_s(2) = elem%i_knot(elem%knot_h)
    iknot_s(3) = elem%i_knot(mod(elem%knot_h,3)+1)

    knots_s = mesh_point(iknot_s) ! sorted knots

    call get_edge_coord(knots_s, r, z, lr, lz)
    kl = get_layer(ktri)

    call common_triangles(knots_s(1), knots_s(3), common_tri)

    Bnflux3 = -sqrt(lr(3)**2+lz(3)**2)*r(3)*(knots_s(1)%b_mod+knots_s(3)%b_mod)/2d0&
         *(mesh_element_rmp(common_tri(1))%bnorm_vac&
         +mesh_element_rmp(common_tri(2))%bnorm_vac)/2d0&
         *1d0/absgradpsi(r(3),z(3))

    ! type 2 triangle
    if (.not. (minval(elem%i_knot) == elem%i_knot(elem%knot_h))) then
       Bnflux3 = -Bnflux3
    end if

    Bnflux(ktri,mod(elem%knot_h,3)+1) = Bnflux3
    Bnp(ktri) = mesh_element_rmp(ktri)%bnorm_vac*sum(knots%b_mod)/3d0&
         *(0d0,1d0)/(n*q(kl))*dqdpsi(kl)*sum(knots%rcoord)/3d0
  end subroutine set_Bnflux3_Bp_notheta
      
  subroutine set_Bnfluxes_test_notheta()
    integer :: kl, kt, nt_loop, it
    type(triangle) :: elem
    
    integer :: nz
    complex(dp), dimension(2*nt_core) :: d, du, x       
    integer, dimension(4*nt_core) :: irow, icol
    complex(dp), dimension(4*nt_core) :: aval

    complex(dp) :: flux_total
    complex(dp) :: flux_pol(2*nt_core)

    nt_loop = nt_core
    
    do kl = 1,nr_max
       flux_total = 0d0 ! total flux through loop
       
       do kt = 1,nt_loop
          it = nt_core+(kl-2)*2*nt_core+kt
          if(kl==1) it = kt
          
          call set_Bnflux3_Bp_notheta(it)
          flux_total = flux_total&
               +Bnflux(it,mod(mesh_element(it)%knot_h,3)+1)&
               +(0d0,1d0)*n*Bnp(it)*mesh_element(it)%det_3/2d0
       end do

       do kt = 1,nt_loop
          it = nt_core+(kl-2)*2*nt_core+kt
          if(kl==1) it = kt

          Bnp(it) = Bnp(it) - flux_total/((0d0,1d0)*n*mesh_element(it)%det_3/2d0*nt_loop)       
       end do

       ! check
       flux_total = 0d0
       do kt = 1,nt_loop
          it = nt_core+(kl-2)*2*nt_core+kt
          if(kl==1) it = kt
          
          flux_total = flux_total&
               +Bnflux(it,mod(mesh_element(it)%knot_h,3)+1)&
               +(0d0,1d0)*n*Bnp(it)*mesh_element(it)%det_3/2d0
       end do
       if (abs(flux_total) > 1d-12) then
          stop 'abs of total flux through triangle loop nonzero'
       end if

       flux_pol = 0d0
       do kt = 1,nt_loop-1
          it = nt_core+(kl-2)*2*nt_core+kt
          if(kl==1) it = kt
          elem = mesh_element(it)
          
          flux_pol(kt+1) = flux_pol(kt) - Bnflux(it,mod(elem%knot_h,3)+1)&
               -(0d0,1d0)*n*Bnp(it)*elem%det_3/2d0 ! counted in counter-clockwise direction
          
          ! type 1 triangle
          if (minval(elem%i_knot) == elem%i_knot(elem%knot_h)) then
             Bnflux(it, elem%knot_h) = flux_pol(kt+1)
             Bnflux(it, mod(elem%knot_h,3)+2) = -flux_pol(kt)
          else ! type 2
             Bnflux(it, elem%knot_h) = -flux_pol(kt)
             Bnflux(it, mod(elem%knot_h,3)+2) = flux_pol(kt+1)
          end if
       end do

       ! last triangle
       kt = nt_loop
       it = nt_core+(kl-2)*2*nt_core+kt
       if(kl==1) it = kt
       elem = mesh_element(it)
       ! type 1 triangle
       if (minval(elem%i_knot) == elem%i_knot(elem%knot_h)) then
          Bnflux(it, elem%knot_h) = flux_pol(1)
          Bnflux(it, mod(elem%knot_h,3)+2) = -flux_pol(kt)
       else ! type 2
          Bnflux(it, elem%knot_h) = -flux_pol(kt)
          Bnflux(it, mod(elem%knot_h,3)+2) = flux_pol(1)
       end if
       
       if(kl==1) nt_loop = 2*nt_core
    end do
  end subroutine set_Bnfluxes_test_notheta
  
  subroutine init_safety_factor
    integer :: kl, kt, nt_loop
    type(triangle) :: elem
    type(knot), dimension(3) :: knots
    real(dp) :: r, z, Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,rinv
    
    q = 0d0
    do kl = 1,nr_max
       do kt=1,nt_loop
          if(kl==1) then
             elem = mesh_element(kt)
          else
             elem = mesh_element(nt_core+(kl-2)*2*nt_core+kt)
          end if
          knots = mesh_point(elem%i_knot)
          r = sum(knots(:)%rcoord)/3d0
          z = sum(knots(:)%zcoord)/3d0
          rinv = sum(1d0/knots(:)%rcoord**2)/3d0
          call field(r,0d0,z,Br,Bp,Bz,dBrdR,&
               dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,&
               dBzdR,dBzdp,dBzdZ)
          q(kl) = q(kl) + Bp/(2d0*pi)*elem%det_3/2d0 !bphicovar*rinv*elem%V_tri
       end do
       q(kl) = q(kl)/(psi(kl+1)-psi(kl))
       if(kl==1) nt_loop = 2*nt_core
    end do
    q(1) = q(2)
    
    dqdpsi(1) = 0d0
    do kl = 2,nr_max-1
       ! q on triangle loop, psi on edge loop
       dqdpsi(kl) = .5d0*(q(kl+1)-q(kl-1))/(psi(kl+1)-psi(kl))
    end do
    dqdpsi(nr_max) = dqdpsi(nr_max-1)
  end subroutine init_safety_factor

  subroutine read_mesh
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
  end subroutine read_mesh

  subroutine read_hpsi
    integer :: k
    real(dp) :: dummy_re, dummy_im
    
    open(1, file='VACFIELD/hpsi_vac.dat')
    !open(1,file='hpsi.dat')
    do k=1, ntri
       read (1,*) dummy_re, dummy_im
       !     mesh_element_rmp(k)%bnorm_vac=cmplx(dummy_re,dummy_im,dp)
       !mesh_element_rmp(k)%bnorm_vac = 1d0 ! check with constant source
       mesh_element_rmp(k)%bnorm_vac = 3.d0*R0*abs(bphicovar) &
            /sum(mesh_point(mesh_element(k)%i_knot(:))%rcoord**2 &
            *mesh_point(mesh_element(k)%i_knot(:))%b_mod) ! check 2
    end do
    close(1)
  end subroutine read_hpsi

  subroutine test_Bnp
    ! for testing with purely radial field, no theta component
    integer :: ktri, kl
    complex(dp) :: Bnp_test
    type(triangle) :: elem
   
    open(1, file='test_bnp.out')
    open(2, file='test_bnflux.out')
    do ktri=1,ntri
       elem = mesh_element(ktri)
       kl = get_layer(ktri)
       if(kl > nr_max) exit
       Bnp_test = (0d0,1d0)/(n*q(kl))*dqdpsi(kl)*mesh_element_rmp(ktri)%bnorm_vac*&
            sum(mesh_point(elem%i_knot)%b_mod)/3d0*sum(mesh_point(elem%i_knot)%rcoord)/3d0
       write(1,*) real(Bnp(ktri)), aimag(Bnp(ktri)), real(Bnp_test), aimag(Bnp_test) 
       write(2,*) real(Bnflux(ktri,mod(elem%knot_h+1,3)+1)),&
            real(Bnflux(ktri,elem%knot_h)), real(Bnflux(ktri,mod(elem%knot_h,3)+1))
    end do
    close(2)
    close(1)
  end subroutine test_Bnp
end module magdif
