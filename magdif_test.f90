program magdif_test

  use from_nrtype
  use constants,     only : pi
  use mesh_mod,      only : ntri, npoint, mesh_point, mesh_element
      
  use magdif,        only : assemble_sparse, assemble_system
  use sparse_mod,    only : sparse_matmul, sparse_solve

  implicit none

  integer, parameter :: nt_core=200
  integer, parameter :: nr_core=75

  print *, "MAGDIF TEST START"
  
  open(1,file='points.dat',form='unformatted')
  read(1) npoint
  print *, npoint
  allocate(mesh_point(npoint))
  read(1) mesh_point
  close(1)
  open(1,file='triangles.dat',form='unformatted')
  read(1) ntri
  print *, ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  close(1)

!  call test_fdm_simple()
!  call test_triangle_strip()
  call test_pressure_profile()

contains
  subroutine test_fdm_simple()
    integer, parameter :: ns = 10
    real(dp) :: ds
    real(dp), dimension(ns) :: s

    integer :: m, n
    complex(dp) :: qmn, pmn
    complex(dp), dimension(ns) :: qn, pn, a, b, c, d, du, x

    integer :: k, nz
    integer, dimension(2*ns) :: irow, icol, pcol
    complex(dp), dimension(2*ns) :: aval
    
    ds = 2d0*pi/ns

    s = (/(k*ds, k=0,ns-1)/)
    n = 1
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

  subroutine test_triangle_strip()

    integer :: k, npmin, npmax, n, nrow
    complex(dp), allocatable, dimension(:) :: a, b, c, q, d, du, rhs, x, y

    integer, allocatable, dimension(:) :: irow, icol, pcol
    complex(dp), allocatable, dimension(:) :: aval
    complex(dp), allocatable, dimension(:,:) :: Amat
    integer :: ncol, nz
    
    npmin = 52
    npmax = 101
    nrow = npmax-npmin+1
    
    open(1,file='test_triangle_strip.out')
    do k=npmin,npmax
       write(1,*) mesh_point(k)%rcoord, mesh_point(k)%zcoord, mesh_point(k)%psi_pol
    end do
    close(1)

    n = 2

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

  subroutine test_pressure_profile()

    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: rold, zold, Brold, Bpold, Bzold ! previous values in loop
    complex(dp), dimension(nt_core) :: qn, pn, a, b, c, d, du, x
    integer :: kp, kl, kt
    integer :: n ! harmonic index of perturbation
    integer :: nz
    integer, dimension(2*nt_core) :: irow, icol, pcol
    complex(dp), dimension(2*nt_core) :: aval
    
    r = 0d0
    p = 0d0
    z = 0d0

    n = 2 ! testing

    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    kl = 10
    open(1,file='test_pressure_rz.out')
    do kp=1, nt_core
       rold = r
       zold = z
       Brold = Br
       Bpold = Bp
       Bzold = Bz
       
       r = mesh_point((kl-1)*nt_core+1 + kp)%rcoord
       p = 0d0
       z = mesh_point((kl-1)*nt_core+1 + kp)%zcoord
       call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

       write(1,*) r,z, (r-rold)**2+(z-zold)**2
       
       c(kp) = 1d0 ! TODO: p0'(psi)
       qn(kp) = cos(atan2(z,r)) ! TODO: deltaBpsictr
       !qn(kp) = 1d0
       
       if(kp==1) cycle ! first point

       a(kp-1)=((Br+Brold)/2d0*(r-rold)+(Bz+Bzold)/2d0*(z-zold))/((r-rold)**2+(z-zold)**2)
       b(kp-1)=(0d0,1d0)*n*(Bp/r+Bpold/rold)/2d0
    end do
    close(1)

    ! once more for the last point
    kp = 1
    rold = r
    zold = z
    Brold = Br
    Bpold = Bp
    Bzold = Bz
       
    r = mesh_point((kl-1)*nt_core+1 + kp)%rcoord
    p = 0d0
    z = mesh_point((kl-1)*nt_core+1 + kp)%zcoord
    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    a(nt_core)=((Br+Brold)/2d0*(r-rold) + (Bz+Bzold)/2d0*(z-zold))/((r-rold)**2+(z-zold)**2)
    b(nt_core)=(0d0,1d0)*n*(Bp/r+Bpold/rold)/2d0

    print *, real(a)
    !print *, aimag(b)
    
    call assemble_system(nt_core, a, b, c, qn, d, du, x)
    call assemble_sparse(nt_core, d, du, nz, irow, icol, aval)
    call sparse_solve(nt_core,nt_core,nz,irow,icol,aval,x)
    
    open(1,file='test_pressure_x.out', recl=1024)
    do kp=1,nt_core
       write(1,*) real(x(kp)), aimag(x(kp))
    end do
    close(1)
  end subroutine test_pressure_profile
end program magdif_test
