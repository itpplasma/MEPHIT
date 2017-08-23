program magdif_test

  use from_nrtype
  use constants,     only : pi
  use mesh_mod,      only : ntri, npoint, mesh_point, mesh_element
      
  use magdif,        only : assemble_system
  use sparse_mod,    only : column_full2pointer, column_pointer2full, full2sparse, remap_rc, &
       sparse2full, sparse_matmul, sparse_solve

  implicit none

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

  call test_linear_system()

contains

  subroutine test_linear_system()

    integer :: k, npmin, npmax, n, nrow
    complex(8), allocatable, dimension(:) :: a, b, c, q, d, du, rhs, xvec
    complex(8) :: alpha

    integer, allocatable, dimension(:) :: irow, icol, pcol
    complex(8), allocatable, dimension(:) :: aval
    complex(8), allocatable, dimension(:,:) :: Amat
    integer :: ncol, nz
    
    npmin = 52
    npmax = 101
    nrow = npmax-npmin+1
    
    open(1,file='points.out')
    do k=npmin,npmax
       write(1,*) mesh_point(k)%rcoord, mesh_point(k)%zcoord
    end do
    close(1)

    n = 2

    allocate(a(nrow), b(nrow), c(nrow), q(nrow))
    allocate(d(nrow), du(nrow-1), rhs(nrow), xvec(nrow))
    
    do k=npmin, npmax
       a(k-npmin+1) = 1d0/sqrt(mesh_point(k)%rcoord**2 + mesh_point(k)%zcoord**2)
    end do

    b = (0d0,1d0)*n
    c = 1d0
    q = 1d0

    call assemble_system(nrow, a, b, c, q, d, du, alpha, rhs)
   
    allocate(irow(2*nrow), icol(2*nrow), aval(2*nrow))
    allocate(Amat(nrow,nrow))
    do k = 1,nrow-1
       ! diagonal entries
       irow(2*k-1) = k
       icol(2*k-1) = k
       aval(2*k-1) = d(k)

       ! off-diagonal upper
       irow(2*k) = k
       icol(2*k) = k+1
       aval(2*k) = du(k)
    end do
    ! off-diagonal last entry
    irow(2*nrow-1) = nrow
    icol(2*nrow-1) = 1
    aval(2*nrow-1) = alpha
    
    ! diagonal last entry
    irow(2*nrow) = nrow
    icol(2*nrow) = nrow
    aval(2*nrow) = d(nrow)

    call remap_rc(2*nrow,nz,irow,icol,aval)
    call sparse_solve(nrow,nrow,nz,irow,icol,aval,rhs)

    ! Test numerical error
    !call sparse_matmul(nrow,nrow,irow,icol,aval,rhs,xvec)
    !print *, xvec
    
    deallocate(a,b,c,q,d,du,rhs,xvec)
    deallocate(irow,icol,aval)
  end subroutine test_linear_system
end program magdif_test
