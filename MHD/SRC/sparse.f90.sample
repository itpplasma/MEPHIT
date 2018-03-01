program magdif_test

  use from_nrtype
  use sparse_mod,    only : column_full2pointer, column_pointer2full, full2sparse, remap_rc, &
       sparse2full, sparse_matmul, sparse_solve

  implicit none

  call test_linear_system()

contains

  subroutine test_linear_system()

    integer, allocatable, dimension(:) :: irow, icol      ! Row and column indexes
    complex(8), allocatable, dimension(:) :: aval         ! Value
    integer :: nz, nz_sq                                  ! Number of nonzero entries
    integer :: k                                          ! Counter
    
    ! Read input file
    open(1,file='matrix.dat')
    read(1) nrow, ncol, nz
    allocate(rhs(nz),xvec(nz))
    allocate(irow(nz),icol(nz),aval(nz))
    do k = 1,nz
       read(1) irow(k), icol(k), aval(k)
    end do
    close(1)

    call remap_rc(nz,nz_sq,irow,icol,aval)                ! Sort and squeeze matrix
    call sparse_solve(nrow,nrow,nz_sq,irow,icol,aval,rhs) ! Solve system

    ! Test numerical error
    call sparse_matmul(nrow,nrow,irow,icol,aval,rhs,xvec)
    print *, rhs
    print *, xvec

    deallocate(rhs,xvec,irow,icol,aval)
  end subroutine test_linear_system
end program magdif_test
