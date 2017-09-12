module magdif
  use from_nrtype
  use sparse_mod, only: remap_rc
  implicit none
  
  contains
  
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

  subroutine solve_full(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_full
  
  subroutine solve_cycl_tridiag(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_cycl_tridiag

  subroutine solve_sparse(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_sparse
end module magdif

