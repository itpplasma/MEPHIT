module magdif
  implicit none
  
  contains
  
  subroutine assemble_system(n, a, b, c, q, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: a(:), b(:), c(:) ! system coefficients
    complex(8), intent(in) :: q(:)             ! right-hand side without mass matrix applied
    complex(8), intent(out) :: d(n)            ! diagonal of stiffness matrix
    complex(8), intent(out) :: du(n-1)         ! superdiagonal of stiffness matrix
    complex(8), intent(out) :: alpha           ! single entry in stiffness matrix to be periodic
    complex(8), intent(out) :: Mq(:)           ! right-hand side with mass matrix applied

    integer :: k

    do k=1,n-1
       d(k) = -a(k) + b(k)/2d0
       du(k) = a(k) + b(k)/2d0
       Mq(k) = (c(k+1)*q(k+1)+c(k)*q(k))/2d0
    end do
    d(n) = -a(n) + b(n)/2d0
    Mq(n) = (c(1)*q(1)+c(n)*q(n))/2d0
    alpha = a(1) + b(1)/2d0
  end subroutine assemble_system

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
