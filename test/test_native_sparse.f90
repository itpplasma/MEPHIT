program test_native_sparse
  use sparse_mod, only: sparse_solve, sparse_solve_method
  implicit none
  integer, parameter :: dp = kind(1d0)
  integer :: method, repeat
  integer :: rows(4) = [1,2,2,3], cols(4) = [1,1,2,3]
  integer :: ptr(4) = [1,3,4,5]
  real(dp) :: vals(4) = [3d0,1d0,4d0,5d0]
  real(dp) :: a(3,3), x(3), rhs(3), block(3,2), expected(3,2)
  complex(dp) :: z(3,3), zx(3), zb(3), zblock(3,2), zexpected(3,2)
  character(16) :: mode

  call get_command_argument(1,mode)
  a = reshape([3d0,1d0,0d0,0d0,4d0,0d0,0d0,0d0,5d0],[3,3])
  x = [1d0,-2d0,.5d0]
  if (trim(mode) == 'singular') then
    a(3,3) = 0d0
    a(1,3) = 1d0
    rhs = matmul(a,x)
    call sparse_solve(a,rhs)
    error stop 'singular native factor unexpectedly succeeded'
  else if (trim(mode) == 'invalid') then
    rhs = matmul(a,x)
    call sparse_solve(3,3,4,rows,[1,3,4,6],vals,rhs)
    error stop 'invalid CSC unexpectedly succeeded'
  end if
  z = cmplx(reshape([4d0,-.5d0,.25d0,1d0,5d0,-1d0,0d0,1d0,6d0], &
                   [3,3]),reshape([1d0,.25d0,-1d0,-2d0,-1d0,2d0,.5d0,.25d0,.5d0],[3,3]),dp)
  zx = cmplx(x,[.25d0,-.5d0,.125d0],dp)
  expected(:,1) = x
  expected(:,2) = -.25d0*x
  zexpected(:,1) = zx
  zexpected(:,2) = cmplx(.5d0,-.25d0,dp)*zx
  if (trim(mode) == 'matched') then
    ! One identical workload also compiles against the original sparse_mod.
    ! Explicit iopt=0 avoids relying on the legacy optional-argument SAVE bug.
    sparse_solve_method=3
    do repeat=1,1000
      rhs=matmul(a,x)
      call sparse_solve(a,rhs,0)
      call check(maxval(abs(rhs-x))<2d-14,'matched real scalar')
      block=matmul(a,expected)
      call sparse_solve(a,block,0)
      call check(maxval(abs(block-expected))<2d-14,'matched real block')
      zb=matmul(z,zx)
      call sparse_solve(z,zb,0)
      call check(maxval(abs(zb-zx))<2d-14,'matched complex scalar')
      zblock=matmul(z,zexpected)
      call sparse_solve(z,zblock,0)
      call check(maxval(abs(zblock-zexpected))<2d-14,'matched complex block')
    end do
    print '(a)', 'PASS matched native/legacy real/complex scalar/block workload'
    stop
  end if
  do method = 2,3
    sparse_solve_method = method
    ! COO factor, canonical CSC reuse (nnz=n+1), explicit release.
    rhs = matmul(a,x)
    call sparse_solve(3,3,4,rows,cols,vals,rhs,1)
    call check(maxval(abs(rhs-matmul(a,x))) < 1d-15,'factor-only changed RHS')
    call sparse_solve(3,3,4,rows,ptr,vals,rhs,2)
    call check(maxval(abs(rhs-x)) < 2d-14,'retained CSC real solve')
    call sparse_solve(3,3,4,rows,ptr,vals,rhs,3)
    ! The omitted iopt must reset after previous explicit lifecycle calls.
    rhs = matmul(a,x)
    call sparse_solve(a,rhs,1)
    rhs = 2d0*matmul(a,x)
    call sparse_solve(2d0*a,rhs)
    call check(maxval(abs(rhs-x)) < 2d-14,'dense real default lifecycle')
    block = matmul(a,expected)
    call sparse_solve(3,3,4,rows,cols,vals,block,1)
    call sparse_solve(3,3,4,rows,ptr,vals,block,2)
    call check(maxval(abs(block-expected)) < 2d-14,'retained real block')
    call sparse_solve(3,3,4,rows,ptr,vals,block,3)
    block = matmul(a,expected)
    call sparse_solve(a,block)
    call check(maxval(abs(block-expected)) < 2d-14,'dense real block')
    zb = matmul(z,zx)
    call sparse_solve(z,zb,1)
    call sparse_solve(z,zb,2)
    call check(maxval(abs(zb-zx)) < 2d-14,'retained complex dense')
    call sparse_solve(z,zb,3)
    zb = 2d0*matmul(z,zx)
    call sparse_solve(2d0*z,zb)
    call check(maxval(abs(zb-zx)) < 2d-14,'default complex lifecycle')
    zblock = matmul(z,zexpected)
    call sparse_solve(z,zblock)
    call check(maxval(abs(zblock-zexpected)) < 2d-14,'complex dense block')
    call complex_sparse(z,zexpected)
  end do
  print '(a)', 'PASS native real/complex COO/CSC/dense/block and factor lifecycles'
contains
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
      print *,label
      error stop 1
    end if
  end subroutine

  subroutine complex_sparse(matrix,truth)
    complex(dp), intent(in) :: matrix(3,3),truth(3,2)
    integer :: ir(9),ic(9),pc(4),i,j,k
    complex(dp) :: values(9),b(3),many(3,2)
    k=0
    do j=1,3
      do i=1,3
        k=k+1; ir(k)=i; ic(k)=j; values(k)=matrix(i,j)
      end do
    end do
    pc=[1,4,7,10]
    b=matmul(matrix,truth(:,1))
    call sparse_solve(3,3,9,ir,ic,values,b,1)
    call sparse_solve(3,3,9,ir,pc,values,b,2)
    call check(maxval(abs(b-truth(:,1))) < 2d-14,'complex COO/CSC retained')
    call sparse_solve(3,3,9,ir,pc,values,b,3)
    many=matmul(matrix,truth)
    call sparse_solve(3,3,9,ir,ic,values,many,1)
    call sparse_solve(3,3,9,ir,pc,values,many,2)
    call check(maxval(abs(many-truth)) < 2d-14,'complex COO/CSC block')
    call sparse_solve(3,3,9,ir,pc,values,many,3)
  end subroutine
end program
