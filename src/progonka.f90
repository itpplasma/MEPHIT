!
  subroutine progonka(isw_f,npoi,x,q,  &
                      a2_in,a2_out,a0, &
                      b2_in,b2_out,b0, &
                      c2_in,c2_out,c0, &
                      d2_in,d2_out,d0, &
                      f,g)
!
! Solves 2-nd order ODE for the unknown $f$
!
!$$
! a_2^{out}(x) \difp{^2}{x^2} f(x) + \difp{^2}{x^2} a_2^{in}(x) f(x) + a_0(x) f(x) = S(x)
!$$
!
! where the source term is given by derivatives of a given function $q(x)$
!
!$$
! S(x) = b_2^{out}(x) \difp{^2}{x^2} q(x) + \difp{^2}{x^2} b_2^{in}(x) q(x) + b_0(x) q(x)
!$$
!
! Boundary conditions at both interval ends are
!
!$$
! a_0 f = b_0 q
!$$
!
! Using the result, computes the following function
!
!$$
! g(x) = c_2^{out}(x) \difp{^2}{x^2} f(x) + \difp{^2}{x^2} c_2^{in}(x) f(x) + c_0(x) f(x)
!      + d_2^{out}(x) \difp{^2}{x^2} q(x) + \difp{^2}{x^2} d_2^{in}(x) q(x) + d_0(x) q(x)
!$$
!
!
! Input:
! isw_f      - (integer) switch for account of f(x) in g(x):
!              0 - skip the first line in the expression for g(x), 1 - compute g(x) as is
! npoi       - (integer) number of grid points
! x(npoi)    - (double precision) array of independent variable, need not to be equidistant
! c_1(npoi)  - (double complex) array of coefficient $C_1$
! c_2(npoi)  - (double complex) array of coefficient $C_2$
! b(npoi)    - (double complex) array of source $B$
!
! Output:
! phi(npoi)   - (double complex) array of the solution $\Phi$
!
  implicit none
!
  integer :: isw_f,npoi,i
!
  double precision :: dxp,dxm,dxt
  double complex :: denom
!
  double precision, dimension(npoi) :: x
  double complex,   dimension(npoi) :: a2_in,a2_out,a0, &
                                       b2_in,b2_out,b0, &
                                       c2_in,c2_out,c0, &
                                       d2_in,d2_out,d0, &
                                       q,f,g
  double complex,   dimension(:),   allocatable :: alp,bet,quelle
  double complex,   dimension(:,:), allocatable :: wsecder,eqmat
!
  allocate(wsecder(-1:1,npoi),eqmat(-1:1,npoi),alp(npoi),bet(npoi),quelle(npoi))
!
  wsecder(:,1)=0.d0
  wsecder(:,npoi)=0.d0
!
! weights wsecder are defined so that the second derivative of function F(x) at the point x(i) is 
! ddF(i)=sum(F(i-1:i+1)*wsecder(:,i))
!
  do i=2,npoi-1
    dxp = x(i+1)-x(i)
    dxm = x(i)-x(i-1)
    dxt = dxp+dxm
    wsecder(-1,i) = 2.d0/(dxm*dxt)
    wsecder( 0,i) = -2.d0/(dxm*dxp)
    wsecder( 1,i) = 2.d0/(dxp*dxt)
!
    quelle(i) = b2_out(i)*sum(q(i-1:i+1)*wsecder(:,i))      &
              + sum(b2_in(i-1:i+1)*q(i-1:i+1)*wsecder(:,i)) &
              + b0(i)*q(i)
!
    eqmat(:,i) = a2_out(i)*wsecder(:,i)      &
               + a2_in(i-1:i+1)*wsecder(:,i)
    eqmat(0,i) = eqmat(0,i) + a0(i)
  enddo
!
  eqmat(0,1)=(1.d0,0.d0)
  eqmat(0,npoi)=(1.d0,0.d0)
!
  alp(2)=(0.d0,0.d0)
  bet(2)=b0(1)*q(1)/a0(1)
!
  do i=2,npoi-1
    denom=eqmat(0,i)+alp(i)*eqmat(-1,i)
    alp(i+1)=-eqmat(1,i)/denom
    bet(i+1)=(quelle(i)-eqmat(-1,i)*bet(i))/denom
  enddo
!
  f(npoi)=b0(npoi)*q(npoi)/a0(npoi)
!
  do i=npoi-1,1,-1
    f(i)=alp(i+1)*f(i+1)+bet(i+1)
  enddo
!
  f=f*dble(isw_f)
!
  do i=2,npoi-1
!
    g(i) = c2_out(i)*sum(f(i-1:i+1)*wsecder(:,i))      &
         + sum(c2_in(i-1:i+1)*f(i-1:i+1)*wsecder(:,i)) &
         + c0(i)*f(i)                                  &
         + d2_out(i)*sum(q(i-1:i+1)*wsecder(:,i))      &
         + sum(d2_in(i-1:i+1)*q(i-1:i+1)*wsecder(:,i)) &
         + d0(i)*q(i)
!
  enddo
!
  deallocate(wsecder,eqmat,alp,bet,quelle)
!
  end subroutine progonka
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine first_deriv(npoi,x,f,df_dx)
!
! Computes first derivative df_dx of function f over variable x
!
! Input:
! npoi        - (integer) number of grid points
! x(npoi)     - (double precision) array of independent variable, need not to be equidistant
! f(npoi)     - (double precision) array of function values
!
! Output:
! df_dx(npoi) - (double precision) array of derivative values
!
!
  implicit none
!
  integer :: npoi,i
!
  double precision :: dxp,dxm,dxt
!
  double precision, dimension(npoi) :: x,f,df_dx
!
  do i=2,npoi-1
    dxp=x(i+1)-x(i)
    dxm=x(i)-x(i-1)
    dxt=dxp+dxm
    df_dx(i)=((f(i+1)-f(i))*dxm**2-(f(i-1)-f(i))*dxp**2)/(dxp*dxm*dxt)
  enddo
!
  df_dx(1)=(df_dx(2)*(x(3)-x(1))+df_dx(3)*(x(1)-x(2)))/(x(3)-x(2))
  df_dx(npoi)=(df_dx(npoi-1)*(x(npoi-2)-x(npoi))+df_dx(npoi-2)*(x(npoi)-x(npoi-1)))/(x(npoi-2)-x(npoi-1))
!
  end subroutine first_deriv
