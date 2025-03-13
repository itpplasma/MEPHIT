module mephit_flr2
  implicit none

  private

  public :: response_current

contains

  subroutine response_current(isw_Phi_m,mpol,ntor,npoi,am_i,z_i,Rtor,     &
    psi,qsaf,bcovar_phi,Phi_0,avR2nabpsi2,      &
    dens_e,temp_e,temp_i,anu_e,anu_i,           &
    bpsi_over_bphi,parcur_over_b0,Phi_m)

    ! All units are Gaussian if not otherwise mentioned

    ! Input ::
    ! isw_Phi_m            - (integer) switch for the potential perturbation in the parallel current:
    !                        0 - set it to zero in the current, 1 - use it as is in the current
    ! mpol                 - (integer) poloidal mode number
    ! ntor                 - (integer) toroidal mode number
    ! npoi                 - (integer) number of radial grid points
    ! am_i                 - (double precision) ion mass number (mass in units of proton mass)
    ! Z_i                  - (double precision) ion charge number (charge in units of elementary charge)
    ! Rtor                 - (double precision) reference major radius
    ! psi(npoi)            - (double precision) unperturbed poloidal flux $\psi$
    ! qsaf(npoi)           - (double precision) safety factor
    ! bcovar_phi(npoi)     - (double precision) co-variant unperturbed toroidal magnetic field component
    ! Phi_0(npoi)          - (double precision) equilibrium electrostatic potential
    ! avR2nabpsi2(npoi)    - (double precision) average over the polodal angle of $R^2 |\nabla \psi|^2$
    ! dens_e(npoi)         - (double precision) electron density
    ! temp_e(npoi)         - (double precision) electron temperature
    ! temp_i(npoi)         - (double precision) ion temperature
    ! anu_e(npoi)          - (double precision) electron collision frequency
    ! anu_i(npoi)          - (double precision) ion collision frequency
    ! bpsi_over_bphi(npoi) - (double complex) psi-component of the perturbation magnetic field divided by 
    !                        contra-variant component of the unperturbed magnetic field
    ! Output:
    ! parcur_over_b0(npoi) - (double complex) parallel current density divided by uperturbed magnetic field module
    ! Phi_m(npoi)          - (double complex) perturbation of the electrostatic potential

    integer,          parameter :: mnmax=3
    double precision, parameter :: pi=3.14159265358979d0
    double precision, parameter :: c=2.9979d10
    double precision, parameter :: e_charge=4.8032d-10
    double precision, parameter :: e_mass=9.1094d-28
    double precision, parameter :: p_mass=1.6726d-24
    double precision, parameter :: ev=1.6022d-12
    double complex,   parameter :: imun=(0.d0,1.d0)

    integer :: isw_Phi_m,mpol,ntor,npoi,i
    double precision :: am_i,z_i,Rtor
    double precision :: e_e,e_i,omega_E,x_1,x_2,v_T,rho2_factor
    double precision :: switch_flr_e,switch_flr_i
    double precision :: switch_cur_e,switch_cur_i
    double complex   :: factor_of_Phi
    double complex   :: F_m,Gtor_m,Htor
    double complex   :: F_me,Gtor_me,Htore
    double complex   :: F_mi,Gtor_mi,Htori
    double precision, dimension(npoi) :: psi,qsaf,bcovar_phi,Phi_0,avR2nabpsi2, &
      dens_e,temp_e,temp_i,anu_e,anu_i
    double complex,   dimension(npoi) :: bpsi_over_bphi,parcur_over_b0,Phi_m

    double complex, dimension(0:mnmax,0:mnmax) :: symbI
    double precision, dimension(:),   allocatable :: dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i
    double complex,   dimension(:),   allocatable :: a2_in,a2_out,a0
    double complex,   dimension(:),   allocatable :: b2_in,b2_out,b0
    double complex,   dimension(:),   allocatable :: c2_in,c2_out,c0
    double complex,   dimension(:),   allocatable :: d2_in,d2_out,d0

    allocate(dens_i(npoi),dPhi_0_dpsi(npoi),derpar(npoi))
    allocate(A1e(npoi),A2e(npoi),A1i(npoi),A2i(npoi))
    allocate(a2_in(npoi),a2_out(npoi),a0(npoi))
    allocate(b2_in(npoi),b2_out(npoi),b0(npoi))
    allocate(c2_in(npoi),c2_out(npoi),c0(npoi))
    allocate(d2_in(npoi),d2_out(npoi),d0(npoi))

    ! For playing around:
    ! switch for FLR effects in electrons (1 - on, 0 - off):
    switch_flr_e=1.d0
    ! switch for FLR effects in ions (1 - on, 0 - off):
    switch_flr_i=1.d0
    ! switch for the electron contribution to the parallel current (1 - on, 0 - off):
    switch_cur_e=1.d0
    ! switch for the ion contribution to the parallel current (1 - on, 0 - off):
    switch_cur_i=1.d0

    ! electron and ion charges and ion density:
    e_e=-e_charge
    e_i=e_charge*z_i
    dens_i=dens_e/z_i

    call first_deriv(npoi,psi,Phi_0,dPhi_0_dpsi)


    ! thermodynamic forces for electrons:

    call first_deriv(npoi,psi,dens_e,derpar)

    A1e=derpar/dens_e+e_e*dPhi_0_dpsi/temp_e

    call first_deriv(npoi,psi,temp_e,derpar)

    A2e=derpar/temp_e
    A1e=A1e-1.5d0*A2e


    ! thermodynamic forces for ions:

    call first_deriv(npoi,psi,dens_i,derpar)

    A1i=derpar/dens_i+e_i*dPhi_0_dpsi/temp_i

    call first_deriv(npoi,psi,temp_i,derpar)

    A2i=derpar/temp_i
    A1i=A1i-1.5d0*A2i


    do i=1,npoi

      factor_of_Phi=imun*(mpol+qsaf(i)*ntor)/(qsaf(i)*dPhi_0_dpsi(i))
      omega_E=ntor*c*dPhi_0_dpsi(i)

      ! functions F_m and Gtor_m:

      ! electron contribution:
      v_T=sqrt(temp_e(i)/e_mass)
      x_1=(mpol+qsaf(i)*ntor)*v_T/(qsaf(i)*Rtor*anu_e(i))
      x_2=-omega_E/anu_e(i)
      rho2_factor=(e_mass*c*v_T/(e_e*bcovar_phi(i)))**2*switch_flr_e

      ! susceptibility functions for electrons:

      call getIfunc(x_1,x_2,symbI)

      F_me = e_e*dens_e(i)*v_T**2/anu_e(i)                               &
        * (symbI(1,1)*(A1e(i)+A2e(i))+0.5d0*symbI(1,3)*A2e(i))
      Gtor_me = 0.5d0*rho2_factor*e_e*dens_e(i)*v_T**2/anu_e(i)          &
        * (symbI(1,1)*(A1e(i)+2.d0*A2e(i))+0.5d0*symbI(1,3)*A2e(i))
      Htore = 0.5d0*rho2_factor*e_e*dens_e(i)*Rtor**2/dPhi_0_dpsi(i)     &
        * (A1e(i)+2.5d0*A2e(i))

      ! ion contribution:
      v_T=sqrt(temp_i(i)/(am_i*p_mass))
      x_1=(mpol+qsaf(i)*ntor)*v_T/(qsaf(i)*Rtor*anu_i(i))
      x_2=-omega_E/anu_i(i)
      rho2_factor=(am_i*p_mass*c*v_T/(e_i*bcovar_phi(i)))**2*switch_flr_i

      ! susceptibility functions for ions:

      call getIfunc(x_1,x_2,symbI)

      F_mi = e_i*dens_i(i)*v_T**2/anu_i(i)                               &
        * (symbI(1,1)*(A1i(i)+A2i(i))+0.5d0*symbI(1,3)*A2i(i))
      Gtor_mi = 0.5d0*rho2_factor*e_i*dens_i(i)*v_T**2/anu_i(i)          &
        * (symbI(1,1)*(A1i(i)+2.d0*A2i(i))+0.5d0*symbI(1,3)*A2i(i))
      Htori = 0.5d0*rho2_factor*e_i*dens_i(i)*Rtor**2/dPhi_0_dpsi(i)     &
        * (A1i(i)+2.5d0*A2i(i))

      F_m = F_me + F_mi
      Gtor_m = Gtor_me + Gtor_mi
      Htor = Htore + Htori

      b0(i) = -4.d0*pi*(mpol+qsaf(i)*ntor)*F_m/(qsaf(i)*omega_E*avR2nabpsi2(i))
      b2_in(i) = -4.d0*pi*(mpol+qsaf(i)*ntor)*Gtor_m/(qsaf(i)*omega_E)
      b2_out(i) = b2_in(i)

      a0(i) = -b0(i)*factor_of_Phi
      a2_in(i) = -b2_in(i)*factor_of_Phi+4.d0*pi*Htor
      a2_out(i) = (1.d0,0.d0)+a2_in(i)

      d0(i) = -(F_me*switch_cur_e+F_mi*switch_cur_i)/bcovar_phi(i)
      d2_in(i) = -(Gtor_me*switch_cur_e+Gtor_mi*switch_cur_i)*avR2nabpsi2(i)/bcovar_phi(i)
      d2_out(i) = d2_in(i)

      c0(i) = d0(i)*factor_of_Phi
      c2_in(i) = d2_in(i)*factor_of_Phi
      c2_out(i) = c2_in(i)
    enddo

    ! Solve Laplace equation:

    call progonka(isw_Phi_m,npoi,       &
      psi,bpsi_over_bphi,   &
      a2_in,a2_out,a0,      &
      b2_in,b2_out,b0,      &
      c2_in,c2_out,c0,      &
      d2_in,d2_out,d0,      &
      Phi_m,parcur_over_b0)

    deallocate(dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i)
    deallocate(a2_in,a2_out,a0)
    deallocate(b2_in,b2_out,b0)
    deallocate(c2_in,c2_out,c0)
    deallocate(d2_in,d2_out,d0)

  end subroutine response_current

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine getIfunc(x1,x2,symbI)
    integer, parameter :: mnmax=3
    integer :: m,n
    double precision :: x1,x2,z
    double complex :: denom
    double complex, dimension(0:mnmax,0:mnmax) :: symbI,Imn

    !  if(.true.) then
    if(.false.) then
      ! collisionless case:
      symbI=(0.d0,0.d0)
      z=x2/(sqrt(2.d0)*x1)

      symbI(0,0)=sqrt(2.d0)*exp(-z**2)/abs(x1)
      symbI(1,0)=symbI(0,0)*x2/x1
      symbI(1,1)=symbI(1,0)*x2/x1

      symbI(2,0)=symbI(0,0)*(x2/x1)**2
      symbI(2,1)=symbI(1,0)*(x2/x1)**2
      symbI(3,0)=symbI(2,1)
      symbI(3,1)=symbI(1,1)*(x2/x1)**2

      symbI(2,2)=symbI(2,0)*(x2/x1)**2
      symbI(3,2)=symbI(2,1)*(x2/x1)**2
      symbI(3,3)=symbI(3,1)*(x2/x1)**2
    else
      ! collisional case:

      call W2_arr(x1,x2,Imn)

      if(.true.) then
        !    if(.false.) then
        ! energy conservation:
        denom=(1.d0,0.d0)-Imn(0,0)+(2.d0,0.d0)*Imn(2,0)-Imn(2,2)
        do m=0,3
          do n=0,3
            symbI(m,n)=Imn(m,n)+(Imn(m,0)-Imn(m,2))*(Imn(n,0)-Imn(n,2))/denom
          enddo
        enddo
      else
        symbI=Imn
      endif
      !write (1234,*) x1,real(Imn(0,0)),dimag(Imn(0,0)),real(Imn(2,0)),dimag(Imn(2,0)), &
      !real(Imn(2,2)),dimag(Imn(2,2)),real(Imn(1,0)-Imn(1,2)),dimag(Imn(1,0)-Imn(1,2))

    endif

  end subroutine getIfunc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !<Calculates array of W2 special functions.

  !------------------------------------------------------------------------------

  subroutine W2_arr (x1_in,x2_in,Imn)

    interface
       subroutine hypergeometric1f1_cont_fract_1_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
         bind(C, name = 'hypergeometric1f1_cont_fract_1_modified_0_ada')
         use iso_c_binding, only: c_double
         real(c_double), intent(in) :: b_re, b_im, z_re, z_im
         real(c_double), intent(out) :: f_re, f_im
       end subroutine hypergeometric1f1_cont_fract_1_modified_0_ada

       subroutine hypergeometric1f1_kummer_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
         bind(C, name = 'hypergeometric1f1_kummer_modified_0_ada')
         use iso_c_binding, only: c_double
         real(c_double), intent(in) :: b_re, b_im, z_re, z_im
         real(c_double), intent(out) :: f_re, f_im
       end subroutine hypergeometric1f1_kummer_modified_0_ada
    end interface

    integer, parameter :: dp = 8
    integer, parameter :: dpc = 8
    real(dp), parameter :: pi    = 3.141592653589793238462643383279502884197_dp
    real(dp) :: x1_in,x2_in;
    complex(dpc), dimension(0:3,0:3) :: Imn
    complex(dpc) :: t1, t2, F11m, x1, x2
    !complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), one = cmplx(1.0d0, 0.0d0, 8);
    double complex, parameter :: I = (0.0d0, 1.0d0), one = (1.0d0, 0.0d0)

    integer :: l,nmax,m,n;

    real(dp) :: F_im, F_re;
    complex(dpc), allocatable, dimension(:,:,:) :: W2;

    nmax=3

    !allocate (W2(0:1, 0:3, -Nmax:Nmax))
    allocate (W2(0:3, 0:3, 0:nmax))
    W2 = (0.0d0, 0.0d0)

    x1 = dcmplx(x1_in,0.d0)
    t1 = x1**2

    do l = 0, nmax ! does not work for 3, 2 and 1!


      x2 = x2_in + I*l

      t2 = - I*x2 + t1;

      !      call Hypergeometric1F1_kummer_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);

      !      F11m = F_re + I*F_im;

      call hypergeometric1f1_cont_fract_1_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);

      F11m = F_re + I*F_im;

      !asymptotic form of expressions: has no fake singularity at kp=0

      W2(0,0,l) =  &
        (-I)*(-x2**2 + 2 + 5*x1**2 -  &
        (3*I)*x2*(1 + x1**2) + (3 + F11m)*x1**4)/ &
        ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
        (x2 + I*(2 + x1**2)))

      W2(0,1,l) =  &
        (-I)*x1*(x2*(-1 + F11m*x1**2) -  &
        I*(2 + 3*x1**2+ x1**4))/ &
        ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
        (x2 + I*(2 + x1**2)))

      W2(0,2,l) =  &
        -((x2 + I)*(2 + 3*x1**2 -  &
        I*x2*(1 - F11m*x1**2) + x1**4))/ &
        ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
        (x2 + I*(2 + x1**2)))

      W2(0,3,l) =  &
        (-I)*x1*(F11m*x2**3- (3 + 2*F11m)*x2 +  &
        I*x2**2*(3*F11m - x1**2) -  &
        I*(6 + (7 + 2*F11m)*x1**2 + x1**4))/ &
        ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
        (x2+ I*(2 + x1**2)))

      W2(1,0,l) =  &
        (-I)*x1*(x2*(-1 + F11m*x1**2) -  &
        I*(2 + 3*x1**2 + x1**4))/ &
        ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
        (x2+ I*(2 + x1**2)))

      W2(1,1,l) =  &
        (-I)*x2*(x2*(-1 + F11m*x1**2) -  &
        I*(2 + 3*x1**2 + x1**4))/ &
        ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
        (x2+ I*(2 + x1**2)))

      W2(1,2,l) =  &
        (-I)*x1*(F11m*x2**3+ I*x2**2*(F11m - x1**2) -  &
        x2*(3 + 2*x1**2) -  &
        I*(2 + 3*x1**2 + x1**4))/ &
        ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
        (x2+ I*(2 + x1**2)))

      W2(1,3,l) =  &
        x2*((-I)*F11m*x2**3+ I*(3 + 2*F11m)*x2 -  &
        6 - (7 + 2*F11m)*x1**2 +  &
        x2**2*(3*F11m - x1**2) - x1**4)/ &
        ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
        (x2+ I*(2 + x1**2)))


    end do

    W2(2,2,0) = I*x1*W2(1,2,1) + 2.d0*W2(1,1,1) - I*x1*W2(1,2,0) + W2(0,2,0)
    W2(2,2,1) = I*x1*W2(1,2,2) + 2.d0*W2(1,1,2) - I*x1*W2(1,2,1) + W2(0,2,1)
    W2(2,3,0) = I*x1*W2(2,2,1) + 2.d0*W2(1,2,1) - I*x1*W2(2,2,0) + 2.0d0*W2(1,2,0)
    W2(2,3,1) = I*x1*W2(1,3,2) + 3.d0*W2(1,2,2) - I*x1*W2(1,3,1) + W2(0,3,1)
    W2(3,3,0) = I*x1*W2(2,3,1) + 3.d0*W2(2,2,1) - I*x1*W2(2,3,0) + 2.0d0*W2(1,3,0)

    Imn = W2(:,:,0)
    do m=0,3
      do n=0,m-1
        Imn(m,n)=Imn(n,m)
      enddo
    enddo

  end subroutine W2_arr

  !------------------------------------------------------------------------------


  subroutine progonka(isw_f,npoi,x,q,  &
    a2_in,a2_out,a0, &
    b2_in,b2_out,b0, &
    c2_in,c2_out,c0, &
    d2_in,d2_out,d0, &
    f,g)

    ! Solves 2-nd order ODE for the unknown $f$

    !$$
    ! a_2^{out}(x) \difp{^2}{x^2} f(x) + \difp{^2}{x^2} a_2^{in}(x) f(x) + a_0(x) f(x) = S(x)
    !$$

    ! where the source term is given by derivatives of a given function $q(x)$

    !$$
    ! S(x) = b_2^{out}(x) \difp{^2}{x^2} q(x) + \difp{^2}{x^2} b_2^{in}(x) q(x) + b_0(x) q(x)
    !$$

    ! Boundary conditions at both interval ends are

    !$$
    ! a_0 f = b_0 q
    !$$

    ! Using the result, computes the following function

    !$$
    ! g(x) = c_2^{out}(x) \difp{^2}{x^2} f(x) + \difp{^2}{x^2} c_2^{in}(x) f(x) + c_0(x) f(x)
    !      + d_2^{out}(x) \difp{^2}{x^2} q(x) + \difp{^2}{x^2} d_2^{in}(x) q(x) + d_0(x) q(x)
    !$$


    ! Input:
    ! isw_f      - (integer) switch for account of f(x) in g(x):
    !              0 - skip the first line in the expression for g(x), 1 - compute g(x) as is
    ! npoi       - (integer) number of grid points
    ! x(npoi)    - (double precision) array of independent variable, need not to be equidistant
    ! c_1(npoi)  - (double complex) array of coefficient $C_1$
    ! c_2(npoi)  - (double complex) array of coefficient $C_2$
    ! b(npoi)    - (double complex) array of source $B$

    ! Output:
    ! phi(npoi)   - (double complex) array of the solution $\Phi$

    integer :: isw_f,npoi,i

    double precision :: dxp,dxm,dxt
    double complex :: denom

    double precision, dimension(npoi) :: x
    double complex,   dimension(npoi) :: a2_in,a2_out,a0, &
      b2_in,b2_out,b0, &
      c2_in,c2_out,c0, &
      d2_in,d2_out,d0, &
      q,f,g
    double complex,   dimension(:),   allocatable :: alp,bet,quelle
    double complex,   dimension(:,:), allocatable :: wsecder,eqmat

    allocate(wsecder(-1:1,npoi),eqmat(-1:1,npoi),alp(npoi),bet(npoi),quelle(npoi))

    wsecder(:,1)=0.d0
    wsecder(:,npoi)=0.d0

    ! weights wsecder are defined so that the second derivative of function F(x) at the point x(i) is 
    ! ddF(i)=sum(F(i-1:i+1)*wsecder(:,i))

    do i=2,npoi-1
      dxp = x(i+1)-x(i)
      dxm = x(i)-x(i-1)
      dxt = dxp+dxm
      wsecder(-1,i) = 2.d0/(dxm*dxt)
      wsecder( 0,i) = -2.d0/(dxm*dxp)
      wsecder( 1,i) = 2.d0/(dxp*dxt)

      quelle(i) = b2_out(i)*sum(q(i-1:i+1)*wsecder(:,i))      &
        + sum(b2_in(i-1:i+1)*q(i-1:i+1)*wsecder(:,i)) &
        + b0(i)*q(i)

      eqmat(:,i) = a2_out(i)*wsecder(:,i)      &
        + a2_in(i-1:i+1)*wsecder(:,i)
      eqmat(0,i) = eqmat(0,i) + a0(i)
    enddo

    eqmat(0,1)=(1.d0,0.d0)
    eqmat(0,npoi)=(1.d0,0.d0)

    alp(2)=(0.d0,0.d0)
    bet(2)=b0(1)*q(1)/a0(1)

    do i=2,npoi-1
      denom=eqmat(0,i)+alp(i)*eqmat(-1,i)
      alp(i+1)=-eqmat(1,i)/denom
      bet(i+1)=(quelle(i)-eqmat(-1,i)*bet(i))/denom
    enddo

    f(npoi)=b0(npoi)*q(npoi)/a0(npoi)

    do i=npoi-1,1,-1
      f(i)=alp(i+1)*f(i+1)+bet(i+1)
    enddo

    f=f*dble(isw_f)

    do i=2,npoi-1

      g(i) = c2_out(i)*sum(f(i-1:i+1)*wsecder(:,i))      &
        + sum(c2_in(i-1:i+1)*f(i-1:i+1)*wsecder(:,i)) &
        + c0(i)*f(i)                                  &
        + d2_out(i)*sum(q(i-1:i+1)*wsecder(:,i))      &
        + sum(d2_in(i-1:i+1)*q(i-1:i+1)*wsecder(:,i)) &
        + d0(i)*q(i)

    enddo

    deallocate(wsecder,eqmat,alp,bet,quelle)

  end subroutine progonka

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine first_deriv(npoi,x,f,df_dx)

    ! Computes first derivative df_dx of function f over variable x

    ! Input:
    ! npoi        - (integer) number of grid points
    ! x(npoi)     - (double precision) array of independent variable, need not to be equidistant
    ! f(npoi)     - (double precision) array of function values

    ! Output:
    ! df_dx(npoi) - (double precision) array of derivative values


    integer :: npoi,i

    double precision :: dxp,dxm,dxt

    double precision, dimension(npoi) :: x,f,df_dx

    do i=2,npoi-1
      dxp=x(i+1)-x(i)
      dxm=x(i)-x(i-1)
      dxt=dxp+dxm
      df_dx(i)=((f(i+1)-f(i))*dxm**2-(f(i-1)-f(i))*dxp**2)/(dxp*dxm*dxt)
    enddo

    df_dx(1)=(df_dx(2)*(x(3)-x(1))+df_dx(3)*(x(1)-x(2)))/(x(3)-x(2))
    df_dx(npoi)=(df_dx(npoi-1)*(x(npoi-2)-x(npoi))+df_dx(npoi-2)*(x(npoi)-x(npoi-1)))/(x(npoi-2)-x(npoi-1))

  end subroutine first_deriv

end module mephit_flr2
