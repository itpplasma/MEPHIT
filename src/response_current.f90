!
  subroutine response_current(isw_Phi_m,mpol,ntor,npoi,am_i,z_i,Rtor,     &
                              psi,qsaf,bcovar_phi,Phi_0,avR2nabpsi2,      &
                              dens_e,temp_e,temp_i,anu_e,anu_i,           &
                              bpsi_over_bphi,parcur_over_b0,Phi_m)
!
! All units are Gaussian if not otherwise mentioned
!
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
! 
  implicit none
!
  integer,          parameter :: mnmax=3
  double precision, parameter :: pi=3.14159265358979d0
  double precision, parameter :: c=2.9979d10
  double precision, parameter :: e_charge=4.8032d-10
  double precision, parameter :: e_mass=9.1094d-28
  double precision, parameter :: p_mass=1.6726d-24
  double precision, parameter :: ev=1.6022d-12
  double complex,   parameter :: imun=(0.d0,1.d0)
!
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
!
  double complex, dimension(0:mnmax,0:mnmax) :: symbI
  double precision, dimension(:),   allocatable :: dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i
  double complex,   dimension(:),   allocatable :: a2_in,a2_out,a0
  double complex,   dimension(:),   allocatable :: b2_in,b2_out,b0
  double complex,   dimension(:),   allocatable :: c2_in,c2_out,c0
  double complex,   dimension(:),   allocatable :: d2_in,d2_out,d0
!
  allocate(dens_i(npoi),dPhi_0_dpsi(npoi),derpar(npoi))
  allocate(A1e(npoi),A2e(npoi),A1i(npoi),A2i(npoi))
  allocate(a2_in(npoi),a2_out(npoi),a0(npoi))
  allocate(b2_in(npoi),b2_out(npoi),b0(npoi))
  allocate(c2_in(npoi),c2_out(npoi),c0(npoi))
  allocate(d2_in(npoi),d2_out(npoi),d0(npoi))
!
! For playing around:
! switch for FLR effects in electrons (1 - on, 0 - off):
  switch_flr_e=1.d0
! switch for FLR effects in ions (1 - on, 0 - off):
  switch_flr_i=1.d0
! switch for the electron contribution to the parallel current (1 - on, 0 - off):
  switch_cur_e=1.d0
! switch for the ion contribution to the parallel current (1 - on, 0 - off):
  switch_cur_i=1.d0
!
! electron and ion charges and ion density:
  e_e=-e_charge
  e_i=e_charge*z_i
  dens_i=dens_e/z_i
!
  call first_deriv(npoi,psi,Phi_0,dPhi_0_dpsi)
!
!
! thermodynamic forces for electrons:
!
  call first_deriv(npoi,psi,dens_e,derpar)
!
  A1e=derpar/dens_e+e_e*dPhi_0_dpsi/temp_e
!
  call first_deriv(npoi,psi,temp_e,derpar)
!
  A2e=derpar/temp_e
  A1e=A1e-1.5d0*A2e
!
!
! thermodynamic forces for ions:
!
  call first_deriv(npoi,psi,dens_i,derpar)
!
  A1i=derpar/dens_i+e_i*dPhi_0_dpsi/temp_i
! 
  call first_deriv(npoi,psi,temp_i,derpar)
!
  A2i=derpar/temp_i
  A1i=A1i-1.5d0*A2i
!
!
  do i=1,npoi
!
    factor_of_Phi=imun*(mpol+qsaf(i)*ntor)/(qsaf(i)*dPhi_0_dpsi(i))
    omega_E=ntor*c*dPhi_0_dpsi(i)
!
! functions F_m and Gtor_m:
!
! electron contribution:
    v_T=sqrt(temp_e(i)/e_mass)
    x_1=(mpol+qsaf(i)*ntor)*v_T/(qsaf(i)*Rtor*anu_e(i))
    x_2=-omega_E/anu_e(i)
    rho2_factor=(e_mass*c*v_T/(e_e*bcovar_phi(i)))**2*switch_flr_e
!
! susceptibility functions for electrons:
!
    call getIfunc(x_1,x_2,symbI)
!
    F_me = e_e*dens_e(i)*v_T**2/anu_e(i)                               &
         * (symbI(1,1)*(A1e(i)+A2e(i))+0.5d0*symbI(1,3)*A2e(i))
    Gtor_me = 0.5d0*rho2_factor*e_e*dens_e(i)*v_T**2/anu_e(i)          &
            * (symbI(1,1)*(A1e(i)+2.d0*A2e(i))+0.5d0*symbI(1,3)*A2e(i))
    Htore = 0.5d0*rho2_factor*e_e*dens_e(i)*Rtor**2/dPhi_0_dpsi(i)     &
          * (A1e(i)+2.5d0*A2e(i))
!
! ion contribution:
    v_T=sqrt(temp_i(i)/(am_i*p_mass))
    x_1=(mpol+qsaf(i)*ntor)*v_T/(qsaf(i)*Rtor*anu_i(i))
    x_2=-omega_E/anu_i(i)
    rho2_factor=(am_i*p_mass*c*v_T/(e_i*bcovar_phi(i)))**2*switch_flr_i
!
! susceptibility functions for ions:
!
    call getIfunc(x_1,x_2,symbI)
!
    F_mi = e_i*dens_i(i)*v_T**2/anu_i(i)                               &
         * (symbI(1,1)*(A1i(i)+A2i(i))+0.5d0*symbI(1,3)*A2i(i))
    Gtor_mi = 0.5d0*rho2_factor*e_i*dens_i(i)*v_T**2/anu_i(i)          &
            * (symbI(1,1)*(A1i(i)+2.d0*A2i(i))+0.5d0*symbI(1,3)*A2i(i))
    Htori = 0.5d0*rho2_factor*e_i*dens_i(i)*Rtor**2/dPhi_0_dpsi(i)     &
          * (A1i(i)+2.5d0*A2i(i))
!
    F_m = F_me + F_mi
    Gtor_m = Gtor_me + Gtor_mi
    Htor = Htore + Htori
!
    b0(i) = -4.d0*pi*(mpol+qsaf(i)*ntor)*F_m/(qsaf(i)*omega_E*avR2nabpsi2(i))
    b2_in(i) = -4.d0*pi*(mpol+qsaf(i)*ntor)*Gtor_m/(qsaf(i)*omega_E)
    b2_out(i) = b2_in(i)
!
    a0(i) = -b0(i)*factor_of_Phi
    a2_in(i) = -b2_in(i)*factor_of_Phi+4.d0*pi*Htor
    a2_out(i) = (1.d0,0.d0)+a2_in(i)
!
    d0(i) = -(F_me*switch_cur_e+F_mi*switch_cur_i)/bcovar_phi(i)
    d2_in(i) = -(Gtor_me*switch_cur_e+Gtor_mi*switch_cur_i)*avR2nabpsi2(i)/bcovar_phi(i)
    d2_out(i) = d2_in(i)
!
    c0(i) = d0(i)*factor_of_Phi
    c2_in(i) = d2_in(i)*factor_of_Phi
    c2_out(i) = c2_in(i)
  enddo
!
! Solve Laplace equation:
!
  call progonka(isw_Phi_m,npoi,       &
                psi,bpsi_over_bphi,   &
                a2_in,a2_out,a0,      &
                b2_in,b2_out,b0,      &
                c2_in,c2_out,c0,      &
                d2_in,d2_out,d0,      &
                Phi_m,parcur_over_b0)
!
  deallocate(dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i)
  deallocate(a2_in,a2_out,a0)
  deallocate(b2_in,b2_out,b0)
  deallocate(c2_in,c2_out,c0)
  deallocate(d2_in,d2_out,d0)
!
  end subroutine response_current
