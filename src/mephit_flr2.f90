module mephit_flr2
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  ! types and associated procedures
  public :: flr2_t, flr2_init, flr2_deinit, flr2_write, flr2_read, &
    flr2_coeff, flr2_response_current

  type :: flr2_t
     private

     integer :: npoi = 0, m_min = 0, m_max = 0
     real(dp), dimension(:), allocatable, public :: psi
     complex(dp), dimension(:, :), allocatable, public :: a2_in, a2_out, a0
     complex(dp), dimension(:, :), allocatable, public :: b2_in, b2_out, b0
     complex(dp), dimension(:, :), allocatable, public :: c2_in, c2_out, c0
     complex(dp), dimension(:, :), allocatable, public :: d2_in, d2_out, d0
  end type flr2_t

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

  integer, parameter :: mnmax = 3

contains

  subroutine flr2_init(flr2, npoi, m_min, m_max)
    type(flr2_t), intent(inout) :: flr2
    integer, intent(in) :: npoi, m_min, m_max

    call flr2_deinit(flr2)
    flr2%npoi = npoi
    flr2%m_min = m_min
    flr2%m_max = m_max
    allocate(flr2%psi(npoi))
    allocate(flr2%a2_in(npoi, m_min:m_max))
    allocate(flr2%a2_out(npoi, m_min:m_max))
    allocate(flr2%a0(npoi, m_min:m_max))
    allocate(flr2%b2_in(npoi, m_min:m_max))
    allocate(flr2%b2_out(npoi, m_min:m_max))
    allocate(flr2%b0(npoi, m_min:m_max))
    allocate(flr2%c2_in(npoi, m_min:m_max))
    allocate(flr2%c2_out(npoi, m_min:m_max))
    allocate(flr2%c0(npoi, m_min:m_max))
    allocate(flr2%d2_in(npoi, m_min:m_max))
    allocate(flr2%d2_out(npoi, m_min:m_max))
    allocate(flr2%d0(npoi, m_min:m_max))
  end subroutine flr2_init

  subroutine flr2_deinit(flr2)
    type(flr2_t), intent(inout) :: flr2

    flr2%npoi = 0
    flr2%m_min = 0
    flr2%m_max = 0
    if (allocated(flr2%psi)) deallocate(flr2%psi)
    if (allocated(flr2%a2_in)) deallocate(flr2%a2_in)
    if (allocated(flr2%a2_out)) deallocate(flr2%a2_out)
    if (allocated(flr2%a0)) deallocate(flr2%a0)
    if (allocated(flr2%b2_in)) deallocate(flr2%b2_in)
    if (allocated(flr2%b2_out)) deallocate(flr2%b2_out)
    if (allocated(flr2%b0)) deallocate(flr2%b0)
    if (allocated(flr2%c2_in)) deallocate(flr2%c2_in)
    if (allocated(flr2%c2_out)) deallocate(flr2%c2_out)
    if (allocated(flr2%c0)) deallocate(flr2%c0)
    if (allocated(flr2%d2_in)) deallocate(flr2%d2_in)
    if (allocated(flr2%d2_out)) deallocate(flr2%d2_out)
    if (allocated(flr2%d0)) deallocate(flr2%d0)
  end subroutine flr2_deinit

  subroutine flr2_write(flr2, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(flr2_t), intent(in) :: flr2
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/npoi', flr2%npoi, &
      comment = 'number of flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_min', flr2%m_min, &
      comment = 'minimum poloidal mode number including sign')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_max', flr2%m_max, &
      comment = 'maximum poloidal mode number including sign')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', &
      flr2%psi, lbound(flr2%psi), ubound(flr2%psi), unit = 'Mx', &
      comment = 'poloidal flux')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', &
      flr2%psi, lbound(flr2%psi), ubound(flr2%psi), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/a2_in', &
      flr2%a2_in, lbound(flr2%a2_in), ubound(flr2%a2_in), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/a2_out', &
      flr2%a2_out, lbound(flr2%a2_out), ubound(flr2%a2_out), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/a0', &
      flr2%a0, lbound(flr2%a0), ubound(flr2%a0), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/b2_in', &
      flr2%b2_in, lbound(flr2%b2_in), ubound(flr2%b2_in), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/b2_out', &
      flr2%b2_out, lbound(flr2%b2_out), ubound(flr2%b2_out), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/b0', &
      flr2%b0, lbound(flr2%b0), ubound(flr2%b0), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/c2_in', &
      flr2%c2_in, lbound(flr2%c2_in), ubound(flr2%c2_in), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/c2_out', &
      flr2%c2_out, lbound(flr2%c2_out), ubound(flr2%c2_out), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/c0', &
      flr2%c0, lbound(flr2%c0), ubound(flr2%c0), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2_in', &
      flr2%d2_in, lbound(flr2%d2_in), ubound(flr2%d2_in), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2_out', &
      flr2%d2_out, lbound(flr2%d2_out), ubound(flr2%d2_out), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d0', &
      flr2%d0, lbound(flr2%d0), ubound(flr2%d0), unit = '?', &
      comment = 'Poisson equation coefficient')
    call h5_close(h5id_root)
  end subroutine flr2_write

  subroutine flr2_read(flr2, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(flr2_t), intent(inout) :: flr2
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer :: npoi, m_min, m_max

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/npoi', npoi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_min', m_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_max', m_max)
    call flr2_init(flr2, npoi, m_min, m_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi', flr2%psi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/a2_in', flr2%a2_in)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/a2_out', flr2%a2_out)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/a0', flr2%a0)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/b2_in', flr2%b2_in)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/b2_out', flr2%b2_out)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/b0', flr2%b0)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/c2_in', flr2%c2_in)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/c2_out', flr2%c2_out)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/c0', flr2%c0)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/d2_in', flr2%d2_in)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/d2_out', flr2%d2_out)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/d0', flr2%d0)
    call h5_close(h5id_root)
  end subroutine flr2_read


  ! All units are Gaussian if not otherwise mentioned

  ! Input/Output ::
  ! flr2                 - (flr2_t) Poisson equation coefficients
  !
  ! Input ::
  ! mpol_min             - (integer) minimum poloidal mode number
  ! mpol_max             - (integer) maximum poloidal mode number
  ! ntor                 - (integer) toroidal mode number
  ! npoi                 - (integer) number of radial grid points
  ! am_i                 - (real) ion mass number (mass in units of proton mass)
  ! Z_i                  - (real) ion charge number (charge in units of elementary charge)
  ! Rtor                 - (real) reference major radius
  ! psi(npoi)            - (real) unperturbed poloidal flux $\psi$
  ! qsaf(npoi)           - (real) safety factor
  ! bcovar_phi(npoi)     - (real) co-variant unperturbed toroidal magnetic field component
  ! Phi_0(npoi)          - (real) equilibrium electrostatic potential
  ! avR2nabpsi2(npoi)    - (real) average over the polodal angle of $R^2 |\nabla \psi|^2$
  ! dens_e(npoi)         - (real) electron density
  ! temp_e(npoi)         - (real) electron temperature
  ! temp_i(npoi)         - (real) ion temperature
  ! anu_e(npoi)          - (real) electron collision frequency
  ! anu_i(npoi)          - (real) ion collision frequency
  subroutine flr2_coeff(flr2, mpol_min, mpol_max, ntor, npoi, am_i, z_i, Rtor, &
    psi, qsaf, bcovar_phi, Phi_0, avR2nabpsi2, &
    dens_e, temp_e, temp_i, anu_e, anu_i)

    use mephit_conf, only: logger
    use mephit_util, only: imun, pi, c => clight, e_charge => elem_charge, e_mass => m_e, p_mass => m_p

    type(flr2_t), intent(inout) :: flr2
    integer, intent(in) :: mpol_min, mpol_max, ntor, npoi
    integer :: i, mpol
    real(dp), intent(in) :: am_i, Z_i, Rtor
    real(dp) :: e_e,e_i,omega_E,x_1,x_2,v_T,rho2_factor
    real(dp) :: switch_flr_e,switch_flr_i
    real(dp) :: switch_cur_e,switch_cur_i
    complex(dp)   :: factor_of_Phi
    complex(dp)   :: F_m,Gtor_m,Htor
    complex(dp)   :: F_me,Gtor_me,Htore
    complex(dp)   :: F_mi,Gtor_mi,Htori
    real(dp), dimension(:), intent(in) :: psi, qsaf, bcovar_phi, Phi_0, avR2nabpsi2, &
      dens_e, temp_e, temp_i, anu_e, anu_i

    complex(dp), dimension(0:mnmax,0:mnmax) :: symbI
    real(dp), dimension(:),   allocatable :: dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i

    character(len=11), dimension(10), parameter :: array_names = [character(len=11) :: &
      'psi', 'qsaf', 'bcovar_phi', 'Phi_0', 'avR2nabpsi2', &
      'dens_e', 'temp_e', 'temp_i', 'anu_e', 'anu_i']
    integer, dimension(10) :: array_sizes
    array_sizes = [ &
      size(psi), size(qsaf), size(bcovar_phi), size(Phi_0), size(avR2nabpsi2), &
      size(dens_e), size(temp_e), size(temp_i), size(anu_e), size(anu_i)]
    do i = 1, 10
      if (npoi /= array_sizes(i)) then
        call logger%msg_arg_size('flr2_coeff', &
          'npoi', 'size(' // trim(array_names(i)) // ')', npoi, array_sizes(i))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do
    if (npoi /= flr2%npoi) then
      call logger%msg_arg_size('flr2_coeff', 'npoi', 'flr2%npoi', npoi, flr2%npoi)
      if (logger%err) call logger%write_msg
      error stop
    end if

    flr2%psi(:) = psi
    allocate(dens_i(npoi),dPhi_0_dpsi(npoi),derpar(npoi))
    allocate(A1e(npoi),A2e(npoi),A1i(npoi),A2i(npoi))

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


    do mpol = mpol_min, mpol_max
      do i = 1, npoi

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

        flr2%b0(i, mpol) = -4.d0 * pi * (mpol + qsaf(i) * ntor) * F_m / (qsaf(i) * omega_E * avR2nabpsi2(i))
        flr2%b2_in(i, mpol) = -4.d0 * pi * (mpol + qsaf(i) * ntor) * Gtor_m / (qsaf(i) * omega_E)
        flr2%b2_out(i, mpol) = flr2%b2_in(i, mpol)

        flr2%a0(i, mpol) = -flr2%b0(i, mpol) * factor_of_Phi
        flr2%a2_in(i, mpol) = -flr2%b2_in(i, mpol) * factor_of_Phi + 4.d0 * pi * Htor
        flr2%a2_out(i, mpol) = (1.d0, 0.d0) + flr2%a2_in(i, mpol)

        flr2%d0(i, mpol) = -(F_me * switch_cur_e + F_mi * switch_cur_i) / bcovar_phi(i)
        flr2%d2_in(i, mpol) = -(Gtor_me * switch_cur_e + Gtor_mi * switch_cur_i) * avR2nabpsi2(i) / bcovar_phi(i)
        flr2%d2_out(i, mpol) = flr2%d2_in(i, mpol)

        flr2%c0(i, mpol) = flr2%d0(i, mpol) * factor_of_Phi
        flr2%c2_in(i, mpol) = flr2%d2_in(i, mpol) * factor_of_Phi
        flr2%c2_out(i, mpol) = flr2%c2_in(i, mpol)
      end do
    end do

    deallocate(dens_i,dPhi_0_dpsi,derpar,A1e,A2e,A1i,A2i)
  end subroutine flr2_coeff


  ! Input/Output ::
  ! flr2                 - (flr2_t) Poisson equation coefficients
  !
  ! Input ::
  ! isw_Phi_m            - (integer) switch for the potential perturbation in the parallel current:
  !                        0 - set it to zero in the current, 1 - use it as is in the current
  ! mpol                 - (integer) poloidal mode number
  ! npoi                 - (integer) number of radial grid points
  ! psi(npoi)            - (real) unperturbed poloidal flux $\psi$
  ! bpsi_over_bphi(npoi) - (complex) psi-component of the perturbation magnetic field divided by
  !                        contra-variant component of the unperturbed magnetic field
  ! Output ::
  ! parcur_over_b0(npoi) - (complex) parallel current density divided by uperturbed magnetic field module
  ! Phi_m(npoi)          - (complex) perturbation of the electrostatic potential
  subroutine flr2_response_current(flr2, isw_Phi_m, mpol, npoi, bpsi_over_bphi, parcur_over_b0, Phi_m)

    use mephit_conf, only: logger

    type(flr2_t), intent(in) :: flr2
    integer, intent(in) :: isw_Phi_m, mpol, npoi
    complex(dp), dimension(:), intent(in) :: bpsi_over_bphi
    complex(dp), dimension(:), intent(out) :: parcur_over_b0, Phi_m

    integer :: i
    character(len=14), dimension(3), parameter :: array_names = &
      [character(len=14) :: 'bpsi_over_bphi', 'parcur_over_b0', 'Phi_m']
    integer, dimension(3) :: array_sizes
    array_sizes = [size(bpsi_over_bphi), size(parcur_over_b0), size(Phi_m)]
    do i = 1, 3
      if (npoi /= array_sizes(i)) then
        call logger%msg_arg_size('flr2_response_current', &
          'npoi', 'size(' // trim(array_names(i)) // ')', npoi, array_sizes(i))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do
    if (npoi /= flr2%npoi) then
      call logger%msg_arg_size('flr2_response_current', 'npoi', 'flr2%npoi', npoi, flr2%npoi)
      if (logger%err) call logger%write_msg
      error stop
    end if

    ! Solve Poisson equation:

    call progonka(isw_Phi_m, npoi, &
      flr2%psi, bpsi_over_bphi, &
      flr2%a2_in(:, mpol), flr2%a2_out(:, mpol), flr2%a0(:, mpol), &
      flr2%b2_in(:, mpol), flr2%b2_out(:, mpol), flr2%b0(:, mpol), &
      flr2%c2_in(:, mpol), flr2%c2_out(:, mpol), flr2%c0(:, mpol), &
      flr2%d2_in(:, mpol), flr2%d2_out(:, mpol), flr2%d0(:, mpol), &
      Phi_m, parcur_over_b0)
  end subroutine flr2_response_current


  subroutine getIfunc(x1, x2, symbI)
    integer :: m, n
    real(dp), intent(in) :: x1, x2
    real(dp) :: z
    complex(dp) :: denom
    complex(dp), dimension(0:mnmax, 0:mnmax), intent(out) :: symbI
    complex(dp), dimension(0:mnmax, 0:mnmax) :: Imn

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
          end do
        end do
      else
        symbI=Imn
      end if

    end if

  end subroutine getIfunc


  !> Calculates array of W2 special functions.
  subroutine W2_arr(x1_in, x2_in, Imn)
    real(dp), intent(in) :: x1_in, x2_in
    complex(dp), dimension(0:3, 0:3), intent(out) :: Imn
    complex(dp) :: t1, t2, F11m, x1, x2
    complex(dp), parameter :: I = (0.0d0, 1.0d0), one = (1.0d0, 0.0d0)

    integer :: l,m,n

    real(dp) :: F_im, F_re
    complex(dp), dimension(0:3, 0:3, 0:mnmax) :: W2

    W2 = (0.0d0, 0.0d0)

    x1 = cmplx(x1_in, 0.d0, dp)
    t1 = x1**2

    do l = 0, mnmax ! does not work for 3, 2 and 1!


      x2 = x2_in + I*l

      t2 = - I*x2 + t1

      ! call Hypergeometric1F1_kummer_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im)

      ! F11m = F_re + I*F_im

      call hypergeometric1f1_cont_fract_1_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im)

      F11m = F_re + I*F_im

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
      end do
    end do

  end subroutine W2_arr


  subroutine progonka(isw_f, npoi, x, q, &
    a2_in, a2_out, a0, &
    b2_in, b2_out, b0, &
    c2_in, c2_out, c0, &
    d2_in, d2_out, d0, &
    f, g)

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
    ! x(npoi)    - (real) array of independent variable, need not to be equidistant
    ! c_1(npoi)  - (complex) array of coefficient $C_1$
    ! c_2(npoi)  - (complex) array of coefficient $C_2$
    ! b(npoi)    - (complex) array of source $B$

    ! Output:
    ! phi(npoi)   - (complex) array of the solution $\Phi$

    use mephit_conf, only: logger

    integer, intent(in) :: isw_f, npoi
    integer :: i

    real(dp) :: dxp,dxm,dxt
    complex(dp) :: denom

    real(dp), dimension(:), intent(in) :: x
    complex(dp), dimension(:), intent(in) :: q, &
      a2_in, a2_out, a0, &
      b2_in, b2_out, b0, &
      c2_in, c2_out, c0, &
      d2_in, d2_out, d0
    complex(dp), dimension(:), intent(out) :: f, g
    complex(dp),   dimension(:),   allocatable :: alp,bet,quelle
    complex(dp),   dimension(:,:), allocatable :: wsecder,eqmat

    character(len=6), dimension(16), parameter :: array_names = [character(len=6) :: &
      'x', 'q', &
      'a2_in', 'a2_out', 'a0', 'b2_in', 'b2_out', 'b0', &
      'c2_in', 'c2_out', 'c0', 'd2_in', 'd2_out', 'd0', &
      'f', 'g']
    integer, dimension(16) :: array_sizes
    array_sizes = [ &
      size(x), size(q), &
      size(a2_in), size(a2_out), size(a0), size(b2_in), size(b2_out), size(b0), &
      size(c2_in), size(c2_out), size(c0), size(d2_in), size(d2_out), size(d0), &
      size(f), size(g)]
    do i = 1, 16
      if (npoi /= array_sizes(i)) then
        call logger%msg_arg_size('progonka', &
          'npoi', 'size(' // trim(array_names(i)) // ')', npoi, array_sizes(i))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do

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
    end do

    eqmat(0,1)=(1.d0,0.d0)
    eqmat(0,npoi)=(1.d0,0.d0)

    alp(2)=(0.d0,0.d0)
    bet(2)=b0(1)*q(1)/a0(1)

    do i=2,npoi-1
      denom=eqmat(0,i)+alp(i)*eqmat(-1,i)
      alp(i+1)=-eqmat(1,i)/denom
      bet(i+1)=(quelle(i)-eqmat(-1,i)*bet(i))/denom
    end do

    f(npoi)=b0(npoi)*q(npoi)/a0(npoi)

    do i=npoi-1,1,-1
      f(i)=alp(i+1)*f(i+1)+bet(i+1)
    end do

    f=f*dble(isw_f)

    do i=2,npoi-1

      g(i) = c2_out(i)*sum(f(i-1:i+1)*wsecder(:,i))      &
        + sum(c2_in(i-1:i+1)*f(i-1:i+1)*wsecder(:,i)) &
        + c0(i)*f(i)                                  &
        + d2_out(i)*sum(q(i-1:i+1)*wsecder(:,i))      &
        + sum(d2_in(i-1:i+1)*q(i-1:i+1)*wsecder(:,i)) &
        + d0(i)*q(i)

    end do

    deallocate(wsecder,eqmat,alp,bet,quelle)

  end subroutine progonka


  subroutine first_deriv(npoi, x, f, df_dx)

    ! Computes first derivative df_dx of function f over variable x

    ! Input:
    ! npoi        - (integer) number of grid points
    ! x(npoi)     - (real) array of independent variable, need not to be equidistant
    ! f(npoi)     - (real) array of function values

    ! Output:
    ! df_dx(npoi) - (real) array of derivative values

    use mephit_conf, only: logger

    integer, intent(in) :: npoi
    integer :: i

    real(dp) :: dxp,dxm,dxt

    real(dp), dimension(:), intent(in) :: x, f
    real(dp), dimension(:), intent(out) :: df_dx

    character(len=5), dimension(3), parameter :: array_names = &
      [character(len=5) :: 'x', 'f', 'df_dx']
    integer, dimension(3) :: array_sizes
    array_sizes = [size(x), size(f), size(df_dx)]
    do i = 1, 3
      if (npoi /= array_sizes(i)) then
        call logger%msg_arg_size('first_deriv', &
          'npoi', 'size(' // trim(array_names(i)) // ')', npoi, array_sizes(i))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do

    do i=2,npoi-1
      dxp=x(i+1)-x(i)
      dxm=x(i)-x(i-1)
      dxt=dxp+dxm
      df_dx(i)=((f(i+1)-f(i))*dxm**2-(f(i-1)-f(i))*dxp**2)/(dxp*dxm*dxt)
    end do

    df_dx(1)=(df_dx(2)*(x(3)-x(1))+df_dx(3)*(x(1)-x(2)))/(x(3)-x(2))
    df_dx(npoi)=(df_dx(npoi-1)*(x(npoi-2)-x(npoi))+df_dx(npoi-2)*(x(npoi)-x(npoi-1)))/(x(npoi-2)-x(npoi-1))

  end subroutine first_deriv

end module mephit_flr2
