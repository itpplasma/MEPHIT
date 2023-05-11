module mephit_equil

  use iso_fortran_env, only: dp => real64
  use mephit_util, only: func1d_t

  implicit none

  private

  ! types and associated procedures
  public :: process_profiles, write_profiles, read_profiles, deinit_profiles

  ! module variables
  public :: m_i, Z_i, dens_e, temp_e, temp_i, E_r, Phi0, dPhi0_dpsi, nu_e, nu_i

  !> mass number of ions
  real(dp) :: m_i = 2d0

  !> charge number of ions
  real(dp) :: Z_i = 1d0

  !> Electron density profile \f$ n \f$ in cm^-3.
  type(func1d_t) :: dens_e

  !> Electron temperature profile \f$ T_{e} \f$ in eV.
  type(func1d_t) :: temp_e

  !> Ion temperature profile \f$ T_{i} \f$ in eV.
  type(func1d_t) :: temp_i

  !> Radial electric field profile \f$ E_{r} \f$ in statV cm^-1.
  type(func1d_t) :: E_r

  !> Electric potential profile \f$ \Phi_{0} \f$ in statV.
  type(func1d_t) :: Phi0

  !> psi derivative of electric potential profile \f$ \Phi_{0} \f$ in statV.
  type(func1d_t) :: dPhi0_dpsi

  !> Electron collision frequency profile \f$ \nu_{e} \f$ in s^-1.
  type(func1d_t) :: nu_e

  !> Ion collision frequency profile \f$ \nu_{e} \f$ in s^-1.
  type(func1d_t) :: nu_i

contains

  subroutine process_profiles
    use mephit_util, only: pi, func1d_init, func1d_deinit, func1d_read_formatted, resample1d
    use mephit_mesh, only: mesh, fs
    type(func1d_t) :: raw, raw_int
    real(dp) :: rsmall(0:mesh%nflux)
    integer :: krad

    rsmall = sqrt(fs%area / pi)
    ! read electron density profile
    call func1d_init(dens_e, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'n.dat')
    call resample1d(raw%x, raw%y, rsmall, dens_e%y, 3)
    dens_e%x(:) = fs%psi
    ! read electron temperature profile
    call func1d_init(temp_e, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'Te.dat')
    call resample1d(raw%x, raw%y, rsmall, temp_e%y, 3)
    temp_e%x(:) = fs%psi
    ! read ion temperature profile
    call func1d_init(temp_i, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'Ti.dat')
    call resample1d(raw%x, raw%y, rsmall, temp_i%y, 3)
    temp_i%x(:) = fs%psi
    ! read radial electric field profile
    call func1d_init(E_r, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'Er.dat')
    call resample1d(raw%x, raw%y, rsmall, E_r%y, 3)
    E_r%x(:) = fs%psi
    ! integrate electric field to yield the electric potential
    call func1d_init(Phi0, 0, mesh%nflux)
    call func1d_init(raw_int, lbound(raw%x, 1), ubound(raw%x, 1))
    raw_int%x(:) = raw%x
    raw_int%y(:) = 0d0
    do krad = lbound(raw%x, 1) + 1, ubound(raw%x, 1)
       raw_int%y(krad) = raw_int%y(krad - 1) + (raw%x(krad) - raw%x(krad - 1)) &
            * 0.5d0 * (raw%y(krad) + raw%y(krad - 1))
    end do
    call resample1d(raw_int%x, raw_int%y, rsmall, Phi0%y, 3)
    Phi0%x(:) = fs%psi
    call func1d_deinit(raw_int)
    call func1d_deinit(raw)
    ! psi derivative of Phi0
    call func1d_init(dPhi0_dpsi, 0, mesh%nflux)
    dPhi0_dpsi%x(:) = fs%psi
    call resample1d(Phi0%x, Phi0%y, dPhi0_dpsi%x, dPhi0_dpsi%y, 3, .true.)
    ! electron collision frequency
    call func1d_init(nu_e, 0, mesh%nflux)
    nu_e%x(:) = fs%psi
    nu_e%y(:) = 7.7d-6 * (1d0 + Z_i) * dens_e%y / temp_e%y ** 1.5d0 &
         * (24.d0 - log(sqrt(dens_e%y) / temp_e%y))
    ! ion collision frequency
    call func1d_init(nu_i, 0, mesh%nflux)
    nu_i%x(:) = fs%psi
    nu_i%y(:) = 1.8d-7 * Z_i ** 3 / sqrt(m_i) * dens_e%y / temp_i%y ** 1.5d0 &
         * (23.d0 - log(Z_i ** 2.5d0 * sqrt(2.d0 * dens_e%y) / temp_i%y ** 1.5d0))
  end subroutine process_profiles

  subroutine write_profiles(file, group)
    use mephit_util, only: func1d_write
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_write(dens_e, file, grp // '/dens_e', &
         'poloidal flux', 'Mx', 'electron density', 'cm^-3')
    call func1d_write(temp_e, file, grp // '/temp_e', &
         'poloidal flux', 'Mx', 'electron temperature', 'eV')
    call func1d_write(temp_i, file, grp // '/temp_i', &
         'poloidal flux', 'Mx', 'ion temperature', 'eV')
    call func1d_write(E_r, file, grp // '/E_r', &
         'poloidal flux', 'Mx', 'radial electric field', 'statV cm^-1')
    call func1d_write(Phi0, file, grp // '/Phi0', &
         'poloidal flux', 'Mx', 'electric potential', 'statV')
    call func1d_write(Phi0, file, grp // '/dPhi0_dpsi', &
         'poloidal flux', 'Mx', 'psi derivative of electric potential', 'statV Mx^-1')
    call func1d_write(nu_e, file, grp // '/nu_e', &
         'poloidal flux', 'Mx', 'electron collision frequency', 's^-1')
    call func1d_write(nu_i, file, grp // '/nu_i', &
         'poloidal flux', 'Mx', 'ion collision frequency', 's^-1')
  end subroutine write_profiles

  subroutine read_profiles(file, group)
    use mephit_util, only: func1d_read
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_read(dens_e, file, grp // '/dens_e')
    call func1d_read(temp_e, file, grp // '/temp_e')
    call func1d_read(temp_i, file, grp // '/temp_i')
    call func1d_read(E_r, file, grp // '/E_r')
    call func1d_read(Phi0, file, grp // '/Phi0')
    call func1d_read(Phi0, file, grp // '/dPhi0_dpsi')
    call func1d_read(nu_e, file, grp // '/nu_e')
    call func1d_read(nu_i, file, grp // '/nu_i')
  end subroutine read_profiles

  subroutine deinit_profiles
    use mephit_util, only: func1d_deinit

    call func1d_deinit(dens_e)
    call func1d_deinit(temp_e)
    call func1d_deinit(temp_i)
    call func1d_deinit(E_r)
    call func1d_deinit(Phi0)
    call func1d_deinit(dPhi0_dpsi)
    call func1d_deinit(nu_e)
    call func1d_deinit(nu_i)
  end subroutine deinit_profiles

end module mephit_equil
