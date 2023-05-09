module mephit_equil

  use iso_fortran_env, only: dp => real64
  use mephit_util, only: func1d_t

  implicit none

  private

  ! types and associated procedures
  public :: process_profiles, write_profiles, read_profiles

  ! module variables
  public :: dens, e_temp, i_temp, E_r, Phi0

  !> Density profile \f$ n \f$ in cm^-3.
  type(func1d_t) :: dens

  !> Electron temperature profile \f$ T_{e} \f$ in eV.
  type(func1d_t) :: e_temp

  !> Ion temperature profile \f$ T_{i} \f$ in eV.
  type(func1d_t) :: i_temp

  !> Radial electric field profile \f$ E_{r} \f$ in statV cm^-1.
  type(func1d_t) :: E_r

  !> Electric potential profile \f$ \Phi_{0} \f$ in statV.
  type(func1d_t) :: Phi0

contains

  subroutine process_profiles
    use mephit_util, only: pi, ev2erg, func1d_init, func1d_deinit, &
         func1d_read_formatted, resample1d
    use mephit_mesh, only: mesh, fs
    type(func1d_t) :: raw, raw_int
    real(dp) :: rsmall(0:mesh%nflux)
    integer :: krad

    rsmall = sqrt(fs%area / pi)
    ! read density profile
    call func1d_init(dens, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'n.dat')
    call resample1d(raw%x, raw%y, rsmall, dens%y, 3)
    dens%x(:) = fs%psi
    ! read and rescale electron temperature profile
    call func1d_init(e_temp, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'Te.dat')
    call resample1d(raw%x, raw%y * ev2erg, rsmall, e_temp%y, 3)
    e_temp%x(:) = fs%psi
    ! read and rescale ion temperature profile
    call func1d_init(i_temp, 0, mesh%nflux)
    call func1d_read_formatted(raw, 'Ti.dat')
    call resample1d(raw%x, raw%y * ev2erg, rsmall, i_temp%y, 3)
    i_temp%x(:) = fs%psi
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
    Phi0%x(:) = rsmall
    call func1d_deinit(raw_int)
    call func1d_deinit(raw)
  end subroutine process_profiles

  subroutine write_profiles(file, group)
    use mephit_util, only: func1d_write
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_write(dens, file, grp // '/dens', &
         'poloidal flux', 'Mx', 'density', 'cm^-3')
    call func1d_write(e_temp, file, grp // '/e_temp', &
         'poloidal flux', 'Mx', 'electron temperature', 'eV')
    call func1d_write(i_temp, file, grp // '/i_temp', &
         'poloidal flux', 'Mx', 'ion temperature', 'eV')
    call func1d_write(E_r, file, grp // '/E_r', &
         'poloidal flux', 'Mx', 'radial electric field', 'statV cm^-1')
    call func1d_write(Phi0, file, grp // '/Phi0', &
         'poloidal flux', 'Mx', 'electric potential', 'statV')
  end subroutine write_profiles

  subroutine read_profiles(file, group)
    use mephit_util, only: func1d_read
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_read(dens, file, grp // '/dens')
    call func1d_read(e_temp, file, grp // '/e_temp')
    call func1d_read(i_temp, file, grp // '/i_temp')
    call func1d_read(E_r, file, grp // '/E_r')
    call func1d_read(Phi0, file, grp // '/Phi0')
  end subroutine read_profiles

end module mephit_equil
