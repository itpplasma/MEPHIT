module mephit_equil

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  ! types and associated procedures
  public :: profile_t, profile_init, profile_deinit, profile_read, read_profiles, &
       profile_export_hdf5, profile_import_hdf5, write_profiles_hdf5, read_profiles_hdf5

  ! module variables
  public :: dens, e_temp, i_temp, E_r, Phi0

  type :: profile_t
     integer :: nrad

     !> Minor radius \f$ r \f$ in cm.
     !>
     !> The minor radius is defined here as the radius of a circle with
     !> area equal to the enclosed flux surface.
     real(dp), dimension(:), allocatable :: rsmall

     !> Profile data.
     real(dp), dimension(:), allocatable :: prof

  end type profile_t

  !> Density profile \f$ n \f$ in cm^-3.
  type(profile_t) :: dens

  !> Electron temperature profile \f$ T_{e} \f$ in eV.
  type(profile_t) :: e_temp

  !> Ion temperature profile \f$ T_{i} \f$ in eV.
  type(profile_t) :: i_temp

  !> Radial electric field profile \f$ E_{r} \f$ in statV cm^-1.
  type(profile_t) :: E_r

  !> Electric potential profile \f$ \Phi_{0} \f$ in statV.
  type(profile_t) :: Phi0

contains

  subroutine profile_init(profile, nrad)
    type(profile_t), intent(inout) :: profile
    integer, intent(in) :: nrad

    call profile_deinit(profile)
    profile%nrad = nrad
    allocate(profile%rsmall(nrad))
    allocate(profile%prof(nrad))
  end subroutine profile_init

  subroutine profile_deinit(profile)
    type(profile_t), intent(inout) :: profile

    profile%nrad = 0
    if (allocated(profile%rsmall)) deallocate(profile%rsmall)
    if (allocated(profile%prof)) deallocate(profile%prof)
  end subroutine profile_deinit

  subroutine profile_read(profile, fname)
    type(profile_t), intent(inout) :: profile
    character(len = *), intent(in) :: fname
    integer :: fid, status, nrad, krad

    open(newunit = fid, file = trim(adjustl(fname)), &
         status = 'old', form = 'formatted', action = 'read')
    nrad = 0
    do
       read(fid, *, iostat = status)
       if (status /= 0) exit
       nrad = nrad + 1
    end do
    rewind fid
    call profile_init(profile, nrad)
    do krad = 1, nrad
       read (fid, *) profile%rsmall(krad), profile%prof(krad)
    end do
    close(fid)
  end subroutine profile_read

  subroutine profile_export_hdf5(profile, file, group, comment, unit)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(profile_t), intent(in) :: profile
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: comment
    character(len = *), intent(in) :: unit
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/nrad', profile%nrad, &
         comment = 'number of radial positions')
    call h5_add(h5id_root, grp // '/rsmall', &
         profile%rsmall, lbound(profile%rsmall), ubound(profile%rsmall), &
         comment = 'minor radius (of equal-area circle)', unit = 'cm')
    call h5_add(h5id_root, grp // '/prof', &
         profile%prof, lbound(profile%prof), ubound(profile%prof), &
         comment = comment, unit = unit)
    call h5_close(h5id_root)
  end subroutine profile_export_hdf5

  subroutine profile_import_hdf5(profile, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(profile_t), intent(inout) :: profile
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer :: nrad

    grp = trim(group)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/nrad', profile%nrad)
    call profile_init(dens, nrad)
    call h5_get(h5id_root, grp // '/rsmall', profile%rsmall)
    call h5_get(h5id_root, grp // '/prof', profile%prof)
    call h5_close(h5id_root)
  end subroutine profile_import_hdf5

  subroutine read_profiles
    use mephit_util, only: ev2erg
    integer :: krad

    call profile_read(dens, 'n.dat')
    call profile_read(e_temp, 'Te.dat')
    call profile_read(i_temp, 'Ti.dat')
    call profile_read(E_r, 'Er.dat')
    ! convert temperatures from electronvolt to erg
    e_temp%prof(:) = e_temp%prof(:) * ev2erg
    i_temp%prof(:) = i_temp%prof(:) * ev2erg
    ! integrate electric field to yield the electric potential
    call profile_init(Phi0, E_r%nrad)
    Phi0%prof(1) = 0d0
    do krad = 2, Phi0%nrad
       Phi0%prof(krad) = Phi0%prof(krad - 1) + (E_r%rsmall(krad) - E_r%rsmall(krad - 1)) &
            & * 0.5d0 * (E_r%prof(krad) + E_r%prof(krad - 1))
    end do
  end subroutine read_profiles

  subroutine write_profiles_hdf5(file, group)
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call profile_export_hdf5(dens, file, grp // '/dens', &
         'density profile', 'cm^-3')
    call profile_export_hdf5(e_temp, file, grp // '/e_temp', &
         'electron temperature profile', 'eV')
    call profile_export_hdf5(i_temp, file, grp // '/i_temp', &
         'ion temperature profile', 'eV')
    call profile_export_hdf5(E_r, file, grp // '/E_r', &
         'radial electric field profile', 'statV cm^-1')
    call profile_export_hdf5(Phi0, file, grp // '/Phi0', &
         'electric potential profile', 'statV')
  end subroutine write_profiles_hdf5

  subroutine read_profiles_hdf5(file, group)
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call profile_import_hdf5(dens, file, grp // '/dens')
    call profile_import_hdf5(e_temp, file, grp // '/e_temp')
    call profile_import_hdf5(i_temp, file, grp // '/i_temp')
    call profile_import_hdf5(E_r, file, grp // '/E_r')
    call profile_import_hdf5(Phi0, file, grp // '/Phi0')
  end subroutine read_profiles_hdf5

end module mephit_equil
