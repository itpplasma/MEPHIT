module mephit_conf
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  ! types and associated procedures
  public :: config_t, config_read, config_export_hdf5
  public :: config_delayed_t
  public :: logger_t

  ! utility procedures
  public :: decorate_filename

  ! module variables
  public :: conf, conf_arr, logger
  public :: longlines, cmplx_fmt, nl_fmt, datafile, &
       runmode_single, runmode_direct, runmode_precon, &
       pres_prof_eps, pres_prof_par, pres_prof_geqdsk, &
       curr_prof_ps, curr_prof_rot, curr_prof_geqdsk, &
       q_prof_flux, q_prof_rot, q_prof_geqdsk, &
       vac_src_nemov, vac_src_gpec, vac_src_fourier

  character(len = *), parameter :: cmplx_fmt = 'es24.16e3, 1x, sp, es24.16e3, s, " i"'
  character(len = *), parameter :: nl_fmt = '"' // new_line('A') // '"'
  character(len = *), parameter :: datafile = 'mephit.h5'

  integer, parameter :: longlines = 1024

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer, parameter :: pres_prof_eps = 0    !< pressure profile from EPS paper
  integer, parameter :: pres_prof_par = 1    !< parabolic pressure profile
  integer, parameter :: pres_prof_geqdsk = 2 !< pressure profile from G EQDSK file

  integer, parameter :: curr_prof_ps = 0     !< only Pfirsch-Schlueter current
  integer, parameter :: curr_prof_rot = 1    !< current via Ampere's law
  integer, parameter :: curr_prof_geqdsk = 2 !< current profile from G EQDSK file

  integer, parameter :: q_prof_flux = 0   !< q profile from flux between flux surfaces
  integer, parameter :: q_prof_rot = 1    !< q profile from rotational transform
  integer, parameter :: q_prof_geqdsk = 2 !< q profile from G EQDSK file

  integer, parameter :: vac_src_nemov = 0   !< vacuum field perturbation from Viktor Nemov's code
  integer, parameter :: vac_src_gpec = 1    !< vacuum field perturbation from GPEC
  integer, parameter :: vac_src_fourier = 2 !< vacuum field perturbation from precomputed Fourier modes

  type :: config_t

     !> Name of the file the configuration is read from.
     character(len = 1024) :: config_file

     !> Verbosity of logging. Possible values are 0 (none), 1 (errors), 2 (warnings),
     !> 3 (info), and 4 (debug, default).
     integer :: log_level = 4

     !> Suppress log messages to stdout. Defaults to false.
     logical :: quiet = .false.

     !> Kind of iterations to run. Possible values are #runmode_single, #runmode_direct,
     !> and #runmode_precon (default).
     integer :: runmode = runmode_precon

     !> Source of pressure profile. Possible values are #pres_prof_eps, #pres_prof_par,
     !> and #pres_prof_geqdsk (default).
     integer :: pres_prof = pres_prof_geqdsk

     !> Source of toroidal equilibrium current. Possible values are #curr_prof_ps,
     !> #curr_prof_rot, and #curr_prof_geqdsk (default).
     integer :: curr_prof = curr_prof_geqdsk

     !> Source of safety factor profile. Possivle values are #q_prof_flux, #q_prof_rot
     !> (default), and #q_prof_geqdsk
     integer :: q_prof = q_prof_rot

     !> Source of vacuum field perturbation. Possible values are #vac_src_nemov (default)
     !> and #vac_src_gpec.
     integer :: vac_src = vac_src_nemov

     !> Generate non-resonant vacuum perturbation for testing. Defaults to false.
     logical :: nonres = .false.

     !> Average over quadrilaterals for non-resonant test case. Defaults to true.
     logical :: quad_avg = .true.

     !> Coil currents for AUG B coils; upper coils come first
     real(dp), dimension(16) :: Ic = 0d0

     !> Maximum number of iterations. Does not apply when #runmode equals
     !> #runmode_single. Defaults to 50.
     integer :: niter = 50

     !> Acceptable relative error to consider L^2 integral in fixed-point
     !> iterations converged. Defaults to 1.0e-12.
     real(dp) :: iter_rel_err = 1d-12

     !> Maxmimum number of iterations in Arnoldi method. Only applies when
     !> #runmode equals #runmode_precon. Defaults to 100.
     integer :: nkrylov = 100

     !> Threshold for eigenvalues' magnitude calculated via Arnoldi method.
     !> Only applies when #runmode equals #runmode_precon. Defaults to 0.5.
     real(dp) :: ritz_threshold = 0.5d0

     !> Acceptable relative error to consider eigenvalues converged in Arnoldi method.
     !> Only applies when #runmode equals #runmode_precon. Defaults to 1.0e-12.
     real(dp) :: ritz_rel_err = 1d-12

     !> Index of toroidal harmonics of perturbation. Defaults to 2.
     integer :: n = 2

     !> Maximum distance between flux surfaces along \f$ \theta = 0 \f$. Defaults
     !> to 0.45 cm.
     real(dp) :: max_Delta_rad = 0.2d0

     !> Number of radial divisions in radially refined regions for parallel current
     !> computations. Defaults to 2048.
     integer :: nrad_Ipar = 2048

     !> Minimum temperature. Does not apply when #pres_prof equals #pres_prof_geqdsk.
     !> Defaults to 20 eV.
     real(dp) :: temp_min = 2d1

     ! Minimum density. Does not apply when #pres_prof equals #pres_prof_geqdsk. Defaults
     !> to 2.0e+11 cm^-3
     real(dp) :: dens_min = 2d11

     !> Temperature at magnetic axis. Does not apply when #pres_prof equals
     !> #pres_prof_geqdsk. Defaults to 3.0e+03 eV.
     real(dp) :: temp_max = 3d3

     !> Density at magnetic axis. Does not apply when #pres_prof equals #pres_prof_geqdsk.
     !> Defaults to 5.0e+13 cm^-3.
     real(dp) :: dens_max = 5d13

     !> Damping factor for resonances. Defaults to 0.
     real(dp) :: damp = 0d0

     !> Single poloidal mode used in comparison with KiLCA code. Defaults to 0 (ASDEX
     !> geometry).
     integer :: kilca_pol_mode = 0

     !> KiLCA HDF5 output file from which to infer coefficients of vacuum perturbation.
     character(len = 1024) :: kilca_vac_output = ''

     !> Scaling factor used for comparison with results from KiLCA code. Defaults to 0
     !> (ASDEX geometry).
     !>
     !> If non-zero, #n and #r0 are scaled by this factor to approximate the cylindrical
     !> topology assumed in KiLCA.
     integer :: kilca_scale_factor = 0

     !> Maximum number of eigenvectors to be dumped in #eigvec_file. Defaults to 10.
     integer :: max_eig_out = 10

     real(dp) :: debug_pol_offset = 0.5d0
     logical :: debug_kilca_geom_theta = .false.

  end type config_t

  type :: config_delayed_t

     integer :: m_min, m_max

     !> Number of unrefined flux surfaces to be replaced by refined ones.
     integer, dimension(:), allocatable :: deletions

     !> Width ratio of neighbouring refined flux surfaces.
     real(dp), dimension(:), allocatable :: refinement

     !> Free parameters setting the magnitudes of sheet currents.
     real(dp), dimension(:), allocatable :: sheet_current_factor

   contains
     procedure :: read => config_delayed_read
     procedure :: export_hdf5 => config_delayed_export_hdf5
     procedure :: import_hdf5 => config_delayed_import_hdf5
     procedure :: deinit => config_delayed_deinit
  end type config_delayed_t

  type :: logger_t
     private
     character(len = 1024), public :: msg = ''
     character(len = 1024), public :: filename = ''
     integer :: fid
     integer, public :: level
     logical, public :: err, warn, info, debug, quiet
   contains
     procedure :: init => logger_init
     procedure :: write_msg => logger_write_msg
     procedure :: msg_arg_size => logger_msg_arg_size
     procedure :: deinit => logger_deinit
  end type logger_t

  type(config_t) :: conf
  type(config_delayed_t) :: conf_arr
  type(logger_t) :: logger

contains

  !> Read configuration namelist.
  subroutine config_read(config, filename)
    type(config_t), intent(inout) :: config
    character(len = *), intent(in) :: filename
    integer :: fid
    namelist /scalars/ config

    open(newunit = fid, file = filename)
    read(fid, nml = scalars)
    close(fid)
    ! override if erroneously set in namelist
    config%config_file = trim(filename)
  end subroutine config_read

  subroutine config_export_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(config_t), intent(in) :: config
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/runmode', config%runmode)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/pres_prof', config%pres_prof)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/curr_prof', config%curr_prof)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/q_prof', config%q_prof)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/vac_src', config%vac_src)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nonres', config%nonres)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/quad_avg', config%quad_avg)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Ic', config%Ic, &
         lbound(config%Ic), ubound(config%Ic), unit = 'c_0 statA', &
         comment = 'coil currents for AUG B coils; upper coils come first')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/niter', config%niter, &
         comment = 'maximum number of iterations')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/iter_rel_err', config%iter_rel_err, &
         comment = 'convergence threshold in L2 integral in fixed-point iteration')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nkrylov', config%nkrylov, &
         comment = 'maximum number of iterations in Arnoldi method')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ritz_threshold', config%ritz_threshold, &
         comment = 'threshold for eigenvalues'' magnitude in Arnoldi method')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ritz_rel_err', config%ritz_rel_err, &
         comment = 'relative error for eigenvalues in Arnoldi method')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n', config%n, &
         comment = 'index of toroidal harmonics of perturbation')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/max_Delta_rad', config%max_Delta_rad, &
         comment = 'maximum distance between flux surfaces along theta = 0', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nrad_Ipar', config%nrad_Ipar, &
         comment = 'radial divisions in radially refined regions for parallel current computations')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/temp_min', config%temp_min, &
         comment = 'minimum temperature', unit = 'eV')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dens_min', config%dens_min, &
         comment = 'minimum density', unit = 'cm^-3')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/temp_max', config%temp_max, &
         comment = 'maximum temperature', unit = 'eV')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dens_max', config%dens_max, &
         comment = 'maximum density', unit = 'cm^-3')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/damp', config%damp, &
         comment = 'damping factor for resonances', unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_pol_mode', config%kilca_pol_mode, &
         comment = 'single poloidal mode used in comparison with KiLCA code')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_scale_factor', config%kilca_scale_factor, &
         comment = 'scaling factor used for comparison with results from KiLCA code')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/max_eig_out', config%max_eig_out, &
         comment = 'maximum number of exported eigenvectors')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/debug_pol_offset', config%debug_pol_offset)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/debug_kilca_geom_theta', config%debug_kilca_geom_theta)
    call h5_close(h5id_root)
  end subroutine config_export_hdf5

  !> Read configuration arrays from namelist.
  subroutine config_delayed_read(config, filename, m_min, m_max)
    class(config_delayed_t), intent(inout) :: config
    character(len = *), intent(in) :: filename
    integer, intent(in) :: m_min, m_max
    integer :: fid
    integer, dimension(m_min:m_max) :: deletions
    real(dp), dimension(m_min:m_max) :: refinement, sheet_current_factor
    namelist /arrays/ deletions, refinement, sheet_current_factor

    config%m_min = m_min
    config%m_max = m_max
    deletions = 0
    refinement = 0d0
    sheet_current_factor = 0d0
    open(newunit = fid, file = filename)
    read(fid, nml = arrays)
    close(fid)
    if (allocated(config%deletions)) deallocate(config%deletions)
    allocate(config%deletions(m_min:m_max))
    config%deletions(:) = deletions
    if (allocated(config%refinement)) deallocate(config%refinement)
    allocate(config%refinement(m_min:m_max))
    config%refinement(:) = refinement
    if (allocated(config%sheet_current_factor)) deallocate(config%sheet_current_factor)
    allocate(config%sheet_current_factor(m_min:m_max))
    config%sheet_current_factor(:) = sheet_current_factor
  end subroutine config_delayed_read

  subroutine config_delayed_export_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(config_delayed_t), intent(in) :: config
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_min', config%m_min, &
         comment = 'minimum poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_max', config%m_max, &
         comment = 'maximum poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/deletions', config%deletions, &
         lbound(config%deletions), ubound(config%deletions), &
         comment = 'number of unrefined flux surfaces to be replaced by refined ones')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/refinement', config%refinement, &
         lbound(config%refinement), ubound(config%refinement), &
         comment = 'width ratio of neighbouring refined flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/sheet_current_factor', config%sheet_current_factor, &
         lbound(config%sheet_current_factor), ubound(config%sheet_current_factor), &
         comment = 'free parameters setting the magnitudes of sheet currents')
    call h5_close(h5id_root)
  end subroutine config_delayed_export_hdf5

  subroutine config_delayed_import_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(config_delayed_t), intent(inout) :: config
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call config_delayed_deinit(config)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_min', config%m_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_max', config%m_max)
    allocate(config%deletions(config%m_min:config%m_max))
    allocate(config%refinement(config%m_min:config%m_max))
    allocate(config%sheet_current_factor(config%m_min:config%m_max))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/deletions', config%deletions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/refinement', config%refinement)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/sheet_current_factor', config%sheet_current_factor)
    call h5_close(h5id_root)
  end subroutine config_delayed_import_hdf5

  subroutine config_delayed_deinit(config)
    class(config_delayed_t), intent(inout) :: config

    config%m_min = 0
    config%m_max = 0
    if (allocated(config%deletions)) deallocate(config%deletions)
    if (allocated(config%refinement)) deallocate(config%refinement)
    if (allocated(config%sheet_current_factor)) deallocate(config%sheet_current_factor)
  end subroutine config_delayed_deinit

  !> Associate logfile and open if necessary.
  subroutine logger_init(this, filename, log_level, quiet)
    use iso_fortran_env, only: output_unit
    class(logger_t), intent(inout) :: this
    character(len = *), intent(in) :: filename
    integer, intent(in) :: log_level
    logical, intent(in) :: quiet

    this%filename = filename
    this%quiet = quiet
    this%level = log_level
    this%err = .false.
    this%warn = .false.
    this%info = .false.
    this%debug = .false.
    if (this%level > 0) this%err = .true.
    if (this%level > 1) this%warn = .true.
    if (this%level > 2) this%info = .true.
    if (this%level > 3) this%debug = .true.
    if (trim(this%filename) == '-') then
       this%fid = output_unit
    else
       open(newunit = this%fid, file = this%filename, status = 'replace', form = 'formatted', action = 'write')
    end if
  end subroutine logger_init

  !> Close logfile if necessary.
  subroutine logger_deinit(this)
    use iso_fortran_env, only: output_unit
    class(logger_t), intent(inout) :: this

    if (this%fid /= output_unit) close(this%fid)
  end subroutine logger_deinit

  !> Generate timestamp
  function timestamp()
    character(len=25) :: timestamp
    character(len=*), parameter :: ISO_8601 = &
         '(ss, i4, 2("-", i2.2), "T", i2.2, 2(":", i2.2), sp, i3.2, ":", ss, i2.2, s)'
    integer :: values(8)
    call date_and_time(values = values)
    write (timestamp, ISO_8601) values(1:3), values(5:7), &
         values(4) / 60, mod(values(4), 60)
  end function timestamp

  !> Write to logfile and flush if necessary.
  subroutine logger_write_msg(this)
    use iso_fortran_env, only: output_unit
    class(logger_t), intent(inout) :: this
    character(len = 28 + len_trim(logger%msg)) :: full_msg
    write (full_msg, '("[", a25, "] ", a)') timestamp(), trim(this%msg)
    write (this%fid, '(a)') full_msg
    if (this%fid /= output_unit .and. .not. this%quiet) write (output_unit, '(a)') full_msg
    flush(this%fid)
  end subroutine logger_write_msg

  subroutine logger_msg_arg_size(this, funcname, name1, name2, value1, value2)
    class(logger_t), intent(inout) :: this
    character(len = *), intent(in) :: funcname, name1, name2
    integer, intent(in) :: value1, value2
    write (this%msg, '("Argument size mismatch in ", a, ": ", a, " = ", i0, ' // &
         '", ", a, " = ", i0)') funcname, name1, value1, name2, value2
  end subroutine logger_msg_arg_size

  !> Generates a new filename from a given template.
  !>
  !> @param in_name undecorated filename
  !> @param prefix string to be prefixed to the file basename
  !> @param postfix string to be postfixed to the file basename
  !>
  !> No attempt is made to check whether the given filename is valid or accessible.
  function decorate_filename(in_name, prefix, postfix) result(out_name)
    character(len = *), intent(in) :: in_name
    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: postfix
    character(len = len_trim(in_name) + len_trim(prefix) + len_trim(postfix)) :: out_name

    integer :: basename_start, basename_end, name_end

    name_end = len_trim(in_name)
    basename_start = scan(in_name, '/', .true.) + 1
    basename_end = basename_start - 2 + scan(in_name(basename_start:), '.', .true.)
    if (basename_end - basename_start == -1) basename_start = basename_start + 1 ! dotfile
    if (basename_end < basename_start) basename_end = name_end
    out_name = in_name(:basename_start-1) // prefix // in_name(basename_start: &
         basename_end) // postfix // in_name(basename_end+1:name_end)
  end function decorate_filename
end module mephit_conf
