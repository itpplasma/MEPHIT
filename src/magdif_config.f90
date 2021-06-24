module magdif_conf
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: magdif_config, magdif_config_read, magdif_config_export_hdf5, magdif_config_delayed, &
       magdif_log, conf, conf_arr, log, runmode_single, runmode_direct, runmode_precon, pres_prof_eps, &
       pres_prof_par, pres_prof_geqdsk, curr_prof_ps, curr_prof_rot, curr_prof_geqdsk, &
       q_prof_flux, q_prof_rot, q_prof_geqdsk, vac_src_nemov, vac_src_gpec, vac_src_fourier, &
       decorate_filename, longlines, cmplx_fmt, nl_fmt, datafile

  character(len = *), parameter :: cmplx_fmt = 'es24.16e3, 1x, sp, es24.16e3, s, " i"'
  character(len = *), parameter :: nl_fmt = '"' // new_line('A') // '"'
  character(len = *), parameter :: datafile = 'magdif.h5'

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

  type :: magdif_config

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

     !> Number of iterations. Does not apply when #runmode equals #runmode_single.
     !> Defaults to 20.
     integer :: niter = 20

     !> Number of Ritz eigenvalues to calculate. Only applies when #runmode equals
     !> #runmode_precon. Defaults to 30.
     integer :: nritz = 30

     !> Threshold for absolute eigenvalues used for Arnoldi preconditioner. Defaults to
     !> 0.5.
     real(dp) :: ritz_threshold = 0.5d0

     !> Index of toroidal harmonics of perturbation. Defaults to 2.
     integer :: n = 2

     !> Number of knots per poloidal loop. Defaults to 300.
     integer :: nkpol = 300

     !> Number of flux surfaces before refinement. Defaults to 100.
     integer :: nflux_unref = 100

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

     !> Error threshold for divergence-freeness of #bnflux and #bnphi. Defaults to 10^-7.
     real(dp) :: rel_err_Bn = 1d-7

     !> Error threshold for divergence-freeness of #jnflux and #jnphi. Defaults to 10^-8.
     real(dp) :: rel_err_currn = 1d-8

     !> Single poloidal mode used in comparison with KiLCA code. Defaults to 0 (ASDEX
     !> geometry).
     integer :: kilca_pol_mode = 0

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

  end type magdif_config

  type :: magdif_config_delayed

     integer :: m_min, m_max

     !> Number of unrefined flux surfaces to be replaced by refined ones.
     integer, dimension(:), allocatable :: deletions

     !> Number of refined flux surfaces.
     integer, dimension(:), allocatable :: additions

     !> Relative size of most refined flux surface.
     real(dp), dimension(:), allocatable :: refinement

     !> Free parameters setting the magnitudes of sheet currents.
     complex(dp), dimension(:), allocatable :: sheet_current_factor

     !> Integration constant for resonant vacuum perturbation in KiLCA comparison.
     complex(dp), dimension(:), allocatable :: kilca_vac_coeff

     !> Data point for resonant vacuum perturbation in KiLCA comparison.
     real(dp), dimension(:), allocatable :: kilca_vac_r

     !> Data point for resonant vacuum perturbation in KiLCA comparison.
     complex(dp), dimension(:), allocatable :: kilca_vac_Bz

   contains
     procedure :: read => magdif_config_delayed_read
     procedure :: export_hdf5 => magdif_config_delayed_export_hdf5
     procedure :: import_hdf5 => magdif_config_delayed_import_hdf5
     final :: magdif_config_delayed_destructor
  end type magdif_config_delayed

  type :: magdif_log
     private
     character(len = 1024), public :: msg = ''
     character(len = 1024), public :: filename = ''
     integer :: fid
     integer, public :: level
     logical, public :: err, warn, info, debug, quiet
   contains
     procedure :: write => magdif_log_write
     procedure :: msg_arg_size => magdif_log_msg_arg_size
     final :: magdif_log_close
  end type magdif_log

  interface magdif_log
     module procedure magdif_log_open
  end interface magdif_log

  type(magdif_config) :: conf
  type(magdif_config_delayed) :: conf_arr
  type(magdif_log) :: log

contains

  !> Read static configuration file for magdif.
  subroutine magdif_config_read(config, filename)
    type(magdif_config), intent(inout) :: config
    character(len = *), intent(in) :: filename
    integer :: fid
    namelist /scalars/ config

    open(newunit = fid, file = filename)
    read(fid, nml = scalars)
    close(fid)
    ! override if erroneously set in namelist
    config%config_file = trim(filename)
  end subroutine magdif_config_read

  subroutine magdif_config_export_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(magdif_config), intent(in) :: config
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
         comment = 'number of iterations')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nritz', config%nritz, &
         comment = 'number of Ritz eigenvalues')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ritz_threshold', config%ritz_threshold, &
         comment = 'threshold for absolute eigenvalues used for Arnoldi preconditioner')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n', config%n, &
         comment = 'index of toroidal harmonics of perturbation')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nkpol', config%nkpol, &
         comment = 'knots per poloidal loop')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nflux_unref', config%nflux_unref, &
         comment = 'number of flux surfaces before refinement')
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rel_err_Bn', config%rel_err_Bn, &
         comment = 'error threshold for divergence-freeness of perturbation fields', unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rel_err_currn', config%rel_err_currn, &
         comment = 'error threshold for divergence-freeness of perturbation currents', unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_pol_mode', config%kilca_pol_mode, &
         comment = 'single poloidal mode used in comparison with KiLCA code')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_scale_factor', config%kilca_scale_factor, &
         comment = 'scaling factor used for comparison with results from KiLCA code')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/max_eig_out', config%max_eig_out, &
         comment = 'maximum number of exported eigenvectors')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/debug_pol_offset', config%debug_pol_offset)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/debug_kilca_geom_theta', config%debug_kilca_geom_theta)
    call h5_close(h5id_root)
  end subroutine magdif_config_export_hdf5

  !> Read dynamic configuration from configuration file for magdif.
  subroutine magdif_config_delayed_read(config, filename, m_min, m_max)
    class(magdif_config_delayed), intent(inout) :: config
    character(len = *), intent(in) :: filename
    integer, intent(in) :: m_min, m_max
    integer :: fid
    integer, dimension(m_min:m_max) :: deletions, additions
    real(dp), dimension(m_min:m_max) :: refinement, kilca_vac_r
    complex(dp), dimension(m_min:m_max) :: kilca_vac_Bz, kilca_vac_coeff, &
         sheet_current_factor
    namelist /arrays/ deletions, additions, refinement, sheet_current_factor, &
         kilca_vac_coeff, kilca_vac_r, kilca_vac_Bz

    config%m_min = m_min
    config%m_max = m_max
    deletions = 0
    additions = 0
    refinement = 0d0
    sheet_current_factor = (0d0, 0d0)
    kilca_vac_coeff = (1d0, 0d0)
    kilca_vac_r = 0d0
    kilca_vac_Bz = (0d0, 0d0)
    open(newunit = fid, file = filename)
    read(fid, nml = arrays)
    close(fid)
    if (allocated(config%deletions)) deallocate(config%deletions)
    allocate(config%deletions(m_min:m_max))
    config%deletions(:) = deletions
    if (allocated(config%additions)) deallocate(config%additions)
    allocate(config%additions(m_min:m_max))
    config%additions(:) = additions
    if (allocated(config%refinement)) deallocate(config%refinement)
    allocate(config%refinement(m_min:m_max))
    config%refinement(:) = refinement
    if (allocated(config%sheet_current_factor)) deallocate(config%sheet_current_factor)
    allocate(config%sheet_current_factor(m_min:m_max))
    config%sheet_current_factor(:) = sheet_current_factor
    if (allocated(config%kilca_vac_coeff)) deallocate(config%kilca_vac_coeff)
    allocate(config%kilca_vac_coeff(m_min:m_max))
    config%kilca_vac_coeff(:) = kilca_vac_coeff
    if (allocated(config%kilca_vac_r)) deallocate(config%kilca_vac_r)
    allocate(config%kilca_vac_r(m_min:m_max))
    config%kilca_vac_r(:) = kilca_vac_r
    if (allocated(config%kilca_vac_Bz)) deallocate(config%kilca_vac_Bz)
    allocate(config%kilca_vac_Bz(m_min:m_max))
    config%kilca_vac_Bz(:) = kilca_vac_Bz
  end subroutine magdif_config_delayed_read

  subroutine magdif_config_delayed_export_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(magdif_config_delayed), intent(in) :: config
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/additions', config%additions, &
         lbound(config%additions), ubound(config%additions), &
         comment = 'number of refined flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/refinement', config%refinement, &
         lbound(config%refinement), ubound(config%refinement), &
         comment = 'relative size of most refined flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/sheet_current_factor', config%sheet_current_factor, &
         lbound(config%sheet_current_factor), ubound(config%sheet_current_factor), &
         comment = 'free parameters setting the magnitudes of sheet currents')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_coeff', config%kilca_vac_coeff, &
         lbound(config%kilca_vac_coeff), ubound(config%kilca_vac_coeff), &
         comment = 'integration constant for resonant vacuum perturbation in KiLCA comparison')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_r', config%kilca_vac_r, &
         lbound(config%kilca_vac_r), ubound(config%kilca_vac_r), &
         comment = 'data point for resonant vacuum perturbation in KiLCA comparison')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_Bz', config%kilca_vac_Bz, &
         lbound(config%kilca_vac_Bz), ubound(config%kilca_vac_Bz), &
         comment = 'data point for resonant vacuum perturbation in KiLCA comparison')
    call h5_close(h5id_root)
  end subroutine magdif_config_delayed_export_hdf5

  subroutine magdif_config_delayed_import_hdf5(config, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(magdif_config_delayed), intent(inout) :: config
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call magdif_config_delayed_destructor(config)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_min', config%m_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_max', config%m_max)
    allocate(config%deletions(config%m_min:config%m_max))
    allocate(config%additions(config%m_min:config%m_max))
    allocate(config%refinement(config%m_min:config%m_max))
    allocate(config%sheet_current_factor(config%m_min:config%m_max))
    allocate(config%kilca_vac_coeff(config%m_min:config%m_max))
    allocate(config%kilca_vac_r(config%m_min:config%m_max))
    allocate(config%kilca_vac_Bz(config%m_min:config%m_max))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/deletions', config%deletions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/additions', config%additions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/refinement', config%refinement)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/sheet_current_factor', config%sheet_current_factor)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_coeff', config%kilca_vac_coeff)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_r', config%kilca_vac_r)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kilca_vac_Bz', config%kilca_vac_Bz)
    call h5_close(h5id_root)
  end subroutine magdif_config_delayed_import_hdf5

  subroutine magdif_config_delayed_destructor(config)
    type(magdif_config_delayed), intent(inout) :: config
    if (allocated(config%deletions)) deallocate(config%deletions)
    if (allocated(config%additions)) deallocate(config%additions)
    if (allocated(config%refinement)) deallocate(config%refinement)
    if (allocated(config%sheet_current_factor)) deallocate(config%sheet_current_factor)
    if (allocated(config%kilca_vac_coeff)) deallocate(config%kilca_vac_coeff)
    if (allocated(config%kilca_vac_r)) deallocate(config%kilca_vac_r)
    if (allocated(config%kilca_vac_Bz)) deallocate(config%kilca_vac_Bz)
  end subroutine magdif_config_delayed_destructor

  !> Associate logfile and open if necessary.
  function magdif_log_open(filename, log_level, quiet) result(log)
    use iso_fortran_env, only: output_unit
    character(len = *), intent(in) :: filename
    integer, intent(in) :: log_level
    logical, intent(in) :: quiet
    type(magdif_log) :: log

    log%filename = filename
    log%quiet = quiet
    log%level = log_level
    log%err = .false.
    log%warn = .false.
    log%info = .false.
    log%debug = .false.
    if (log%level > 0) log%err = .true.
    if (log%level > 1) log%warn = .true.
    if (log%level > 2) log%info = .true.
    if (log%level > 3) log%debug = .true.
    if (trim(log%filename) == '-') then
       log%fid = output_unit
    else
       open(newunit = log%fid, file = log%filename, status = 'replace', form = 'formatted', action = 'write')
    end if
  end function magdif_log_open

  !> Close logfile if necessary.
  subroutine magdif_log_close(log)
    use iso_fortran_env, only: output_unit
    type(magdif_log), intent(inout) :: log
    if (log%fid /= output_unit) close(log%fid)
  end subroutine magdif_log_close

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
  subroutine magdif_log_write(log)
    use iso_fortran_env, only: output_unit
    class(magdif_log), intent(inout) :: log
    character(len = 28 + len_trim(log%msg)) :: full_msg
    write (full_msg, '("[", a25, "] ", a)') timestamp(), trim(log%msg)
    write (log%fid, '(a)') full_msg
    if (log%fid /= output_unit .and. .not. log%quiet) write (output_unit, '(a)') full_msg
    flush(log%fid)
  end subroutine magdif_log_write

  subroutine magdif_log_msg_arg_size(log, funcname, name1, name2, value1, value2)
    class(magdif_log), intent(inout) :: log
    character(len = *), intent(in) :: funcname, name1, name2
    integer, intent(in) :: value1, value2
    write (log%msg, '("Argument size mismatch in ", a, ": ", a, " = ", i0, ' // &
         '", ", a, " = ", i0)') funcname, name1, value1, name2, value2
  end subroutine magdif_log_msg_arg_size

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
end module magdif_conf
