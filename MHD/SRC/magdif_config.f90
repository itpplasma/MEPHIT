module magdif_config
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90

  implicit none

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer, parameter :: logfile = 6             !< log to stdout, TODO: make this configurable

  integer :: log_level
  integer :: runmode
  logical :: log_err, log_warn, log_info, log_debug ! specify log levels
  logical :: nonres = .false.  !< use non-resonant test case

  character(len=1024) :: point_file       !< input data file for mesh points
  character(len=1024) :: tri_file         !< input data file for triangles and edges
  character(len=1024) :: Bnflux_file      !< input data file for magnetic field perturbation
  character(len=1024) :: Bnflux_vac_file  !< input data file for magnetic field perturbation
  character(len=1024) :: Bn_sum_file      !< output data file for accumulated magnetic field perturbation
  character(len=1024) :: hpsi_file        !< input data file for \f$ h_{n}^{\psi} \f$
  character(len=1024) :: config_file      !< input config file for namelist settings
  character(len=1024) :: presn_file       !< output data file for pressure perturbation
  character(len=1024) :: currn_file       !< output data file for current perturbation

  integer  :: niter = 20
  integer  :: n               !< harmonic index of perturbation
  integer  :: nkpol           !< number of knots per poloidal loop
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density

  namelist / settings / log_level, runmode, nonres, point_file, tri_file, Bnflux_file, &
       Bnflux_vac_file, Bn_sum_file, hpsi_file, presn_file, currn_file, niter, n, nkpol, &
       nflux, ti0, di0
  !< namelist for input parameters

contains

  !> Read configuration file for magdif
  !! @param config_filename file name of config file
  subroutine read_config(config_filename)
    character(len = *) :: config_filename

    open(1, file = config_filename)
    read(1, nml = settings)
    close(1)

    log_err = .false.
    log_warn = .false.
    log_info = .false.
    log_debug = .false.
    if (log_level > 0) log_err = .true.
    if (log_level > 1) log_warn = .true.
    if (log_level > 2) log_info = .true.
    if (log_level > 3) log_debug = .true.
  end subroutine read_config
end module magdif_config
