module magdif_config
  use from_nrtype, only: dp                                     ! PRELOAD/SRC/from_nrtype.f90
  
  implicit none

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer log_level
  integer runmode
  logical :: log_err, log_warn, log_info, log_debug ! specify log levels
  logical :: nonres = .false.  !< use non-resonant test case

  character(len=1024) :: point_file   !< input data file for mesh points
  character(len=1024) :: tri_file     !< input data file for triangles and edges
  character(len=1024) :: Bnflux_file  !< input data file for magnetic field perturbation
  character(len=1024) :: hpsi_file    !< input data file for \f$ h_{n}^{\psi} \f$
  character(len=1024) :: config_file  !< input config file for namelist settings
  character(len=1024) :: presn_file   !< output data file for pressure perturbation
  character(len=1024) :: currn_file   !< output data file for current perturbation

  integer  :: n               !< harmonic index of perturbation
  integer  :: nkpol           !< number of knots per poloidal loop
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density

  namelist / settings / log_level, runmode, nonres, point_file, tri_file, Bnflux_file,&
       hpsi_file, presn_file, currn_file, n, nkpol, nflux, ti0, di0
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

program magdif_test

  use magdif_config, only: config_file, read_config, runmode, runmode_single,&
       runmode_direct, runmode_precon
  use magdif, only: magdif_init, magdif_cleanup, &
       compute_presn, compute_j0phi, compute_currn, read_bnflux, Bnflux_file, currn_file

  implicit none

  if (command_argument_count() == 1) then
     call get_command_argument(1, config_file)
  else
     stop 'Error: expected path to config file as first parameter'
  endif

  call read_config(config_file)
  call magdif_init
  if (runmode == runmode_direct) call magdif_single
  if (runmode == runmode_direct) call magdif_direct
  call magdif_cleanup

contains

  subroutine magdif_single
    call compute_presn
    call compute_j0phi
    call compute_currn
  end subroutine magdif_single
   
  subroutine magdif_direct
    integer, parameter :: niter = 1 ! TODO: make this configurable
    integer :: kiter
    integer :: stat

    ! TODO: make this configurable over tmpdir parameter
    currn_file = '/tmp/currents.dat'
    Bnflux_file = '/tmp/bnflux.dat'
    do kiter=1,niter
       call compute_presn       ! compute pressure based on previous perturbation field
       call compute_j0phi       ! compute currents based on previous perturbation field
       call compute_currn 
       call compute_field(stat) ! use field code to generate new field from currents
       call read_bnflux         ! read new bnflux from field code
    end do   
  end subroutine magdif_direct

  subroutine compute_field(stat)
    integer, intent(out) :: stat
    
    call execute_command_line (&
         "PATH=/temp/ert/local/bin:$PATH /temp/ert/local/bin/FreeFem++ "//&
         "../FEM/maxwell.edp ../PRELOAD/inputformaxwell_ext.msh /tmp/currents.dat 2 > /tmp/freefem.out 2>&1; cd ..",&
         exitstat=stat)
    
  end subroutine compute_field
end program magdif_test
