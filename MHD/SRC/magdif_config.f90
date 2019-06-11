module magdif_config
  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90
  use arnoldi_mod, only: tol ! RUN/SRC/arnoldi.f90

  implicit none

  integer, parameter :: longlines = 1024

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer, parameter :: logfile = 6             !< log to stdout, TODO: make this configurable

  integer :: log_level
  integer :: runmode
  logical :: log_err, log_warn, log_info, log_debug ! specify log levels
  logical :: nonres = .false.  !< use non-resonant test case
  logical :: quad_avg = .true. !< average over quadrilaterals for non-resonant test case

  integer  :: niter = 20      !< number of iterations
  integer  :: nritz = 20      !< number of Ritz eigenvalues
  integer  :: n               !< harmonic index of perturbation
  integer  :: nkpol           !< number of knots per poloidal loop
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density
  real(dp) :: damp            !< damping factor for resonances

  !> Unformatted input data file containing mesh points.
  !>
  !> This file is generated in PRELOAD and is read in by magdif::read_mesh().
  character(len = 1024) :: point_file

  !> Unformatted input data file containing data on triangles and edges.
  !>
  !> This file is generated in PRELOAD and is read in by magdif::read_mesh().
  character(len = 1024) :: tri_file

  !> Formatted output data file for magnetic field perturbation (vacuum).
  !>
  !> This file is generated by magdif::init().
  character(len = 1024) :: Bn_vacout_file

  !> Formatted input data file for magnetic field perturbation (vacuum).
  !>
  !> This file is generated in PRELOAD and is read in by magdif::init() via
  !> magdif::read_bn().
  character(len = 1024) :: Bn_vac_file

  !> Formatted data file for accumulated magnetic field perturbation.
  !>
  !> This file is generated in every iteration in magdif::magdif_direct() and
  !> contains the final solution as degrees of freedom, i.e. the same format as
  !> #bn_vac_file.
  character(len = 1024) :: Bn_file

  !> Formatted output data file for magnetic field perturbation (plasma response).
  !>
  !> This file is generated in every iteration by magdif::compute_bn() and read in
  !> afterwards via magdif::read_bn().
  character(len = 1024) :: Bn_diff_file

  !> Formatted output data file containing flux variables.
  !>
  !> This file is generated by magdif::magdif_init(). Its columns contain, in that order,
  !> flux \f$ \psi \f$, safety factor \$f q \f$, density, temperature and pressure.
  character(len = 1024) :: fluxvar_file

  !> output data file for calculated magdif::j0phi
  character(len = 1024) :: j0phi_file

  !> output data file for pressure perturbation
  character(len = 1024) :: presn_file

  !> output data file for current perturbation
  character(len = 1024) :: currn_file

  !> output data file for eigenvectors in magdif::magdif_precon().
  character(len = 1024) :: eigvec_file

  !> namelist for input parameters
  namelist / settings / log_level, runmode, nonres, quad_avg, niter, nritz, n, &
       nkpol, nflux, ti0, di0, damp, point_file, tri_file, Bn_vacout_file, &
       Bn_vac_file, Bn_file, Bn_diff_file, fluxvar_file, j0phi_file, presn_file, &
       currn_file, eigvec_file, tol

contains

  !> Read configuration file for magdif.
  !>
  !> @param config_file file name of config file
  subroutine read_config(config_file)
    character(len = *) :: config_file

    open(1, file = config_file, recl = longlines)
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

  !> Generates a new filename from a given template.
  !>
  !> @param in_name undecorated filename
  !> @param prefix string to be prefixed to the file basename
  !> @param iter optional postfix indicating an iteration step
  !>
  !> No attempt is made to check whether the given filename is valid or accessible.
  !> Without an \p iter and an empty prefix, \p in_name is returned.
  function decorate_filename(in_name, prefix, iter) result(out_name)
    character(len = *), intent(in) :: in_name
    character(len = *), intent(in) :: prefix
    integer, intent(in), optional :: iter
    character(len = len_trim(in_name) + len_trim(prefix) + 4) :: out_name

    integer :: basename_start, basename_end, name_end

    name_end = len_trim(in_name)
    basename_start = scan(in_name, '/', .true.) + 1
    basename_end = basename_start - 2 + scan(in_name(basename_start:), '.', .true.)
    if (basename_end - basename_start == -1) basename_start = basename_start + 1 ! dotfile
    if (basename_end < basename_start) basename_end = name_end
    if (present(iter)) then
       write (out_name, "(a, '_', i0.3, a)") in_name(:basename_start-1) // prefix // &
            in_name(basename_start:basename_end), iter, in_name(basename_end+1:name_end)
    else
       out_name = in_name(:basename_start-1) // prefix // in_name(basename_start:name_end)
    end if
  end function decorate_filename
end module magdif_config
