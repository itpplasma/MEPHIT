module magdif_config
  use from_nrtype, only: dp  ! PRELOAD/SRC/from_nrtype.f90
  use arnoldi_mod, only: tol ! RUN/SRC/arnoldi.f90
  use for_macrostep, only: t_min, d_min  ! PRELOAD/SRC/orbit_mod.f90

  implicit none

  character(len = 1024) :: config_file
  character(len = 1024) :: log_file = ''
  character(len = 1024) :: log_msg = ''

  character(len = *), parameter :: cmplx_fmt = 'es23.16, 1x, sp, es23.16, " i"'
  character(len = *), parameter :: nl_fmt = '"' // new_line('A') // '"'

  integer, parameter :: longlines = 1024

  integer, parameter :: runmode_single = 0 !< single iteration mode
  integer, parameter :: runmode_direct = 1 !< direct iteration mode
  integer, parameter :: runmode_precon = 2 !< preconditioned iteration mode

  integer, parameter :: pres_prof_eps = 0  !< pressure profile from EPS paper
  integer, parameter :: pres_prof_par = 1  !< parabolic pressure profile
  integer, parameter :: pres_prof_efit = 2 !< pressure profile from EFIT file

  integer, parameter :: curr_prof_ps = 0   !< only Pfirsch-Schlueter current
  integer, parameter :: curr_prof_rot = 1  !< current via Ampere's law
  integer, parameter :: curr_prof_efit = 2 !< current profile from EFIT file

  integer, parameter :: q_prof_flux = 0 !< q profile from flux between flux surfaces
  integer, parameter :: q_prof_efit = 2 !< q profile from EFIT file

  integer :: log = 6  !< log to stdout

  integer :: log_level
  integer :: runmode
  integer :: pres_prof = pres_prof_par  !< pressure profile
  integer :: curr_prof = curr_prof_ps   !< current profile
  integer :: q_prof = q_prof_flux       !< q profile
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

  !> Distance of magnetic axis from center \f$ R_{0} \f$ in cm.
  real(dp) :: R0 = 172.74467899999999d0

  !> Free parameters setting the magnitudes of sheet currents.
  complex(dp), dimension(:), allocatable :: sheet_current_factor

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

  !> output data file for eigenvectors in magdif::magdif_iterated().
  character(len = 1024) :: eigvec_file

  !> output data file for convergence estimation in magdif::magdif_iterated().
  character(len = 1024) :: conv_file = 'convergence.dat'

  !> output data file for poloidal_modes in magdif::magdif_iterated().
  character(len = 1024) :: kilca_pol_mode_file

  !> error threshold for divergence-freeness of #bnflux and #bnphi
  real(dp) :: rel_err_Bn = 1d-7

  !> error threshold for divergence-freeness of #jnflux and #jnphi
  real(dp) :: rel_err_currn = 1d-8

  !> Integration constant for resonant vacuum perturbation in KiLCA comparison.
  complex(dp) :: kilca_vac_coeff

  !> single poloidal mode used in comparison with KiLCA code
  integer :: kilca_pol_mode

  !> Scaling factor used for comparison with results from KiLCA code.
  !>
  !> If non-zero, #n and #r0 are scaled by this factor to approximate the cylindrical
  !> topology assumed in KiLCA.
  integer :: kilca_scale_factor

  !> maximum number of eigenvectors to be dumped in #eigvec_file
  integer :: max_eig_out = 10

  !> namelists for input parameters
  namelist /settings/ log_level, runmode, pres_prof, nonres, quad_avg, niter, nritz, &
       tol, n, nkpol, nflux, ti0, di0, t_min, d_min, damp, R0, point_file, tri_file, &
       Bn_vacout_file, Bn_vac_file, Bn_file, Bn_diff_file, fluxvar_file, j0phi_file, &
       presn_file, currn_file, eigvec_file, rel_err_Bn, rel_err_currn, kilca_pol_mode, &
       kilca_vac_coeff, kilca_scale_factor, kilca_pol_mode_file, max_eig_out, curr_prof, &
       q_prof, conv_file
  namelist /delayed/ sheet_current_factor

contains

  !> Read #settings configuration file #config_file for magdif.
  subroutine read_config
    open(1, file = config_file)
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

  !> Read #delayed from configuration file #config_file for magdif.
  subroutine read_delayed_config
    open(1, file = config_file)
    read(1, nml = delayed)
    close(1)
  end subroutine read_delayed_config

  !> Associate logfile and open if necessary.
  subroutine log_open
    use iso_fortran_env, only: output_unit
    if (trim(log_file) == '-') then
       log = output_unit
    elseif (trim(log_file) == '') then
       open(newunit = log, file = trim(config_file) // '.log', status = 'replace')
    else
       open(newunit = log, file = log_file, status = 'replace')
    end if
  end subroutine log_open

  !> Close logfile if necessary.
  subroutine log_close
    use iso_fortran_env, only: output_unit
    if (log /= output_unit) close(log)
  end subroutine log_close

  !> Write to logfile and flush if necessary.
  subroutine log_write
    use iso_fortran_env, only: output_unit
    write (log, '(a)') trim(log_msg)
    if (log /= output_unit) flush(log) ! or: write (output_unit, '(a)') trim(log_msg)
  end subroutine log_write

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
end module magdif_config
