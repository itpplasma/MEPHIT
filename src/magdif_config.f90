module magdif_config
  use arnoldi_mod, only: tol ! arnoldi.f90
  use for_macrostep, only: t_min, d_min  ! orbit_mod.f90
  use iso_fortran_env, only: output_unit, dp => real64

  implicit none

  character(len = 1024) :: bin_dir = '.'
  character(len = 1024) :: config_file = 'magdif.inp'
  character(len = 1024) :: log_file = ''
  character(len = 1024) :: log_msg = ''

  character(len = *), parameter :: cmplx_fmt = 'es23.16, 1x, sp, es23.16, s, " i"'
  character(len = *), parameter :: nl_fmt = '"' // new_line('A') // '"'

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

  integer, private :: log = output_unit  !< log file id
  logical :: quiet = .false.  !< suppress log messages to stdout

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
  integer  :: nflux_unref     !< number of flux surfaces before refinement
  integer  :: nflux           !< number of flux surfaces
  real(dp) :: ti0             !< interpolation step for temperature
  real(dp) :: di0             !< interpolation step for density
  real(dp) :: damp            !< damping factor for resonances

  integer, dimension(:), allocatable :: deletions
  integer, dimension(:), allocatable :: additions
  real(dp), dimension(:), allocatable :: refinement

  !> Free parameters setting the magnitudes of sheet currents.
  complex(dp), dimension(:), allocatable :: sheet_current_factor

  !> Unformatted input data file containing mesh data.
  !>
  !> This file is generated in magdif_mesh_mod::generate_mesh() and is read in by
  !> magdif::read_mesh().
  character(len = 1024) :: meshdata_file

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

  !> error threshold for divergence-freeness of #bnflux and #bnphi
  real(dp) :: rel_err_Bn = 1d-7

  !> error threshold for divergence-freeness of #jnflux and #jnphi
  real(dp) :: rel_err_currn = 1d-8

  !> Integration constant for resonant vacuum perturbation in KiLCA comparison.
  complex(dp), dimension(:), allocatable :: kilca_vac_coeff

  !> Data point for resonant vacuum perturbation in KiLCA comparison.
  real(dp), dimension(:), allocatable :: kilca_vac_r

  !> Data point for resonant vacuum perturbation in KiLCA comparison.
  complex(dp), dimension(:), allocatable :: kilca_vac_Bz

  !> single poloidal mode used in comparison with KiLCA code
  integer :: kilca_pol_mode = 0

  !> Scaling factor used for comparison with results from KiLCA code.
  !>
  !> If non-zero, #n and #r0 are scaled by this factor to approximate the cylindrical
  !> topology assumed in KiLCA.
  integer :: kilca_scale_factor = 0

  !> maximum number of eigenvectors to be dumped in #eigvec_file
  integer :: max_eig_out = 10

  !> namelists for input parameters
  namelist /settings/ log_level, runmode, pres_prof, nonres, quad_avg, niter, nritz, &
       tol, n, nkpol, nflux_unref, ti0, di0, t_min, d_min, damp, meshdata_file, &
       Bn_vacout_file, Bn_vac_file, Bn_file, Bn_diff_file, fluxvar_file, j0phi_file, &
       presn_file, currn_file, eigvec_file, rel_err_Bn, rel_err_currn, kilca_pol_mode, &
       kilca_scale_factor, max_eig_out, curr_prof, q_prof, conv_file, &
       quiet
  namelist /delayed/ refinement, deletions, additions, sheet_current_factor, &
       kilca_vac_coeff, kilca_vac_r, kilca_vac_Bz

contains

  !> Read #settings configuration file #config_file for magdif.
  subroutine read_config
    integer :: fid

    open(newunit = fid, file = config_file)
    read(fid, nml = settings)
    close(fid)

    log_file = '-'  ! always write to stdout for now
    log_err = .false.
    log_warn = .false.
    log_info = .false.
    log_debug = .false.
    if (log_level > 0) log_err = .true.
    if (log_level > 1) log_warn = .true.
    if (log_level > 2) log_info = .true.
    if (log_level > 3) log_debug = .true.

    if (kilca_scale_factor /= 0) then
       n = n * kilca_scale_factor
    end if
  end subroutine read_config

  !> Read #delayed from configuration file #config_file for magdif.
  subroutine read_delayed_config(m_res_min, m_res_max)
    integer, intent(in) :: m_res_min, m_res_max
    integer :: fid
    if (allocated(refinement)) deallocate(refinement)
    allocate(refinement(m_res_min:m_res_max))
    refinement = 0d0
    if (allocated(deletions)) deallocate(deletions)
    allocate(deletions(m_res_min:m_res_max))
    deletions = 0
    if (allocated(additions)) deallocate(additions)
    allocate(additions(m_res_min:m_res_max))
    additions = 0
    if (allocated(sheet_current_factor)) deallocate(sheet_current_factor)
    allocate(sheet_current_factor(m_res_min:m_res_max))
    sheet_current_factor = (0d0, 0d0)
    if (allocated(kilca_vac_coeff)) deallocate(kilca_vac_coeff)
    allocate(kilca_vac_coeff(m_res_min:m_res_max))
    kilca_vac_coeff = (1d0, 0d0)
    if (allocated(kilca_vac_r)) deallocate(kilca_vac_r)
    allocate(kilca_vac_r(m_res_min:m_res_max))
    kilca_vac_r = 0d0
    if (allocated(kilca_vac_Bz)) deallocate(kilca_vac_Bz)
    allocate(kilca_vac_Bz(m_res_min:m_res_max))
    kilca_vac_Bz = (0d0, 0d0)
    open(newunit = fid, file = config_file)
    read(fid, nml = delayed)
    close(fid)
  end subroutine read_delayed_config

  !> Associate logfile and open if necessary.
  subroutine log_open
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
    if (log /= output_unit) close(log)
  end subroutine log_close

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
  subroutine log_write
    character(len = 28 + len_trim(log_msg)) :: full_msg
    write (full_msg, '("[", a25, "] ", a)') timestamp(), trim(log_msg)
    write (log, '(a)') full_msg
    if (log /= output_unit .and. .not. quiet) write (output_unit, '(a)') full_msg
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
