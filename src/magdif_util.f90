module magdif_util

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: clight, imun, initialize_globals, get_equil_filenames, interp_psi_pol, &
       ring_centered_avg_coord, assemble_sparse, linspace, straight_cyl2bent_cyl, &
       bent_cyl2straight_cyl, binsearch, interleave, &
       gauss_legendre_unit_interval, heapsort_complex, complex_abs_asc

  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

  type, public :: sign_convention
    integer :: exp_Bpol, sgn_cyl, sgn_dpsi, sgn_Btor, sgn_Itor, &
         sgn_F, sgn_q, sgn_Bpol, sgn_pol, index
  end type sign_convention

  type, public :: g_eqdsk
    type(sign_convention) :: cocos
    character(len = 1024) :: fname
    character(len = 10) :: text(6)
    integer :: nw, nh, nbbbs, limitr
    real(dp) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, &
         current
    real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime, qpsi, &
         rbbbs, zbbbs, rlim, zlim, psi_eqd, R_eqd, Z_eqd
    real(dp), dimension(:, :), allocatable :: psirz
  contains
    procedure :: read => g_eqdsk_read
    procedure :: check_consistency => g_eqdsk_check_consistency
    procedure :: classify => g_eqdsk_classify
    procedure :: standardise => g_eqdsk_standardise
    procedure :: scale => g_eqdsk_scale
    procedure :: write => g_eqdsk_write
    procedure :: grad_shafranov_normalization => g_eqdsk_grad_shafranov_normalization
    final :: g_eqdsk_destructor
  end type g_eqdsk

  character(len = *), parameter :: geqdsk_2000 = '(6a8,3i4)'
  character(len = *), parameter :: geqdsk_2020 = '(5e16.9)'
  character(len = *), parameter :: geqdsk_2022 = '(2i5)'

  character(len = *), parameter :: &
       incons_fmt = '("Signs of ", a, " and ", a, " are inconsistent.")', &
       invert_fmt = '("Inverting sign of ", a, "...")', &
       unscaled_fmt = '("Flux in ", a, " is not normalized by 2 pi.")', &
       rescale_fmt = '("Rescaling ", a, "...")'

  type, public :: flux_func
    private
    integer :: n_lag, n_var
    real(dp), dimension(:), allocatable :: indep_var
  contains
    procedure :: init => flux_func_init
    procedure :: interp => flux_func_interp
    final :: flux_func_destructor
  end type flux_func

  !> Structure containing flux functions evaluated at a specific flux surface, indicated
  !> by a common array index. For details see flux_func_cache_init().
  type, public :: flux_func_cache
     private
     integer :: nflux

     !> Magnetic flux surface label \f$ \psi \f$ in maxwell.
     !>
     !> \f$ \psi \f$ is the disc poloidal flux divided by \f$ 2 \pi \f$. Its sign is
     !> positive and its magnitude is growing in the radially inward direction.
     real(dp), dimension(:), allocatable, public :: psi

     !> Minor radius \f$ r \f$ in centimeter.
     real(dp), dimension(:), allocatable, public :: rad

     real(dp), dimension(:), allocatable, public :: F

     !> Unperturbed pressure \f$ p_{0} \f$ in barye.
     real(dp), dimension(:), allocatable, public :: p

     real(dp), dimension(:), allocatable, public :: FdF_dpsi

     !> Derivative of unperturbed pressure w.r.t. flux surface label,
     !> \f$ p_{0}'(\psi) \f$, in barye per maxwell.
     real(dp), dimension(:), allocatable, public :: dp_dpsi

     !> Safety factor \f$ q \f$ (dimensionless).
     real(dp), dimension(:), allocatable, public :: q
   contains
     procedure :: init => flux_func_cache_init
     final :: flux_func_cache_destructor
  end type flux_func_cache

  type, public :: neumaier_accumulator_real
     private
     real(dp) :: sum, c, t
   contains
     procedure :: init => neumaier_accumulator_real_init
     procedure :: add => neumaier_accumulator_real_add
     procedure :: get_sum => neumaier_accumulator_real_get_sum
  end type neumaier_accumulator_real

  type, public :: neumaier_accumulator_complex
     private
     type(neumaier_accumulator_real) :: real_part, imag_part
   contains
     procedure :: init => neumaier_accumulator_complex_init
     procedure :: add => neumaier_accumulator_complex_add
     procedure :: get_sum => neumaier_accumulator_complex_get_sum
  end type neumaier_accumulator_complex

  interface interleave
     procedure interleave_vv, interleave_vs, interleave_sv, interleave_ss
  end interface interleave

contains

  ! better future solution: put this in a separate subroutine in field_divB0.f90
  subroutine get_equil_filenames(gfile, convexfile)
    character(len = *), intent(out) :: gfile, convexfile
    integer :: fid
    open(newunit = fid, file = 'field_divB0.inp', status = 'old')
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *) gfile        ! equilibrium file
    read (fid, *)
    read (fid, *) convexfile   ! convex file for stretchcoords
    close(fid)
  end subroutine get_equil_filenames

  subroutine initialize_globals(r, z)
    real(dp), intent(in) :: r, z
    real(dp) :: dum

    call field(r, 0d0, z, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
  end subroutine initialize_globals

  function interp_psi_pol(r, z) result(psi_pol)
    use field_eq_mod, only: psif, psib

    real(dp), intent(in) :: r, z
    real(dp) :: psi_pol

    real(dp) :: dum

    call field(r, 0d0, z, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
    ! field_divB0.f90 adds psib (SIBRY, i.e., flux at the boundary) to interpolated psi[f]
    psi_pol = psif - psib
  end function interp_psi_pol

  !> Computes the "weighted" centroid for a triangle so that it is approximately
  !> equidistant between the enclosing flux surfaces, independent of triangle orientation.
  !>
  !> @param elem the triangle for which the centroid is to be computed
  !> @param r radial cylindrical coordinate of the centroid
  !> @param z axial cylindrical coordinate of the centroid
  !>
  !> Depending on the orientation of the triangle (see also \p orient of
  !> get_labeled_edges()), two knots lie on the inner flux surface and one on the outer
  !> one, or vice versa. A simple arithmetic mean of the three knots' coordinates would
  !> place the centroid closer to the inner flux surface for one orientation and closer
  !> to the outer one for the other. To counteract this, the "lonely" knot is counted
  !> twice in the averaging procedure, i.e. with double weighting.
  pure subroutine ring_centered_avg_coord(elem, r, z)
    use mesh_mod, only: triangle, knot, mesh_point
    type(triangle), intent(in) :: elem
    real(dp), intent(out) :: r, z
    type(knot), dimension(3) :: knots

    knots = mesh_point(elem%i_knot)
    r = (sum(knots%rcoord) + knots(elem%knot_h)%rcoord) * 0.25d0
    z = (sum(knots%zcoord) + knots(elem%knot_h)%zcoord) * 0.25d0
  end subroutine ring_centered_avg_coord

  !> Assembles a sparse matrix in coordinate list (COO) representation for use with
  !> sparse_mod::sparse_solve().
  !>
  !> @param nrow number \f$ n \f$ of rows in the matrix
  !> @param d \f$ n \f$ diagnonal elements of stiffness matrix \f$ A \f$
  !> @param du \f$ n-1 \f$ superdiagonal elements of stiffness matrix \f$ A \f$ and
  !> \f$ A_{n, 1} \f$ (lower left corner)
  !> @param nz number of non-zero entries (2*nrow)
  !> @param irow nz row indices of non-zero entries
  !> @param icol nz column indices of non-zero entries
  !> @param aval nz values of non-zero entries
  !>
  !> The input is a stiffness matrix \f$ K \f$ with non-zero entries on the main diagonal,
  !> the upper diagonal and, due to periodicity, in the lower left corner. This shape
  !> results from the problems in compute_presn() and compute_currn().
  subroutine assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    use magdif_conf, only: log
    integer, intent(in)  :: nrow
    complex(dp), intent(in)  :: d(1:)  !nrow
    complex(dp), intent(in)  :: du(1:)  !nrow
    integer, intent(out) :: nz
    integer, intent(out) :: irow(1:), icol(1:)  !2*nrow
    complex(dp), intent(out) :: aval(1:)  !2*nrow

    integer :: k

    nz = 2*nrow

    if (nrow /= size(d)) then
       call log%msg_arg_size('assemble_sparse', 'nrow', 'size(d)', nrow, size(d))
       if (log%err) call log%write
       error stop
    end if
    if (nrow /= size(du)) then
       call log%msg_arg_size('assemble_sparse', 'nrow', 'size(du)', nrow, size(du))
       if (log%err) call log%write
       error stop
    end if
    if (nz /= size(irow)) then
       call log%msg_arg_size('assemble_sparse', 'nz', 'size(irow)', nz, size(irow))
       if (log%err) call log%write
       error stop
    end if
    if (nz /= size(icol)) then
       call log%msg_arg_size('assemble_sparse', 'nz', 'size(icol)', nz, size(icol))
       if (log%err) call log%write
       error stop
    end if
    if (nz /= size(aval)) then
       call log%msg_arg_size('assemble_sparse', 'nz', 'size(aval)', nz, size(aval))
       if (log%err) call log%write
       error stop
    end if

    irow(1) = 1
    icol(1) = 1
    aval(1) = d(1)

    irow(2) = nrow
    icol(2) = 1
    aval(2) = du(nrow)

    do k = 2, nrow
       ! off-diagonal
       irow(2*k-1) = k-1
       icol(2*k-1) = k
       aval(2*k-1) = du(k-1)

       ! diagonal
       irow(2*k) = k
       icol(2*k) = k
       aval(2*k) = d(k)
    end do
  end subroutine assemble_sparse

  function linspace(lo, hi, cnt, excl_lo, excl_hi)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: cnt, excl_lo, excl_hi
    real(dp) :: linspace(cnt)
    real(dp) :: step
    integer :: k

    step = (hi - lo) / dble(cnt - 1 + excl_lo + excl_hi)
    linspace = lo + [(k * step, k = excl_lo, cnt - 1 + excl_lo)]
  end function linspace


  !> Binary search to find index \p k of ordered array \p x so that \p xi lies between
  !> x(k-1) and x(k). The lower bound of \p x is given by \p lb.
  !>
  !> When \p x is strictly monotonically increasing, \f$ x_{k-1} < \xi < x_{k} \f$.
  !> When \p x is strictly monotonically decreasing, \f$ x_{k-1} > \xi > x_{k} \f$.
  !> It is not checked whether \p x is ordered; only the first and last value are used
  !> to determine wether the array values are increasing or decrasing.
  subroutine binsearch(x, lb, xi, k)
    real(dp), intent(in) :: x(lb:)
    integer, intent(in) :: lb
    real(dp), intent(in) :: xi
    integer, intent(out) :: k
    integer :: k_min, k_max, iter

    k_min = lbound(x, 1)
    k_max = ubound(x, 1)
    if (x(k_min) < x(k_max)) then
       do iter = 1, size(x) - 1
          k = (k_max - k_min) / 2 + k_min
          if (x(k) > xi) then
             k_max = k
          else
             k_min = k
          end if
          if (k_max == k_min + 1) exit
       end do
    else
       do iter = 1, size(x) - 1
          k = (k_max - k_min) / 2 + k_min
          if (x(k) < xi) then
             k_max = k
          else
             k_min = k
          end if
          if (k_max == k_min + 1) exit
       end do
    end if
    k = k_max
  end subroutine binsearch

  subroutine gauss_legendre_unit_interval(order, points, weights)
    use magdif_conf, only: log
    use fgsl, only: fgsl_size_t, fgsl_double, fgsl_int, fgsl_success, &
         fgsl_integration_glfixed_point, fgsl_integration_glfixed_table, &
         fgsl_integration_glfixed_table_alloc, fgsl_integration_glfixed_table_free
    integer, intent(in) :: order
    real(fgsl_double), dimension(:), intent(out) :: points, weights
    type(fgsl_integration_glfixed_table) :: table
    integer(fgsl_size_t) :: k
    integer(fgsl_int) :: err
    if (order /= size(points)) then
       call log%msg_arg_size('gauss_legendre_unit_interval', 'order', 'size(points)', &
            order, size(points))
       if (log%err) call log%write
       error stop
    end if
    if (order /= size(weights)) then
       call log%msg_arg_size('gauss_legendre_unit_interval', 'order', 'size(weights)', &
            order, size(weights))
       if (log%err) call log%write
       error stop
    end if
    table = fgsl_integration_glfixed_table_alloc(int(order, fgsl_size_t))
    do k = 1, int(order, fgsl_size_t)
       err = fgsl_integration_glfixed_point(0d0, 1d0, k-1, points(k), weights(k), table)
       if (err /= fgsl_success) then
          write (log%msg, '("fgsl_integration_glfixed_point returned error ", i0)') err
          if (log%err) call log%write
          error stop
       end if
    end do
    call fgsl_integration_glfixed_table_free(table)
  end subroutine gauss_legendre_unit_interval


  !> Transform components of a vector \f$ \vec{v} \f$ from straight cylinder coordinates
  !> \f$ (r, \theta, z) \f$ to bent cylinder coordinates \f$ (R, \varphi, Z) \f$.
  !>
  !> @param comp_rad physical component \f$ v_{r} \f$
  !> @param comp_pol physical component \f$ v_{(\theta)} \f$
  !> @param comp_tor physical component \f$ v_{z} \f$
  !> @param theta geometrical poloidal angle \f$ theta \f$ (coinciding with symmetry flux
  !> coordinates' poloidal angle in this geometry)
  !> @param comp_R physical component \f$ v_{R} \f$
  !> @param comp_phi physical component \f$ v_{(\varphi)} \f$
  !> @param comp_Z physical component \f$ v_{Z} \f$
  subroutine straight_cyl2bent_cyl(comp_rad, comp_pol, comp_tor, theta, &
       comp_R, comp_phi, comp_Z)
    complex(dp), intent(in) :: comp_rad, comp_pol, comp_tor
    real(dp), intent(in) :: theta
    complex(dp), intent(out) :: comp_R, comp_phi, comp_Z

    comp_R = comp_rad * cos(theta) - comp_pol * sin(theta)
    comp_phi = comp_tor ! * (1d0 + r / R_0 * cos(theta))  ! exact version
    comp_Z = comp_rad * sin(theta) + comp_pol * cos(theta)
  end subroutine straight_cyl2bent_cyl

  !> Transform components of a vector \f$ \vec{v} \f$ from bent cylinder coordinates
  !> \f$ (R, \varphi, Z) \f$ to straight cylinder coordinates \f$ (r, \theta, z) \f$.
  !>
  !> @param comp_R physical component \f$ v_{R} \f$
  !> @param comp_phi physical component \f$ v_{(\varphi)} \f$
  !> @param comp_Z physical component \f$ v_{Z} \f$
  !> @param theta geometrical poloidal angle \f$ theta \f$ (coinciding with symmetry flux
  !> coordinates' poloidal angle in this geometry)
  !> @param comp_rad physical component \f$ v_{r} \f$
  !> @param comp_pol physical component \f$ v_{(\theta)} \f$
  !> @param comp_tor physical component \f$ v_{z} \f$
  subroutine bent_cyl2straight_cyl(comp_R, comp_phi, comp_Z, theta, &
       comp_rad, comp_pol, comp_tor)
    complex(dp), intent(in) :: comp_R, comp_phi, comp_Z
    real(dp), intent(in) :: theta
    complex(dp), intent(out) :: comp_rad, comp_pol, comp_tor

    comp_rad = comp_Z * sin(theta) + comp_R * cos(theta)
    comp_pol = comp_Z * cos(theta) - comp_R * sin(theta)
    comp_tor = comp_phi ! / (1d0 + r / R_0 * cos(theta))  ! exact version
  end subroutine bent_cyl2straight_cyl

  subroutine g_eqdsk_read(this, fname)
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    call g_eqdsk_destructor(this)
    this%fname = fname
    open(newunit = fid, file = this%fname, status = 'old')
    read (fid, geqdsk_2000) (this%text(kw), kw = 1, 6), idum, this%nw, this%nh
    allocate(this%fpol(this%nw))
    allocate(this%pres(this%nw))
    allocate(this%ffprim(this%nw))
    allocate(this%pprime(this%nw))
    allocate(this%qpsi(this%nw))
    allocate(this%psirz(this%nw, this%nh))
    read (fid, geqdsk_2020) this%rdim, this%zdim, this%rcentr, this%rleft, this%zmid
    read (fid, geqdsk_2020) this%rmaxis, this%zmaxis, this%simag, this%sibry, this%bcentr
    read (fid, geqdsk_2020) this%current, (xdum, kw = 1, 4)
    read (fid, geqdsk_2020) (xdum, kw = 1, 5)
    read (fid, geqdsk_2020) (this%fpol(kw), kw = 1, this%nw)
    read (fid, geqdsk_2020) (this%pres(kw), kw = 1, this%nw)
    read (fid, geqdsk_2020) (this%ffprim(kw), kw = 1, this%nw)
    read (fid, geqdsk_2020) (this%pprime(kw), kw = 1, this%nw)
    read (fid, geqdsk_2020) ((this%psirz(kw, kh), kw = 1, this%nw), kh = 1, this%nh)
    read (fid, geqdsk_2020) (this%qpsi(kw), kw = 1, this%nw)
    read (fid, geqdsk_2022) this%nbbbs, this%limitr
    allocate(this%rbbbs(this%nbbbs))
    allocate(this%zbbbs(this%nbbbs))
    allocate(this%rlim(this%limitr))
    allocate(this%zlim(this%limitr))
    read (fid, geqdsk_2020) (this%rbbbs(kw), this%zbbbs(kw), kw = 1, this%nbbbs)
    read (fid, geqdsk_2020) (this%rlim(kw), this%zlim(kw), kw = 1, this%limitr)
    close(fid)
    ! convert meter to centimeter
    this%rdim = this%rdim * 1d2
    this%zdim = this%zdim * 1d2
    this%rcentr = this%rcentr * 1d2
    this%rleft = this%rleft * 1d2
    this%zmid = this%zmid * 1d2
    this%rmaxis = this%rmaxis * 1d2
    this%zmaxis = this%zmaxis * 1d2
    ! convert weber to maxwell
    this%simag = this%simag * 1d8
    this%sibry = this%sibry * 1d8
    ! convert tesla to gauss
    this%bcentr = this%bcentr * 1d4
    ! convert ampere to statampere
    this%current = this%current * 1d-1 * clight
    ! convert tesla meter to gauss centimeter
    this%fpol = this%fpol * 1d6
    ! convert pascal to barye
    this%pres = this%pres * 1d1
    ! convert tesla to gauss
    this%ffprim = this%ffprim * 1d4
    ! convert pascal per weber to barye per maxwell
    this%pprime = this%pprime * 1d-7
    ! convert weber to maxwell
    this%psirz = this%psirz * 1d8
    ! convert meter to centimeter
    this%rbbbs = this%rbbbs * 1d2
    this%zbbbs = this%zbbbs * 1d2
    this%rlim = this%rlim * 1d2
    this%zlim = this%zlim * 1d2
    ! initialize equidistant values
    allocate(this%psi_eqd(this%nw))
    this%psi_eqd(:) = linspace(this%simag, this%sibry, this%nw, 0, 0)
    allocate(this%R_eqd(this%nw))
    this%R_eqd(:) = this%rleft + linspace(0d0, this%rdim, this%nw, 0, 0)
    allocate(this%Z_eqd(this%nh))
    this%Z_eqd(:) = this%zmid + 0.5d0 * linspace(-this%zdim, this%zdim, this%nh, 0, 0)
  end subroutine g_eqdsk_read

  subroutine g_eqdsk_check_consistency(this)
    use constants, only: pi  ! src/orbit_mod.f90
    use magdif_conf, only: log
    class(g_eqdsk), intent(inout) :: this
    integer, parameter :: ignore = 3
    type(flux_func) :: psi_interpolator
    real(dp) :: deriv_eqd(this%nw), factor
    logical :: mask(this%nw)
    integer :: k

    if (this%cocos%sgn_Btor /= this%cocos%sgn_F) then
       write (log%msg, incons_fmt) 'FPOL', 'BCENTR'
       if (log%warn) call log%write
       ! only modify array if signs are consistent
       if (this%cocos%sgn_F /= 0) then
          write (log%msg, invert_fmt) 'FPOL'
          if (log%info) call log%write
          this%fpol = -this%fpol
          this%cocos%sgn_F = -this%cocos%sgn_F
       end if
    end if
    call psi_interpolator%init(4, this%psi_eqd)
    deriv_eqd(:) = [(psi_interpolator%interp(this%pres, this%psi_eqd(k), .true.), &
         k = 1, this%nw)]
    ! ignore a few possibly unreliable values on both ends of the interval
    mask = .true.
    mask(1:1+ignore) = .false.
    mask(this%nw-ignore:this%nw) = .false.
    factor = sum(deriv_eqd / this%pprime, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2d0 * pi)) then
       write (log%msg, unscaled_fmt) 'PPRIME'
       if (log%warn) call log%write
       write (log%msg, rescale_fmt) 'PPRIME'
       if (log%info) call log%write
       this%pprime = this%pprime * (2d0 * pi)
    end if
    if (factor < 0d0) then
       write (log%msg, incons_fmt) 'PPRIME', 'PRES/SIBRY-SIMAG'
       if (log%warn) call log%write
       write (log%msg, invert_fmt) 'PPRIME'
       if (log%info) call log%write
       this%pprime = -this%pprime
    end if
    deriv_eqd(:) = [(psi_interpolator%interp(this%fpol, this%psi_eqd(k), .true.), &
         k = 1, this%nw)] * this%fpol
    factor = sum(deriv_eqd / this%ffprim, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2d0 * pi)) then
       write (log%msg, unscaled_fmt) 'FFPRIM'
       if (log%warn) call log%write
       write (log%msg, rescale_fmt) 'FFPRIM'
       if (log%info) call log%write
       this%ffprim = this%ffprim * (2d0 * pi)
    end if
    if (factor < 0d0) then
       write (log%msg, incons_fmt) 'FFPRIM', 'FPOL/SIBRY-SIMAG'
       if (log%warn) call log%write
       write (log%msg, invert_fmt) 'FFPRIM'
       if (log%info) call log%write
       this%ffprim = -this%ffprim
    end if
  end subroutine g_eqdsk_check_consistency

  !> Estimates terms of Grad-Shafranov equation to determine sign_convention::exp_bpol.
  function g_eqdsk_grad_shafranov_normalization(this) result(gs_factor)
    use constants, only: pi  ! src/orbit_mod.f90
    use magdif_conf, only: log
    class(g_eqdsk), intent(inout) :: this
    real(dp) :: gs_factor
    type(flux_func) :: psi_interpolator
    real(dp) :: Delta_R, Delta_Z, psi, gs_rhs, gs_lhs
    integer :: kw, kh

    ! find evaluation point to the upper right of the magnetic axis
    call binsearch(this%R_eqd, lbound(this%R_eqd, 1), this%rmaxis, kw)
    call binsearch(this%Z_eqd, lbound(this%Z_eqd, 1), this%zmaxis, kh)
    kw = kw + 2
    kh = kh + 2
    psi = this%psirz(kw, kh)
    ! evaluate flux functions on RHS of GS equation
    call psi_interpolator%init(4, this%psi_eqd)
    gs_rhs = -4d0 * pi * this%R_eqd(kw) ** 2 * psi_interpolator%interp(this%pprime, psi) &
         - psi_interpolator%interp(this%ffprim, psi)
    ! approximate differential operator on LHS of GS equation
    Delta_R = this%rdim / dble(this%nw - 1)
    Delta_Z = this%zdim / dble(this%nh - 1)
    gs_lhs = (this%psirz(kw, kh + 1) - 2d0 * psi + this%psirz(kw, kh - 1)) / Delta_Z ** 2 &
         + (this%psirz(kw + 1, kh) - 2d0 * psi + this%psirz(kw - 1, kh)) / Delta_R ** 2 &
         - (this%psirz(kw + 1, kh) - this%psirz(kw - 1, kh)) / (2d0 * Delta_R * this%R_eqd(kw))
    gs_factor = gs_lhs / gs_rhs
    write (log%msg, '("Grad-Shafranov equation LHS / RHS: ", f19.16)') gs_factor
    if (log%info) call log%write
  end function g_eqdsk_grad_shafranov_normalization

  function sign_array(array, name, most)
    use magdif_conf, only: log
    real(dp), intent(in), dimension(:) :: array
    character(len = *), intent(in) :: name
    logical, intent(in), optional :: most
    logical :: strict
    integer :: sign_array, pos, neg, all

    strict = .true.
    if (present(most)) then
       strict = .not. most
    end if
    pos = count(array > 0d0)
    neg = count(array < 0d0)
    all = size(array)
    if ((strict .and. all == pos) .or. (.not. strict .and. pos > neg)) then
       sign_array = +1
    elseif ((strict .and. all == neg) .or. (.not. strict .and. neg > pos)) then
       sign_array = -1
    else
       sign_array = 0
       write (log%msg, '("Sign of ", a, " is inconsistent.")') trim(name)
       if (log%warn) call log%write
    end if
  end function sign_array

  subroutine g_eqdsk_classify(this)
    use constants, only: pi  ! src/orbit_mod.f90
    use magdif_conf, only: log
    class(g_eqdsk), intent(inout) :: this

    this%cocos = sign_convention(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    this%cocos%sgn_cyl = +1  ! specified by G EQDSK format and assumed in field_divB0.f90
    this%cocos%sgn_dpsi = sign_array([this%sibry - this%simag], 'SIBRY-SIMAG')
    this%cocos%sgn_Btor = sign_array([this%bcentr], 'BCENTR')
    this%cocos%sgn_Itor = sign_array([this%current], 'CURRENT')
    this%cocos%sgn_F = sign_array(this%fpol, 'FPOL')
    ! q is often unreliable in gfiles and its sign is only necessary for exact
    ! classification, not for calculations where the poloidal direction is fixed anyway
    this%cocos%sgn_q = sign_array(this%qpsi, 'QPSI', .true.)
    this%cocos%sgn_Bpol = this%cocos%sgn_dpsi * this%cocos%sgn_Itor
    this%cocos%sgn_pol = this%cocos%sgn_q * this%cocos%sgn_Itor * this%cocos%sgn_Btor
    if (this%cocos%sgn_Bpol == +1) then
       if (this%cocos%sgn_pol == +1) then
          this%cocos%index = 1
       elseif (this%cocos%sgn_pol == -1) then
          this%cocos%index = 5
       end if
    elseif (this%cocos%sgn_Bpol == -1) then
       if (this%cocos%sgn_pol == +1) then
          this%cocos%index = 7
       elseif (this%cocos%sgn_pol == -1) then
          this%cocos%index = 3
       end if
    end if
    if (this%cocos%index == 0) then
       log%msg = 'COCOS index could not be determined for ' // trim(this%fname)
       if (log%warn) call log%write
    end if
    call this%check_consistency
    if (abs(this%grad_shafranov_normalization()) >= 2d0 * pi) then
       this%cocos%index = this%cocos%index + 10
       this%cocos%exp_Bpol = 1
    else
       ! specified by G EQDSK format and assumed in field_divB0.f90
       this%cocos%exp_Bpol = 0
    end if
  end subroutine g_eqdsk_classify

  subroutine g_eqdsk_standardise(this)
    use constants, only: pi  ! src/orbit_mod.f90
    use magdif_conf, only: log
    class(g_eqdsk), intent(inout) :: this

    if (this%cocos%sgn_Bpol == +1) then
       write (log%msg, invert_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (log%info) call log%write
       this%simag = -this%simag
       this%sibry = -this%sibry
       this%psirz = -this%psirz
       this%psi_eqd = -this%psi_eqd
       this%cocos%sgn_dpsi = -this%cocos%sgn_dpsi
       this%pprime = -this%pprime
       this%ffprim = -this%ffprim
       this%cocos%sgn_Bpol = -this%cocos%sgn_Bpol
    end if
    if (this%cocos%sgn_pol == +1) then
       write (log%msg, invert_fmt) 'QPSI'
       if (log%info) call log%write
       this%qpsi = -this%qpsi
       this%cocos%sgn_q = -this%cocos%sgn_q
       this%cocos%sgn_pol = -this%cocos%sgn_pol
    end if
    if (this%cocos%exp_Bpol == 1) then
       write (log%msg, unscaled_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (log%warn) call log%write
       write (log%msg, rescale_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (log%info) call log%write
       this%simag = this%simag / (2d0 * pi)
       this%sibry = this%sibry / (2d0 * pi)
       this%psirz = this%psirz / (2d0 * pi)
       this%psi_eqd = this%psi_eqd / (2d0 * pi)
       this%pprime = this%pprime * (2d0 * pi)
       this%ffprim = this%ffprim * (2d0 * pi)
       this%cocos%exp_Bpol = 0
    end if
    this%cocos%index = 3
  end subroutine g_eqdsk_standardise

  subroutine g_eqdsk_scale(this, gamma)
    class(g_eqdsk), intent(inout) :: this
    integer, intent(in) :: gamma
    real(dp) :: r_shift

    r_shift = this%rcentr * (gamma - 1)
    this%rcentr = this%rcentr * gamma
    this%rleft = this%rleft + r_shift
    this%rmaxis = this%rmaxis + r_shift
    this%simag = this%simag * gamma
    this%sibry = this%sibry * gamma
    this%fpol = this%fpol * gamma
    this%ffprim = this%ffprim * gamma
    this%pprime = this%pprime / gamma
    this%psirz = this%psirz * gamma
    this%qpsi = this%qpsi / gamma
    this%rbbbs = this%rbbbs + r_shift
    this%rlim = this%rlim + r_shift
    ! auxiliary values
    this%psi_eqd = this%psi_eqd * gamma
    this%R_eqd = this%R_eqd + r_shift
  end subroutine g_eqdsk_scale

  subroutine g_eqdsk_write(this, fname)
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    idum = 0
    xdum = 0d0
    this%fname = fname
    open(newunit = fid, file = this%fname, status = 'replace')
    write (fid, geqdsk_2000) (this%text(kw), kw = 1, 6), idum, this%nw, this%nh
    write (fid, geqdsk_2020) this%rdim * 1d-2, this%zdim * 1d-2, this%rcentr * 1d-2, &
         this%rleft * 1d-2, this%zmid * 1d-2
    write (fid, geqdsk_2020) this%rmaxis * 1d-2, this%zmaxis * 1d-2, this%simag * 1d-8, &
         this%sibry * 1d-8, this%bcentr * 1d-4
    write (fid, geqdsk_2020) this%current * 1d1 / clight, this%simag * 1d-8, xdum, &
         this%rmaxis * 1d-2, xdum
    write (fid, geqdsk_2020) this%zmaxis * 1d-2, xdum, this%sibry * 1d-8, xdum, xdum
    write (fid, geqdsk_2020) (this%fpol(kw) * 1d-6, kw = 1, this%nw)
    write (fid, geqdsk_2020) (this%pres(kw) * 1d-1, kw = 1, this%nw)
    write (fid, geqdsk_2020) (this%ffprim(kw) * 1d-4, kw = 1, this%nw)
    write (fid, geqdsk_2020) (this%pprime(kw) * 1d7, kw = 1, this%nw)
    write (fid, geqdsk_2020) ((this%psirz(kw, kh) * 1d-8, kw = 1, this%nw), &
         kh = 1, this%nh)
    write (fid, geqdsk_2020) (this%qpsi(kw), kw = 1, this%nw)
    write (fid, geqdsk_2022) this%nbbbs, this%limitr
    write (fid, geqdsk_2020) (this%rbbbs(kw) * 1d-2, this%zbbbs(kw) * 1d-2, &
         kw = 1, this%nbbbs)
    write (fid, geqdsk_2020) (this%rlim(kw) * 1d-2, this%zlim(kw) * 1d-2, &
         kw = 1, this%limitr)
    close(fid)
  end subroutine g_eqdsk_write

  subroutine g_eqdsk_destructor(this)
    type(g_eqdsk) :: this

    if (allocated(this%fpol)) deallocate(this%fpol)
    if (allocated(this%pres)) deallocate(this%pres)
    if (allocated(this%ffprim)) deallocate(this%ffprim)
    if (allocated(this%pprime)) deallocate(this%pprime)
    if (allocated(this%qpsi)) deallocate(this%qpsi)
    if (allocated(this%rbbbs)) deallocate(this%rbbbs)
    if (allocated(this%zbbbs)) deallocate(this%zbbbs)
    if (allocated(this%rlim)) deallocate(this%rlim)
    if (allocated(this%zlim)) deallocate(this%zlim)
    if (allocated(this%psi_eqd)) deallocate(this%psi_eqd)
    if (allocated(this%R_eqd)) deallocate(this%R_eqd)
    if (allocated(this%Z_eqd)) deallocate(this%Z_eqd)
    if (allocated(this%psirz)) deallocate(this%psirz)
  end subroutine g_eqdsk_destructor

  subroutine flux_func_init(this, n_lag, indep_var)
    use magdif_conf, only: log
    class(flux_func), intent(inout) :: this
    integer, intent(in) :: n_lag
    real(dp), intent(in), dimension(:) :: indep_var

    if (n_lag >= size(indep_var)) then
       write (log%msg, '("Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points", i0)') n_lag, size(indep_var)
       if (log%err) call log%write
       return
    end if
    call flux_func_destructor(this)
    this%n_lag = n_lag
    this%n_var = size(indep_var)
    allocate(this%indep_var(this%n_var))
    this%indep_var = indep_var
  end subroutine flux_func_init

  function flux_func_interp(this, sample, position, deriv) result(interp)
    use magdif_conf, only: log
    class(flux_func) :: this
    real(dp), intent(in) :: sample(:)
    real(dp), intent(in) :: position
    logical, intent(in), optional :: deriv
    real(dp) :: interp
    real(dp) :: lag_coeff(0:1, this%n_lag)
    integer :: k, kder

    kder = 0
    if (present(deriv)) then
       if (deriv) then
          kder = 1
       end if
    end if
    if (this%n_var /= size(sample)) then
       call log%msg_arg_size('flux_func_interp', 'this%n_var', 'size(sample)', &
            this%n_var, size(sample))
       if (log%err) call log%write
       error stop
    end if
    call binsearch(this%indep_var, lbound(this%indep_var, 1), position, k)
    ! ensure that polynomial sample points do not go below lower bound of 1
    if (k < 1 + this%n_lag / 2) then
       k = 1 + this%n_lag / 2
    ! ensure that polynomial sample points do not go above upper bound of this%n_var
    elseif (k > this%n_var - this%n_lag / 2 + 1) then
       k = this%n_var - this%n_lag / 2 + 1
    end if
    call plag_coeff(this%n_lag, 1, position, &
         this%indep_var(k - this%n_lag / 2:k + this%n_lag / 2 - 1), lag_coeff)
    interp = sum(sample(k - this%n_lag / 2:k + this%n_lag / 2 - 1) * lag_coeff(kder, :))
  end function flux_func_interp

  subroutine flux_func_destructor(this)
    type(flux_func), intent(inout) :: this

    if (allocated(this%indep_var)) deallocate(this%indep_var)
  end subroutine flux_func_destructor

  !> Set up arrays of cached values of flux functions.
  !>
  !> @nflux number of flux surfaces
  !> @half_step values are taken at flux surfaces (false) or between flux surfaces (true)
  !>
  !> For full-grid quantities, values are taken on flux surfaces with indices running
  !> from 0 to \p nflux, i.e. from the magnetic axis to the separatrix. An exception is
  !> made for \psi, where the index runs up to \p nflux +1. This value is extrapolated for
  !> finite differences in magdif::compute_presn() and magdif::compute_bn_nonres().
  !> For half-grid quantities, values are taken between two flux surfaces with indices
  !> running from 1 to \p nflux, i.e. from the triangle strip surrounding the magnetic
  !> axis to the triangle strip just inside the separatrix.
  subroutine flux_func_cache_init(this, nflux, half_step)
    class(flux_func_cache), intent(inout) :: this
    integer, intent(in) :: nflux
    logical, intent(in) :: half_step

    call flux_func_cache_destructor(this)
    if (half_step) then
       allocate(this%psi(nflux))
       allocate(this%rad(nflux))
       allocate(this%F(nflux))
       allocate(this%p(nflux))
       allocate(this%FdF_dpsi(nflux))
       allocate(this%dp_dpsi(nflux))
       allocate(this%q(nflux))
    else
       allocate(this%psi(0:nflux))
       allocate(this%rad(0:nflux))
       allocate(this%F(0:nflux))
       allocate(this%p(0:nflux))
       allocate(this%FdF_dpsi(0:nflux))
       allocate(this%dp_dpsi(0:nflux))
       allocate(this%q(0:nflux))
    end if
    this%psi = 0d0
    this%rad = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
  end subroutine flux_func_cache_init

  subroutine flux_func_cache_destructor(this)
    type(flux_func_cache), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%F)) deallocate(this%F)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%FdF_dpsi)) deallocate(this%FdF_dpsi)
    if (allocated(this%dp_dpsi)) deallocate(this%dp_dpsi)
    if (allocated(this%q)) deallocate(this%q)
  end subroutine flux_func_cache_destructor

  subroutine neumaier_accumulator_real_init(this)
    class(neumaier_accumulator_real), intent(inout) :: this
    this%sum = 0d0
    this%c = 0d0
    this%t = 0d0
  end subroutine neumaier_accumulator_real_init

  subroutine neumaier_accumulator_real_add(this, val)
    class(neumaier_accumulator_real), intent(inout) :: this
    real(dp) :: val
    this%t = this%sum + val
    if (abs(this%sum) >= abs(val)) then
       this%c = this%c + ((this%sum - this%t) + val)
    else
       this%c = this%c + ((val - this%t) + this%sum)
    end if
    this%sum = this%t
  end subroutine neumaier_accumulator_real_add

  function neumaier_accumulator_real_get_sum(this) result(acc)
    class(neumaier_accumulator_real), intent(inout) :: this
    real(dp) :: acc
    acc = this%sum + this%c
  end function neumaier_accumulator_real_get_sum

  subroutine neumaier_accumulator_complex_init(this)
    class(neumaier_accumulator_complex), intent(inout) :: this
    call this%real_part%init
    call this%imag_part%init
  end subroutine neumaier_accumulator_complex_init

  subroutine neumaier_accumulator_complex_add(this, val)
    class(neumaier_accumulator_complex), intent(inout) :: this
    complex(dp) :: val
    call this%real_part%add(real(val))
    call this%imag_part%add(aimag(val))
  end subroutine neumaier_accumulator_complex_add

  function neumaier_accumulator_complex_get_sum(this) result(acc)
    class(neumaier_accumulator_complex), intent(inout) :: this
    complex(dp) :: acc
    acc = cmplx(this%real_part%get_sum(), this%imag_part%get_sum(), dp)
  end function neumaier_accumulator_complex_get_sum

  subroutine heapsort_complex(array, comparison)
    complex(dp), intent(inout) :: array(0:)
    interface
       function comparison(val1, val2)
         import :: dp
         complex(dp), intent(in) :: val1, val2
         logical :: comparison
       end function comparison
    end interface
    integer :: n, k, child, root
    complex(dp) :: temp
    n = size(array)
    do k = (n - 2) / 2, 0, -1
       call siftdown(k, n);
    end do
    do k = n - 1, 1, -1
       temp = array(0)
       array(0) = array(k)
       array(k) = temp;
       call siftdown(0, k)
    end do
  contains
    subroutine siftdown(start, bottom)
      integer, intent(in) :: start, bottom
      root = start
      do while (root * 2 + 1 < bottom)
         child = root * 2 + 1
         if (child + 1 < bottom) then
            if (comparison(array(child), array(child+1))) then
               child = child + 1
            end if
         end if
         if (comparison(array(root), array(child))) then
            temp = array(child)
            array(child) = array (root)
            array(root) = temp
            root = child
         else
            return
         end if
      end do
    end subroutine siftdown
  end subroutine heapsort_complex

  function complex_abs_asc(val1, val2)
    complex(dp), intent(in) :: val1, val2
    logical :: complex_abs_asc
    complex_abs_asc = abs(val1) < abs(val2)
  end function complex_abs_asc

  pure function interleave_vv(first_v, second_v, num) result(merged)
    integer, intent(in) :: first_v(:), second_v(:), num
    integer :: merged(num)
    integer :: k
    merged = merge(first_v, second_v, [([.true., .false.], k = 1, num / 2)])
  end function interleave_vv

  pure function interleave_vs(first_v, second_s, num) result(merged)
    integer, intent(in) :: first_v(:), second_s, num
    integer :: merged(num)
    integer :: k
    merged = merge(first_v, [(second_s, k = 1, num)], &
         [([.true., .false.], k = 1, num / 2)])
  end function interleave_vs

  pure function interleave_sv(first_s, second_v, num) result(merged)
    integer, intent(in) :: first_s, second_v(:), num
    integer :: merged(num)
    integer :: k
    merged = merge([(first_s, k = 1, num)], second_v, &
         [([.true., .false.], k = 1, num / 2)])
  end function interleave_sv

  pure function interleave_ss(first_s, second_s, num) result(merged)
    integer, intent(in) :: first_s, second_s, num
    integer :: merged(num)
    integer :: k
    merged = [([first_s, second_s], k = 1, num / 2)]
  end function interleave_ss

end module magdif_util
