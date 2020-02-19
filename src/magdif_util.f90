module magdif_util

  use from_nrtype, only: dp

  implicit none

  private

  public :: clight, imun, interp_psi_pol, ring_centered_avg_coord, assemble_sparse

  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

  type, public :: sign_convention
    integer :: exp_Bpol, sgn_cyl, sgn_dpsi, sgn_Btor, sgn_Itor, &
         sgn_F, sgn_dp_dpsi, sgn_q, sgn_Bpol, sgn_pol, index
  end type sign_convention

  type, public :: g_eqdsk
    type(sign_convention) :: cocos
    character(len = 1024) :: fname
    character(len = 10) :: text(6)
    integer :: nw, nh, nbbbs, limitr
    real(dp) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, &
         current
    real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime, qpsi, &
         rbbbs, zbbbs, rlim, zlim
    real(dp), dimension(:, :), allocatable :: psirz
  contains
    procedure :: read => g_eqdsk_read
    procedure :: classify => g_eqdsk_classify
    procedure :: standardise => g_eqdsk_standardise
    procedure :: scale => g_eqdsk_scale
    procedure :: write => g_eqdsk_write
    final :: g_eqdsk_destructor
  end type g_eqdsk

  character(len = *), parameter :: geqdsk_2000 = '(6a8,3i4)'
  character(len = *), parameter :: geqdsk_2020 = '(5e16.9)'
  character(len = *), parameter :: geqdsk_2022 = '(2i5)'

  type, public :: flux_func
    private
    integer :: n_lag, nw, nflux
    real(dp), dimension(:, :), allocatable :: lag_coeff, lag_coeff_half
    integer, dimension(:), allocatable :: interval, interval_half
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

contains

  function interp_psi_pol(r, z) result(psi_pol)
    use field_eq_mod, only: psif

    real(dp), intent(in) :: r, z
    real(dp) :: psi_pol

    real(dp) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
    psi_pol = psif
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
    integer, intent(in)  :: nrow
    complex(dp), intent(in)  :: d(nrow)
    complex(dp), intent(in)  :: du(nrow)
    integer, intent(out) :: nz
    integer, intent(out) :: irow(2*nrow), icol(2*nrow)
    complex(dp), intent(out) :: aval(2*nrow)

    integer :: k

    irow(1) = 1
    icol(1) = 1
    aval(1) = d(1)

    irow(2) = nrow
    icol(2) = 1
    aval(2) = du(nrow)

    do k = 2,nrow
       ! off-diagonal
       irow(2*k-1) = k-1
       icol(2*k-1) = k
       aval(2*k-1) = du(k-1)

       ! diagonal
       irow(2*k) = k
       icol(2*k) = k
       aval(2*k) = d(k)
    end do

    nz = 2*nrow
  end subroutine assemble_sparse

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
  end subroutine g_eqdsk_read

  subroutine g_eqdsk_classify(this)
    use magdif_config, only: log_msg, log_write, log_warn
    class(g_eqdsk), intent(inout) :: this

    this%cocos = sign_convention(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    this%cocos%exp_Bpol = 0  ! specified by G EQDSK format and assumed in field_divB0.f90
    this%cocos%sgn_cyl = +1  ! specified by G EQDSK format and assumed in field_divB0.f90
    this%cocos%sgn_dpsi = sign_array([this%sibry - this%simag], 'SIMAG/SIBRY', this%fname)
    this%cocos%sgn_Btor = sign_array([this%bcentr], 'BCENTR', this%fname)
    this%cocos%sgn_Itor = sign_array([this%current], 'CURRENT', this%fname)
    this%cocos%sgn_F = sign_array(this%fpol, 'FPOL', this%fname)
    this%cocos%sgn_dp_dpsi = sign_array(this%pprime, 'PPRIME', this%fname)
    this%cocos%sgn_q = sign_array(this%qpsi, 'QPSI', this%fname)
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
       log_msg = 'COCOS index could not be determined for ' // trim(this%fname)
       if (log_warn) call log_write
    end if

  contains
    function sign_array(array, name, fname) result(sign)
      use magdif_config, only: log_msg, log_write, log_warn
      real(dp), intent(in), dimension(:) :: array
      character(len = *), intent(in) :: name, fname
      integer :: sign
      if (all(array > 0d0)) then
         sign = +1
      elseif (all(array < 0d0)) then
         sign = -1
      else
         sign = 0
         write (log_msg, '("Sign of ", a, " is inconsistent in ", a)') &
              trim(name), trim(fname)
         if (log_warn) call log_write
      end if
    end function sign_array
  end subroutine g_eqdsk_classify

  subroutine g_eqdsk_standardise(this)
    use magdif_config, only: log_msg, log_write, log_info, log_warn
    class(g_eqdsk), intent(inout) :: this
    character(len = *), parameter :: invert_fmt = '("Inverting sign of ", a, "...")', &
         incons_fmt = '("Signs of ", a, " and ", a, " are inconsistent")'

    if (this%cocos%sgn_Btor /= this%cocos%sgn_F) then
       write (log_msg, incons_fmt) 'FPOL', 'BCENTR'
       if (log_warn) call log_write
       write (log_msg, invert_fmt) 'FPOL'
       if (log_info) call log_write
       this%fpol = -this%fpol
       this%cocos%sgn_F = -this%cocos%sgn_F
    end if
    if (this%cocos%sgn_dp_dpsi /= -this%cocos%sgn_dpsi) then
       write (log_msg, incons_fmt) 'PPRIME', 'SIMAG/SIBRY'
       if (log_warn) call log_write
       write (log_msg, invert_fmt) 'PPRIME'
       if (log_info) call log_write
       this%pprime = -this%pprime
       this%cocos%sgn_dp_dpsi = -this%cocos%sgn_dp_dpsi
    end if
    if (this%cocos%sgn_Bpol /= -1) then
       write (log_msg, invert_fmt) 'SIMAG/SIBRY, PSIRZ, PPRIME, FFPRIM'
       this%simag = -this%simag
       this%sibry = -this%sibry
       this%psirz = -this%psirz
       this%cocos%sgn_dpsi = -this%cocos%sgn_dpsi
       this%pprime = -this%pprime
       this%cocos%sgn_dp_dpsi = -this%cocos%sgn_dp_dpsi
       this%ffprim = -this%ffprim
       this%cocos%sgn_Bpol = -this%cocos%sgn_Bpol
    end if
    if (this%cocos%sgn_pol /= -1) then
       write (log_msg, invert_fmt) 'QPSI'
       this%qpsi = -this%qpsi
       this%cocos%sgn_q = -this%cocos%sgn_q
       this%cocos%sgn_pol = -this%cocos%sgn_pol
    end if
    this%cocos%index = 3
  end subroutine g_eqdsk_standardise

  subroutine g_eqdsk_scale(this, gamma)
    class(g_eqdsk), intent(inout) :: this
    integer, intent(in) :: gamma
    real(dp) :: unscaled_rcentr, r_shift

    unscaled_rcentr = this%rcentr
    r_shift = this%rcentr * (gamma - 1)
    this%rcentr = this%rcentr * gamma
    this%rleft = this%rcentr - 0.5d0 * this%rdim
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
    if (allocated(this%psirz)) deallocate(this%psirz)
  end subroutine g_eqdsk_destructor

  subroutine flux_func_init(this, n_lag, nw, nflux, psi, psi_half)
    use magdif_config, only: log_msg, log_write, log_err
    class(flux_func), intent(inout) :: this
    integer, intent(in) :: n_lag
    integer, intent(in) :: nw
    integer, intent(in) :: nflux
    real(dp), intent(in), dimension(0:nflux) :: psi
    real(dp), intent(in), dimension(1:nflux) :: psi_half
    real(dp) :: coeff(n_lag), psi_eqd(nw)
    integer :: k, kf

    if (n_lag >= nw) then
       write (log_msg, '("Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points", i0)') n_lag, nw
       if (log_err) call log_write
       return
    end if
    call flux_func_destructor(this)
    this%n_lag = n_lag
    this%nw = nw
    this%nflux = nflux
    psi_eqd = psi(0) + [(dble(k-1) / dble(nw-1) * (psi(nflux) - psi(0)), k = 1, nw)]
    allocate(this%lag_coeff(n_lag, nflux-1))
    allocate(this%lag_coeff_half(n_lag, nflux))
    allocate(this%interval(nflux-1))
    allocate(this%interval_half(nflux))
    this%lag_coeff = 0d0
    this%lag_coeff_half = 0d0
    this%interval = 0
    this%interval_half = 0
    do kf = 1, nflux-1
       if (psi_eqd(1) < psi_eqd(nw)) then
          call binsrc(psi_eqd, 1, nw, psi(kf), k)
       else
          ! binsrc expects strictly monotonically increasing arrays, so we reverse psi_eqd
          call binsrc(psi_eqd(nw:1:-1), 1, nw, psi(kf), k)
          k = nw - (k - 1) + 1
       end if
       if (k < n_lag / 2 + 1) then
          k = n_lag / 2 + 1
       elseif (k > nw - n_lag / 2 + 1) then
          k = nw - n_lag / 2 + 1
       end if
       this%interval(kf) = k
       call plag_coeff(n_lag, 0, psi(kf), psi_eqd(k-n_lag/2:k+n_lag/2-1), coeff)
       this%lag_coeff(:, kf) = coeff
    end do
    do kf = 1, nflux
       if (psi_eqd(1) < psi_eqd(nw)) then
          call binsrc(psi_eqd, 1, nw, psi_half(kf), k)
       else
          ! binsrc expects strictly monotonically increasing arrays, so we reverse psi_eqd
          call binsrc(psi_eqd(nw:1:-1), 1, nw, psi_half(kf), k)
          k = nw - (k - 1) + 1
       end if
       if (k < n_lag / 2 + 1) then
          k = n_lag / 2 + 1
       elseif (k > nw - n_lag / 2 + 1) then
          k = nw - n_lag / 2 + 1
       end if
       this%interval_half(kf) = k
       call plag_coeff(n_lag, 0, psi_half(kf), psi_eqd(k-n_lag/2:k+n_lag/2-1), coeff)
       this%lag_coeff_half(:, kf) = coeff
    end do
  end subroutine flux_func_init

  function flux_func_interp(this, sample_eqd, kf, half_step) result(interp)
    class(flux_func) :: this
    real(dp), intent(in) :: sample_eqd(this%nw)
    integer, intent(in) :: kf
    logical, intent(in) :: half_step
    real(dp) :: interp
    integer :: k

    if (half_step) then
       if (0 < kf .and. kf <= this%nflux) then
          k = this%interval_half(kf)
          interp = sum(sample_eqd(k - this%n_lag / 2:k + this%n_lag / 2 - 1) * &
               this%lag_coeff_half(:, kf))
       else
          interp = 0d0
       end if
    else
       if (kf == 0) then
          interp = sample_eqd(1)
       elseif (kf == this%nflux) then
          interp = sample_eqd(this%nw)
       elseif (0 < kf .and. kf < this%nflux) then
          k = this%interval(kf)
          interp = sum(sample_eqd(k - this%n_lag / 2:k + this%n_lag / 2 - 1) * &
               this%lag_coeff(:, kf))
       else
          interp = 0d0
       end if
    end if
  end function flux_func_interp

  subroutine flux_func_destructor(this)
    type(flux_func), intent(inout) :: this

    if (allocated(this%lag_coeff)) deallocate(this%lag_coeff)
    if (allocated(this%lag_coeff_half)) deallocate(this%lag_coeff_half)
    if (allocated(this%interval)) deallocate(this%interval)
    if (allocated(this%interval_half)) deallocate(this%interval_half)
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
       allocate(this%F(nflux))
       allocate(this%p(nflux))
       allocate(this%FdF_dpsi(nflux))
       allocate(this%dp_dpsi(nflux))
       allocate(this%q(nflux))
    else
       allocate(this%psi(0:nflux))
       allocate(this%F(0:nflux))
       allocate(this%p(0:nflux))
       allocate(this%FdF_dpsi(0:nflux))
       allocate(this%dp_dpsi(0:nflux))
       allocate(this%q(0:nflux))
    end if
    this%psi = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
  end subroutine flux_func_cache_init

  subroutine flux_func_cache_destructor(this)
    type(flux_func_cache), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%F)) deallocate(this%F)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%FdF_dpsi)) deallocate(this%FdF_dpsi)
    if (allocated(this%dp_dpsi)) deallocate(this%dp_dpsi)
    if (allocated(this%q)) deallocate(this%q)
  end subroutine flux_func_cache_destructor

end module magdif_util
