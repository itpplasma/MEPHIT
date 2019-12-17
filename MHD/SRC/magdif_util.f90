module magdif_util

  use from_nrtype, only: dp

  implicit none

  private

  public :: ring_centered_avg_coord, assemble_sparse

  type, public :: g_eqdsk
    character(len = 1024) :: fname
    integer :: nw 
    real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime, qpsi
  contains
    procedure :: read => g_eqdsk_read
    final :: g_eqdsk_destructor
  end type g_eqdsk

  character(len = *), parameter :: geqdsk_2000 = '(6a8,3i4)'
  character(len = *), parameter :: geqdsk_2020 = '(5e16.9)'

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
    character(len = 1024), intent(in) :: fname
    character(len = 10) :: text(6)
    integer :: fid, k, nh, idum
    real(dp) :: xdum

    call g_eqdsk_destructor(this)
    this%fname = fname
    open(newunit = fid, file = this%fname)
    read(fid, geqdsk_2000) (text(k), k = 1, 6), idum, this%nw, nh
    allocate(this%fpol(this%nw))
    allocate(this%pres(this%nw))
    allocate(this%ffprim(this%nw))
    allocate(this%pprime(this%nw))
    allocate(this%qpsi(this%nw))
    read(fid, geqdsk_2020) (xdum, k = 1, 20)
    read(fid, geqdsk_2020) (this%fpol(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%pres(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%ffprim(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%pprime(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (xdum, k = 1, this%nw * nh)
    read(fid, geqdsk_2020) (this%qpsi(k), k = 1, this%nw)
    close(fid)
    ! convert tesla meter to gauss centimeter
    this%fpol = this%fpol * 1e6
    ! convert pascal to barye
    this%pres = this%pres * 1e1
    ! convert tesla to gauss
    this%ffprim = this%ffprim * 1e4
    ! convert pascal per weber to barye per maxwell
    this%pprime = this%pprime * 1e-7
  end subroutine g_eqdsk_read

  subroutine g_eqdsk_destructor(this)
    type(g_eqdsk) :: this

    if (allocated(this%fpol)) deallocate(this%fpol)
    if (allocated(this%pres)) deallocate(this%pres)
    if (allocated(this%ffprim)) deallocate(this%ffprim)
    if (allocated(this%pprime)) deallocate(this%pprime)
    if (allocated(this%qpsi)) deallocate(this%qpsi)
  end subroutine g_eqdsk_destructor

  subroutine flux_func_init(this, n_lag, nw, nflux, psi, psi_half)
    use magdif_config, only: logfile, log_err
    class(flux_func), intent(inout) :: this
    integer, intent(in) :: n_lag
    integer, intent(in) :: nw
    integer, intent(in) :: nflux
    real(dp), intent(in), dimension(0:nflux) :: psi
    real(dp), intent(in), dimension(1:nflux) :: psi_half
    real(dp) :: coeff(n_lag), psi_eqd(nw)
    integer :: k, kf

    if (n_lag >= nw) then
       if (log_err) write (logfile, *) 'Lagrange polynomial order', n_lag, &
            'must be lower than number of sample points', nw
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
       ! binsrc expects strictly monotonically increasing arrays, so we reverse psi_eqd
       call binsrc(psi_eqd(nw:1:-1), 1, nw, psi(kf), k)
       k = nw - (k - 1) + 1
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
       ! binsrc expects strictly monotonically increasing arrays, so we reverse psi_eqd
       call binsrc(psi_eqd(nw:1:-1), 1, nw, psi_half(kf), k)
       k = nw - (k - 1) + 1
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
