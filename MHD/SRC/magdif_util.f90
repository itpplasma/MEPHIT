module magdif_util

  use from_nrtype, only: dp

  implicit none

  private

  type, public :: g_eqdsk
    character(len = 1024) :: fname
    integer :: nw 
    real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime
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

contains

  subroutine g_eqdsk_read(this, fname)
    class(g_eqdsk), intent(inout) :: this
    character(len = 1024), intent(in) :: fname
    character(len = 10) :: text(6)
    integer :: fid, k, idum
    real(dp) :: xdum

    call g_eqdsk_destructor(this)
    this%fname = fname
    open(newunit = fid, file = this%fname)
    read(fid, geqdsk_2000) (text(k), k = 1, 6), idum, this%nw, idum
    allocate(this%fpol(this%nw))
    allocate(this%pres(this%nw))
    allocate(this%ffprim(this%nw))
    allocate(this%pprime(this%nw))
    read(fid, geqdsk_2020) (xdum, k = 1, 20)
    read(fid, geqdsk_2020) (this%fpol(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%pres(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%ffprim(k), k = 1, this%nw)
    read(fid, geqdsk_2020) (this%pprime(k), k = 1, this%nw)
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
  end subroutine g_eqdsk_destructor

  subroutine flux_func_init(this, n_lag, nw, nflux, psi)
    use magdif_config, only: logfile, log_err
    class(flux_func), intent(inout) :: this
    integer, intent(in) :: n_lag
    integer, intent(in) :: nw
    integer, intent(in) :: nflux
    real(dp), intent(in), dimension(0:nflux) :: psi
    real(dp) :: coeff(n_lag), psi_eqd(nw), psi_half(nflux)
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
    ! use linear interpolation for half-grid steps for now
    psi_half = 0.5d0 * (psi(0:nflux-1) + psi(1:nflux))
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

end module magdif_util
