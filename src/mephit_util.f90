module mephit_util

  use iso_fortran_env, only: dp => real64
  use iso_c_binding, only: c_ptr, c_null_ptr, c_double_complex

  implicit none

  private

  ! types and associated procedures
  public :: fft_t
  public :: sign_convention, g_eqdsk
  public :: func1d_t, func1d_init, func1d_deinit, func1d_write, func1d_read, func1d_read_formatted
  public :: neumaier_accumulator_real, neumaier_accumulator_complex

  ! utility procedures
  public :: init_field, deinit_field, interp_psi_pol, interp1d, resample1d, &
       pos_angle, linspace, straight_cyl2bent_cyl, bent_cyl2straight_cyl, zd_cross, dd_cross, &
       binsearch, interleave, heapsort_real, heapsort_complex, complex_abs_asc, complex_abs_desc, &
       arnoldi_break, hessenberg_eigvals, hessenberg_eigvecs, &
       gauss_legendre_unit_interval, C_F_string

  ! module variables
  public :: pi, clight, ev2erg, imun

  real(dp), parameter :: pi = 4d0 * atan(1d0)
  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  real(dp), parameter :: ev2erg = 1.6021766d-12      !< convert eV to erg
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

  ! types and interfaces

  type :: fft_t
     integer :: N = 0
     complex(c_double_complex), dimension(:), pointer :: samples => null()
     complex(c_double_complex), dimension(:), pointer, private :: modes => null()
     type(c_ptr), private :: plan_p = c_null_ptr, samples_p = c_null_ptr, modes_p = c_null_ptr
   contains
     procedure :: init => fft_init
     procedure :: apply => fft_apply
     procedure :: deinit => fft_deinit
  end type fft_t

  type :: sign_convention
    integer :: exp_Bpol, sgn_cyl, sgn_dpsi, sgn_Btor, sgn_Itor, &
         sgn_F, sgn_q, sgn_Bpol, sgn_pol, index
  end type sign_convention

  type :: g_eqdsk
    type(sign_convention) :: cocos
    character(len = 1024) :: fname
    character(len = 48) :: header
    integer :: nw, nh, nbbbs, limitr
    real(dp) :: rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, &
         current
    real(dp), dimension(:), allocatable :: fpol, pres, ffprim, pprime, qpsi, &
         rbbbs, zbbbs, rlim, zlim, psi_eqd, R_eqd, Z_eqd, fprime
    real(dp), dimension(:, :), allocatable :: psirz
  contains
    procedure :: read => g_eqdsk_read
    procedure :: check_consistency => g_eqdsk_check_consistency
    procedure :: classify => g_eqdsk_classify
    procedure :: standardise => g_eqdsk_standardise
    procedure :: scale => g_eqdsk_scale
    procedure :: write => g_eqdsk_write
    procedure :: grad_shafranov_normalization => g_eqdsk_grad_shafranov_normalization
    procedure :: import_hdf5 => g_eqdsk_import_hdf5
    procedure :: export_hdf5 => g_eqdsk_export_hdf5
    procedure :: deinit => g_eqdsk_deinit
  end type g_eqdsk

  character(len = *), parameter :: geqdsk_2000 = '(a48, 3i4)'
  character(len = *), parameter :: geqdsk_2020 = '(5es16.9)'
  character(len = *), parameter :: geqdsk_2022 = '(2i5)'

  character(len = *), parameter :: &
       incons_fmt = '("Signs of ", a, " and ", a, " are inconsistent.")', &
       invert_fmt = '("Inverting sign of ", a, "...")', &
       unscaled_fmt = '("Flux in ", a, " is not normalized by 2 pi.")', &
       rescale_fmt = '("Rescaling ", a, "...")'

  !> 1D function sample
  type :: func1d_t
     real(dp), dimension(:), allocatable :: x
     real(dp), dimension(:), allocatable :: y
  end type func1d_t

  type :: neumaier_accumulator_real
     private
     real(dp) :: sum, c, t
   contains
     procedure :: init => neumaier_accumulator_real_init
     procedure :: add => neumaier_accumulator_real_add
     procedure :: get_sum => neumaier_accumulator_real_get_sum
  end type neumaier_accumulator_real

  type :: neumaier_accumulator_complex
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

  !> Set module variables for initialization of subroutine field
  subroutine init_field(equil)
    use field_mod, only: icall, ipert, iequil
    use field_eq_mod, only: skip_read, icall_eq, nwindow_r, nwindow_z, &
         nrad, nzet, psi_axis, psi_sep, btf, rtf, splfpol, rad, zet, psi, psi0
    type(g_eqdsk), intent(in) :: equil
    real(dp) :: dum

    ! compute equiibrium field
    ipert = 0
    iequil = 1
    ! default values - TODO: options in config_t
    nwindow_r = 0
    nwindow_z = 0
    ! don't let subroutine field read from input file
    icall = 1
    ! let subroutine field_eq do only initialization
    icall_eq = -1
    ! use preprocessed gfile data
    skip_read = .true.
    nrad = equil%nw
    nzet = equil%nh
    allocate(rad(nrad), zet(nzet), psi0(nrad, nzet), psi(nrad, nzet), splfpol(0:5, nrad))
    psi_axis = equil%simag * 1d-8
    psi_sep = equil%sibry * 1d-8
    btf = equil%bcentr * 1d-4
    rtf = equil%rcentr * 1d-2
    splfpol(0, :) = equil%fpol * 1d-6
    psi(:, :) = equil%psirz * 1d-8
    rad(:) = equil%R_eqd * 1d-2
    zet(:) = equil%Z_eqd * 1d-2
    call field_eq(0d0, 0d0, 0d0, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
  end subroutine init_field

  subroutine deinit_field()
    use field_eq_mod, only: ima, imi, jma, jmi, splpsi, splfpol, rad, zet, psi, psi0, ipoint

    if (allocated(ima)) deallocate(ima)
    if (allocated(imi)) deallocate(imi)
    if (allocated(jma)) deallocate(jma)
    if (allocated(jmi)) deallocate(jmi)
    if (allocated(splpsi)) deallocate(splpsi)
    if (allocated(splfpol)) deallocate(splfpol)
    if (allocated(rad)) deallocate(rad)
    if (allocated(zet)) deallocate(zet)
    if (allocated(psi)) deallocate(psi)
    if (allocated(psi0)) deallocate(psi0)
    if (allocated(ipoint)) deallocate(ipoint)
  end subroutine deinit_field

  function interp_psi_pol(r, z) result(psi_pol)
    use field_eq_mod, only: psif, psib

    real(dp), intent(in) :: r, z
    real(dp) :: psi_pol

    real(dp) :: dum

    call field(r, 0d0, z, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
    ! field_divB0.f90 adds psib (SIBRY, i.e., flux at the boundary) to interpolated psi[f]
    psi_pol = psif - psib
  end function interp_psi_pol

  pure elemental function pos_angle(atan_angle)
    real(dp), intent(in) :: atan_angle
    real(dp) :: pos_angle

    pos_angle = atan_angle
    if (pos_angle < 0d0) then
       pos_angle = pos_angle + 2d0 * pi
    end if
    if (pos_angle > 2d0 * pi) then
       pos_angle = pos_angle - 2d0 * pi
    end if
  end function pos_angle

  function linspace(lo, hi, cnt, excl_lo, excl_hi)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: cnt, excl_lo, excl_hi
    real(dp) :: linspace(cnt)
    real(dp) :: step
    integer :: k

    step = (hi - lo) / dble(cnt - 1 + excl_lo + excl_hi)
    linspace = lo + [(k * step, k = excl_lo, cnt - 1 + excl_lo)]
    if (excl_hi == 0) linspace(cnt) = hi
  end function linspace


  !> Binary search to find index \p k of ordered array \p x so that \p xi lies between
  !> x(k-1) and x(k). The lower bound of \p x is given by \p lb.
  !>
  !> When \p x is strictly monotonically increasing, \f$ x_{k-1} < \xi < x_{k} \f$.
  !> When \p x is strictly monotonically decreasing, \f$ x_{k-1} > \xi > x_{k} \f$.
  !> It is not checked whether \p x is ordered; only the first and last value are used
  !> to determine wether the array values are increasing or decrasing.
  subroutine binsearch(x, lb, xi, k)
    integer, intent(in) :: lb
    real(dp), intent(in) :: x(lb:)
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
    use mephit_conf, only: logger
    use fgsl, only: fgsl_size_t, fgsl_double, fgsl_int, fgsl_success, &
         fgsl_integration_glfixed_point, fgsl_integration_glfixed_table, &
         fgsl_integration_glfixed_table_alloc, fgsl_integration_glfixed_table_free
    integer, intent(in) :: order
    real(fgsl_double), dimension(:), intent(out) :: points, weights
    type(fgsl_integration_glfixed_table) :: table
    integer(fgsl_size_t) :: k
    integer(fgsl_int) :: err
    if (order /= size(points)) then
       call logger%msg_arg_size('gauss_legendre_unit_interval', 'order', 'size(points)', &
            order, size(points))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (order /= size(weights)) then
       call logger%msg_arg_size('gauss_legendre_unit_interval', 'order', 'size(weights)', &
            order, size(weights))
       if (logger%err) call logger%write_msg
       error stop
    end if
    table = fgsl_integration_glfixed_table_alloc(int(order, fgsl_size_t))
    do k = 1, int(order, fgsl_size_t)
       err = fgsl_integration_glfixed_point(0d0, 1d0, k-1, points(k), weights(k), table)
       if (err /= fgsl_success) then
          write (logger%msg, '("fgsl_integration_glfixed_point returned error ", i0)') err
          if (logger%err) call logger%write_msg
          error stop
       end if
    end do
    call fgsl_integration_glfixed_table_free(table)
  end subroutine gauss_legendre_unit_interval

  subroutine fft_init(fft, N)
    use iso_c_binding, only: c_size_t, c_f_pointer
    use fftw3, only: fftw_alloc_complex, fftw_plan_dft_1d, &
         FFTW_FORWARD, FFTW_PATIENT, FFTW_DESTROY_INPUT
    class(fft_t), intent(inout) :: fft
    integer, intent(in) :: N

    call fft_deinit(fft)
    fft%N = N
    fft%samples_p = fftw_alloc_complex(int(N, c_size_t))
    call c_f_pointer(fft%samples_p, fft%samples, [N])
    fft%modes_p = fftw_alloc_complex(int(N, c_size_t))
    call c_f_pointer(fft%modes_p, fft%modes, [N])
    fft%plan_p = fftw_plan_dft_1d(N, fft%samples, fft%modes, &
         FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))
  end subroutine fft_init

  subroutine fft_apply(fft, m_min, m_max, modes)
    use iso_c_binding, only: c_associated
    use fftw3, only: fftw_execute_dft
    use mephit_conf, only: logger
    class(fft_t), intent(inout) :: fft
    integer, intent(in) :: m_min, m_max
    complex(dp), dimension(m_min:), intent(inout) :: modes

    if (ubound(modes, 1) /= m_max) then
       call logger%msg_arg_size('fft_apply', &
            'ubound(modes, 1)', 'm_max', &
            ubound(modes, 1), m_max)
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (.not. (c_associated(fft%plan_p) .and. &
         c_associated(fft%modes_p) .and. &
         c_associated(fft%samples_p))) then
       logger%msg = 'Attempt to dereference null pointer in fft_apply'
       if (logger%err) call logger%write_msg
       error stop
    end if
    call fftw_execute_dft(fft%plan_p, fft%samples, fft%modes)
    if (m_min >= 0 .and. m_max >= 0) then
       modes(:) = fft%modes(m_min+1:m_max+1) / dble(fft%N)
    else if (m_min < 0 .and. m_max < 0) then
       modes(:) = fft%modes(fft%N+m_min+1:fft%N+m_max+1) / dble(fft%N)
    else
       modes(m_min:-1) = fft%modes(fft%N+m_min+1:fft%N) / dble(fft%N)
       modes(0:m_max) = fft%modes(:m_max+1) / dble(fft%N)
    end if
  end subroutine fft_apply

  subroutine fft_deinit(fft)
    use iso_c_binding, only: c_associated, c_null_ptr
    use fftw3, only: fftw_destroy_plan, fftw_free
    class(fft_t), intent(inout) :: fft

    fft%N = 0
    fft%samples => null()
    fft%modes => null()
    if (c_associated(fft%plan_p)) then
       call fftw_destroy_plan(fft%plan_p)
       fft%plan_p = c_null_ptr
    end if
    if (c_associated(fft%modes_p)) then
       call fftw_free(fft%modes_p)
       fft%modes_p = c_null_ptr
    end if
    if (c_associated(fft%samples_p)) then
       call fftw_free(fft%samples_p)
       fft%samples_p = c_null_ptr
    end if
  end subroutine fft_deinit

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

  pure function zd_cross(z, d)
    complex(dp), intent(in) :: z(3)
    real(dp), intent(in) :: d(3)
    complex(dp) :: zd_cross(3)

    zd_cross = z([2, 3, 1]) * d([3, 1, 2]) - z([3, 1, 2]) * d([2, 3, 1])
  end function zd_cross

  pure function dd_cross(d1, d2)
    real(dp), intent(in) :: d1(3), d2(3)
    real(dp) :: dd_cross(3)

    dd_cross = d1([2, 3, 1]) * d2([3, 1, 2]) - d1([3, 1, 2]) * d2([2, 3, 1])
  end function dd_cross

  subroutine g_eqdsk_read(this, fname)
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    call g_eqdsk_deinit(this)
    this%fname = fname
    open(newunit = fid, file = this%fname, status = 'old', form = 'formatted', action = 'read')
    read (fid, geqdsk_2000) this%header, idum, this%nw, this%nh
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
    ! cache repeatedly used values
    allocate(this%fprime(this%nw))
    this%fprime(:) = this%ffprim / this%fpol
  end subroutine g_eqdsk_read

  subroutine g_eqdsk_check_consistency(this)
    use mephit_conf, only: logger
    class(g_eqdsk), intent(inout) :: this
    integer, parameter :: ignore = 3
    real(dp) :: deriv_eqd(this%nw), factor
    logical :: mask(this%nw)

    if (this%cocos%sgn_Btor /= this%cocos%sgn_F) then
       write (logger%msg, incons_fmt) 'FPOL', 'BCENTR'
       if (logger%warn) call logger%write_msg
       ! only modify array if signs are consistent
       if (this%cocos%sgn_F /= 0) then
          write (logger%msg, invert_fmt) 'FPOL'
          if (logger%info) call logger%write_msg
          this%fpol = -this%fpol
          this%cocos%sgn_F = -this%cocos%sgn_F
       end if
    end if
    call resample1d(this%psi_eqd, this%pres, this%psi_eqd, deriv_eqd, 3, .true.)
    ! ignore a few possibly unreliable values on both ends of the interval
    mask = .true.
    mask(1:1+ignore) = .false.
    mask(this%nw-ignore:this%nw) = .false.
    factor = sum(deriv_eqd / this%pprime, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2d0 * pi)) then
       write (logger%msg, unscaled_fmt) 'PPRIME'
       if (logger%warn) call logger%write_msg
       write (logger%msg, rescale_fmt) 'PPRIME'
       if (logger%info) call logger%write_msg
       this%pprime = this%pprime * (2d0 * pi)
    end if
    if (factor < 0d0) then
       write (logger%msg, incons_fmt) 'PPRIME', 'PRES/SIBRY-SIMAG'
       if (logger%warn) call logger%write_msg
       write (logger%msg, invert_fmt) 'PPRIME'
       if (logger%info) call logger%write_msg
       this%pprime = -this%pprime
    end if
    call resample1d(this%psi_eqd, this%fpol, this%psi_eqd, deriv_eqd, 3, .true.)
    deriv_eqd(:) = deriv_eqd * this%fpol
    factor = sum(deriv_eqd / this%ffprim, mask = mask) / dble(count(mask))
    if (abs(factor) >= sqrt(2d0 * pi)) then
       write (logger%msg, unscaled_fmt) 'FFPRIM'
       if (logger%warn) call logger%write_msg
       write (logger%msg, rescale_fmt) 'FFPRIM'
       if (logger%info) call logger%write_msg
       this%ffprim = this%ffprim * (2d0 * pi)
    end if
    if (factor < 0d0) then
       write (logger%msg, incons_fmt) 'FFPRIM', 'FPOL/SIBRY-SIMAG'
       if (logger%warn) call logger%write_msg
       write (logger%msg, invert_fmt) 'FFPRIM'
       if (logger%info) call logger%write_msg
       this%ffprim = -this%ffprim
    end if
  end subroutine g_eqdsk_check_consistency

  !> Estimates terms of Grad-Shafranov equation to determine sign_convention::exp_bpol.
  function g_eqdsk_grad_shafranov_normalization(this) result(gs_factor)
    use mephit_conf, only: logger
    class(g_eqdsk), intent(inout) :: this
    real(dp) :: gs_factor
    real(dp) :: Delta_R, Delta_Z, psi, gs_rhs, gs_lhs, dp_dpsi, FdF_dpsi
    integer :: kw, kh

    ! find evaluation point to the upper right of the magnetic axis
    call binsearch(this%R_eqd, lbound(this%R_eqd, 1), this%rmaxis, kw)
    call binsearch(this%Z_eqd, lbound(this%Z_eqd, 1), this%zmaxis, kh)
    kw = kw + 2
    kh = kh + 2
    psi = this%psirz(kw, kh)
    ! evaluate flux functions on RHS of GS equation
    dp_dpsi = interp1d(this%psi_eqd, this%pprime, psi, 3)
    FdF_dpsi = interp1d(this%psi_eqd, this%ffprim, psi, 3)
    gs_rhs = -4d0 * pi * this%R_eqd(kw) ** 2 * dp_dpsi - FdF_dpsi
    ! approximate differential operator on LHS of GS equation
    Delta_R = this%rdim / dble(this%nw - 1)
    Delta_Z = this%zdim / dble(this%nh - 1)
    gs_lhs = (this%psirz(kw, kh + 1) - 2d0 * psi + this%psirz(kw, kh - 1)) / Delta_Z ** 2 &
         + (this%psirz(kw + 1, kh) - 2d0 * psi + this%psirz(kw - 1, kh)) / Delta_R ** 2 &
         - (this%psirz(kw + 1, kh) - this%psirz(kw - 1, kh)) / (2d0 * Delta_R * this%R_eqd(kw))
    gs_factor = gs_lhs / gs_rhs
    write (logger%msg, '("Grad-Shafranov equation LHS / RHS: ", es24.16e3)') gs_factor
    if (logger%info) call logger%write_msg
  end function g_eqdsk_grad_shafranov_normalization

  function sign_array(array, name, most)
    use mephit_conf, only: logger
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
       write (logger%msg, '("Sign of ", a, " is inconsistent.")') trim(name)
       if (logger%warn) call logger%write_msg
    end if
  end function sign_array

  subroutine g_eqdsk_classify(this)
    use mephit_conf, only: logger
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
       logger%msg = 'COCOS index could not be determined for ' // trim(this%fname)
       if (logger%warn) call logger%write_msg
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
    use mephit_conf, only: logger
    class(g_eqdsk), intent(inout) :: this

    if (this%cocos%sgn_Bpol == +1) then
       write (logger%msg, invert_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (logger%info) call logger%write_msg
       this%simag = -this%simag
       this%sibry = -this%sibry
       this%psirz(:, :) = -this%psirz
       this%psi_eqd(:) = -this%psi_eqd
       this%cocos%sgn_dpsi = -this%cocos%sgn_dpsi
       this%pprime(:) = -this%pprime
       this%ffprim(:) = -this%ffprim
       this%fprime(:) = -this%fprime
       this%cocos%sgn_Bpol = -this%cocos%sgn_Bpol
    end if
    if (this%cocos%sgn_pol == +1) then
       write (logger%msg, invert_fmt) 'QPSI'
       if (logger%info) call logger%write_msg
       this%qpsi(:) = -this%qpsi
       this%cocos%sgn_q = -this%cocos%sgn_q
       this%cocos%sgn_pol = -this%cocos%sgn_pol
    end if
    if (this%cocos%exp_Bpol == 1) then
       write (logger%msg, unscaled_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (logger%warn) call logger%write_msg
       write (logger%msg, rescale_fmt) 'SIMAG, SIBRY, PSIRZ, PPRIME, FFPRIM'
       if (logger%info) call logger%write_msg
       this%simag = this%simag / (2d0 * pi)
       this%sibry = this%sibry / (2d0 * pi)
       this%psirz(:, :) = this%psirz / (2d0 * pi)
       this%psi_eqd(:) = this%psi_eqd / (2d0 * pi)
       this%pprime(:) = this%pprime * (2d0 * pi)
       this%ffprim(:) = this%ffprim * (2d0 * pi)
       this%fprime(:) = this%fprime * (2d0 * pi)
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
    this%fpol(:) = this%fpol * gamma
    this%ffprim(:) = this%ffprim * gamma
    this%pprime(:) = this%pprime / gamma
    this%psirz(:, :) = this%psirz * gamma
    this%qpsi(:) = this%qpsi / gamma
    this%rbbbs(:) = this%rbbbs + r_shift
    this%rlim(:) = this%rlim + r_shift
    ! auxiliary values
    this%psi_eqd(:) = this%psi_eqd * gamma
    this%R_eqd(:) = this%R_eqd + r_shift
  end subroutine g_eqdsk_scale

  subroutine g_eqdsk_write(this, fname)
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: fname
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    idum = 0
    xdum = 0d0
    this%fname = fname
    open(newunit = fid, file = this%fname, status = 'replace', form = 'formatted', action = 'write')
    write (fid, geqdsk_2000) this%header, idum, this%nw, this%nh
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

  subroutine g_eqdsk_import_hdf5(this, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call g_eqdsk_deinit(this)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/exp_Bpol', this%cocos%exp_Bpol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_cyl', this%cocos%sgn_cyl)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_dpsi', this%cocos%sgn_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Btor', this%cocos%sgn_Btor)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Itor', this%cocos%sgn_Itor)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_F', this%cocos%sgn_F)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_q', this%cocos%sgn_q)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Bpol', this%cocos%sgn_Bpol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_pol', this%cocos%sgn_pol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cocos/index', this%cocos%index)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/fname', this%fname)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/header', this%header)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nw', this%nw)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nh', this%nh)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nbbbs', this%nbbbs)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/limitr', this%limitr)
    allocate(this%fpol(this%nw), this%pres(this%nw), this%ffprim(this%nw), this%pprime(this%nw), &
         this%qpsi(this%nw), this%rbbbs(this%nbbbs), this%zbbbs(this%nbbbs), &
         this%rlim(this%limitr), this%zlim(this%limitr), this%psi_eqd(this%nw), &
         this%R_eqd(this%nw), this%Z_eqd(this%nh), this%fprime(this%nw), &
         this%psirz(this%nw, this%nh))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rdim', this%rdim)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/zdim', this%zdim)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rcentr', this%rcentr)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rleft', this%rleft)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/zmid', this%zmid)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rmaxis', this%rmaxis)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/zmaxis', this%zmaxis)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/simag', this%simag)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/sibry', this%sibry)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/bcentr', this%bcentr)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/current', this%current)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/fpol', this%fpol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/pres', this%pres)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/ffprim', this%ffprim)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/pprime', this%pprime)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/qpsi', this%qpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rbbbs', this%rbbbs)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/zbbbs', this%zbbbs)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rlim', this%rlim)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/zlim', this%zlim)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi_eqd', this%psi_eqd)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_eqd', this%R_eqd)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_eqd', this%Z_eqd)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/fprime', this%fprime)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psirz', this%psirz)
    call h5_close(h5id_root)
  end subroutine g_eqdsk_import_hdf5

  subroutine g_eqdsk_export_hdf5(this, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(g_eqdsk), intent(in) :: this
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/cocos/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/exp_Bpol', this%cocos%exp_Bpol)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_cyl', this%cocos%sgn_cyl)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_dpsi', this%cocos%sgn_dpsi)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Btor', this%cocos%sgn_Btor)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Itor', this%cocos%sgn_Itor)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_F', this%cocos%sgn_F)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_q', this%cocos%sgn_q)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_Bpol', this%cocos%sgn_Bpol)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/sgn_pol', this%cocos%sgn_pol)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cocos/index', this%cocos%index)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/fname', this%fname, &
         comment = 'original GEQDSK filename')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/header', this%header)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nw', this%nw)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nh', this%nh)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nbbbs', this%nbbbs)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/limitr', this%limitr)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rdim', this%rdim, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/zdim', this%zdim, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rcentr', this%rcentr, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rleft', this%rleft, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/zmid', this%zmid, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rmaxis', this%rmaxis, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/zmaxis', this%zmaxis, unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/simag', this%simag, unit = 'Mx')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/sibry', this%sibry, unit = 'Mx')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/bcentr', this%bcentr, unit = 'G')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/current', this%current, unit = 'statA')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/fpol', this%fpol, &
         lbound(this%fpol), ubound(this%fpol), unit = 'G cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/pres', this%pres, &
         lbound(this%pres), ubound(this%pres), unit = 'dyn cm^-2')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ffprim', this%ffprim, &
         lbound(this%ffprim), ubound(this%ffprim), unit = 'G')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/pprime', this%pprime, &
         lbound(this%pprime), ubound(this%pprime), unit = 'dyn cm^-2 Mx^-1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/qpsi', this%qpsi, &
         lbound(this%qpsi), ubound(this%qpsi), unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rbbbs', this%rbbbs, &
         lbound(this%rbbbs), ubound(this%rbbbs), unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/zbbbs', this%zbbbs, &
         lbound(this%zbbbs), ubound(this%zbbbs), unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rlim', this%rlim, &
         lbound(this%rlim), ubound(this%rlim), unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/zlim', this%zlim, &
         lbound(this%zlim), ubound(this%zlim), unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi_eqd', this%psi_eqd, &
         lbound(this%psi_eqd), ubound(this%psi_eqd), unit = 'Mx', comment = 'equidistant psi grid values')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_eqd', this%R_eqd, &
         lbound(this%R_eqd), ubound(this%R_eqd), unit = 'cm', comment = 'equidistant R grid values')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_eqd', this%Z_eqd, &
         lbound(this%Z_eqd), ubound(this%Z_eqd), unit = 'cm', comment = 'equidistant Z grid values')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/fprime', this%fprime, &
         lbound(this%fprime), ubound(this%fprime), unit = 'cm^-1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psirz', this%psirz, &
         lbound(this%psirz), ubound(this%psirz), unit = 'Mx')
    call h5_close(h5id_root)
  end subroutine g_eqdsk_export_hdf5

  subroutine g_eqdsk_deinit(this)
    class(g_eqdsk) :: this

    this%fname = ''
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
    if (allocated(this%fprime)) deallocate(this%fprime)
    if (allocated(this%psirz)) deallocate(this%psirz)
  end subroutine g_eqdsk_deinit

  subroutine func1d_init(func1d, lb, ub)
    type(func1d_t), intent(inout) :: func1d
    integer, intent(in) :: lb
    integer, intent(in) :: ub

    call func1d_deinit(func1d)
    allocate(func1d%x(lb:ub))
    allocate(func1d%y(lb:ub))
  end subroutine func1d_init

  subroutine func1d_deinit(func1d)
    type(func1d_t), intent(inout) :: func1d

    if (allocated(func1d%x)) deallocate(func1d%x)
    if (allocated(func1d%y)) deallocate(func1d%y)
  end subroutine func1d_deinit

  subroutine func1d_write(func1d, file, group, x_comment, x_unit, y_comment, y_unit)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(func1d_t), intent(in) :: func1d
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: x_comment
    character(len = *), intent(in) :: x_unit
    character(len = *), intent(in) :: y_comment
    character(len = *), intent(in) :: y_unit
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/x', &
         func1d%x, lbound(func1d%x), ubound(func1d%x), &
         comment = trim(adjustl(x_comment)), unit = trim(adjustl(x_unit)))
    call h5_add(h5id_root, grp // '/y', &
         func1d%y, lbound(func1d%y), ubound(func1d%y), &
         comment = trim(adjustl(y_comment)), unit = trim(adjustl(y_unit)))
    call h5_close(h5id_root)
  end subroutine func1d_write

  subroutine func1d_read(func1d, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get_bounds, h5_get, h5_close
    type(func1d_t), intent(inout) :: func1d
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer :: lb, ub

    grp = trim(group)
    call h5_open(file, h5id_root)
    call h5_get_bounds(h5id_root, grp // '/x', lb, ub)
    call func1d_init(func1d, lb, ub)
    call h5_get(h5id_root, grp // '/x', func1d%x)
    call h5_get(h5id_root, grp // '/y', func1d%y)
    call h5_close(h5id_root)
  end subroutine func1d_read

  subroutine func1d_read_formatted(func1d, fname)
    type(func1d_t), intent(inout) :: func1d
    character(len = *), intent(in) :: fname
    integer :: fid, status, n, k

    open(newunit = fid, file = trim(adjustl(fname)), &
         status = 'old', form = 'formatted', action = 'read')
    n = 0
    do
       read(fid, *, iostat = status)
       if (status /= 0) exit
       n = n + 1
    end do
    rewind fid
    call func1d_init(func1d, 1, n)
    do k = 1, n
       read (fid, *) func1d%x(k), func1d%y(k)
    end do
    close(fid)
  end subroutine func1d_read_formatted

  !> 1D piecewise Lagrange polynomial interpolator
  function interp1d(sample_x, sample_y, resampled_x, n_lag, deriv) result(resampled_y)
    use mephit_conf, only: logger
    real(dp), intent(in) :: sample_x(:)
    real(dp), intent(in) :: sample_y(:)
    real(dp), intent(in) :: resampled_x
    real(dp) :: resampled_y
    integer, intent(in) :: n_lag
    logical, intent(in), optional :: deriv
    real(dp) :: lag_coeff(0:1, n_lag + 1)
    integer :: kder, k, k_lo, k_hi

    if (size(sample_x) /= size(sample_y)) then
       call logger%msg_arg_size('interp1d', 'size(sample_x)', &
            'size(sample_y)', size(sample_x), size(sample_y))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (n_lag >= size(sample_x)) then
       write (logger%msg, '("Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points ", i0)') &
            n_lag, size(sample_x)
       if (logger%err) call logger%write_msg
       error stop
    end if
    kder = 0
    if (present(deriv)) then
       if (deriv) then
          kder = 1
       end if
    end if
    call binsearch(sample_x, lbound(sample_x, 1), resampled_x, k)
    k_lo = k - (n_lag + 1) / 2
    k_hi = k_lo + n_lag
    ! ensure that polynomial sample points remain within the bounds of sample_x
    if (k_lo < lbound(sample_x, 1)) then
       k_lo = lbound(sample_x, 1)
       k_hi = k_lo + n_lag
    elseif (k_hi > ubound(sample_x, 1)) then
       k_hi = ubound(sample_x, 1)
       k_lo = k_hi - n_lag
    end if
    call plag_coeff(n_lag + 1, 1, resampled_x, sample_x(k_lo:k_hi), lag_coeff)
    resampled_y = sum(sample_y(k_lo:k_hi) * lag_coeff(kder, :))
  end function interp1d

  !> 1D piecewise Lagrange polynomial resampler
  subroutine resample1d(sample_x, sample_y, resampled_x, resampled_y, n_lag, deriv)
    use mephit_conf, only: logger
    real(dp), intent(in) :: sample_x(:)
    real(dp), intent(in) :: sample_y(:)
    real(dp), intent(in) :: resampled_x(:)
    real(dp), intent(out) :: resampled_y(:)
    integer, intent(in) :: n_lag
    logical, intent(in), optional :: deriv
    real(dp) :: lag_coeff(0:1, n_lag + 1)
    integer :: kder, k, k_lo, k_hi, k_re

    if (size(sample_x) /= size(sample_y)) then
       call logger%msg_arg_size('resample1d', 'size(sample_x)', &
            'size(sample_y)', size(sample_x), size(sample_y))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (size(resampled_x) /= size(resampled_y)) then
       call logger%msg_arg_size('resample1d', 'size(resampled_x)', &
            'size(resampled_y)', size(resampled_x), size(resampled_y))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (n_lag >= size(sample_x)) then
       write (logger%msg, '("Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points ", i0)') &
            n_lag, size(sample_x)
       if (logger%err) call logger%write_msg
       error stop
    end if
    kder = 0
    if (present(deriv)) then
       if (deriv) then
          kder = 1
       end if
    end if
    do k_re = lbound(resampled_x, 1), ubound(resampled_x, 1)
       call binsearch(sample_x, lbound(sample_x, 1), resampled_x(k_re), k)
       k_lo = k - (n_lag + 1) / 2
       k_hi = k_lo + n_lag
       ! ensure that polynomial sample points remain within the bounds of sample_x
       if (k_lo < lbound(sample_x, 1)) then
          k_lo = lbound(sample_x, 1)
          k_hi = k_lo + n_lag
       elseif (k_hi > ubound(sample_x, 1)) then
          k_hi = ubound(sample_x, 1)
          k_lo = k_hi - n_lag
       end if
       call plag_coeff(n_lag + 1, 1, resampled_x(k_re), sample_x(k_lo:k_hi), lag_coeff)
       resampled_y(k_re) = sum(sample_y(k_lo:k_hi) * lag_coeff(kder, :))
    end do
  end subroutine resample1d

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

  subroutine heapsort_real(array, ascending, permutation)
    use mephit_conf, only: logger
    real(dp), intent(inout) :: array(0:)
    logical, intent(in) :: ascending
    integer, intent(out), optional :: permutation(0:)
    integer :: n, k, child, root, itemp
    real(dp) :: temp

    n = size(array)
    if (present(permutation)) then
       if (n /= size(permutation)) then
          call logger%msg_arg_size('heapsort_real', 'size(array)', 'size(permutation)', &
               n, size(permutation))
          if (logger%err) call logger%write_msg
          error stop
       end if
       permutation = [(k, k = 1, n)]
    end if
    do k = (n - 2) / 2, 0, -1
       call siftdown(k, n)
    end do
    do k = n - 1, 1, -1
       temp = array(0)
       array(0) = array(k)
       array(k) = temp
       if (present(permutation)) then
          itemp = permutation(0)
          permutation(0) = permutation(k)
          permutation(k) = itemp
       end if
       call siftdown(0, k)
    end do
  contains
    subroutine siftdown(start, bottom)
      integer, intent(in) :: start, bottom
      root = start
      do while (root * 2 + 1 < bottom)
         child = root * 2 + 1
         if (child + 1 < bottom) then
            if ((ascending .and. array(child) < array(child+1)) .or. &
                 (.not. ascending .and. array(child) > array(child+1))) then
               child = child + 1
            end if
         end if
         if ((ascending .and. array(root) < array(child)) .or. &
              (.not. ascending .and. array(root) > array(child))) then
            temp = array(child)
            array(child) = array(root)
            array(root) = temp
            if (present(permutation)) then
               itemp = permutation(child)
               permutation(child) = permutation(root)
               permutation(root) = itemp
            end if
            root = child
         else
            return
         end if
      end do
    end subroutine siftdown
  end subroutine heapsort_real

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
       call siftdown(k, n)
    end do
    do k = n - 1, 1, -1
       temp = array(0)
       array(0) = array(k)
       array(k) = temp
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
            array(child) = array(root)
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

  function complex_abs_desc(val1, val2)
    complex(dp), intent(in) :: val1, val2
    logical :: complex_abs_desc
    complex_abs_desc = abs(val1) > abs(val2)
  end function complex_abs_desc

  !> Find eigenvalues above given magnitude using Arnoldi iterations.
  !>
  !> @ndim dimension of the linear operator
  !> @nkrylov maximum dimension of the Krylov subspace
  !> @nmin minimum dimension of the Krylov subspace
  !> @threshold find eigenvalues with absolute value abive this threshold
  !> @tol consider eigenvalues converged when this relative error is reached
  !> @next_iteration subroutine yielding matrix-vector product for given input vector
  !> @ierr error flag
  !> @nritz number of eigenvalues fulfilling \p threshold condition
  !> @eigvals sorted eigenvalues fulfilling \p threshold condition
  !> @eigvecs eigenvectors associated with \p eigvals
  !>
  !> If the requested eigenvalues are converged, \p ierr is set to 0. If convergence
  !> is not reached after \p nkrylov iterations, \p ierr is set to the number of
  !> converged eigenvalues and \p eigvals and \p eigvecs are not allocated. If the
  !> LAPACK subroutines fail, \p eigvecs and possibly \p eigvals are not allocated
  !> and the `info` parameter is propagated via \p ierr.
  subroutine arnoldi_break(ndim, nkrylov, nmin, threshold, tol, next_iteration, &
       ierr, nritz, eigvals, eigvecs)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    integer, intent(in) :: ndim, nkrylov, nmin
    real(dp), intent(in) :: threshold, tol
    interface
       subroutine next_iteration(old_val, new_val)
         import :: dp
         complex(dp), intent(in) :: old_val(:)
         complex(dp), intent(out) :: new_val(:)
       end subroutine next_iteration
    end interface
    integer, intent(out) :: ierr, nritz
    complex(dp), intent(out), allocatable :: eigvals(:)
    complex(dp), intent(out), allocatable, optional :: eigvecs(:, :)
    real(dp), parameter :: nearzero = 1d0 / huge(1d0)
    complex(dp), allocatable :: fold(:), fnew(:), fzero(:), &
         qvecs(:, :), hmat(:, :), ritzvals(:), ritzvecs(:, :), progression(:, :)
    integer :: j, k, n
    logical, allocatable :: selection(:), converged(:)

    ierr = 0
    nritz = 0
    if (ndim < 1) then
       ierr = -1
       print '("arnoldi_break: ndim = ", i0, " < 1")', ndim
       return
    end if
    if (nkrylov < 1) then
       ierr = -2
       print '("arnoldi_break: nkrylov = ", i0, " < 1")', nkrylov
       return
    end if
    if (nkrylov > ndim) then
       ierr = -2
       print '("arnoldi_break: nkrylov = ", i0, " > ndim = ", i0)', nkrylov, ndim
       return
    end if
    if (nmin < 1) then
       ierr = -3
       print '("arnoldi_break: nmin = ", i0, " < 1")', nmin
       return
    end if
    if (nmin > nkrylov) then
       ierr = -3
       print '("arnoldi_break: nmin = ", i0, " > nkrylov = ", i0)', nmin, nkrylov
       return
    end if
    allocate(fold(ndim), fnew(ndim), fzero(ndim), qvecs(ndim, nkrylov), hmat(nkrylov, nkrylov), &
         ritzvals(nkrylov), progression(nkrylov, nkrylov), selection(nkrylov), converged(nkrylov))
    ! initialize
    fold = (0d0, 0d0)
    print '("Iteration 1 of ", i0)', nkrylov
    call next_iteration(fold, fnew)
    fzero(:) = fnew
    qvecs(:, 1) = fnew / sqrt(sum(conjg(fnew) * fnew))
    hmat = (0d0, 0d0)
    ritzvals = (0d0, 0d0)
    selection = .false.
    converged = .false.
    progression = cmplx(ieee_value(0d0, ieee_quiet_nan), ieee_value(0d0, ieee_quiet_nan), dp)
    ! Arnoldi iterations
    n = nkrylov
    do k = 2, nkrylov
       print '("Iteration ", i0, " of ", i0)', k, nkrylov
       fold(:) = qvecs(:, k-1)
       call next_iteration(fold, fnew)
       qvecs(:, k) = fnew - fzero
       do j = 1, k-1
          hmat(j, k-1) = sum(conjg(qvecs(:, j)) * qvecs(:, k))
          qvecs(:, k) = qvecs(:, k) - hmat(j, k-1) * qvecs(:, j)
       end do
       hmat(k, k-1) = sqrt(sum(conjg(qvecs(:, k)) * qvecs(:, k)))
       if (abs(hmat(k, k-1)) < nearzero) then
          n = k
          exit
       end if
       qvecs(:, k) = qvecs(:, k) / hmat(k, k-1)
       ! calculate Ritz values
       call hessenberg_eigvals(hmat(:k, :k), ritzvals(:k), ierr)
       if (ierr /= 0) return
       progression(:k, k) = ritzvals(:k)
       selection(:k) = abs(ritzvals(:k)) >= threshold
       nritz = count(selection)
       converged(:k) = abs(progression(:k, k) - progression(:k, k - 1)) / &
            abs(progression(:k, k - 1)) < tol .and. k >= nmin
       if (all(pack(converged, selection))) then
          n = k
          exit
       end if
    end do
    if (.not. all(pack(converged, selection))) then
       ierr = count(converged)
       print '("arnoldi_break: only ", i0, " eigenvalues of ", i0, " converged")', &
            ierr, nritz
       return
    end if
    ! sort eigenvalues in descending order and optionall compute eigenvectors in this order
    call heapsort_complex(ritzvals(:n), complex_abs_desc)
    allocate(eigvals(nritz))
    eigvals(:) = pack(ritzvals, selection)
    if (present(eigvecs)) then
       allocate(ritzvecs(n, nritz), eigvecs(ndim, nritz))
       call hessenberg_eigvecs(hmat(:n, :n), ritzvals(:n), selection(:n), ritzvecs, ierr)
       if (ierr /= 0) return
       eigvecs = matmul(qvecs(:, :n), ritzvecs)
       deallocate(ritzvecs)
    end if
    deallocate(fold, fnew, fzero, qvecs, hmat, ritzvals, progression, selection, converged)
  end subroutine arnoldi_break

  !> Compute eigenvalues of square Hessenberg matrix.
  !>
  !> @hmat square Hessenberg matrix
  !> @eigvals calculated eigenvalues
  !> @ierr error flag
  !>
  !> The `info` parameter from ZHSEQR is propagated via \p ierr on failure.
  !> Note that no check is performed to verify that \p hmat is indeed a
  !> Hessenberg matrix.
  subroutine hessenberg_eigvals(hmat, eigvals, ierr)
    complex(dp), intent(in) :: hmat(:, :)
    complex(dp), intent(out) :: eigvals(:)
    integer, intent(out) :: ierr
    integer :: ndim, lwork
    complex(dp), allocatable :: hmat_work(:, :), zdum(:, :), work(:)

    ndim = size(hmat, 2)
    if (size(hmat, 1) /= ndim) then
       print '("hessenberg_eigvals: hmat has shape (", i0, ", ", i0, "), ' // &
            'but expected (", i0, ", ", i0, ")")', shape(hmat), ndim, ndim
       error stop
    end if
    if (size(eigvals, 1) /= ndim) then
       print '("hessenberg_eigvals: eigvals has shape (", i0, "), ' // &
            'but expected (", i0, ")")', shape(hmat), ndim
       error stop
    end if
    allocate(hmat_work(ndim, ndim), zdum(1, ndim))
    hmat_work(:, :) = hmat
    allocate(work(1))
    lwork = -1
    call zhseqr('E', 'N', ndim, 1, ndim, hmat_work, ndim, eigvals, &
         zdum, 1, work, lwork, ierr)
    if (ierr /= 0) then
       if (ierr < 0) then
          print '("ZHSEQR: illegal value in argument #", i0)', -ierr
       else
          print '("ZHSEQR: only ", i0, " of ", i0, " eigenvalues converged")', &
               ierr, ndim
       end if
       return
    end if
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call zhseqr('E', 'N', ndim, 1, ndim, hmat_work, ndim, eigvals, &
         zdum, 1, work, lwork, ierr)
    deallocate(work, hmat_work, zdum)
    if (ierr /= 0) then
       if (ierr < 0) then
          print '("ZHSEQR: illegal value in argument #", i0)', -ierr
       else
          print '("ZHSEQR: only ", i0, " of ", i0, " eigenvalues converged")', &
               ierr, ndim
       end if
    end if
  end subroutine hessenberg_eigvals

  !> Compute right eigenvectors of square Hessenberg matrix.
  !>
  !> @hmat square Hessenberg matrix
  !> @eigvals calculated eigenvalues
  !> @mask selection of eigenvalues for which to compute eigenvectors
  !> @eigvecs calculated eigenvectors; second dimension is given by count(mask)
  !> @ierr error flag
  !>
  !> The `info` parameter from ZHSEIN is propagated via \p ierr on failure.
  !> Note that no check is performed to verify that \p hmat is indeed a
  !> Hessenberg matrix. \p eigvals are assumed to be calculated by ZHSEQR
  !> via hessenberg_eigvals().
  subroutine hessenberg_eigvecs(hmat, eigvals, mask, eigvecs, ierr)
    complex(dp), intent(in) :: hmat(:, :)
    complex(dp), intent(in) :: eigvals(:)
    logical, intent(in) :: mask(:)
    complex(dp), intent(out) :: eigvecs(:, :)
    integer, intent(out) :: ierr
    integer :: ndim, neff
    integer, allocatable :: idum(:), ifailr(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: eigvals_work(:), zdum(:, :), work(:, :)

    ndim = size(hmat, 2)
    if (size(hmat, 1) /= ndim) then
       print '("hessenberg_eigvecs: hmat has shape (", i0, ", ", i0, "), ' // &
            'but expected (", i0, ", ", i0, ")")', shape(hmat), ndim, ndim
       error stop
    end if
    if (size(eigvals, 1) /= ndim) then
       print '("hessenberg_eigvecs: eigvals has shape (", i0, "), ' // &
            'but expected (", i0, ")")', shape(eigvals), ndim
       error stop
    end if
    if (size(mask, 1) /= ndim) then
       print '("hessenberg_eigvecs: mask has shape (", i0, "), ' // &
            'but expected (", i0, ")")', shape(mask), ndim
       error stop
    end if
    if (size(eigvecs, 1) /= ndim .or. size(eigvecs, 2) /= count(mask)) then
       print '("hessenberg_eigvecs: eigvecs has shape (", i0, ", ", i0, "), ' // &
            'but expected (", i0, ", ", i0, ")")', shape(eigvecs), ndim, count(mask)
       error stop
    end if
    allocate(eigvals_work(ndim), work(ndim, ndim), rwork(ndim), &
         zdum(1, count(mask)), idum(count(mask)), ifailr(count(mask)))
    eigvals_work(:) = eigvals
    call zhsein('R', 'Q', 'N', mask, ndim, hmat, ndim, eigvals_work, &
         zdum, 1, eigvecs, ndim, count(mask), neff, work, rwork, &
         idum, ifailr, ierr)
    if (ierr /= 0) then
       if (ierr < 0) then
          print '("ZHSEIN: illegal value in argument #", i0)', -ierr
       else
          print '("ZHSEIN: only ", i0, " of ", i0, " eigenvectors converged")', &
               ierr, ndim
       end if
    end if
    deallocate(eigvals_work, work, rwork, zdum, idum, ifailr)
  end subroutine hessenberg_eigvecs

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

  subroutine C_F_string(C, F)
    use iso_c_binding, only: c_ptr, c_char, c_associated, c_f_pointer, c_null_char
    type(c_ptr), intent(in), value :: C
    character(len = *), intent(out) :: F
    character(kind = c_char, len = 1), dimension(:), pointer :: p
    integer :: k

    if (c_associated(C)) then
       call c_f_pointer(C, p, [huge(0)])
       k = 1
       do while (p(k) /= c_null_char .and. k <= len(F))
          F(k:k) = p(k)
          k = k + 1
       end do
       if (k < len(f)) then
          F(k:) = ' '
       end if
    else
       F = ''
    end if
  end subroutine C_F_string

end module mephit_util
