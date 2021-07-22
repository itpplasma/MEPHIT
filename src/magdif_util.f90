module magdif_util

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: pi, clight, imun, get_field_filenames, init_field, deinit_field, interp_psi_pol, &
       pos_angle, linspace, straight_cyl2bent_cyl, bent_cyl2straight_cyl, binsearch, interleave, &
       gauss_legendre_unit_interval, heapsort_complex, complex_abs_asc, C_F_string

  real(dp), parameter :: pi = 4d0 * atan(1d0)
  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

  type, public :: sign_convention
    integer :: exp_Bpol, sgn_cyl, sgn_dpsi, sgn_Btor, sgn_Itor, &
         sgn_F, sgn_q, sgn_Bpol, sgn_pol, index
  end type sign_convention

  type, public :: g_eqdsk
    type(sign_convention) :: cocos
    character(len = 1024) :: fname
    character(len = 1024) :: convexfile
    character(len = 48) :: header
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

  !> 1D piecewise Lagrange polynomial interpolator
  type, public :: interp1d
     private
     integer :: n_lag, n_var
     real(dp), dimension(:), allocatable :: indep_var
   contains
     procedure :: init => interp1d_init
     procedure :: eval => interp1d_eval
     procedure :: deinit => interp1d_deinit
  end type interp1d

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
  subroutine get_field_filenames(gfile, pfile, convexfile)
    character(len = *), intent(out) :: gfile, pfile, convexfile
    integer :: fid
    open(newunit = fid, file = 'field_divB0.inp', status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *)
    read (fid, *) gfile        ! equilibrium file
    read (fid, *) pfile        ! coil file
    read (fid, *) convexfile   ! convex file for stretchcoords
    close(fid)
  end subroutine get_field_filenames

  !> Set module variables for initialization of subroutine field
  subroutine init_field(equil)
    use input_files, only: convexfile
    use field_mod, only: icall, ipert, iequil
    use field_eq_mod, only: skip_read, icall_eq, nwindow_r, nwindow_z, &
         nrad, nzet, psi_axis, psi_sep, btf, rtf, splfpol, rad, zet, psi, psi0
    type(g_eqdsk), intent(in) :: equil
    real(dp) :: dum

    ! use supplied filename
    convexfile = equil%convexfile
    ! compute equiibrium field
    ipert = 0
    iequil = 1
    ! default values - TODO: options in magdif_conf
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

  subroutine g_eqdsk_read(this, fname, convexfile)
    class(g_eqdsk), intent(inout) :: this
    character(len = *), intent(in) :: fname, convexfile
    integer :: fid, kw, kh, idum
    real(dp) :: xdum

    call g_eqdsk_deinit(this)
    this%fname = fname
    this%convexfile = convexfile
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
  end subroutine g_eqdsk_read

  subroutine g_eqdsk_check_consistency(this)
    use constants, only: pi  ! src/orbit_mod.f90
    use magdif_conf, only: log
    class(g_eqdsk), intent(inout) :: this
    integer, parameter :: ignore = 3
    type(interp1d) :: psi_interpolator
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
    deriv_eqd(:) = [(psi_interpolator%eval(this%pres, this%psi_eqd(k), .true.), &
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
    deriv_eqd(:) = [(psi_interpolator%eval(this%fpol, this%psi_eqd(k), .true.), &
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
    type(interp1d) :: psi_interpolator
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
    gs_rhs = -4d0 * pi * this%R_eqd(kw) ** 2 * psi_interpolator%eval(this%pprime, psi) &
         - psi_interpolator%eval(this%ffprim, psi)
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
       this%psirz(:, :) = -this%psirz
       this%psi_eqd(:) = -this%psi_eqd
       this%cocos%sgn_dpsi = -this%cocos%sgn_dpsi
       this%pprime(:) = -this%pprime
       this%ffprim(:) = -this%ffprim
       this%cocos%sgn_Bpol = -this%cocos%sgn_Bpol
    end if
    if (this%cocos%sgn_pol == +1) then
       write (log%msg, invert_fmt) 'QPSI'
       if (log%info) call log%write
       this%qpsi(:) = -this%qpsi
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
       this%psirz(:, :) = this%psirz / (2d0 * pi)
       this%psi_eqd(:) = this%psi_eqd / (2d0 * pi)
       this%pprime(:) = this%pprime * (2d0 * pi)
       this%ffprim(:) = this%ffprim * (2d0 * pi)
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/convexfile', this%convexfile)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/header', this%header)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nw', this%nw)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nh', this%nh)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nbbbs', this%nbbbs)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/limitr', this%limitr)
    allocate(this%fpol(this%nw), this%pres(this%nw), this%ffprim(this%nw), this%pprime(this%nw), &
         this%qpsi(this%nw), this%rbbbs(this%nbbbs), this%zbbbs(this%nbbbs), &
         this%rlim(this%limitr), this%zlim(this%limitr), this%psi_eqd(this%nw), &
         this%R_eqd(this%nw), this%Z_eqd(this%nh), this%psirz(this%nw, this%nh))
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/convexfile', this%convexfile, &
         comment = 'associated convexfile for subroutine stretch_coords')
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
         lbound(this%ffprim), ubound(this%ffprim), unit = 'cm^-1')
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psirz', this%psirz, &
         lbound(this%psirz), ubound(this%psirz), unit = 'Mx')
    call h5_close(h5id_root)
  end subroutine g_eqdsk_export_hdf5

  subroutine g_eqdsk_deinit(this)
    class(g_eqdsk) :: this

    this%fname = ''
    this%convexfile = ''
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
  end subroutine g_eqdsk_deinit

  subroutine interp1d_init(this, n_lag, indep_var)
    use magdif_conf, only: log
    class(interp1d), intent(inout) :: this
    integer, intent(in) :: n_lag
    real(dp), intent(in), dimension(:) :: indep_var

    if (n_lag >= size(indep_var)) then
       write (log%msg, '("Lagrange polynomial order ", i0, ' // &
            '" must be lower than number of sample points", i0)') n_lag, size(indep_var)
       if (log%err) call log%write
       return
    end if
    call interp1d_deinit(this)
    this%n_lag = n_lag
    this%n_var = size(indep_var)
    allocate(this%indep_var(this%n_var))
    this%indep_var(:) = indep_var
  end subroutine interp1d_init

  function interp1d_eval(this, sample, position, deriv) result(interp)
    use magdif_conf, only: log
    class(interp1d) :: this
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
       call log%msg_arg_size('interp1d_eval', 'this%n_var', 'size(sample)', &
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
  end function interp1d_eval

  subroutine interp1d_deinit(this)
    class(interp1d), intent(inout) :: this

    if (allocated(this%indep_var)) deallocate(this%indep_var)
  end subroutine interp1d_deinit

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

end module magdif_util
