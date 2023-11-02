module mephit_util

  use iso_fortran_env, only: dp => real64
  use iso_c_binding, only: c_ptr, c_null_ptr, c_double_complex

  implicit none

  private

  ! types and associated procedures
  public :: fft_t
  public :: func1d_t, func1d_init, func1d_deinit, func1d_write, func1d_read, func1d_read_formatted
  public :: neumaier_accumulator_real, neumaier_accumulator_complex

  ! utility procedures
  public :: init_field, deinit_field, generate_symfluxcoord, load_symfluxcoord, save_symfluxcoord, &
    geqdsk_scale, geqdsk_export_hdf5, geqdsk_import_hdf5, &
    straight_cyl2bent_cyl, bent_cyl2straight_cyl, zd_cross, dd_cross, &
    interp_psi_pol, interp1d, resample1d, pos_angle, linspace, &
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
    use geqdsk_tools, only: geqdsk_t
    type(geqdsk_t), intent(in) :: equil
    real(dp) :: dum

    ! compute equilibrium field
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

  !> calculates points on a fine grid in the core region by integrating along field lines
  !> and spline interpolates them for use with magdata_in_symfluxcoord_ext
  subroutine generate_symfluxcoord
    use magdata_in_symfluxcoor_mod, only: unload_magdata_in_symfluxcoord, load, &
      nspl, nlabel, ntheta, twopi, h_theta, &
      rmn, rmx, zmn, zmx, raxis, zaxis, h_theta, psipol_max, psitor_max, &
      rbeg, rsmall, qsaf, psisurf, phitor, circumf, R_st, Z_st, bmod_st, sqgnorm_st
    real(dp), dimension(:, :), allocatable :: R_transp, Z_transp, bmod_transp, sqgnorm_transp
    integer :: fid, nstep, nsurfmax, k

    open(newunit = fid, file = 'preload_for_SYNCH.inp', status = 'old', form = 'formatted', action = 'read')
    ! number of integration steps
    read (fid, *) nstep
    ! grid size over radial variable
    read (fid, *) nlabel
    ! grid size over poloidal angle
    read (fid, *) ntheta
    ! number of starting points between the magnetic axis and right box boundary when searching for the separatrix
    read (fid, *) nsurfmax
    close(fid)
    call unload_magdata_in_symfluxcoord
    allocate(rbeg(nlabel), rsmall(nlabel), qsaf(nlabel), psisurf(0:nlabel), phitor(0:nlabel), circumf(nlabel))
    allocate(R_transp(nlabel, ntheta), Z_transp(nlabel, ntheta), bmod_transp(nlabel, ntheta), sqgnorm_transp(nlabel, ntheta))
    call field_line_integration_for_SYNCH(nstep, nsurfmax, nlabel, ntheta, &
      rmn, rmx, zmn, zmx, raxis, zaxis, rbeg, rsmall, qsaf, psisurf(1:), phitor(1:), circumf, &
      R_transp, Z_transp, bmod_transp, sqgnorm_transp)
    psisurf(0) = 0.0d0
    phitor(0) = 0.0d0
    psipol_max = psisurf(nlabel)
    psitor_max = phitor(nlabel)
    psisurf = psisurf / psipol_max
    phitor = phitor / psitor_max
    h_theta = twopi / ntheta
    allocate(R_st(0:nspl, 0:ntheta, nlabel))
    allocate(Z_st(0:nspl, 0:ntheta, nlabel))
    allocate(bmod_st(0:nspl, 0:ntheta, nlabel))
    allocate(sqgnorm_st(0:nspl, 0:ntheta, nlabel))
    R_st(0, 1:, :) = transpose(R_transp)
    R_st(0, 0, :) = R_st(0, ntheta, :)
    do k = 1, nlabel
      call spl_per(nspl, ntheta + 1, h_theta, R_st(:, :, k))
    end do
    Z_st(0, 1:, :) = transpose(Z_transp)
    Z_st(0, 0, :) = Z_st(0, ntheta, :)
    do k = 1, nlabel
      call spl_per(nspl, ntheta + 1, h_theta, Z_st(:, :, k))
    end do
    bmod_st(0, 1:, :) = transpose(bmod_transp)
    bmod_st(0, 0, :) = bmod_st(0, ntheta, :)
    do k = 1, nlabel
      call spl_per(nspl, ntheta + 1, h_theta, bmod_st(:, :, k))
    end do
    sqgnorm_st(0, 1:, :) = transpose(sqgnorm_transp)
    sqgnorm_st(0, 0, :) = sqgnorm_st(0, ntheta, :)
    do k = 1, nlabel
      call spl_per(nspl, ntheta + 1, h_theta, sqgnorm_st(:, :, k))
    end do
    load = .false.  ! just in case it is used externally
    deallocate(R_transp, Z_transp, bmod_transp, sqgnorm_transp)
  end subroutine generate_symfluxcoord

  subroutine save_symfluxcoord(file, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdata_in_symfluxcoor_mod, only: nlabel, ntheta, h_theta, &
      rmn, rmx, zmn, zmx, raxis, zaxis, h_theta, psipol_max, psitor_max, &
      rbeg, rsmall, qsaf, psisurf, phitor, circumf, R_st, Z_st, bmod_st, sqgnorm_st
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/nlabel', nlabel, comment = 'number of flux surfaces')
    call h5_add(h5id_root, grp // '/ntheta', ntheta, comment = 'number of poloidal divisions')
    call h5_add(h5id_root, grp // '/rmn', rmn, unit = 'cm', &
      comment = 'minimum R of computational domain')
    call h5_add(h5id_root, grp // '/rmx', rmx, unit = 'cm', &
      comment = 'maximum R of computational domain')
    call h5_add(h5id_root, grp // '/zmn', zmn, unit = 'cm', &
      comment = 'minimum Z of computational domain')
    call h5_add(h5id_root, grp // '/zmx', zmx, unit = 'cm', &
      comment = 'maximum Z of computational domain')
    call h5_add(h5id_root, grp // '/raxis', raxis, unit = 'cm', &
      comment = 'R coordinate of magnetic axis')
    call h5_add(h5id_root, grp // '/zaxis', zaxis, unit = 'cm', &
      comment = 'Z coordinate of magnetic axis')
    call h5_add(h5id_root, grp // '/h_theta', h_theta, unit = 'rad', &
      comment = 'angle between poloidal divisions')
    call h5_add(h5id_root, grp // '/psipol_max', psipol_max, unit = 'Mx rad^-1', &
      comment = 'poloidal flux over 2 pi at separatrix')
    call h5_add(h5id_root, grp // '/psitor_max', psitor_max, unit = 'Mx rad^-1', &
      comment = 'toroidal flux over 2 pi at separatrix')
    call h5_add(h5id_root, grp // '/rbeg', rbeg, lbound(rbeg), ubound(rbeg), &
      unit = 'cm', comment = 'small radius (outboard from O point, or towards X point)')
    call h5_add(h5id_root, grp // '/rsmall', rsmall, lbound(rsmall), ubound(rsmall), &
      unit = 'cm', comment = 'equivalent radius of poloidal cross-section area')
    call h5_add(h5id_root, grp // '/qsaf', qsaf, lbound(qsaf), ubound(qsaf), &
      unit = '1', comment = 'safety factor')
    call h5_add(h5id_root, grp // '/psisurf', psisurf, lbound(psisurf), ubound(psisurf), &
      unit = '1', comment = 'normalized poloidal flux')
    call h5_add(h5id_root, grp // '/phitor', phitor, lbound(phitor), ubound(phitor), &
      unit = '1', comment = 'normalized toroidal flux')
    call h5_add(h5id_root, grp // '/circumf', circumf, lbound(circumf), ubound(circumf), &
      unit = 'cm', comment = 'poloidal cross-section circumference')
    call h5_add(h5id_root, grp // '/R_st', R_st, lbound(R_st), ubound(R_st), &
      comment = 'spline of R in s and theta')
    call h5_add(h5id_root, grp // '/Z_st', Z_st, lbound(Z_st), ubound(Z_st), &
      comment = 'spline of Z in s and theta')
    call h5_add(h5id_root, grp // '/bmod_st', bmod_st, lbound(bmod_st), ubound(bmod_st), &
      comment = 'spline of magnetic field modulus in s and theta')
    call h5_add(h5id_root, grp // '/sqgnorm_st', sqgnorm_st, lbound(sqgnorm_st), ubound(sqgnorm_st), &
      comment = 'spline of the (s, theta, phi) Jacobian in s and theta')
    call h5_close(h5id_root)
  end subroutine save_symfluxcoord

  !> loads points that are calculated in preload_for_SYNCH into module variables
  subroutine load_symfluxcoord(file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use magdata_in_symfluxcoor_mod, only: unload_magdata_in_symfluxcoord, load, &
      nlabel, ntheta, nspl, twopi, h_theta, &
      rmn, rmx, zmn, zmx, raxis, zaxis, h_theta, psipol_max, psitor_max, &
      rbeg, rsmall, qsaf, psisurf, phitor, circumf, R_st, Z_st, bmod_st, sqgnorm_st
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call unload_magdata_in_symfluxcoord
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/nlabel', nlabel)
    call h5_get(h5id_root, grp // '/ntheta', ntheta)
    call h5_get(h5id_root, grp // '/rmn', rmn)
    call h5_get(h5id_root, grp // '/rmx', rmx)
    call h5_get(h5id_root, grp // '/zmn', zmn)
    call h5_get(h5id_root, grp // '/zmx', zmx)
    call h5_get(h5id_root, grp // '/raxis', raxis)
    call h5_get(h5id_root, grp // '/zaxis', zaxis)
    call h5_get(h5id_root, grp // '/h_theta', h_theta)
    call h5_get(h5id_root, grp // '/psipol_max', psipol_max)
    call h5_get(h5id_root, grp // '/psitor_max', psitor_max)
    allocate(rbeg(nlabel), rsmall(nlabel), qsaf(nlabel), psisurf(0:nlabel), phitor(0:nlabel), circumf(nlabel))
    call h5_get(h5id_root, grp // '/rbeg', rbeg)
    call h5_get(h5id_root, grp // '/rsmall', rsmall)
    call h5_get(h5id_root, grp // '/qsaf', qsaf)
    call h5_get(h5id_root, grp // '/psisurf', psisurf)
    call h5_get(h5id_root, grp // '/phitor', phitor)
    call h5_get(h5id_root, grp // '/circumf', circumf)
    allocate(R_st(0:nspl, 0:ntheta, nlabel))
    allocate(Z_st(0:nspl, 0:ntheta, nlabel))
    allocate(bmod_st(0:nspl, 0:ntheta, nlabel))
    allocate(sqgnorm_st(0:nspl, 0:ntheta, nlabel))
    call h5_get(h5id_root, grp // '/R_st', R_st)
    call h5_get(h5id_root, grp // '/Z_st', Z_st)
    call h5_get(h5id_root, grp // '/bmod_st', bmod_st)
    call h5_get(h5id_root, grp // '/sqgnorm_st', sqgnorm_st)
    call h5_close(h5id_root)
    h_theta = twopi / ntheta
    load = .false.  ! just in case it is used externally
  end subroutine load_symfluxcoord

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

  subroutine geqdsk_scale(geqdsk, alpha)
    use geqdsk_tools, only: geqdsk_t
    type(geqdsk_t), intent(inout) :: geqdsk
    integer, intent(in) :: alpha
    real(dp) :: r_shift

    r_shift = geqdsk%rcentr * (alpha - 1)
    geqdsk%rcentr = geqdsk%rcentr * alpha
    geqdsk%rleft = geqdsk%rleft + r_shift
    geqdsk%rmaxis = geqdsk%rmaxis + r_shift
    geqdsk%simag = geqdsk%simag * alpha
    geqdsk%sibry = geqdsk%sibry * alpha
    geqdsk%fpol(:) = geqdsk%fpol * alpha
    geqdsk%ffprim(:) = geqdsk%ffprim * alpha
    geqdsk%pprime(:) = geqdsk%pprime / alpha
    geqdsk%psirz(:, :) = geqdsk%psirz * alpha
    geqdsk%qpsi(:) = geqdsk%qpsi / alpha
    geqdsk%rbbbs(:) = geqdsk%rbbbs + r_shift
    geqdsk%rlim(:) = geqdsk%rlim + r_shift
    ! auxiliary values
    geqdsk%psi_eqd(:) = geqdsk%psi_eqd * alpha
    geqdsk%R_eqd(:) = geqdsk%R_eqd + r_shift
  end subroutine geqdsk_scale

  subroutine geqdsk_import_hdf5(geqdsk, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use geqdsk_tools, only: geqdsk_t, geqdsk_deinit
    type(geqdsk_t), intent(inout) :: geqdsk
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call geqdsk_deinit(geqdsk)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/cocos/exp_Bpol', geqdsk%cocos%exp_Bpol)
    call h5_get(h5id_root, grp // '/cocos/sgn_cyl', geqdsk%cocos%sgn_cyl)
    call h5_get(h5id_root, grp // '/cocos/sgn_dpsi', geqdsk%cocos%sgn_dpsi)
    call h5_get(h5id_root, grp // '/cocos/sgn_Btor', geqdsk%cocos%sgn_Btor)
    call h5_get(h5id_root, grp // '/cocos/sgn_Itor', geqdsk%cocos%sgn_Itor)
    call h5_get(h5id_root, grp // '/cocos/sgn_F', geqdsk%cocos%sgn_F)
    call h5_get(h5id_root, grp // '/cocos/sgn_q', geqdsk%cocos%sgn_q)
    call h5_get(h5id_root, grp // '/cocos/sgn_Bpol', geqdsk%cocos%sgn_Bpol)
    call h5_get(h5id_root, grp // '/cocos/sgn_pol', geqdsk%cocos%sgn_pol)
    call h5_get(h5id_root, grp // '/cocos/index', geqdsk%cocos%index)
    call h5_get(h5id_root, grp // '/fname', geqdsk%fname)
    call h5_get(h5id_root, grp // '/header', geqdsk%header)
    call h5_get(h5id_root, grp // '/nw', geqdsk%nw)
    call h5_get(h5id_root, grp // '/nh', geqdsk%nh)
    call h5_get(h5id_root, grp // '/nbbbs', geqdsk%nbbbs)
    call h5_get(h5id_root, grp // '/limitr', geqdsk%limitr)
    allocate(geqdsk%fpol(geqdsk%nw))
    allocate(geqdsk%pres(geqdsk%nw))
    allocate(geqdsk%ffprim(geqdsk%nw))
    allocate(geqdsk%pprime(geqdsk%nw))
    allocate(geqdsk%psirz(geqdsk%nw, geqdsk%nh))
    allocate(geqdsk%qpsi(geqdsk%nw))
    allocate(geqdsk%rbbbs(geqdsk%nbbbs))
    allocate(geqdsk%zbbbs(geqdsk%nbbbs))
    allocate(geqdsk%rlim(geqdsk%limitr))
    allocate(geqdsk%zlim(geqdsk%limitr))
    allocate(geqdsk%psi_eqd(geqdsk%nw))
    allocate(geqdsk%R_eqd(geqdsk%nw))
    allocate(geqdsk%Z_eqd(geqdsk%nh))
    allocate(geqdsk%fprime(geqdsk%nw))
    call h5_get(h5id_root, grp // '/rdim', geqdsk%rdim)
    call h5_get(h5id_root, grp // '/zdim', geqdsk%zdim)
    call h5_get(h5id_root, grp // '/rcentr', geqdsk%rcentr)
    call h5_get(h5id_root, grp // '/rleft', geqdsk%rleft)
    call h5_get(h5id_root, grp // '/zmid', geqdsk%zmid)
    call h5_get(h5id_root, grp // '/rmaxis', geqdsk%rmaxis)
    call h5_get(h5id_root, grp // '/zmaxis', geqdsk%zmaxis)
    call h5_get(h5id_root, grp // '/simag', geqdsk%simag)
    call h5_get(h5id_root, grp // '/sibry', geqdsk%sibry)
    call h5_get(h5id_root, grp // '/bcentr', geqdsk%bcentr)
    call h5_get(h5id_root, grp // '/current', geqdsk%current)
    call h5_get(h5id_root, grp // '/fpol', geqdsk%fpol)
    call h5_get(h5id_root, grp // '/pres', geqdsk%pres)
    call h5_get(h5id_root, grp // '/ffprim', geqdsk%ffprim)
    call h5_get(h5id_root, grp // '/pprime', geqdsk%pprime)
    call h5_get(h5id_root, grp // '/psirz', geqdsk%psirz)
    call h5_get(h5id_root, grp // '/qpsi', geqdsk%qpsi)
    call h5_get(h5id_root, grp // '/rbbbs', geqdsk%rbbbs)
    call h5_get(h5id_root, grp // '/zbbbs', geqdsk%zbbbs)
    call h5_get(h5id_root, grp // '/rlim', geqdsk%rlim)
    call h5_get(h5id_root, grp // '/zlim', geqdsk%zlim)
    call h5_get(h5id_root, grp // '/psi_eqd', geqdsk%psi_eqd)
    call h5_get(h5id_root, grp // '/R_eqd', geqdsk%R_eqd)
    call h5_get(h5id_root, grp // '/Z_eqd', geqdsk%Z_eqd)
    call h5_get(h5id_root, grp // '/fprime', geqdsk%fprime)
    call h5_close(h5id_root)
  end subroutine geqdsk_import_hdf5

  subroutine geqdsk_export_hdf5(geqdsk, file, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use geqdsk_tools, only: geqdsk_t
    type(geqdsk_t), intent(in) :: geqdsk
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/cocos/')
    call h5_add(h5id_root, grp // '/cocos/exp_Bpol', geqdsk%cocos%exp_Bpol)
    call h5_add(h5id_root, grp // '/cocos/sgn_cyl', geqdsk%cocos%sgn_cyl)
    call h5_add(h5id_root, grp // '/cocos/sgn_dpsi', geqdsk%cocos%sgn_dpsi)
    call h5_add(h5id_root, grp // '/cocos/sgn_Btor', geqdsk%cocos%sgn_Btor)
    call h5_add(h5id_root, grp // '/cocos/sgn_Itor', geqdsk%cocos%sgn_Itor)
    call h5_add(h5id_root, grp // '/cocos/sgn_F', geqdsk%cocos%sgn_F)
    call h5_add(h5id_root, grp // '/cocos/sgn_q', geqdsk%cocos%sgn_q)
    call h5_add(h5id_root, grp // '/cocos/sgn_Bpol', geqdsk%cocos%sgn_Bpol)
    call h5_add(h5id_root, grp // '/cocos/sgn_pol', geqdsk%cocos%sgn_pol)
    call h5_add(h5id_root, grp // '/cocos/index', geqdsk%cocos%index)
    call h5_add(h5id_root, grp // '/fname', geqdsk%fname, &
      comment = 'original GEQDSK filename')
    call h5_add(h5id_root, grp // '/header', geqdsk%header)
    call h5_add(h5id_root, grp // '/nw', geqdsk%nw)
    call h5_add(h5id_root, grp // '/nh', geqdsk%nh)
    call h5_add(h5id_root, grp // '/nbbbs', geqdsk%nbbbs)
    call h5_add(h5id_root, grp // '/limitr', geqdsk%limitr)
    call h5_add(h5id_root, grp // '/rdim', geqdsk%rdim, unit = 'cm')
    call h5_add(h5id_root, grp // '/zdim', geqdsk%zdim, unit = 'cm')
    call h5_add(h5id_root, grp // '/rcentr', geqdsk%rcentr, unit = 'cm')
    call h5_add(h5id_root, grp // '/rleft', geqdsk%rleft, unit = 'cm')
    call h5_add(h5id_root, grp // '/zmid', geqdsk%zmid, unit = 'cm')
    call h5_add(h5id_root, grp // '/rmaxis', geqdsk%rmaxis, unit = 'cm')
    call h5_add(h5id_root, grp // '/zmaxis', geqdsk%zmaxis, unit = 'cm')
    call h5_add(h5id_root, grp // '/simag', geqdsk%simag, unit = 'Mx')
    call h5_add(h5id_root, grp // '/sibry', geqdsk%sibry, unit = 'Mx')
    call h5_add(h5id_root, grp // '/bcentr', geqdsk%bcentr, unit = 'G')
    call h5_add(h5id_root, grp // '/current', geqdsk%current, unit = 'statA')
    call h5_add(h5id_root, grp // '/fpol', geqdsk%fpol, &
      lbound(geqdsk%fpol), ubound(geqdsk%fpol), unit = 'G cm')
    call h5_add(h5id_root, grp // '/pres', geqdsk%pres, &
      lbound(geqdsk%pres), ubound(geqdsk%pres), unit = 'dyn cm^-2')
    call h5_add(h5id_root, grp // '/ffprim', geqdsk%ffprim, &
      lbound(geqdsk%ffprim), ubound(geqdsk%ffprim), unit = 'G')
    call h5_add(h5id_root, grp // '/pprime', geqdsk%pprime, &
      lbound(geqdsk%pprime), ubound(geqdsk%pprime), unit = 'dyn cm^-2 Mx^-1')
    call h5_add(h5id_root, grp // '/psirz', geqdsk%psirz, &
      lbound(geqdsk%psirz), ubound(geqdsk%psirz), unit = 'Mx')
    call h5_add(h5id_root, grp // '/qpsi', geqdsk%qpsi, &
      lbound(geqdsk%qpsi), ubound(geqdsk%qpsi), unit = '1')
    call h5_add(h5id_root, grp // '/rbbbs', geqdsk%rbbbs, &
      lbound(geqdsk%rbbbs), ubound(geqdsk%rbbbs), unit = 'cm')
    call h5_add(h5id_root, grp // '/zbbbs', geqdsk%zbbbs, &
      lbound(geqdsk%zbbbs), ubound(geqdsk%zbbbs), unit = 'cm')
    call h5_add(h5id_root, grp // '/rlim', geqdsk%rlim, &
      lbound(geqdsk%rlim), ubound(geqdsk%rlim), unit = 'cm')
    call h5_add(h5id_root, grp // '/zlim', geqdsk%zlim, &
      lbound(geqdsk%zlim), ubound(geqdsk%zlim), unit = 'cm')
    call h5_add(h5id_root, grp // '/psi_eqd', geqdsk%psi_eqd, &
      lbound(geqdsk%psi_eqd), ubound(geqdsk%psi_eqd), unit = 'Mx', comment = 'equidistant psi grid values')
    call h5_add(h5id_root, grp // '/R_eqd', geqdsk%R_eqd, &
      lbound(geqdsk%R_eqd), ubound(geqdsk%R_eqd), unit = 'cm', comment = 'equidistant R grid values')
    call h5_add(h5id_root, grp // '/Z_eqd', geqdsk%Z_eqd, &
      lbound(geqdsk%Z_eqd), ubound(geqdsk%Z_eqd), unit = 'cm', comment = 'equidistant Z grid values')
    call h5_add(h5id_root, grp // '/fprime', geqdsk%fprime, &
      lbound(geqdsk%fprime), ubound(geqdsk%fprime), unit = 'cm^-1')
    call h5_close(h5id_root)
  end subroutine geqdsk_export_hdf5

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
      if (ierr < 0) return
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
    end if
    ! sort eigenvalues in descending order and optionall compute eigenvectors in this order
    call heapsort_complex(ritzvals(:n), complex_abs_desc)
    allocate(eigvals(nritz))
    eigvals(:) = pack(ritzvals, selection)
    if (present(eigvecs)) then
      allocate(ritzvecs(n, nritz), eigvecs(ndim, nritz))
      call hessenberg_eigvecs(hmat(:n, :n), ritzvals(:n), selection(:n), ritzvecs, ierr)
      if (ierr < 0) return
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
      ierr = -1
      return
    end if
    if (size(eigvals, 1) /= ndim) then
      print '("hessenberg_eigvals: eigvals has shape (", i0, "), ' // &
        'but expected (", i0, ")")', shape(hmat), ndim
      ierr = -2
      return
    end if
    allocate(hmat_work(ndim, ndim), zdum(1, ndim))
    hmat_work(:, :) = hmat
    allocate(work(1))
    lwork = -1
    call zhseqr('E', 'N', ndim, 1, ndim, hmat_work, ndim, eigvals, &
      zdum, 1, work, lwork, ierr)
    if (ierr < 0) then
      print '("ZHSEQR: illegal value in argument #", i0)', -ierr
      return
    elseif (ierr > 0) then
      print '("ZHSEQR: only ", i0, " of ", i0, " eigenvalues converged")', &
        ierr, ndim
    end if
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call zhseqr('E', 'N', ndim, 1, ndim, hmat_work, ndim, eigvals, &
      zdum, 1, work, lwork, ierr)
    deallocate(work, hmat_work, zdum)
    if (ierr < 0) then
      print '("ZHSEQR: illegal value in argument #", i0)', -ierr
    elseif (ierr > 0) then
      print '("ZHSEQR: only ", i0, " of ", i0, " eigenvalues converged")', &
        ierr, ndim
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
      ierr = -1
      return
    end if
    if (size(eigvals, 1) /= ndim) then
      print '("hessenberg_eigvecs: eigvals has shape (", i0, "), ' // &
        'but expected (", i0, ")")', shape(eigvals), ndim
      ierr = -2
      return
    end if
    if (size(mask, 1) /= ndim) then
      print '("hessenberg_eigvecs: mask has shape (", i0, "), ' // &
        'but expected (", i0, ")")', shape(mask), ndim
      ierr = -3
      return
    end if
    if (size(eigvecs, 1) /= ndim .or. size(eigvecs, 2) /= count(mask)) then
      print '("hessenberg_eigvecs: eigvecs has shape (", i0, ", ", i0, "), ' // &
        'but expected (", i0, ", ", i0, ")")', shape(eigvecs), ndim, count(mask)
      ierr = -4
      return
    end if
    allocate(eigvals_work(ndim), work(ndim, ndim), rwork(ndim), &
      zdum(1, count(mask)), idum(count(mask)), ifailr(count(mask)))
    eigvals_work(:) = eigvals
    call zhsein('R', 'Q', 'N', mask, ndim, hmat, ndim, eigvals_work, &
      zdum, 1, eigvecs, ndim, count(mask), neff, work, rwork, &
      idum, ifailr, ierr)
    if (ierr < 0) then
      print '("ZHSEIN: illegal value in argument #", i0)', -ierr
    elseif (ierr > 0) then
      print '("ZHSEIN: only ", i0, " of ", i0, " eigenvectors converged")', &
        ierr, ndim
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
