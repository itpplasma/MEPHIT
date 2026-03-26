module mephit_kim
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: kim_t, kim_mode_t, kim_init, kim_deinit, compute_kim_current

  type :: kim_mode_t
    integer :: mpol = 0
    integer :: npts = 0
    integer :: npts_rg = 0
    real(dp), allocatable :: xl_grid(:)
    real(dp), allocatable :: rg_grid(:)
    real(dp), allocatable :: mass_matrix(:,:)
    real(dp), allocatable :: laplace(:,:)
    complex(dp), allocatable :: K_rho_phi(:,:)
    complex(dp), allocatable :: K_rho_B(:,:)
    complex(dp), allocatable :: K_j_phi(:,:)
    complex(dp), allocatable :: K_j_B(:,:)
    real(dp), allocatable :: Er0(:)
    real(dp), allocatable :: B0(:)
    real(dp), allocatable :: kp(:)
    real(dp) :: B_ref = 0d0
    ! Pre-computed quantities
    real(dp), allocatable :: mass_matrix_inv(:,:)
    real(dp), allocatable :: Er0_xl(:)
    real(dp), allocatable :: B0_xl(:)
    real(dp), allocatable :: kp_xl(:)
    ! Sparse system matrix A = L + 4*pi*K_rho_phi (with BC applied)
    complex(dp), allocatable :: A_val(:)
    integer, allocatable :: A_row(:)
    integer, allocatable :: A_col(:)
    integer :: A_nnz = 0
    ! Precomputed coordinate mapping (MEPHIT flux surfaces -> KIM grid)
    real(dp), allocatable :: r_eff(:)       ! effective radius on MEPHIT grid [cm]
    real(dp), allocatable :: dr_dpsi_pol(:) ! Jacobian dr/dpsi_pol on MEPHIT grid
  end type kim_mode_t

  type :: kim_t
    type(kim_mode_t), allocatable :: modes(:)
    integer :: m_min = 0, m_max = 0
  end type kim_t

  real(dp), parameter :: pi = 4d0 * atan(1d0)

contains

  ! === FILE I/O (TEXT FORMAT) ===

  subroutine read_matrix_text(filename, mat, nrow, ncol)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: mat(:,:)
    integer, intent(out) :: nrow, ncol
    integer :: funit, ios, j
    real(dp), allocatable :: row_buf(:)
    character(len=65536) :: line

    ! First pass: determine dimensions
    open(newunit=funit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_matrix_text: cannot open ' // trim(filename)
    nrow = 0
    ncol = 0
    do
      read(funit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      nrow = nrow + 1
      if (nrow == 1) then
        ! Count columns from first row
        ncol = count_columns(line)
      end if
    end do
    close(funit)

    allocate(mat(nrow, ncol))
    allocate(row_buf(ncol))

    ! Second pass: read data
    open(newunit=funit, file=filename, status='old', action='read')
    do j = 1, nrow
      read(funit, *) row_buf
      mat(j, :) = row_buf
    end do
    close(funit)
    deallocate(row_buf)
  end subroutine read_matrix_text

  function count_columns(line) result(ncol)
    character(len=*), intent(in) :: line
    integer :: ncol
    integer :: ios, pos
    real(dp) :: tmp
    character(len=65536) :: buf

    buf = adjustl(line)
    ncol = 0
    do
      if (len_trim(buf) == 0) exit
      read(buf, *, iostat=ios) tmp
      if (ios /= 0) exit
      ncol = ncol + 1
      ! Skip past the number we just read
      buf = adjustl(buf)
      pos = scan(buf, ' ' // achar(9))
      if (pos == 0) exit
      buf = adjustl(buf(pos:))
    end do
  end function count_columns

  subroutine read_grid_text(filename, grid, npts)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: grid(:)
    integer, intent(out) :: npts
    integer :: funit, ios, k
    real(dp) :: x1, x2

    ! Count lines
    open(newunit=funit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_grid_text: cannot open ' // trim(filename)
    npts = 0
    do
      read(funit, *, iostat=ios) x1, x2
      if (ios /= 0) exit
      npts = npts + 1
    end do
    close(funit)

    allocate(grid(npts))
    open(newunit=funit, file=filename, status='old', action='read')
    do k = 1, npts
      read(funit, *) grid(k), x2
    end do
    close(funit)
  end subroutine read_grid_text

  subroutine read_profile_text(filename, grid, profile, npts)
    ! Reads 4-column format: x Re(y) Im(y) |y| -- takes real part only
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: grid(:), profile(:)
    integer, intent(out) :: npts
    integer :: funit, ios, k
    real(dp) :: x, re_y, im_y, abs_y

    open(newunit=funit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_profile_text: cannot open ' // trim(filename)
    npts = 0
    do
      read(funit, *, iostat=ios) x, re_y, im_y, abs_y
      if (ios /= 0) exit
      npts = npts + 1
    end do
    close(funit)

    allocate(grid(npts), profile(npts))
    open(newunit=funit, file=filename, status='old', action='read')
    do k = 1, npts
      read(funit, *) grid(k), profile(k), im_y, abs_y
    end do
    close(funit)
  end subroutine read_profile_text

  subroutine read_scalar_text(filename, val)
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: val
    integer :: funit, ios

    open(newunit=funit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'read_scalar_text: cannot open ' // trim(filename)
    read(funit, *) val
    close(funit)
  end subroutine read_scalar_text

  subroutine read_complex_kernel_text(base_name, kernel, n)
    character(len=*), intent(in) :: base_name
    complex(dp), allocatable, intent(out) :: kernel(:,:)
    integer, intent(out) :: n
    real(dp), allocatable :: mat_re(:,:), mat_im(:,:)
    integer :: nrow_re, ncol_re, nrow_im, ncol_im

    call read_matrix_text(trim(base_name) // '_re.dat', mat_re, nrow_re, ncol_re)
    call read_matrix_text(trim(base_name) // '_im.dat', mat_im, nrow_im, ncol_im)
    if (nrow_re /= nrow_im .or. ncol_re /= ncol_im) &
      error stop 'read_complex_kernel_text: dimension mismatch between re and im'
    n = nrow_re
    allocate(kernel(n, n))
    kernel = cmplx(mat_re, mat_im, dp)
    deallocate(mat_re, mat_im)
  end subroutine read_complex_kernel_text

  subroutine kim_read_mode_text(mode, dir_path)
    type(kim_mode_t), intent(inout) :: mode
    character(len=*), intent(in) :: dir_path
    integer :: n_kern
    real(dp), allocatable :: prof_grid(:)
    integer :: nprof, nrow, ncol

    ! Grid
    call read_grid_text(trim(dir_path) // '/grid/xl_xb.dat', mode%xl_grid, mode%npts)
    call read_grid_text(trim(dir_path) // '/grid/rg_xb.dat', mode%rg_grid, mode%npts_rg)

    ! Kernels
    call read_complex_kernel_text(trim(dir_path) // '/kernel/K_rho_phi', mode%K_rho_phi, n_kern)
    if (n_kern /= mode%npts) error stop 'kim_read_mode_text: K_rho_phi size /= grid size'
    call read_complex_kernel_text(trim(dir_path) // '/kernel/K_rho_B', mode%K_rho_B, n_kern)
    call read_complex_kernel_text(trim(dir_path) // '/kernel/K_j_phi', mode%K_j_phi, n_kern)
    call read_complex_kernel_text(trim(dir_path) // '/kernel/K_j_B', mode%K_j_B, n_kern)

    ! Mass matrix and Laplace
    call read_matrix_text(trim(dir_path) // '/kernel/mass_matrix.dat', mode%mass_matrix, nrow, ncol)
    if (nrow /= mode%npts) error stop 'kim_read_mode_text: mass_matrix size /= grid size'
    call read_matrix_text(trim(dir_path) // '/kernel/Laplace_in_FEM.dat', mode%laplace, nrow, ncol)

    ! B_ref
    call read_scalar_text(trim(dir_path) // '/setup/B_ref.dat', mode%B_ref)

    ! Equilibrium profiles on rg_grid
    call read_profile_text(trim(dir_path) // '/backs/Er.dat', prof_grid, mode%Er0, nprof)
    if (nprof /= mode%npts_rg) error stop 'kim_read_mode_text: Er profile size /= rg_grid size'
    deallocate(prof_grid)
    call read_profile_text(trim(dir_path) // '/backs/B0.dat', prof_grid, mode%B0, nprof)
    deallocate(prof_grid)
    call read_profile_text(trim(dir_path) // '/backs/kp.dat', prof_grid, mode%kp, nprof)
    deallocate(prof_grid)
  end subroutine kim_read_mode_text

  subroutine kim_read_mode_hdf5(mode, h5file)
    use hdf5_tools, only: HID_T, h5_open, h5_close, h5_get, h5_get_bounds
    type(kim_mode_t), intent(inout) :: mode
    character(len=*), intent(in) :: h5file
    integer(HID_T) :: h5id
    integer :: lb1, ub1, n
    real(dp), allocatable :: re_part(:,:), im_part(:,:)

    call h5_open(h5file, h5id)

    ! Grid: query bounds, allocate, then read
    call h5_get_bounds(h5id, 'grid/xl_xb', lb1, ub1)
    mode%npts = ub1 - lb1 + 1
    allocate(mode%xl_grid(mode%npts))
    call h5_get(h5id, 'grid/xl_xb', mode%xl_grid)

    call h5_get_bounds(h5id, 'grid/rg_xb', lb1, ub1)
    mode%npts_rg = ub1 - lb1 + 1
    allocate(mode%rg_grid(mode%npts_rg))
    call h5_get(h5id, 'grid/rg_xb', mode%rg_grid)

    ! Kernels: KIM stores real and imaginary parts as separate datasets
    n = mode%npts
    allocate(re_part(n, n), im_part(n, n))

    call h5_get(h5id, 'kernel/K_rho_phi_re', re_part)
    call h5_get(h5id, 'kernel/K_rho_phi_im', im_part)
    allocate(mode%K_rho_phi(n, n))
    mode%K_rho_phi = cmplx(re_part, im_part, dp)

    call h5_get(h5id, 'kernel/K_rho_B_re', re_part)
    call h5_get(h5id, 'kernel/K_rho_B_im', im_part)
    allocate(mode%K_rho_B(n, n))
    mode%K_rho_B = cmplx(re_part, im_part, dp)

    call h5_get(h5id, 'kernel/K_j_phi_re', re_part)
    call h5_get(h5id, 'kernel/K_j_phi_im', im_part)
    allocate(mode%K_j_phi(n, n))
    mode%K_j_phi = cmplx(re_part, im_part, dp)

    call h5_get(h5id, 'kernel/K_j_B_re', re_part)
    call h5_get(h5id, 'kernel/K_j_B_im', im_part)
    allocate(mode%K_j_B(n, n))
    mode%K_j_B = cmplx(re_part, im_part, dp)

    deallocate(re_part, im_part)

    ! Mass matrix (KIM stores as 'mass_matrix.dat') and Laplace
    allocate(mode%mass_matrix(n, n))
    call h5_get(h5id, 'kernel/mass_matrix.dat', mode%mass_matrix)

    allocate(mode%laplace(n, n))
    call h5_get(h5id, 'kernel/Laplace_in_FEM', mode%laplace)

    ! B_ref
    call h5_get(h5id, 'setup/B_ref', mode%B_ref)

    ! Equilibrium profiles: KIM stores E0r and kp on xl_grid (not rg_grid),
    ! and B0 on a theta grid (not radial). We populate *_xl arrays directly
    ! so kim_init can skip the rg->xl interpolation step.
    allocate(mode%Er0_xl(n))
    call h5_get(h5id, 'backs/E0r', mode%Er0_xl)

    ! B0 is on the theta grid in KIM output, not radial.
    ! Use constant |B_ref| approximation (valid for narrow radial domains).
    ! Must use abs() since B_ref carries the field sign, but B0 here is the
    ! field magnitude used in j_par/B0 and the Poisson boundary conditions.
    allocate(mode%B0_xl(n))
    mode%B0_xl = abs(mode%B_ref)

    allocate(mode%kp_xl(n))
    call h5_get(h5id, 'backs/kp', mode%kp_xl)

    call h5_close(h5id)
  end subroutine kim_read_mode_hdf5

  ! === INTERPOLATION ===

  subroutine interp_linear(x_src, y_src, x_dst, y_dst)
    ! Linear interpolation from (x_src, y_src) to (x_dst, y_dst)
    ! with constant extrapolation (clamping) outside the source domain.
    real(dp), intent(in) :: x_src(:), y_src(:)
    real(dp), intent(in) :: x_dst(:)
    real(dp), intent(out) :: y_dst(:)
    integer :: i, j, n_src, n_dst
    real(dp) :: t

    n_src = size(x_src)
    n_dst = size(x_dst)
    j = 1
    do i = 1, n_dst
      if (x_dst(i) <= x_src(1)) then
        y_dst(i) = y_src(1)
      else if (x_dst(i) >= x_src(n_src)) then
        y_dst(i) = y_src(n_src)
      else
        ! Find bracketing interval
        do while (j < n_src - 1 .and. x_src(j + 1) < x_dst(i))
          j = j + 1
        end do
        t = (x_dst(i) - x_src(j)) / (x_src(j + 1) - x_src(j))
        y_dst(i) = (1d0 - t) * y_src(j) + t * y_src(j + 1)
      end if
    end do
  end subroutine interp_linear

  subroutine interp_linear_zeropad(x_src, y_src, x_dst, y_dst)
    ! Linear interpolation with zero-padding outside the source domain.
    ! The KIM domain should be wide enough that the physical quantity
    ! (current, potential) naturally decays near the boundaries.
    real(dp), intent(in) :: x_src(:), y_src(:)
    real(dp), intent(in) :: x_dst(:)
    real(dp), intent(out) :: y_dst(:)
    integer :: i, j, n_src, n_dst
    real(dp) :: t

    n_src = size(x_src)
    n_dst = size(x_dst)
    j = 1
    do i = 1, n_dst
      if (x_dst(i) <= x_src(1) .or. x_dst(i) >= x_src(n_src)) then
        y_dst(i) = 0d0
      else
        ! Find bracketing interval
        do while (j < n_src - 1 .and. x_src(j + 1) < x_dst(i))
          j = j + 1
        end do
        t = (x_dst(i) - x_src(j)) / (x_src(j + 1) - x_src(j))
        y_dst(i) = (1d0 - t) * y_src(j) + t * y_src(j + 1)
      end if
    end do
  end subroutine interp_linear_zeropad

  subroutine interp_linear_complex(x_src, y_src, x_dst, y_dst)
    real(dp), intent(in) :: x_src(:), x_dst(:)
    complex(dp), intent(in) :: y_src(:)
    complex(dp), intent(out) :: y_dst(:)
    real(dp), allocatable :: re_src(:), im_src(:), re_dst(:), im_dst(:)
    integer :: n_src, n_dst

    n_src = size(x_src)
    n_dst = size(x_dst)
    allocate(re_src(n_src), im_src(n_src), re_dst(n_dst), im_dst(n_dst))
    re_src = real(y_src, dp)
    im_src = aimag(y_src)
    call interp_linear(x_src, re_src, x_dst, re_dst)
    call interp_linear(x_src, im_src, x_dst, im_dst)
    y_dst = cmplx(re_dst, im_dst, dp)
    deallocate(re_src, im_src, re_dst, im_dst)
  end subroutine interp_linear_complex

  subroutine interp_linear_complex_zeropad(x_src, y_src, x_dst, y_dst)
    ! Complex interpolation with zero-padding outside the source domain.
    real(dp), intent(in) :: x_src(:), x_dst(:)
    complex(dp), intent(in) :: y_src(:)
    complex(dp), intent(out) :: y_dst(:)
    real(dp), allocatable :: re_src(:), im_src(:), re_dst(:), im_dst(:)
    integer :: n_src, n_dst

    n_src = size(x_src)
    n_dst = size(x_dst)
    allocate(re_src(n_src), im_src(n_src), re_dst(n_dst), im_dst(n_dst))
    re_src = real(y_src, dp)
    im_src = aimag(y_src)
    call interp_linear_zeropad(x_src, re_src, x_dst, re_dst)
    call interp_linear_zeropad(x_src, im_src, x_dst, im_dst)
    y_dst = cmplx(re_dst, im_dst, dp)
    deallocate(re_src, im_src, re_dst, im_dst)
  end subroutine interp_linear_complex_zeropad

  ! === DENSE TO SPARSE CONVERSION ===

  subroutine dense_to_csc(A, n, threshold, val, row_idx, col_ptr, nnz)
    ! Convert dense complex matrix to CSC (1-based indexing)
    complex(dp), intent(in) :: A(:,:)
    integer, intent(in) :: n
    real(dp), intent(in) :: threshold
    complex(dp), allocatable, intent(out) :: val(:)
    integer, allocatable, intent(out) :: row_idx(:), col_ptr(:)
    integer, intent(out) :: nnz
    integer :: i, j, k

    ! Count nonzeros
    nnz = 0
    do j = 1, n
      do i = 1, n
        if (abs(A(i, j)) > threshold) nnz = nnz + 1
      end do
    end do

    allocate(val(nnz), row_idx(nnz), col_ptr(n + 1))

    k = 1
    do j = 1, n
      col_ptr(j) = k
      do i = 1, n
        if (abs(A(i, j)) > threshold) then
          val(k) = A(i, j)
          row_idx(k) = i
          k = k + 1
        end if
      end do
    end do
    col_ptr(n + 1) = k
  end subroutine dense_to_csc

  ! === MATRIX INVERSION ===

  subroutine invert_real_matrix(A, A_inv, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n, n)
    real(dp), intent(out) :: A_inv(n, n)
    real(dp) :: work_mat(n, n)
    integer :: ipiv(n), info
    real(dp), allocatable :: work(:)
    integer :: lwork

    ! External LAPACK routines
    external :: dgetrf, dgetri

    work_mat = A
    ! LU factorization
    call dgetrf(n, n, work_mat, n, ipiv, info)
    if (info /= 0) error stop 'invert_real_matrix: dgetrf failed'
    ! Workspace query
    lwork = -1
    allocate(work(1))
    call dgetri(n, work_mat, n, ipiv, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    ! Inversion
    call dgetri(n, work_mat, n, ipiv, work, lwork, info)
    if (info /= 0) error stop 'invert_real_matrix: dgetri failed'
    A_inv = work_mat
    deallocate(work)
  end subroutine invert_real_matrix

  ! === INITIALIZATION ===

  subroutine kim_init(kim, m_res_min, m_res_max, nflux, psi_pol, q_prof)
    use mephit_conf, only: conf, logger
    type(kim_t), intent(inout) :: kim
    integer, intent(in) :: m_res_min, m_res_max, nflux
    real(dp), intent(in) :: psi_pol(0:)  ! poloidal flux on MEPHIT flux surfaces
    real(dp), intent(in) :: q_prof(0:)   ! safety factor on MEPHIT flux surfaces

    integer :: m, kf
    character(len=2048) :: dir_path, h5_path
    character(len=16) :: m_str
    logical :: dir_exists, h5_exists
    complex(dp), allocatable :: A_dense(:,:)
    integer :: n
    real(dp) :: psi_tor(0:nflux)

    kim%m_min = m_res_min
    kim%m_max = m_res_max
    allocate(kim%modes(m_res_min:m_res_max))

    ! Compute toroidal flux on MEPHIT grid (mode-independent)
    psi_tor(0) = 0d0
    do kf = 1, nflux
      psi_tor(kf) = psi_tor(kf - 1) + &
        0.5d0 * (q_prof(kf - 1) + q_prof(kf)) * (psi_pol(kf) - psi_pol(kf - 1))
    end do

    do m = m_res_min, m_res_max
      write(m_str, '(i0)') m

      ! Determine file format: check for HDF5 file first, then text directory.
      ! Use inquire(file=...) for portability (no GNU directory extension).
      h5_path = trim(conf%kim_data_dir) // '/kim_m' // trim(m_str) // '.h5'
      dir_path = trim(conf%kim_data_dir) // '/kim_m' // trim(m_str)
      inquire(file=trim(h5_path), exist=h5_exists)
      ! Check for text directory by probing for a known file within it
      inquire(file=trim(dir_path) // '/grid/xl_xb.dat', exist=dir_exists)

      kim%modes(m)%mpol = m
      if (dir_exists) then
        if (logger%info) then
          write(logger%msg, '("KIM: reading text data for mode m = ", i0)') m
          call logger%write_msg
        end if
        call kim_read_mode_text(kim%modes(m), dir_path)
      else if (h5_exists) then
        if (logger%info) then
          write(logger%msg, '("KIM: reading HDF5 data for mode m = ", i0)') m
          call logger%write_msg
        end if
        call kim_read_mode_hdf5(kim%modes(m), h5_path)
      else
        write(logger%msg, '("KIM: no data for mode m = ", i0, &
          &" in ", a, " - skipping (zero current)")') m, trim(conf%kim_data_dir)
        if (logger%warn) call logger%write_msg
        kim%modes(m)%npts = 0
      end if

      ! Skip processing for modes with no data (npts == 0)
      if (kim%modes(m)%npts == 0) cycle

      n = kim%modes(m)%npts

      ! Interpolate profiles from rg_grid to xl_grid (text format).
      ! For HDF5, Er0_xl/B0_xl/kp_xl are already populated by the reader.
      if (.not. allocated(kim%modes(m)%Er0_xl)) then
        allocate(kim%modes(m)%Er0_xl(n))
        allocate(kim%modes(m)%B0_xl(n))
        allocate(kim%modes(m)%kp_xl(n))
        call interp_linear(kim%modes(m)%rg_grid, kim%modes(m)%Er0, &
          kim%modes(m)%xl_grid, kim%modes(m)%Er0_xl)
        call interp_linear(kim%modes(m)%rg_grid, kim%modes(m)%B0, &
          kim%modes(m)%xl_grid, kim%modes(m)%B0_xl)
        call interp_linear(kim%modes(m)%rg_grid, kim%modes(m)%kp, &
          kim%modes(m)%xl_grid, kim%modes(m)%kp_xl)
      end if

      ! Invert mass matrix
      allocate(kim%modes(m)%mass_matrix_inv(n, n))
      call invert_real_matrix(kim%modes(m)%mass_matrix, kim%modes(m)%mass_matrix_inv, n)

      ! Assemble system matrix A = L + 4*pi * K_rho_phi
      allocate(A_dense(n, n))
      A_dense = cmplx(kim%modes(m)%laplace, 0d0, dp) + 4d0 * pi * kim%modes(m)%K_rho_phi

      ! Apply BC: boundary rows become identity
      A_dense(1, :) = (0d0, 0d0)
      A_dense(1, 1) = (1d0, 0d0)
      A_dense(n, :) = (0d0, 0d0)
      A_dense(n, n) = (1d0, 0d0)

      ! Convert to sparse CSC
      call dense_to_csc(A_dense, n, 1d-14, &
        kim%modes(m)%A_val, kim%modes(m)%A_row, kim%modes(m)%A_col, kim%modes(m)%A_nnz)
      deallocate(A_dense)

      ! Precompute coordinate mapping: r_eff and dr/dpsi_pol on MEPHIT grid
      allocate(kim%modes(m)%r_eff(0:nflux))
      allocate(kim%modes(m)%dr_dpsi_pol(0:nflux))
      kim%modes(m)%r_eff(0) = 0d0
      kim%modes(m)%dr_dpsi_pol(0) = 0d0
      do kf = 1, nflux
        kim%modes(m)%r_eff(kf) = sqrt(2d0 * abs(psi_tor(kf)) / abs(kim%modes(m)%B_ref))
        if (kim%modes(m)%r_eff(kf) > 0d0) then
          kim%modes(m)%dr_dpsi_pol(kf) = q_prof(kf) / (abs(kim%modes(m)%B_ref) * kim%modes(m)%r_eff(kf))
        else
          kim%modes(m)%dr_dpsi_pol(kf) = 0d0
        end if
      end do

      if (logger%info) then
        write(logger%msg, '("KIM mode m = ", i0, ": grid npts = ", i0, &
          &", sparse nnz = ", i0)') m, n, kim%modes(m)%A_nnz
        call logger%write_msg
      end if
    end do
  end subroutine kim_init

  ! === POISSON SOLVER ===

  subroutine kim_solve_poisson(mode, Br_kim, Phi_kim)
    use sparse_mod, only: sparse_solve_suitesparseComplex_b1
    type(kim_mode_t), intent(in) :: mode
    complex(dp), intent(in) :: Br_kim(:)
    complex(dp), intent(out) :: Phi_kim(:)

    complex(dp) :: rhs(mode%npts)
    complex(dp) :: Phi_aligned_lo, Phi_aligned_hi
    complex(dp), parameter :: imun = (0d0, 1d0)
    integer :: n

    n = mode%npts

    ! Assemble RHS: b = -4*pi * K_rho_B @ Br_kim
    rhs = -4d0 * pi * matmul(mode%K_rho_B, Br_kim)

    ! Boundary conditions: zero misalignment electric field
    ! Phi_aligned = i * E_r * B_r / (B0 * k_parallel)
    ! (KIM convention: com_unit = imaginary unit)
    if (abs(mode%B0_xl(1) * mode%kp_xl(1)) > 1d-30) then
      Phi_aligned_lo = imun * mode%Er0_xl(1) * Br_kim(1) / &
        (mode%B0_xl(1) * mode%kp_xl(1))
    else
      Phi_aligned_lo = (0d0, 0d0)
    end if
    if (abs(mode%B0_xl(n) * mode%kp_xl(n)) > 1d-30) then
      Phi_aligned_hi = imun * mode%Er0_xl(n) * Br_kim(n) / &
        (mode%B0_xl(n) * mode%kp_xl(n))
    else
      Phi_aligned_hi = (0d0, 0d0)
    end if
    rhs(1) = -Phi_aligned_lo
    rhs(n) = -Phi_aligned_hi

    ! Solve A * Phi = rhs using UMFPACK
    Phi_kim = rhs
    call sparse_solve_suitesparseComplex_b1(n, n, mode%A_nnz, &
      mode%A_row, mode%A_col, mode%A_val, Phi_kim, 0)
  end subroutine kim_solve_poisson

  ! === CURRENT DENSITY ===

  subroutine kim_compute_current(mode, Phi_kim, Br_kim, jpar_over_B0)
    type(kim_mode_t), intent(in) :: mode
    complex(dp), intent(in) :: Phi_kim(:), Br_kim(:)
    complex(dp), intent(out) :: jpar_over_B0(:)

    complex(dp) :: jpar(mode%npts)
    integer :: i

    ! j_par = M^{-1} @ (K_j_phi @ Phi + K_j_B @ Br)
    jpar = matmul(cmplx(mode%mass_matrix_inv, 0d0, dp), &
      matmul(mode%K_j_phi, Phi_kim) + matmul(mode%K_j_B, Br_kim))

    ! j_par / B0
    do i = 1, mode%npts
      jpar_over_B0(i) = jpar(i) / mode%B0_xl(i)
    end do
  end subroutine kim_compute_current

  ! === TOP-LEVEL: CALLED FROM mephit_iter.F90 ===

  subroutine compute_kim_current(Bn, jnpar_B0_mhd, jmnpar_over_Bmod, kim, Phi_mn, Phi_aligned_mn)
    use mephit_conf, only: conf, transition_func_none
    use mephit_util, only: imun, resample1d
    use mephit_mesh, only: equil, mesh, fs, fs_half, dPhi0_dpsi, cache
    use mephit_pert, only: polmodes_t, polmodes_init, polmodes_deinit, &
      L1_poloidal_modes, &
      vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, RT0_poloidal_modes, &
      RT0_t, L1_t

    type(RT0_t), intent(in) :: Bn
    type(L1_t), intent(in) :: jnpar_B0_mhd
    type(polmodes_t), intent(inout) :: jmnpar_over_Bmod
    type(kim_t), intent(in) :: kim
    complex(dp), intent(out) :: Phi_mn(0:, mesh%m_res_min:)
    complex(dp), intent(out) :: Phi_aligned_mn(0:, mesh%m_res_min:)

    integer :: m, m_res, kf
    complex(dp) :: Bmnpsi_over_B0phi(1:mesh%nflux)
    real(dp) :: Bmnrho_over_B0theta(1:mesh%nflux), rtemp(1:mesh%nflux)
    type(vec_polmodes_t) :: Bmn
    type(polmodes_t) :: jmnpar_over_Bmod_mhd

    ! KIM grid quantities
    complex(dp), allocatable :: Br_kim(:), Phi_kim(:), jpar_B0_kim(:)

    jmnpar_over_Bmod%coeff(:, :) = (0d0, 0d0)
    Phi_mn(:, :) = (0d0, 0d0)
    Phi_aligned_mn(:, :) = (0d0, 0d0)

    call vec_polmodes_init(Bmn, conf%m_max, mesh%nflux)
    call RT0_poloidal_modes(Bn, Bmn)

    ! Get MHD baseline current for transition function (same as FLR2 path)
    if (conf%transition_func /= transition_func_none) then
      call polmodes_init(jmnpar_over_Bmod_mhd, conf%m_max, mesh%nflux)
      call L1_poloidal_modes(jnpar_B0_mhd, jmnpar_over_Bmod_mhd)
    end if

    do m = mesh%m_res_min, mesh%m_res_max
      ! Skip modes with no KIM data (zero current contribution)
      if (kim%modes(m)%npts == 0) cycle

      m_res = -equil%cocos%sgn_q * m

      ! Extract Bmnpsi_over_B0phi (same as in compute_flr2_current)
      Bmnrho_over_B0theta(:) = Bmn%coeff_rad(m_res, :)%Re
      call resample1d(fs_half%psi, Bmnrho_over_B0theta / fs_half%q * &
        equil%cocos%sgn_dpsi, fs%psi(1:), rtemp, 3)
      Bmnpsi_over_B0phi%Re = rtemp
      Bmnrho_over_B0theta(:) = Bmn%coeff_rad(m_res, :)%Im
      call resample1d(fs_half%psi, Bmnrho_over_B0theta / fs_half%q * &
        equil%cocos%sgn_dpsi, fs%psi(1:), rtemp, 3)
      Bmnpsi_over_B0phi%Im = rtemp

      allocate(Br_kim(kim%modes(m)%npts), Phi_kim(kim%modes(m)%npts), &
        jpar_B0_kim(kim%modes(m)%npts))

      ! Transform B^psi to B^r and interpolate to KIM grid
      ! B^r on MEPHIT = dr/dpsi_pol * Bmnpsi_over_B0phi * B0^phi
      ! where B0^phi = F/R, using F/R_O as cylindrical approximation
      block
        complex(dp) :: Br_mephit(1:mesh%nflux)

        do kf = 1, mesh%nflux
          Br_mephit(kf) = kim%modes(m)%dr_dpsi_pol(kf) * Bmnpsi_over_B0phi(kf) * &
            fs%F(kf) / mesh%R_O
        end do

        ! Interpolate from MEPHIT's r_eff to KIM's xl_grid
        call interp_linear_complex(kim%modes(m)%r_eff(1:), Br_mephit, &
          kim%modes(m)%xl_grid, Br_kim)
      end block

      ! Solve Poisson equation
      call kim_solve_poisson(kim%modes(m), Br_kim, Phi_kim)

      ! Compute parallel current density
      call kim_compute_current(kim%modes(m), Phi_kim, Br_kim, jpar_B0_kim)

      ! Interpolate results back to MEPHIT's flux surface grid.
      ! Use zero-padding: KIM current/potential is zero outside the
      ! restricted radial domain around the resonance.
      block
        complex(dp) :: Phi_back(1:mesh%nflux), jpar_back(1:mesh%nflux)

        call interp_linear_complex_zeropad(kim%modes(m)%xl_grid, Phi_kim, &
          kim%modes(m)%r_eff(1:), Phi_back)
        call interp_linear_complex_zeropad(kim%modes(m)%xl_grid, jpar_B0_kim, &
          kim%modes(m)%r_eff(1:), jpar_back)

        Phi_mn(1:, m) = Phi_back
        jmnpar_over_Bmod%coeff(m_res, 1:) = jpar_back
      end block

      ! Compute Phi_aligned for diagnostics (same formula as FLR2 path)
      Phi_aligned_mn(1:, m) = imun * Bmnpsi_over_B0phi * &
        fs%q(1:) * dPhi0_dpsi%y(1:) / (m_res + mesh%n * fs%q(1:))

      ! Apply transition function: cross-fade between KIM and MHD current.
      ! This matches the FLR2 path (mephit_iter.F90 compute_flr2_current).
      ! Without this, the MHD current would be double-counted at the resonance
      ! because compute_currn adds resonant_jmnpar_over_Bmod to the MHD baseline.
      if (conf%transition_func /= transition_func_none) then
        jmnpar_over_Bmod%coeff(m_res, :) = cache%shielding(m)%cross_fade * &
          (jmnpar_over_Bmod%coeff(m_res, :) - jmnpar_over_Bmod_mhd%coeff(m_res, :))
      end if

      deallocate(Br_kim, Phi_kim, jpar_B0_kim)
    end do

    call vec_polmodes_deinit(Bmn)
    if (conf%transition_func /= transition_func_none) then
      call polmodes_deinit(jmnpar_over_Bmod_mhd)
    end if
  end subroutine compute_kim_current

  ! === CLEANUP ===

  subroutine kim_mode_deinit(mode)
    type(kim_mode_t), intent(inout) :: mode

    if (allocated(mode%xl_grid)) deallocate(mode%xl_grid)
    if (allocated(mode%rg_grid)) deallocate(mode%rg_grid)
    if (allocated(mode%mass_matrix)) deallocate(mode%mass_matrix)
    if (allocated(mode%laplace)) deallocate(mode%laplace)
    if (allocated(mode%K_rho_phi)) deallocate(mode%K_rho_phi)
    if (allocated(mode%K_rho_B)) deallocate(mode%K_rho_B)
    if (allocated(mode%K_j_phi)) deallocate(mode%K_j_phi)
    if (allocated(mode%K_j_B)) deallocate(mode%K_j_B)
    if (allocated(mode%Er0)) deallocate(mode%Er0)
    if (allocated(mode%B0)) deallocate(mode%B0)
    if (allocated(mode%kp)) deallocate(mode%kp)
    if (allocated(mode%mass_matrix_inv)) deallocate(mode%mass_matrix_inv)
    if (allocated(mode%Er0_xl)) deallocate(mode%Er0_xl)
    if (allocated(mode%B0_xl)) deallocate(mode%B0_xl)
    if (allocated(mode%kp_xl)) deallocate(mode%kp_xl)
    if (allocated(mode%A_val)) deallocate(mode%A_val)
    if (allocated(mode%A_row)) deallocate(mode%A_row)
    if (allocated(mode%A_col)) deallocate(mode%A_col)
    if (allocated(mode%r_eff)) deallocate(mode%r_eff)
    if (allocated(mode%dr_dpsi_pol)) deallocate(mode%dr_dpsi_pol)
    mode%npts = 0
    mode%npts_rg = 0
    mode%A_nnz = 0
  end subroutine kim_mode_deinit

  subroutine kim_deinit(kim)
    type(kim_t), intent(inout) :: kim
    integer :: m

    if (allocated(kim%modes)) then
      do m = kim%m_min, kim%m_max
        call kim_mode_deinit(kim%modes(m))
      end do
      deallocate(kim%modes)
    end if
    kim%m_min = 0
    kim%m_max = 0
  end subroutine kim_deinit

end module mephit_kim
