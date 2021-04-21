module magdif
  use iso_fortran_env, only: dp => real64
  use magdif_pert, only: L1_t, RT0_t

  implicit none

  private

  public :: magdif_init, magdif_cleanup, magdif_iterate, magdif_postprocess, freefem_pipe

  character(len = 1024) :: freefem_pipe = 'maxwell.dat'

  !> Pressure perturbation \f$ p_{n} \f$ in dyn cm^-1.
  type(L1_t) :: pn

  !> Perturbation field in units of Gauss.
  type(RT0_t) :: Bn

  !> Vacuum perturbation field in units of Gauss.
  type(RT0_t) :: Bnvac

  !> Plasma perturbation field in units of Gauss.
  type(RT0_t) :: Bnplas

  !> Current density perturbation in units of statampere cm^-2.
  type(RT0_t) :: jn

contains

  !> Initialize magdif module
  subroutine magdif_init
    use magdata_in_symfluxcoor_mod, only: load_magdata_in_symfluxcoord
    use magdif_conf, only: conf, conf_arr, log, magdif_log, datafile
    use magdif_util, only: get_field_filenames, init_field
    use magdif_mesh, only: equil, mesh, read_mesh_cache, fluxvar, flux_func_cache_check, &
         check_curr0, check_safety_factor
    use magdif_pert, only: RT0_check_div_free, RT0_check_redundant_edges, &
         RT0_init, RT0_read, L1_init
    character(len = 1024) :: gfile, pfile, convexfile
    integer :: dum

    log = magdif_log('-', conf%log_level, conf%quiet)

    call get_field_filenames(gfile, pfile, convexfile)
    call equil%read(gfile)
    call equil%classify
    if (equil%cocos%index /= 3) then
       write (log%msg, '("GEQDSK file ", a, " is not conforming to COCOS 3")') trim(gfile)
       if (log%err) call log%write
       error stop
    end if
    call init_field(gfile, pfile, convexfile)

    ! read in preprocessed data
    call read_mesh_cache
    call load_magdata_in_symfluxcoord
    ! TODO: save previously processed config parameters to HDF5 and load here
    call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
    ! TODO: cache Lagrange polynomials instead
    call fluxvar%init(4, equil%psi_eqd)

    ! check preprocessed data
    call flux_func_cache_check
    call check_curr0
    call check_safety_factor

    ! initialize perturbation
    call L1_init(pn, mesh%npoint)
    call RT0_init(Bn, mesh%ntri)
    call RT0_init(Bnvac, mesh%ntri)
    call RT0_init(Bnplas, Bn%ntri)
    call RT0_init(jn, mesh%ntri)

    call RT0_read(Bnvac, datafile, 'Bnvac')
    call RT0_check_redundant_edges(Bnvac, 'Bnvac')
    call RT0_check_div_free(Bnvac, mesh%n, conf%rel_err_Bn, 'Bnvac')

    ! pass effective toroidal mode number to FreeFem++
    call send_flag_to_freefem(mesh%n, freefem_pipe)
    call receive_flag_from_freefem(dum, freefem_pipe)

    log%msg = 'magdif initialized'
    if (log%info) call log%write
  end subroutine magdif_init

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    use magdif_conf, only: log
    use magdif_mesh, only: B0R, B0phi, B0Z, B0R_Omega, B0phi_Omega, B0Z_Omega, B0flux, &
         j0phi
    use magdif_pert, only: L1_deinit, RT0_deinit
    integer :: dum

    ! tell FreeFem++ to stop processing
    call send_flag_to_freefem(-3, freefem_pipe)
    call receive_flag_from_freefem(dum, freefem_pipe)
    if (allocated(B0r)) deallocate(B0r)
    if (allocated(B0phi)) deallocate(B0phi)
    if (allocated(B0z)) deallocate(B0z)
    if (allocated(B0r_Omega)) deallocate(B0r_Omega)
    if (allocated(B0phi_Omega)) deallocate(B0phi_Omega)
    if (allocated(B0z_Omega)) deallocate(B0z_Omega)
    if (allocated(B0flux)) deallocate(B0flux)
    if (allocated(j0phi)) deallocate(j0phi)
    call L1_deinit(pn)
    call RT0_deinit(Bn)
    call RT0_deinit(Bnvac)
    call RT0_deinit(Bnplas)
    call RT0_deinit(jn)
    log%msg = 'magdif cleanup finished'
    if (log%info) call log%write
  end subroutine magdif_cleanup

  subroutine magdif_single
    ! compute pressure based on previous perturbation field
    call compute_presn
    ! compute currents based on previous perturbation field
    call compute_currn
    ! use field code to generate new field from currents
    call compute_Bn
  end subroutine magdif_single

  subroutine magdif_iterate
    use arnoldi_mod, only: ieigen, ngrow, tol, eigvecs  ! arnoldi.f90
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdif_mesh, only: mesh
    use magdif_pert, only: RT0_init, RT0_deinit, RT0_write, L1_write, &
         RT0_check_div_free, RT0_check_redundant_edges, vec_polmodes_t, vec_polmodes_init, &
         vec_polmodes_deinit, vec_polmodes_write, RT0_poloidal_modes
    use magdif_conf, only: conf, log, runmode_precon, runmode_single, datafile, cmplx_fmt, &
         decorate_filename

    logical :: preconditioned
    integer :: kiter, niter, ndim, i, j, info
    integer(HID_T) :: h5id_root
    real(dp), allocatable :: L2int(:)
    type(RT0_t) :: Bn_diff
    complex(dp) :: packed_Bn(mesh%nedge), packed_Bn_prev(mesh%nedge)
    complex(dp) :: eigvals(conf%nritz)
    complex(dp), allocatable :: Lr(:,:), Yr(:,:)
    integer, allocatable :: ipiv(:)
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"
    integer, parameter :: m_max = 24
    type(vec_polmodes_t) :: jmn

    ! system dimension: number of non-redundant edges in core plasma
    ndim = mesh%nedge
    ! runmodes
    preconditioned = runmode_precon == conf%runmode
    if (runmode_single == conf%runmode) then
       niter = 0
    else
       niter = conf%niter
    end if
    if (preconditioned) then
       tol = conf%ritz_threshold
       ! calculate eigenvectors
       ieigen = 1
       call arnoldi(ndim, conf%nritz, eigvals, next_iteration_arnoldi)
       call h5_open_rw(datafile, h5id_root)
       call h5_create_parent_groups(h5id_root, 'iter/')
       call h5_add(h5id_root, 'iter/eigvals', eigvals, lbound(eigvals), ubound(eigvals), &
         comment = 'iteration eigenvalues')
       call h5_close(h5id_root)
       if (log%info) then
          write (log%msg, '("Arnoldi method yields ", i0, " Ritz eigenvalues > ", f0.2)') &
               ngrow, tol
          call log%write
          do i = 1, ngrow
             write (log%msg, '("lambda ", i0, ": ", ' // cmplx_fmt // ')') i, eigvals(i)
             call log%write
          end do
       end if
       if (ngrow > 0) then
          do i = 1, min(ngrow, conf%max_eig_out)
             write (postfix, postfix_fmt) i
             call unpack_dof(Bn, eigvecs(:, i))
             call RT0_write(Bn, datafile, 'iter/eigvec' // postfix, &
                  'iteration eigenvector', 'G', 1)
          end do
          allocate(Lr(ngrow, ngrow), Yr(ngrow, ngrow))
          Yr = (0d0, 0d0)
          do i = 1, ngrow
             Yr(i, i) = (1d0, 0d0)
             do j = 1, ngrow
                Lr(i, j) = sum(conjg(eigvecs(:, i)) * eigvecs(:, j)) * (eigvals(j) - (1d0, 0d0))
             end do
          end do
          allocate(ipiv(ngrow))
          call zgesv(ngrow, ngrow, Lr, ngrow, ipiv, Yr, ngrow, info)
          if (allocated(ipiv)) deallocate(ipiv)
          if (info == 0) then
             log%msg = 'Successfully inverted matrix for preconditioner'
             if (log%info) call log%write
          else
             write (log%msg, '("Matrix inversion for preconditioner failed: ' // &
                  'zgesv returns error ", i0)') info
             if (log%err) call log%write
             error stop
          end if
          do i = 1, ngrow
             Lr(i, :) = eigvals(i) * Yr(i, :)
          end do
          if (allocated(Yr)) deallocate(Yr)
       else
          preconditioned = .false.
       end if
    end if

    call vec_polmodes_init(jmn, m_max, mesh%nflux)
    call RT0_init(Bn_diff, mesh%ntri)
    allocate(L2int(0:niter))
    call compute_L2int(Bnvac, L2int(0))
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, 'iter/')
    call h5_add(h5id_root, 'iter/L2int_Bnvac', L2int(0), &
         comment = 'L2 integral of magnetic field (vacuum)', unit = 'Mx')
    call h5_close(h5id_root)
    Bn%DOF(:, :) = Bnvac%DOF
    Bn%comp_phi(:) = Bnvac%comp_phi
    call pack_dof(Bnvac, packed_Bn_prev)
    if (preconditioned) then
       packed_Bn_prev(:) = packed_Bn_prev - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
            matmul(transpose(conjg(eigvecs(:, 1:ngrow))), packed_Bn_prev)))
    end if
    do kiter = 0, niter
       write (log%msg, '("Iteration ", i2, " of ", i2)') kiter, niter
       if (log%info) call log%write
       write (postfix, postfix_fmt) kiter

       call next_iteration(ndim, packed_Bn_prev, packed_Bn)
       if (preconditioned) then
          packed_Bn(:) = packed_Bn - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
               matmul(transpose(conjg(eigvecs(:, 1:ngrow))), packed_Bn - packed_Bn_prev)))
          call unpack_dof(Bn, packed_Bn)
          call RT0_check_redundant_edges(Bn, 'Bn')
          call RT0_check_div_free(Bn, mesh%n, conf%rel_err_Bn, 'Bn')
       end if

       call unpack_dof(Bn_diff, packed_Bn - packed_Bn_prev)
       call compute_L2int(Bn_diff, L2int(kiter))
       if (kiter <= 1) then
          call L1_write(pn, datafile, 'iter/pn' // postfix, &
               'pressure (after iteration)', 'dyn cm^-2')
          call RT0_write(jn, datafile, 'iter/jn' // postfix, &
               'current density (after iteration)', 'statA cm^-2', 1)
          call RT0_write(Bn, datafile, 'iter/Bn' // postfix, &
               'magnetic field (after iteration)', 'G', 1)
          call RT0_write(Bn_diff, datafile, 'iter/Bn_diff' // postfix, &
               'magnetic field (difference between iterations)', 'G', 1)
          call RT0_poloidal_modes(jn, jmn)
          call vec_polmodes_write(jmn, datafile, 'postprocess/jmn' // postfix, &
               'current density (after iteration)', 'statA cm^-2')
       end if

       call pack_dof(Bn, packed_Bn_prev)
    end do
    Bnplas%DOF(:, :) = Bn%DOF - Bnvac%DOF
    Bnplas%comp_phi(:) = Bn%comp_phi - Bnvac%comp_phi
    call h5_open_rw(datafile, h5id_root)
    call h5_add(h5id_root, 'iter/L2int_Bn_diff', L2int, lbound(L2int), ubound(L2int), &
         comment = 'L2 integral of magnetic field (difference between iterations)', unit = 'Mx')
    call h5_close(h5id_root)
    call L1_write(pn, datafile, 'iter/pn', &
         'pressure (full perturbation)', 'dyn cm^-2')
    call RT0_write(Bn, datafile, 'iter/Bn', &
         'magnetic field (full perturbation)', 'G', 2)
    call RT0_write(Bnplas, datafile, 'iter/Bnplas', &
         'magnetic field (plasma response)', 'G', 2)
    call RT0_write(jn, datafile, 'iter/jn', &
         'current density (full perturbation)', 'statA cm^-2', 1)

    call RT0_deinit(Bn_diff)
    call vec_polmodes_deinit(jmn)
    if (allocated(Lr)) deallocate(Lr)
    if (allocated(L2int)) deallocate(L2int)

  contains

    pure subroutine pack_dof(elem, packed)
      type(RT0_t), intent(in) :: elem
      complex(dp), intent(out) :: packed(mesh%nedge)
      integer :: kedge
      do kedge = 1, mesh%nedge
         packed(kedge) = elem%DOF(mesh%edge_map2ke(1, kedge), mesh%edge_map2ktri(1, kedge))
      end do
    end subroutine pack_dof

    pure subroutine unpack_dof(elem, packed)
      use magdif_util, only: imun
      type(RT0_t), intent(inout) :: elem
      complex(dp), intent(in) :: packed(mesh%nedge)
      integer :: kedge, ktri
      do kedge = 1, mesh%nedge
         elem%DOF(mesh%edge_map2ke(1, kedge), mesh%edge_map2ktri(1, kedge)) = &
              packed(kedge)
         if (mesh%edge_map2ktri(2, kedge) > 0) then
            elem%DOF(mesh%edge_map2ke(2, kedge), mesh%edge_map2ktri(2, kedge)) = &
                 -packed(kedge)
         end if
      end do
      do ktri = 1, mesh%ntri
         elem%comp_phi(ktri) = sum(elem%DOF(:, ktri)) * imun / mesh%n / mesh%area(ktri)
      end do
    end subroutine unpack_dof

    ! computes B_(n+1) = K*B_n + B_vac ... different from kin2d.f90
    subroutine next_iteration(n, xold, xnew)
      integer, intent(in) :: n
      complex(dp), intent(in) :: xold(:)
      complex(dp), intent(out) :: xnew(:)
      if (n /= size(xold)) then
         call log%msg_arg_size('next_iteration', 'n', 'size(xold)', n, size(xold))
         if (log%err) call log%write
         error stop
      end if
      if (n /= size(xnew)) then
         call log%msg_arg_size('next_iteration', 'n', 'size(xnew)', n, size(xnew))
         if (log%err) call log%write
         error stop
      end if
      call unpack_dof(Bn, xold)
      call magdif_single
      Bn%DOF(:, :) = Bn%DOF + Bnvac%DOF
      Bn%comp_phi(:) = Bn%comp_phi + Bnvac%comp_phi
      call pack_dof(Bn, xnew)
    end subroutine next_iteration

    ! computes B_(n+1) = K*(B_n + B_vac) ... as in kin2d.f90
    ! next_iteration in arnoldi_mod is still declared external and has no interface,
    ! so we use explicit-shape arrays here
    subroutine next_iteration_arnoldi(n, xold, xnew)
      integer, intent(in) :: n
      complex(dp), intent(in) :: xold(n)
      complex(dp), intent(out) :: xnew(n)
      call unpack_dof(Bn, xold)
      Bn%DOF(:, :) = Bn%DOF + Bnvac%DOF
      Bn%comp_phi(:) = Bn%comp_phi + Bnvac%comp_phi
      call magdif_single
      call pack_dof(Bn, xnew)
    end subroutine next_iteration_arnoldi
  end subroutine magdif_iterate

  subroutine send_flag_to_freefem(flag, namedpipe)
    use iso_c_binding, only: c_long
    integer, intent(in) :: flag
    character(len = *), intent(in) :: namedpipe
    integer(c_long) :: long_flag
    integer :: fid

    long_flag = int(flag, c_long)
    open(newunit = fid, file = namedpipe, status = 'old', access = 'stream', &
         form = 'unformatted', action = 'write')
    write (fid) long_flag
    close(fid)
  end subroutine send_flag_to_freefem

  subroutine receive_flag_from_freefem(flag, namedpipe)
    use iso_c_binding, only: c_long
    integer, intent(out) :: flag
    character(len = *), intent(in) :: namedpipe
    integer(c_long) :: long_flag
    integer :: fid

    open(newunit = fid, file = namedpipe, status = 'old', access = 'stream', &
         form = 'unformatted', action = 'read')
    read (fid) long_flag
    close(fid)
    flag = int(long_flag)
  end subroutine receive_flag_from_freefem

  subroutine send_RT0_to_freefem(elem, outfile)
    use iso_c_binding, only: c_long
    use magdif_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
    character(len = *), intent(in) :: outfile
    integer(c_long) :: length
    integer :: ktri, fid

    length = 8 * mesh%ntri  ! (Re, Im) of RT0 DoFs + toroidal component
    ! status = 'old' for writing to named pipe
    open(newunit = fid, file = outfile, access = 'stream', status = 'old', &
         action = 'write', form = 'unformatted')
    write (fid) length
    do ktri = 1, mesh%ntri
       write (fid) elem%DOF(:, ktri), elem%comp_phi(ktri) * mesh%area(ktri)
    end do
    close(fid)
  end subroutine send_RT0_to_freefem

  subroutine receive_RT0_from_freefem(elem, infile)
    use iso_c_binding, only: c_long
    use magdif_conf, only: log
    use magdif_mesh, only: mesh
    type(RT0_t), intent(inout) :: elem
    character(len = *), intent(in) :: infile
    integer(c_long) :: length
    integer :: ktri, fid

    open(newunit = fid, file = infile, access = 'stream', status = 'old', &
         action = 'read', form = 'unformatted')
    read (fid) length
    if (length < 8 * mesh%ntri) then
       ! (Re, Im) of RT0 DoFs + toroidal component
       write (log%msg, '("Pipe ", a, " only contains ", i0, " real values, ' // &
            'expected ", i0, ".")') infile, length, 8 * mesh%ntri
       if (log%err) call log%write
       error stop
    end if
    do ktri = 1, mesh%ntri
       read (fid) elem%DOF(:, ktri), elem%comp_phi(ktri)
    end do
    elem%comp_phi = elem%comp_phi / mesh%area
    close(fid)
  end subroutine receive_RT0_from_freefem

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi via an external program. No data
  !> is read yet; this is done by read_bn().
  !>
  !> Currently, this subroutine Calls FreeFem++ via shell script maxwell.sh in the current
  !> directory and handles the exit code. For further information see maxwell.sh and the
  !> script called therein.
  subroutine compute_Bn
    use iso_c_binding, only: c_long
    use magdif_conf, only: conf, decorate_filename
    use magdif_mesh, only: mesh
    use magdif_pert, only: RT0_check_redundant_edges, RT0_check_div_free
    integer :: dummy

    call send_flag_to_freefem(-1, freefem_pipe)
    call receive_flag_from_freefem(dummy, freefem_pipe)
    call send_RT0_to_freefem(jn, freefem_pipe)
    call receive_RT0_from_freefem(Bn, freefem_pipe)
    call RT0_check_redundant_edges(Bn, 'Bn')
    call RT0_check_div_free(Bn, mesh%n, conf%rel_err_Bn, 'Bn')
  end subroutine compute_Bn

  subroutine compute_L2int(elem, integral)
    use iso_c_binding, only: c_long
    use magdif_conf, only: log
    type(RT0_t), intent(in) :: elem
    real(dp), intent(out) :: integral
    integer(c_long) :: length
    integer :: fid, dummy

    call send_flag_to_freefem(-2, freefem_pipe)
    call receive_flag_from_freefem(dummy, freefem_pipe)
    call send_RT0_to_freefem(elem, freefem_pipe)
    open(newunit = fid, file = freefem_pipe, access = 'stream', status = 'old', &
         action = 'read', form = 'unformatted')
    read (fid) length
    if (length /= 1) then
       close(fid)
       write (log%msg, '("Expected 1 double value from FreeFem++, but got ", i0)') length
       if (log%err) call log%write
       error stop
    end if
    read (fid) integral
    close(fid)
  end subroutine compute_L2int

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

  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  subroutine compute_presn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use magdif_conf, only: conf, log
    use magdif_util, only: imun
    use magdif_mesh, only: fs, mesh, B0R, B0phi, B0Z
    real(dp) :: r
    real(dp) :: lr, lz  ! edge vector components
    complex(dp), dimension(conf%nkpol) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(conf%nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, ktri, kp
    integer :: nz
    integer, dimension(2 * conf%nkpol) :: irow, icol
    complex(dp), dimension(2 * conf%nkpol) :: aval
    integer :: base, tip
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux
       inner: do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          if (.not. mesh%orient(ktri)) cycle inner
          ! use midpoint of edge f
          base = mesh%lf(1, ktri)
          tip = mesh%lf(2, ktri)
          R = (mesh%node_R(base) + mesh%node_R(tip)) * 0.5d0
          lR = mesh%node_R(tip) - mesh%node_R(base)
          lZ = mesh%node_Z(tip) - mesh%node_Z(base)

          kp = base - mesh%kp_low(kf)
          a(kp) = (B0R(mesh%ef(ktri), ktri) * lR + B0Z(mesh%ef(ktri), ktri) * lZ) / &
               (lR ** 2 + lZ ** 2)
          x(kp) = -fs%dp_dpsi(kf) * a(kp) * Bn%DOF(mesh%ef(ktri), ktri)
          b(kp) = imun * (mesh%n + imun * conf%damp) * B0phi(mesh%ef(ktri), ktri) / R
       end do inner

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       call assemble_sparse(conf%nkpol, d, du, nz, irow, icol, aval)
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(conf%nkpol, conf%nkpol, nz, irow, icol, aval, x)
       call sparse_matmul(conf%nkpol, conf%nkpol, irow, icol, aval, x, resid)
       resid = resid - inhom
       where (abs(inhom) >= small)
          rel_err = abs(resid) / abs(inhom)
       elsewhere
          rel_err = 0d0
       end where
       max_rel_err = max(max_rel_err, maxval(rel_err))
       avg_rel_err = avg_rel_err + sum(rel_err)

       if (kf == 1) then ! first point on axis - average over enclosing flux surface
          pn%DOF(1) = sum(x) / size(x)
       end if
       do kp = 1, mesh%kp_max(kf)
          pn%DOF(mesh%kp_low(kf) + kp) = x(kp)
       end do
    end do

    avg_rel_err = avg_rel_err / sum(mesh%kp_max(1:mesh%nflux))
    write (log%msg, '("compute_presn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (log%debug) call log%write
    if (allocated(resid)) deallocate(resid)
  end subroutine compute_presn

  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #magdif_conf::currn_file.
  subroutine compute_currn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use magdif_conf, only: conf, log, decorate_filename
    use magdif_mesh, only: fs, mesh, B0phi, B0flux, j0phi
    use magdif_pert, only: RT0_check_div_free, RT0_check_redundant_edges
    use magdif_util, only: imun, clight
    complex(dp), dimension(2 * conf%nkpol) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(2 * conf%nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, ktri, ke
    integer :: nz
    integer, dimension(4 * conf%nkpol) :: irow, icol
    complex(dp), dimension(4 * conf%nkpol) :: aval
    complex(dp) :: Bnphi_Gamma
    real(dp) :: R, area
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux ! loop through flux surfaces
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          area = mesh%area(ktri)
          ! first term on source side: flux through edge f
          ke = mesh%ef(ktri)
          R = sum(mesh%node_R(mesh%lf(:, ktri))) * 0.5d0
          jn%DOF(ke, ktri) = j0phi(ke, ktri) / B0phi(ke, ktri) * Bn%DOF(ke, ktri) + &
               clight * R / B0phi(ke, ktri) * (pn%DOF(mesh%lf(2, ktri)) - pn%DOF(mesh%lf(1, ktri)))
          x(kt) = -jn%DOF(ke, ktri)
          ! diagonal matrix element - edge i
          ke = mesh%ei(ktri)
          d(kt) = -1d0 - imun * (mesh%n + imun * conf%damp) * area * 0.5d0 * &
               B0phi(ke, ktri) / B0flux(ke, ktri)
          ! additional term from edge i on source side
          R = sum(mesh%node_R(mesh%li(:, ktri))) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bn%comp_phi(ktri) + Bn%comp_phi(mesh%adj_tri(ke, ktri)))
          x(kt) = x(kt) - imun * mesh%n * area * 0.5d0 * (clight * R / B0flux(ke, ktri) * &
               (Bnphi_Gamma / B0phi(ke, ktri) * (fs%p(kf) - fs%p(kf-1)) - &
               (pn%DOF(mesh%li(2, ktri)) - pn%DOF(mesh%li(1, ktri)))) + j0phi(ke, ktri) * &
               (Bnphi_Gamma / B0phi(ke, ktri) - Bn%DOF(ke, ktri) / B0flux(ke, ktri)))
          ! superdiagonal matrix element - edge o
          ke = mesh%eo(ktri)
          du(kt) = 1d0 + imun * (mesh%n + imun * conf%damp) * area * 0.5d0 * &
               B0phi(ke, ktri) / B0flux(ke, ktri)
          ! additional term from edge o on source side
          R = sum(mesh%node_R(mesh%lo(:, ktri))) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bn%comp_phi(ktri) + Bn%comp_phi(mesh%adj_tri(ke, ktri)))
          x(kt) = x(kt) - imun * mesh%n * area * 0.5d0 * (clight * r / B0flux(ke, ktri) * &
               (Bnphi_Gamma / B0phi(ke, ktri) * (fs%p(kf-1) - fs%p(kf)) - &
               (pn%DOF(mesh%lo(2, ktri)) - pn%DOF(mesh%lo(1, ktri)))) + j0phi(ke, ktri) * &
               (Bnphi_Gamma / B0phi(ke, ktri) - Bn%DOF(ke, ktri) / B0flux(ke, ktri)))
       end do
       associate (ndim => mesh%kt_max(kf))
         call assemble_sparse(ndim, d(:ndim), du(:ndim), nz, &
              irow(:2*ndim), icol(:2*ndim), aval(:2*ndim))
         inhom = x  ! remember inhomogeneity before x is overwritten with the solution
         call sparse_solve(ndim, ndim, nz, irow(:nz), icol(:nz), aval(:nz), x(:ndim))
         call sparse_matmul(ndim, ndim, irow(:nz), icol(:nz), aval(:nz), x(:ndim), resid)
         resid = resid - inhom(:ndim)
         where (abs(inhom(:ndim)) >= small)
            rel_err(:ndim) = abs(resid(:ndim)) / abs(inhom(:ndim))
         elsewhere
            rel_err(:ndim) = 0d0
         end where
         max_rel_err = max(max_rel_err, maxval(rel_err(:ndim)))
         avg_rel_err = avg_rel_err + sum(rel_err(:ndim))
       end associate
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          jn%DOF(mesh%ei(ktri), ktri) = -x(kt)
          jn%DOF(mesh%eo(ktri), ktri) = x(mod(kt, mesh%kt_max(kf)) + 1)
          jn%comp_phi(ktri) = sum(jn%DOF(:, ktri)) * imun / mesh%n / mesh%area(ktri)
       end do
    end do
    avg_rel_err = avg_rel_err / sum(mesh%kt_max(1:mesh%nflux))
    write (log%msg, '("compute_currn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (log%debug) call log%write
    if (allocated(resid)) deallocate(resid)

    call add_sheet_current
    call RT0_check_redundant_edges(jn, 'jn')
    call RT0_check_div_free(jn, mesh%n, conf%rel_err_currn, 'jn')
  end subroutine compute_currn

  subroutine add_sheet_current
    use magdif_conf, only: conf_arr
    use magdif_util, only: imun
    use magdif_mesh, only: mesh, B0flux
    integer :: kf, kt, ktri, ktri_adj
    complex(dp) :: pn_half

    do kf = 1, mesh%nflux
       if (mesh%m_res(kf) > 0) then
          if (abs(conf_arr%sheet_current_factor(mesh%m_res(kf))) > 0d0) then
             do kt = 1, mesh%kt_max(kf)
                ktri = mesh%kt_low(kf) + kt
                ! add sheet current on edge i
                if (mod(kt, 2) == 0) then
                   ! edge i is diagonal
                   if (mesh%orient(ktri)) then
                      pn_half = 0.5d0 * sum(pn%DOF(mesh%lf(:, ktri)))
                   else
                      ktri_adj = mesh%adj_tri(mesh%ei(ktri), ktri)
                      pn_half = 0.5d0 * sum(pn%DOF(mesh%lf(:, ktri_adj)))
                   end if
                else
                   pn_half = pn%DOF(mesh%li(2, ktri))
                end if
                jn%DOF(mesh%ei(ktri), ktri) = jn%DOF(mesh%ei(ktri), ktri) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * &
                     B0flux(mesh%ei(ktri), ktri) * pn_half
                ! add sheet current on edge o
                if (mod(kt, 2) == 1) then
                   ! edge o is diagonal
                   if (mesh%orient(ktri)) then
                      pn_half = 0.5d0 * sum(pn%DOF(mesh%lf(:, ktri)))
                   else
                      ktri_adj = mesh%adj_tri(mesh%eo(ktri), ktri)
                      pn_half = 0.5d0 * sum(pn%DOF(mesh%lf(:, ktri_adj)))
                   end if
                else
                   pn_half = pn%DOF(mesh%lo(1, ktri))
                end if
                jn%DOF(mesh%eo(ktri), ktri) = jn%DOF(mesh%eo(ktri), ktri) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * &
                     B0flux(mesh%eo(ktri), ktri) * pn_half
                ! adjust toroidal current density
                jn%comp_phi(ktri) = sum(jn%DOF(:, ktri)) * imun / mesh%n / mesh%area(ktri)
             end do
          end if
       end if
    end do
  end subroutine add_sheet_current

  subroutine magdif_postprocess
    use magdata_in_symfluxcoor_mod, only: ntheta
    use magdif_conf, only: conf, datafile
    use magdif_mesh, only: mesh, coord_cache_ext, coord_cache_ext_init, coord_cache_ext_deinit, &
         compute_sample_Ipar
    use magdif_pert, only: vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, &
         vec_polmodes_write, RT0_poloidal_modes
    integer, parameter :: m_max = 24
    integer :: k
    type(vec_polmodes_t) :: vec_polmodes
    type(coord_cache_ext) :: sample_Ipar

    ! poloidal modes
    call vec_polmodes_init(vec_polmodes, m_max, mesh%nflux)
    call RT0_poloidal_modes(Bn, vec_polmodes)
    call vec_polmodes_write(vec_polmodes, datafile, 'postprocess/Bmn', &
         'poloidal modes of magnetic field (full perturbation)', 'G')
    call RT0_poloidal_modes(Bnvac, vec_polmodes)
    call vec_polmodes_write(vec_polmodes, datafile, 'postprocess/Bmn_vac', &
         'poloidal modes of magnetic field (vacuum perturbation)', 'G')
    call RT0_poloidal_modes(Bnplas, vec_polmodes)
    call vec_polmodes_write(vec_polmodes, datafile, 'postprocess/Bmn_plas', &
         'poloidal modes of magnetic field (plasma response)', 'G')
    if (conf%kilca_scale_factor /= 0) then
       call check_furth(jn, vec_polmodes)
    end if
    call RT0_poloidal_modes(jn, vec_polmodes)
    call vec_polmodes_write(vec_polmodes, datafile, 'postprocess/jmn', &
         'poloidal modes of current density perturbation', 'statA cm^-2')
    call vec_polmodes_deinit(vec_polmodes)
    ! parallel currents
    call coord_cache_ext_init(sample_Ipar, conf%nrad_Ipar, ntheta)
    do k = lbound(mesh%res_modes, 1), ubound(mesh%res_modes, 1)
       call compute_sample_Ipar(sample_Ipar, mesh%res_modes(k))
       call write_Ipar_symfluxcoord(sample_Ipar)
       if (conf%kilca_scale_factor /= 0) then
          call write_Ipar(sample_Ipar)
       end if
    end do
    call coord_cache_ext_deinit(sample_Ipar)
  end subroutine magdif_postprocess

  subroutine check_furth(jn, Bmn_plas)
    use magdif_conf, only: conf, longlines
    use magdif_util, only: imun, clight
    use magdif_mesh, only: equil, fs_half, mesh, sample_polmodes
    use magdif_pert, only: RT0_t, vec_polmodes_t
    type(RT0_t), intent(in) :: jn
    type(vec_polmodes_t), intent(in) :: Bmn_plas
    integer :: kf, kt, ktri, ktri_eff, kilca_m_res, fid_furth
    complex(dp) :: sheet_flux
    real(dp) :: k_z, k_theta

    kilca_m_res = -equil%cocos%sgn_q * abs(conf%kilca_pol_mode)
    k_z = mesh%n / mesh%R_O
    open(newunit = fid_furth, file = 'check_furth.dat', recl = longlines, status = 'replace')
    do kf = 1, mesh%nflux
       sheet_flux = 0d0
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ktri_eff = sample_polmodes%ktri(ktri)
          sheet_flux = sheet_flux + mesh%area(ktri_eff) * jn%comp_phi(ktri_eff) * &
               exp(-imun * kilca_m_res * sample_polmodes%theta(ktri))
       end do
       k_theta = kilca_m_res / fs_half%rad(kf)
       sheet_flux = -2d0 * imun / clight / k_theta * sheet_flux
       write (fid_furth, '(7(1x, es24.16e3))') fs_half%rad(kf), k_z, k_theta, &
            Bmn_plas%coeff_rad(-kilca_m_res, kf), sheet_flux
    end do
    close(fid_furth)
  end subroutine check_furth


  !> calculate parallel current (density) on a finer grid
  subroutine write_Ipar_symfluxcoord(s)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, longlines
    use magdif_util, only: imun, clight
    use magdif_mesh, only: coord_cache_ext
    use magdif_pert, only: RT0_interp
    type(coord_cache_ext), intent(in) :: s  ! shorten names
    character(len = *), parameter :: fname_fmt = '("currn_par_", i0, ".dat")'
    character(len = 1024) :: fname
    integer :: krad, kpol, k, fid_jpar
    real(dp) :: B0_psi, dB0R_dpsi, dB0phi_dpsi, dB0Z_dpsi, &
         B0_dB0_dpsi, dhphi2_dpsi, B0_theta, dB0R_dtheta, dB0phi_dtheta, dB0Z_dtheta, &
         B0_dB0_dtheta, dhphi2_dtheta, common_term, dB0theta_dpsi, dB0psi_dtheta, &
         dhphihtheta_dpsi, dhphihpsi_dtheta
    complex(dp) :: jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
         dBnR_dR, dBnR_dZ, dBnZ_dR, dBnZ_dZ, &
         Bn_psi, dBnpsi_dpsi, Bn_theta, Delta_mn, part_int, bndry
    real(dp), dimension(s%nrad) :: rad, psi, I_char
    complex(dp), dimension(s%nrad) :: jmn_par_neg, jmn_par_pos, &
         part_int_neg, part_int_pos, bndry_neg, bndry_pos, Delta_mn_neg, Delta_mn_pos

    psi(:) = s%psi(::s%npol)
    rad(:) = s%rad(::s%npol)
    jmn_par_neg = (0d0, 0d0)
    jmn_par_pos = (0d0, 0d0)
    part_int_neg = (0d0, 0d0)
    part_int_pos = (0d0, 0d0)
    bndry_neg = (0d0, 0d0)
    bndry_pos = (0d0, 0d0)
    I_char = 0d0
    Delta_mn_neg = (0d0, 0d0)
    Delta_mn_pos = (0d0, 0d0)
    write (fname, fname_fmt) s%m
    open(newunit = fid_jpar, file = fname, status = 'replace', recl = longlines)
    do krad = 1, s%nrad
       do kpol = 1, s%npol
          k = (krad - 1) * s%npol + kpol
          associate (ktri => s%ktri(k), R => s%R(k), Z => s%Z(k), theta => s%theta(k), m => s%m, &
               sqrt_g => s%sqrt_g(k), q => s%q(k), dq_dpsi => s%dq_dpsi(k), B0_2 => s%B0_2(k), &
               B0_R => s%B0_R(k),  dB0R_dR => s%dB0R_dR(k), dB0R_dZ => s%dB0R_dZ(k), &
               B0_phi => s%B0_phi(k), dB0phi_dR => s%dB0phi_dR(k), dB0phi_dZ => s%dB0phi_dZ(k), &
               B0_Z => s%B0_Z(k), dB0Z_dR => s%dB0Z_dR(k), dB0Z_dZ => s%dB0Z_dZ(k), &
               dR_dtheta => s%dR_dtheta(k), dR_dpsi => s%dR_dpsi(k), d2R_dpsi_dtheta => s%d2R_dpsi_dtheta(k), &
               dZ_dtheta => s%dZ_dtheta(k), dZ_dpsi => s%dZ_dpsi(k), d2Z_dpsi_dtheta => s%d2Z_dpsi_dtheta(k))
            call RT0_interp(ktri, jn, R, Z, jn_R, jn_Z, jn_phi)
            ! include h^phi in current density
            jn_par = (jn_R * B0_R + jn_Z * B0_Z + jn_phi * B0_phi) * &
                 B0_phi / B0_2 * sqrt_g / R
            jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * exp(imun * abs(m) * theta)
            jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * exp(-imun * abs(m) * theta)
            ! comparison with indirect calculation (Boozer and Nuehrenberg 2006)
            call RT0_interp(ktri, Bn, R, Z, Bn_R, Bn_Z, Bn_phi, &
                 dBnR_dR, dBnR_dZ, dBnZ_dR, dBnZ_dZ)
            dB0phi_dpsi = (dB0phi_dR * dR_dpsi + dB0phi_dZ * dZ_dpsi) / R - &
                 B0_phi * dR_dpsi / R ** 2
            Bn_psi = R * (Bn_R * B0_Z - Bn_Z * B0_R)
            dBnpsi_dpsi = dR_dpsi * (Bn_R * B0_Z - Bn_Z * B0_R) + R * ( &
                 Bn_R * (dB0Z_dR * dR_dpsi + dB0Z_dZ * dZ_dpsi) - &
                 Bn_Z * (dB0R_dR * dR_dpsi + dB0R_dZ * dZ_dpsi) + &
                 (dBnR_dR * dR_dpsi + dBnR_dZ * dZ_dpsi) * B0_Z - &
                 (dBnZ_dR * dR_dpsi + dBnZ_dZ * dZ_dpsi) * B0_R)
            Delta_mn = dq_dpsi / q * ((Bn_psi + dBnpsi_dpsi) / B0_phi * R - &
                 Bn_psi * dB0phi_dpsi / B0_phi ** 2 * R ** 2)
            Delta_mn_neg(krad) = Delta_mn_neg(krad) + Delta_mn * exp(imun * abs(m) * theta)
            Delta_mn_pos(krad) = Delta_mn_pos(krad) + Delta_mn * exp(-imun * abs(m) * theta)
            I_char(krad) = I_char(krad) + B0_2 / &
                 ((B0_R ** 2 + B0_Z ** 2) * q * q * R * B0_phi)
            ! comparison with indirect calculation (Ampere's law and integration by parts)
            B0_psi = B0_R * dR_dpsi + B0_Z * dZ_dpsi  ! covariant component
            dB0R_dpsi = dB0R_dR * dR_dpsi + dB0R_dZ * dZ_dpsi
            dB0phi_dpsi = dB0phi_dR * dR_dpsi + dB0phi_dZ * dZ_dpsi
            dB0Z_dpsi = dB0Z_dR * dR_dpsi + dB0Z_dZ * dZ_dpsi
            B0_dB0_dpsi = B0_R * dB0R_dpsi + B0_phi * dB0phi_dpsi + B0_Z * dB0Z_dpsi
            dhphi2_dpsi = 2d0 * B0_phi * (dB0phi_dpsi / B0_2 - &
                 B0_phi * B0_dB0_dpsi / B0_2 ** 2)
            B0_theta = B0_R * dR_dtheta + B0_Z * dZ_dtheta  ! covariant component
            dB0R_dtheta = dB0R_dR * dR_dtheta + dB0R_dZ * dZ_dtheta
            dB0phi_dtheta = dB0phi_dR * dR_dtheta + dB0phi_dZ * dZ_dtheta
            dB0Z_dtheta = dB0Z_dR * dR_dtheta + dB0Z_dZ * dZ_dtheta
            B0_dB0_dtheta = B0_R * dB0R_dtheta + B0_phi * dB0phi_dtheta + B0_Z * dB0Z_dtheta
            dhphi2_dtheta = 2d0 * B0_phi * (dB0phi_dtheta / B0_2 - &
                 B0_phi * B0_dB0_dtheta / B0_2 ** 2)
            common_term = B0_R * d2R_dpsi_dtheta + B0_Z * d2Z_dpsi_dtheta + &
                 dB0R_dR * dR_dpsi * dR_dtheta + dB0Z_dZ * dZ_dpsi * dZ_dtheta
            dB0theta_dpsi = common_term + &
                 dB0R_dZ * dZ_dpsi * dR_dtheta + dB0Z_dR * dR_dpsi * dZ_dtheta
            dB0psi_dtheta = common_term + &
                 dB0R_dZ * dZ_dtheta * dR_dpsi + dB0Z_dR * dR_dtheta * dZ_dpsi
            dhphihtheta_dpsi = ((dB0phi_dpsi / R - B0_phi * dR_dpsi / R ** 2) * B0_theta + &
                 B0_phi / R * dB0theta_dpsi) / B0_2 - 2d0 * B0_phi / R * B0_theta * &
                 B0_dB0_dpsi / B0_2 ** 2
            dhphihpsi_dtheta = ((dB0phi_dtheta / R - B0_phi * dR_dtheta / R ** 2) * B0_psi + &
                 B0_phi / R * dB0psi_dtheta) / B0_2 - 2d0 * B0_phi / R * B0_psi * &
                 B0_dB0_dtheta / B0_2 ** 2
            Bn_psi = Bn_R * dR_dpsi + Bn_Z * dZ_dpsi
            Bn_theta = Bn_R * dR_dtheta + Bn_Z * dZ_dtheta
            part_int = dhphi2_dpsi * Bn_theta - dhphi2_dtheta * Bn_psi + Bn_phi / R * &
                 (dhphihtheta_dpsi - dhphihpsi_dtheta) + imun * conf%n * B0_phi / R * &
                 (B0_psi * Bn_theta - B0_theta * Bn_psi) / B0_2
            part_int_neg(krad) = part_int_neg(krad) + (part_int - imun * abs(m) * &
                 B0_phi * (B0_phi * Bn_psi - B0_psi * Bn_phi)) * &
                 exp(imun * abs(m) * theta)
            part_int_pos(krad) = part_int_pos(krad) + (part_int + imun * abs(m) * &
                 B0_phi * (B0_phi * Bn_psi - B0_psi * Bn_phi)) * &
                 exp(-imun * abs(m) * theta)
            bndry = B0_phi * (Bn_phi * B0_theta - Bn_theta * B0_phi) / B0_2
            bndry_neg(krad) = bndry_neg(krad) + bndry * exp(imun * abs(m) * theta)
            bndry_pos(krad) = bndry_pos(krad) + bndry * exp(-imun * abs(m) * theta)
          end associate
       end do
       jmn_par_neg(krad) = jmn_par_neg(krad) / dble(s%npol)
       jmn_par_pos(krad) = jmn_par_pos(krad) / dble(s%npol)
       part_int_neg(krad) = part_int_neg(krad) / dble(s%npol) * 0.25d0 * clight / pi
       part_int_pos(krad) = part_int_pos(krad) / dble(s%npol) * 0.25d0 * clight / pi
       bndry_neg(krad) = bndry_neg(krad) / dble(s%npol) * 0.25d0 * clight / pi
       bndry_pos(krad) = bndry_pos(krad) / dble(s%npol) * 0.25d0 * clight / pi
       Delta_mn_neg(krad) = Delta_mn_neg(krad) / dble(s%npol)
       Delta_mn_pos(krad) = Delta_mn_pos(krad) / dble(s%npol)
       I_char(krad) = dble(s%npol) * 0.5d0 * clight / I_char(krad)
       write (fid_jpar, '(19(1x, es24.16e3))') psi(krad), rad(krad), &
            jmn_par_neg(krad), jmn_par_pos(krad), &
            part_int_neg(krad), part_int_pos(krad), bndry_neg(krad), bndry_pos(krad), &
            I_char(krad), Delta_mn_neg(krad), Delta_mn_pos(krad)
    end do
    close(fid_jpar)
  end subroutine write_Ipar_symfluxcoord


  !> calculate parallel current (density) on a finer grid
  subroutine write_Ipar(s)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: longlines, decorate_filename
    use magdif_util, only: imun, clight, bent_cyl2straight_cyl
    use magdif_mesh, only: mesh, coord_cache_ext
    use magdif_pert, only: RT0_interp
    type(coord_cache_ext), intent(in) :: s  ! shorten names
    integer :: krad, kpol, k, fid_jpar
    real(dp) :: dB0R_drad, dB0phi_drad, dB0Z_drad, B0_dB0_drad, &
         B0_theta, dB0theta_drad, dhz2_drad, dradhthetahz_drad
    complex(dp) :: jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
         Bn_rad, Bn_pol, Bn_tor, part_int, bndry
    real(dp), dimension(s%nrad) :: rad, psi
    complex(dp), dimension(s%nrad) :: jmn_par_neg, jmn_par_pos, &
         part_int_neg, part_int_pos, bndry_neg, bndry_pos

    psi(:) = s%psi(::s%npol)
    rad(:) = s%rad(::s%npol)
    jmn_par_neg = (0d0, 0d0)
    jmn_par_pos = (0d0, 0d0)
    part_int_neg = (0d0, 0d0)
    part_int_pos = (0d0, 0d0)
    bndry_neg = (0d0, 0d0)
    bndry_pos = (0d0, 0d0)
    open(newunit = fid_jpar, status = 'replace', recl = longlines, file = 'currn_par.dat')
    do krad = 1, s%nrad
       do kpol = 1, s%npol
          k = (krad - 1) * s%npol + kpol
          associate (ktri => s%ktri(k), R => s%R(k), Z => s%Z(k), &
               theta => s%theta(k), m => s%m, B0_2 => s%B0_2(k), &
               B0_R => s%B0_R(k), dB0R_dR => s%dB0R_dR(k), dB0R_dZ => s%dB0R_dZ(k), &
               B0_phi => s%B0_phi(k), dB0phi_dR => s%dB0phi_dR(k), dB0phi_dZ => s%dB0phi_dZ(k), &
               B0_Z => s%B0_Z(k), dB0Z_dR => s%dB0Z_dR(k), dB0Z_dZ => s%dB0Z_dZ(k))
            call RT0_interp(ktri, jn, R, Z, jn_R, jn_Z, jn_phi)
            ! include h^z in current density
            jn_par = (jn_R * B0_R + jn_Z * B0_Z + jn_phi * B0_phi) * B0_phi / B0_2
            jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * exp(imun * abs(m) * theta)
            jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * exp(-imun * abs(m) * theta)
            ! comparison with indirect calculation (Ampere's law and integration by parts)
            dB0R_drad = dB0R_dR * cos(theta) + dB0R_dZ * sin(theta)
            dB0phi_drad = dB0phi_dR * cos(theta) + dB0phi_dZ * sin(theta)
            dB0Z_drad = dB0Z_dR * cos(theta) + dB0Z_dZ * sin(theta)
            B0_dB0_drad = B0_R * dB0R_drad + B0_phi * dB0phi_drad + B0_Z * dB0Z_drad
            B0_theta = B0_Z * cos(theta) - B0_R * sin(theta)
            dB0theta_drad = dB0Z_dR * cos(theta) ** 2 - dB0R_dZ * sin(theta) ** 2 + &
                 (dB0Z_dZ - dB0R_dR) * sin(theta) * cos(theta)
            dhz2_drad = 2d0 * B0_phi * (B0_2 * dB0phi_drad - B0_phi * B0_dB0_drad) / &
                 B0_2 ** 2
            dradhthetahz_drad = B0_theta * B0_phi / B0_2 + rad(krad) * &
                 ((B0_theta * dB0phi_drad + dB0theta_drad * B0_phi) / B0_2 - &
                 2d0 * B0_theta * B0_phi * B0_dB0_drad / B0_2 ** 2)
            call RT0_interp(ktri, Bn, R, Z, Bn_R, Bn_Z, Bn_phi)
            call bent_cyl2straight_cyl(Bn_R, Bn_phi, Bn_Z, theta, &
                 Bn_rad, Bn_pol, Bn_tor)
            part_int = -rad(krad) * Bn_pol * dhz2_drad + Bn_tor * dradhthetahz_drad
            part_int_neg(krad) = part_int_neg(krad) + (part_int + imun * B0_phi * &
                 (dble(mesh%n) / mesh%R_O * rad(krad) * B0_theta + &
                 abs(m) * B0_phi) * Bn_rad / B0_2) * exp(imun * abs(m) * theta)
            part_int_pos(krad) = part_int_pos(krad) + (part_int + imun * B0_phi * &
                 (dble(mesh%n) / mesh%R_O * rad(krad) * B0_theta - &
                 abs(m) * B0_phi) * Bn_rad / B0_2) * exp(-imun * abs(m) * theta)
            bndry = B0_phi * rad(krad) * (B0_phi * Bn_pol - B0_theta * Bn_tor) / B0_2
            bndry_neg(krad) = bndry_neg(krad) + bndry * exp(imun * abs(m) * theta)
            bndry_pos(krad) = bndry_pos(krad) + bndry * exp(-imun * abs(m) * theta)
          end associate
       end do
       jmn_par_neg(krad) = jmn_par_neg(krad) / dble(s%npol)
       jmn_par_pos(krad) = jmn_par_pos(krad) / dble(s%npol)
       part_int_neg(krad) = part_int_neg(krad) / dble(s%npol) * 0.25d0 * clight / pi
       part_int_pos(krad) = part_int_pos(krad) / dble(s%npol) * 0.25d0 * clight / pi
       bndry_neg(krad) = bndry_neg(krad) / dble(s%npol) * 0.25d0 * clight / pi
       bndry_pos(krad) = bndry_pos(krad) / dble(s%npol) * 0.25d0 * clight / pi
       write (fid_jpar, '(14(1x, es24.16e3))') psi(krad), rad(krad), &
            jmn_par_neg(krad), jmn_par_pos(krad), &
            part_int_neg(krad), part_int_pos(krad), bndry_neg(krad), bndry_pos(krad)
    end do
    close(fid_jpar)
  end subroutine write_Ipar
end module magdif
