module magdif
  use iso_fortran_env, only: dp => real64
  use magdif_pert, only: L1_t, RT0_t

  implicit none

  private

  public :: magdif_run, magdif_init, magdif_cleanup, magdif_iterate, magdif_postprocess

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

  interface
     subroutine FEM_init(tormode, runmode) bind(C, name = 'FEM_init')
       use iso_c_binding, only: c_int
       integer(c_int), intent(in), value :: tormode, runmode
     end subroutine FEM_init

     subroutine FEM_extend_mesh() bind(C, name = 'FEM_extend_mesh')
     end subroutine FEM_extend_mesh

     subroutine FEM_compute_Bn(shape, Jn, Bn) bind(C, name = 'FEM_compute_Bn')
       use iso_c_binding, only: c_int, c_double_complex
       integer(c_int), intent(in) :: shape(2)
       complex(c_double_complex), intent(in) :: Jn(1:shape(1), 1:shape(2))
       complex(c_double_complex), intent(out) :: Bn(1:shape(1), 1:shape(2))
     end subroutine FEM_compute_Bn

     subroutine FEM_compute_L2int(shape, elem, L2int) bind(C, name = 'FEM_compute_L2int')
       use iso_c_binding, only: c_int, c_double_complex, c_double
       integer(c_int), intent(in) :: shape(2)
       complex(c_double_complex), intent(in) :: elem(1:shape(1), 1:shape(2))
       real(c_double), intent(out) :: L2int
     end subroutine FEM_compute_L2int

     subroutine FEM_deinit() bind(C, name = 'FEM_deinit')
     end subroutine FEM_deinit
  end interface

contains

  subroutine magdif_run(runmode, config) bind(C, name = 'magdif_run')
    use iso_c_binding, only: c_int, c_ptr
    use magdata_in_symfluxcoor_mod, only: load_magdata_in_symfluxcoord, unload_magdata_in_symfluxcoord
    use magdif_util, only: C_F_string, get_field_filenames, init_field, deinit_field
    use magdif_conf, only: conf, magdif_config_read, magdif_config_export_hdf5, conf_arr, magdif_log, log, datafile
    use magdif_mesh, only: equil, fs, fs_half, mesh, generate_mesh, write_mesh_cache, read_mesh_cache, &
         magdif_mesh_deinit, sample_polmodes, coord_cache_deinit, psi_interpolator, psi_fine_interpolator, &
         B0R_edge, B0phi_edge, B0Z_edge, B0R_Omega, B0phi_Omega, B0Z_Omega, B0_flux, j0phi_edge
    use magdif_pert, only: generate_vacfield
    use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
    integer(c_int), intent(in), value :: runmode
    type(c_ptr), intent(in), value :: config
    integer(c_int) :: runmode_flags
    character(len = 1024) :: config_filename, gfile, pfile, convexfile
    logical :: meshing, analysis, iterations

    meshing = iand(runmode, ishft(1, 0)) /= 0
    iterations = iand(runmode, ishft(1, 1)) /= 0
    analysis = iand(runmode, ishft(1, 2)) /= 0
    runmode_flags = runmode
    if (.not. (meshing .or. iterations .or. analysis)) then
       meshing = .true.
       iterations = .true.
       analysis = .true.
       runmode_flags = ior(ior(ishft(1, 0), ishft(1, 1)), ishft(1, 2))
    end if
    call C_F_string(config, config_filename)
    call magdif_config_read(conf, config_filename)
    log = magdif_log('-', conf%log_level, conf%quiet)
    call h5_init
    h5overwrite = .true.
    call magdif_config_export_hdf5(conf, datafile, 'config')
    if (meshing) then
       ! initialize equilibrium field
       call get_field_filenames(gfile, pfile, convexfile)
       call equil%read(trim(gfile), trim(convexfile))
       call equil%classify
       call equil%standardise
       if (conf%kilca_scale_factor /= 0) then
          call equil%scale(conf%kilca_scale_factor)
       end if
       call equil%export_hdf5(datafile, 'equil')
       call init_field(equil)
       ! generate mesh and vacuum field
       call generate_mesh
       call write_mesh_cache
       call generate_vacfield
       ! pass effective toroidal mode number and runmode to FreeFem++
       call FEM_init(mesh%n, runmode_flags)
       call FEM_extend_mesh
    else
       ! initialize equilibrium field
       call equil%import_hdf5(datafile, 'equil')
       call init_field(equil)
       ! read in preprocessed data
       call read_mesh_cache
       call load_magdata_in_symfluxcoord
       ! reload config parameters here in case they changed since the meshing phase
       call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
       call conf_arr%export_hdf5(datafile, 'config')
       ! pass effective toroidal mode number and runmode to FreeFem++
       call FEM_init(mesh%n, runmode)
    end if
    if (iterations .or. analysis) then
       call magdif_init
       if (iterations) then
          call magdif_iterate
          call FEM_deinit
       else
          if (analysis) then
             call magdif_read
          end if
          ! FEM_deinit is not needed because scripts/maxwell_daemon.edp exits
          ! if iterations are not requested
       end if
       if (analysis) then
          call magdif_postprocess
       end if
       call magdif_cleanup
    end if
    if (allocated(B0R_edge)) deallocate(B0R_edge)
    if (allocated(B0phi_edge)) deallocate(B0phi_edge)
    if (allocated(B0Z_edge)) deallocate(B0Z_edge)
    if (allocated(B0r_Omega)) deallocate(B0r_Omega)
    if (allocated(B0phi_Omega)) deallocate(B0phi_Omega)
    if (allocated(B0z_Omega)) deallocate(B0z_Omega)
    if (allocated(B0_flux)) deallocate(B0_flux)
    if (allocated(j0phi_edge)) deallocate(j0phi_edge)
    call psi_interpolator%deinit
    call psi_fine_interpolator%deinit
    call coord_cache_deinit(sample_polmodes)
    call fs%deinit
    call fs_half%deinit
    call magdif_mesh_deinit(mesh)
    call unload_magdata_in_symfluxcoord
    call deinit_field
    call equil%deinit
    call conf_arr%deinit
    call log%deinit
    call h5_deinit
  end subroutine magdif_run

  !> Initialize magdif module
  subroutine magdif_init
    use magdif_conf, only: conf, log, datafile
    use magdif_mesh, only: mesh
    use magdif_pert, only: RT0_check_redundant_edges, RT0_check_div_free, &
         RT0_init, RT0_read, L1_init

    ! initialize perturbation
    call L1_init(pn, mesh%npoint)
    call RT0_init(Bn, mesh%ntri)
    call RT0_init(Bnvac, mesh%ntri)
    call RT0_init(Bnplas, mesh%ntri)
    call RT0_init(jn, mesh%ntri)
    call RT0_read(Bnvac, datafile, 'Bnvac')
    call RT0_check_redundant_edges(Bnvac, 'Bnvac')
    call RT0_check_div_free(Bnvac, mesh%n, conf%rel_err_Bn, 'Bnvac')

    log%msg = 'magdif initialized'
    if (log%info) call log%write
  end subroutine magdif_init

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    use magdif_conf, only: log
    use magdif_pert, only: L1_deinit, RT0_deinit

    call L1_deinit(pn)
    call RT0_deinit(Bn)
    call RT0_deinit(Bnvac)
    call RT0_deinit(Bnplas)
    call RT0_deinit(jn)
    log%msg = 'magdif cleanup finished'
    if (log%info) call log%write
  end subroutine magdif_cleanup

  subroutine magdif_read
    use magdif_conf, only: datafile
    use magdif_pert, only: L1_read, RT0_read

    call L1_read(pn, datafile, 'iter/pn')
    call RT0_read(Bn, datafile, 'iter/Bn')
    call RT0_read(Bnplas, datafile, 'iter/Bnplas')
    call RT0_read(jn, datafile, 'iter/jn')
  end subroutine magdif_read

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
          deallocate(ipiv)
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
          deallocate(Yr)
       else
          preconditioned = .false.
       end if
    end if

    call vec_polmodes_init(jmn, m_max, mesh%nflux)
    call RT0_init(Bn_diff, mesh%ntri)
    allocate(L2int(0:niter))
    call FEM_compute_L2int(shape(Bnvac%DOF), Bnvac%DOF, L2int(0))
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
       call FEM_compute_L2int(shape(Bn_diff%DOF), Bn_diff%DOF, L2int(kiter))
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
    if (allocated(Lr)) deallocate(Lr)
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
    ! TODO: refactor arnoldi_mod?
    deallocate(L2int, eigvecs)

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
      integer :: kedge
      do kedge = 1, mesh%nedge
         elem%DOF(mesh%edge_map2ke(1, kedge), mesh%edge_map2ktri(1, kedge)) = &
              packed(kedge)
         if (mesh%edge_map2ktri(2, kedge) > 0) then
            elem%DOF(mesh%edge_map2ke(2, kedge), mesh%edge_map2ktri(2, kedge)) = &
                 -packed(kedge)
         end if
      end do
      elem%comp_phi(:) = sum(elem%DOF, 1) * imun / mesh%n / mesh%area
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

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi.
  !>
  !> This subroutine calls a C function that pipes the data to/from FreeFem.
  subroutine compute_Bn
    use magdif_util, only: imun
    use magdif_mesh, only: mesh
    use magdif_pert, only: RT0_check_redundant_edges, RT0_check_div_free

    call FEM_compute_Bn(shape(jn%DOF), jn%DOF, Bn%DOF)
    call RT0_check_redundant_edges(Bn, 'Bn')
    Bn%comp_phi(:) = imun / mesh%n * sum(Bn%DOF, 1) / mesh%area
  end subroutine compute_Bn

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
    use magdif_mesh, only: fs, mesh, B0R_edge, B0phi_edge, B0Z_edge
    real(dp) :: r
    real(dp) :: lr, lz  ! edge vector components
    complex(dp), dimension(conf%nkpol) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(conf%nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kp, kedge, ktri
    integer :: nz
    integer, dimension(2 * conf%nkpol) :: irow, icol
    complex(dp), dimension(2 * conf%nkpol) :: aval
    integer :: base, tip
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ! use midpoint of poloidal edge
          base = mesh%edge_node(1, kedge)
          tip = mesh%edge_node(2, kedge)
          R = (mesh%node_R(base) + mesh%node_R(tip)) * 0.5d0
          lR = mesh%node_R(tip) - mesh%node_R(base)
          lZ = mesh%node_Z(tip) - mesh%node_Z(base)
          ktri = mesh%edge_tri(1, kedge)
          a(kp) = (B0R_edge(kedge) * lR + B0Z_edge(kedge) * lZ) / (lR ** 2 + lZ ** 2)
          x(kp) = -fs%dp_dpsi(kf) * a(kp) * Bn%DOF(mesh%ef(ktri), ktri)
          b(kp) = imun * (mesh%n + imun * conf%damp) * B0phi_edge(kedge) / R
       end do

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       call assemble_sparse(conf%nkpol, d, du, nz, irow, icol, aval)
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(conf%nkpol, conf%nkpol, nz, irow, icol, aval, x)
       call sparse_matmul(conf%nkpol, conf%nkpol, irow, icol, aval, x, resid)
       resid(:) = resid - inhom
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
    use magdif_mesh, only: fs, mesh, B0phi_edge, B0_flux, j0phi_edge
    use magdif_pert, only: RT0_check_div_free, RT0_check_redundant_edges
    use magdif_util, only: imun, clight
    complex(dp), dimension(2 * conf%nkpol) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(2 * conf%nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kp, kt, ktri, ke, kedge
    integer :: nz
    integer, dimension(4 * conf%nkpol) :: irow, icol
    complex(dp), dimension(4 * conf%nkpol) :: aval
    complex(dp) :: Bnphi_edge, Delta_pn
    real(dp) :: R, Delta_p0
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ktri = mesh%edge_tri(1, kedge)
          ke = mesh%edge_map2ke(1, kedge)
          R = sum(mesh%node_R(mesh%edge_node(:, kedge))) * 0.5d0
          jn%DOF(ke, ktri) = j0phi_edge(kedge) / B0phi_edge(kedge) * Bn%DOF(ke, ktri) + &
               clight * R / B0phi_edge(kedge) * (pn%DOF(mesh%edge_node(2, kedge)) - pn%DOF(mesh%edge_node(1, kedge)))
          if (mesh%edge_tri(2, kedge) > 0) then
             jn%DOF(mesh%edge_map2ke(2, kedge), mesh%edge_tri(2, kedge)) = -jn%DOF(ke, ktri)
          end if
       end do
    end do
    do kf = 1, mesh%nflux
       Delta_p0 = fs%p(kf) - fs%p(kf-1)
       x = (0d0, 0d0)
       do kt = 1, mesh%kt_max(kf)
          kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
          R = sum(mesh%node_R(mesh%edge_node(:, kedge))) * 0.5d0
          Bnphi_edge = 0.5d0 * sum(Bn%comp_phi(mesh%edge_tri(:, kedge)))
          Delta_pn = pn%DOF(mesh%edge_node(2, kedge)) - pn%DOF(mesh%edge_node(1, kedge))
          ! first term on source side: flux through edge f
          ktri = mesh%edge_tri(2, kedge)
          ke = kt
          x(ke) = x(ke) - jn%DOF(mesh%ef(ktri), ktri)
          ! diagonal matrix element - edge i
          d(ke) = -1d0 - imun * (mesh%n + imun * conf%damp) * &
               mesh%area(ktri) * 0.5d0 * B0phi_edge(kedge) / B0_flux(kedge)
          ! additional term from edge i on source side
          x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * R / B0_flux(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) * Delta_p0 - Delta_pn) + j0phi_edge(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) - Bn%DOF(mesh%ei(ktri), ktri) / B0_flux(kedge)))
          ! superdiagonal matrix element - edge o
          ktri = mesh%edge_tri(1, kedge)
          ke = mod(kt + mesh%kt_max(kf) - 2, mesh%kt_max(kf)) + 1
          du(ke) = 1d0 + imun * (mesh%n + imun * conf%damp) * &
               mesh%area(ktri) * 0.5d0 * B0phi_edge(kedge) / (-B0_flux(kedge))
          ! additional term from edge o on source side
          x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * R / (-B0_flux(kedge)) * &
               (Bnphi_edge / B0phi_edge(kedge) * (-Delta_p0) - (-Delta_pn)) + j0phi_edge(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) - Bn%DOF(mesh%eo(ktri), ktri) / (-B0_flux(kedge))))
       end do
       associate (ndim => mesh%kt_max(kf))
         call assemble_sparse(ndim, d(:ndim), du(:ndim), nz, &
              irow(:2*ndim), icol(:2*ndim), aval(:2*ndim))
         inhom = x  ! remember inhomogeneity before x is overwritten with the solution
         call sparse_solve(ndim, ndim, nz, irow(:nz), icol(:nz), aval(:nz), x(:ndim))
         call sparse_matmul(ndim, ndim, irow(:nz), icol(:nz), aval(:nz), x(:ndim), resid)
         resid(:) = resid - inhom(:ndim)
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
    use magdif_mesh, only: mesh, B0_flux
    integer :: kf, kt, ktri, ktri_adj, kedge
    complex(dp) :: pn_half

    do kf = 1, mesh%nflux
       if (mesh%m_res(kf) > 0) then
          if (abs(conf_arr%sheet_current_factor(mesh%m_res(kf))) > 0d0) then
             do kt = 1, mesh%kt_max(kf)
                kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
                ktri = mesh%edge_tri(2, kedge)
                ktri_adj = mesh%edge_tri(1, kedge)
                if (mesh%orient(ktri)) then
                   pn_half = pn%DOF(mesh%edge_node(2, kedge))
                else
                   ! edge i is diagonal
                   pn_half = 0.5d0 * sum(pn%DOF(mesh%lf(:, ktri_adj)))
                end if
                jn%DOF(mesh%ei(ktri), ktri) = jn%DOF(mesh%ei(ktri), ktri) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * B0_flux(kedge) * pn_half
                jn%DOF(mesh%eo(ktri_adj), ktri_adj) = jn%DOF(mesh%eo(ktri_adj), ktri_adj) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * (-B0_flux(kedge)) * pn_half
             end do
             do kt = 1, mesh%kt_max(kf)
                ktri = mesh%kt_low(kf) + kt
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
    character(len = 24) :: dataset
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
       write (dataset, '("postprocess/Imn_par_", i0)') mesh%res_modes(k)
       call write_Ipar_symfluxcoord(sample_Ipar, datafile, dataset)
       if (conf%kilca_scale_factor /= 0) then
          call write_Ipar(sample_Ipar, datafile, 'postprocess/Imn_par_KiLCA')
       end if
    end do
    call coord_cache_ext_deinit(sample_Ipar)
  end subroutine magdif_postprocess

  subroutine check_furth(jn, Bmn_plas)
    use magdif_conf, only: conf, datafile
    use magdif_util, only: imun, clight
    use magdif_mesh, only: equil, fs_half, mesh, sample_polmodes
    use magdif_pert, only: RT0_t, vec_polmodes_t
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(RT0_t), intent(in) :: jn
    type(vec_polmodes_t), intent(in) :: Bmn_plas
    character(len = *), parameter :: dataset = 'debug_furth'
    integer(HID_T) :: h5id_root
    integer :: kf, kt, ktri, ktri_eff, kilca_m_res
    complex(dp) :: sheet_flux(mesh%nflux)
    real(dp) :: k_z, k_theta(mesh%nflux)

    kilca_m_res = -equil%cocos%sgn_q * abs(conf%kilca_pol_mode)
    k_z = mesh%n / mesh%R_O
    k_theta(:) = kilca_m_res / fs_half%rad
    sheet_flux(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ktri_eff = sample_polmodes%ktri(ktri)
          sheet_flux(kf) = sheet_flux(kf) + mesh%area(ktri_eff) * jn%comp_phi(ktri_eff) * &
               exp(-imun * kilca_m_res * sample_polmodes%theta(ktri))
       end do
    end do
    sheet_flux(:) = -2d0 * imun / clight / k_theta * sheet_flux
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/k_z', k_z, unit = 'cm^-1')
    call h5_add(h5id_root, dataset // '/k_theta', k_theta, &
         lbound(k_theta), ubound(k_theta), unit = 'cm^-1')
    call h5_add(h5id_root, dataset // '/rad', fs_half%rad, &
         lbound(fs_half%rad), ubound(fs_half%rad), unit = 'cm')
    call h5_add(h5id_root, dataset // '/Bmn_plas_rad', Bmn_plas%coeff_rad(-kilca_m_res, :), &
         [1], [mesh%nflux], unit = 'G')
    call h5_add(h5id_root, dataset // '/sheet_flux', sheet_flux, &
         lbound(sheet_flux), ubound(sheet_flux), unit = 'statA')
    call h5_close(h5id_root)
  end subroutine check_furth


  !> calculate parallel current (density) on a finer grid
  subroutine write_Ipar_symfluxcoord(s, file, dataset)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf
    use magdif_util, only: imun, clight
    use magdif_mesh, only: coord_cache_ext
    use magdif_pert, only: RT0_interp
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(coord_cache_ext), intent(in) :: s  ! shorten names
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root
    integer :: krad, kpol, k
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
    end do
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', psi, &
         lbound(psi), ubound(psi))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad', rad, &
         lbound(rad), ubound(rad))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/jmn_par_pos', jmn_par_pos, &
         lbound(jmn_par_pos), ubound(jmn_par_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/jmn_par_neg', jmn_par_neg, &
         lbound(jmn_par_neg), ubound(jmn_par_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_pos', part_int_pos, &
         lbound(part_int_pos), ubound(part_int_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_neg', part_int_neg, &
         lbound(part_int_neg), ubound(part_int_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_pos', part_int_pos, &
         lbound(part_int_pos), ubound(part_int_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_neg', part_int_neg, &
         lbound(part_int_neg), ubound(part_int_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/bndry_pos', bndry_pos, &
         lbound(bndry_pos), ubound(bndry_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/bndry_neg', bndry_neg, &
         lbound(bndry_neg), ubound(bndry_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/I_char', I_char, &
         lbound(I_char), ubound(I_char))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Delta_mn_pos', Delta_mn_pos, &
         lbound(Delta_mn_pos), ubound(Delta_mn_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Delta_mn_neg', Delta_mn_neg, &
         lbound(Delta_mn_neg), ubound(Delta_mn_neg))
  end subroutine write_Ipar_symfluxcoord


  !> calculate parallel current (density) on a finer grid
  subroutine write_Ipar(s, file, dataset)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_util, only: imun, clight, bent_cyl2straight_cyl
    use magdif_mesh, only: mesh, coord_cache_ext
    use magdif_pert, only: RT0_interp
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(coord_cache_ext), intent(in) :: s  ! shorten names
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root
    integer :: krad, kpol, k
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
    end do
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', psi, &
         lbound(psi), ubound(psi))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad', rad, &
         lbound(rad), ubound(rad))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/jmn_par_pos', jmn_par_pos, &
         lbound(jmn_par_pos), ubound(jmn_par_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/jmn_par_neg', jmn_par_neg, &
         lbound(jmn_par_neg), ubound(jmn_par_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_pos', part_int_pos, &
         lbound(part_int_pos), ubound(part_int_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_neg', part_int_neg, &
         lbound(part_int_neg), ubound(part_int_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_pos', part_int_pos, &
         lbound(part_int_pos), ubound(part_int_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/part_int_neg', part_int_neg, &
         lbound(part_int_neg), ubound(part_int_neg))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/bndry_pos', bndry_pos, &
         lbound(bndry_pos), ubound(bndry_pos))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/bndry_neg', bndry_neg, &
         lbound(bndry_neg), ubound(bndry_neg))
  end subroutine write_Ipar
end module magdif
