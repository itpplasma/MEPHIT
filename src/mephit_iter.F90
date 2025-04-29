module mephit_iter
  use iso_fortran_env, only: dp => real64
  use mephit_pert, only: L1_t, RT0_t

  implicit none

  private

  public :: mephit_run, mephit_deinit, perteq_read

  type :: perteq_t
    !> Pressure perturbation \f$ p_{n} \f$ in dyn cm^-1.
    type(L1_t) :: pn

    !> Perturbation field in units of Gauss.
    type(RT0_t) :: Bn

    !> Plasma perturbation field in units of Gauss.
    type(RT0_t) :: Bnplas

    !> Current density perturbation in units of statampere cm^-2.
    type(RT0_t) :: jn

    !> Parallel current density perturbation in units of statampere cm^-2 G^-1.
    type(L1_t) :: jnpar_B0

    !> Vector potential components for GORILLA
    type(L1_t) :: AnR, AnZ

    !> Poloidal modes of electric potential perturbation in units of statV.
    complex(dp), allocatable :: Phi_mn(:, :)

    !> Poloidal modes of aligned electric potential perturbation in units of statV.
    complex(dp), allocatable :: Phi_aligned_mn(:, :)
  end type perteq_t

  type :: precond_t
    integer :: nritz = 0
    complex(dp), allocatable :: Lr(:, :), eigvecs(:, :)
  end type precond_t

  type :: FDM_t
    integer :: nnz = 0
    integer, allocatable :: irow(:), icol(:)
    complex(dp), allocatable :: aval_MDE(:), aval_MDE_damped(:), aval_jnperp(:)
  end type FDM_t

  ! interfaces
  abstract interface
    subroutine real_vector_field(R, Z, vector) bind(C)
      use iso_c_binding, only: c_double
      real(c_double), intent(in), value :: R
      real(c_double), intent(in), value :: Z
      real(c_double), intent(out) :: vector(3)
    end subroutine real_vector_field

    subroutine complex_scalar_field(R, Z, scalar) bind(C)
      use iso_c_binding, only: c_double, c_double_complex
      real(c_double), intent(in), value :: R
      real(c_double), intent(in), value :: Z
      complex(c_double_complex), intent(out) :: scalar
    end subroutine complex_scalar_field
  end interface

  interface
    subroutine FEM_init(tormode, nedge, npoint, runmode) bind(C, name = 'FEM_init')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in), value :: tormode, nedge, npoint, runmode
    end subroutine FEM_init

    subroutine FEM_extend_mesh() bind(C, name = 'FEM_extend_mesh')
    end subroutine FEM_extend_mesh

    subroutine FEM_compute_magfn(nedge, npoint, Jn, Bn, AnR, AnZ) bind(C, name = 'FEM_compute_magfn')
      use iso_c_binding, only: c_int, c_double_complex
      integer(c_int), intent(in), value :: nedge
      integer(c_int), intent(in), value :: npoint
      complex(c_double_complex), intent(in) :: Jn(1:nedge)
      complex(c_double_complex), intent(out) :: Bn(1:nedge)
      complex(c_double_complex), intent(out) :: AnR(1:npoint)
      complex(c_double_complex), intent(out) :: AnZ(1:npoint)
    end subroutine FEM_compute_magfn

    subroutine FEM_compute_L2int(nedge, elem, L2int) bind(C, name = 'FEM_compute_L2int')
      use iso_c_binding, only: c_int, c_double_complex, c_double
      integer(c_int), intent(in), value :: nedge
      complex(c_double_complex), intent(in) :: elem(1:nedge)
      real(c_double), intent(out) :: L2int
    end subroutine FEM_compute_L2int

    subroutine FEM_deinit() bind(C, name = 'FEM_deinit')
    end subroutine FEM_deinit

    function FEM_test(mesh_file, tor_mode, n_dof, dof, unit_B0, MDE_inhom) &
      bind(C, name = 'FEM_test')
      use iso_c_binding, only: c_int, c_char, c_double_complex, c_funptr
      character(c_char), intent(in) :: mesh_file(*)
      integer(c_int), intent(in), value :: tor_mode
      integer(c_int), intent(in), value :: n_dof
      complex(c_double_complex), intent(inout) :: dof(1:n_dof)
      type(c_funptr), intent(in), value :: unit_B0
      type(c_funptr), intent(in), value :: MDE_inhom
      integer(c_int) :: FEM_test
    end function FEM_test
  end interface

contains

  subroutine mephit_run(runmode, config, suffix) bind(C, name = 'mephit_run')
    use iso_c_binding, only: c_int, c_ptr
    use input_files, only: gfile
    use field_sub, only : read_field_input
    use geqdsk_tools, only: geqdsk_read, geqdsk_classify, geqdsk_standardise
    use hdf5_tools, only: h5_init, h5overwrite
    use mephit_util, only: C_F_string, init_field, geqdsk_scale, geqdsk_export_hdf5, geqdsk_import_hdf5, &
      save_symfluxcoord, load_symfluxcoord
    use mephit_conf, only: conf, config_read, config_export_hdf5, conf_arr, logger, &
      datafile, basename_suffix, decorate_filename
    use mephit_mesh, only: equil, mesh, generate_mesh, mesh_write, mesh_read, write_cache, read_cache, &
      resample_profiles, write_profiles_hdf5, read_profiles_hdf5
    use mephit_pert, only: generate_vacfield, vac, vac_init, vac_write, vac_read
    use mephit_flr2, only: flr2_t, flr2_write, flr2_read, flr2_deinit
    integer(c_int), intent(in), value :: runmode
    type(c_ptr), intent(in), value :: config
    type(c_ptr), intent(in), value :: suffix
    character(len = 1024) :: config_filename
    integer(c_int) :: runmode_flags
    logical :: meshing, preconditioner, iterations
    type(fdm_t) :: fdm
    type(flr2_t) :: flr2
    type(perteq_t) :: perteq
    type(precond_t) :: precond

    meshing = iand(runmode, ishft(1, 0)) /= 0
    preconditioner = iand(runmode, ishft(1, 1)) /= 0
    iterations = iand(runmode, ishft(1, 2)) /= 0
    runmode_flags = runmode
    if (.not. (meshing .or. preconditioner .or. iterations)) then
      meshing = .true.
      preconditioner = .true.
      iterations = .true.
      runmode_flags = ior(ior(ishft(1, 0), ishft(1, 1)), ishft(1, 2))
    end if
    call C_F_string(suffix, basename_suffix)
    datafile = decorate_filename(datafile, '', basename_suffix)
    call C_F_string(config, config_filename)
    call config_read(conf, config_filename)
    call logger%init('-', conf%log_level, conf%quiet)
    call h5_init
    h5overwrite = .true.
    call config_export_hdf5(conf, datafile, 'config')
    if (meshing) then
      ! initialize equilibrium field
      call read_field_input
      call geqdsk_read(equil, trim(gfile))
      call geqdsk_classify(equil)
      call geqdsk_standardise(equil)
      if (conf%kilca_scale_factor /= 0) then
        call geqdsk_scale(equil, conf%kilca_scale_factor)
      end if
      call geqdsk_export_hdf5(equil, datafile, 'equil')
      call init_field(equil)
      ! generate mesh and vacuum field
      call generate_mesh
      call save_symfluxcoord(datafile, 'symfluxcoord')
      call mesh_write(mesh, datafile, 'mesh')
      call write_cache
      call resample_profiles
      call write_profiles_hdf5(datafile, 'equil/profiles')
      if (conf%resonance_sweep > 0) then
        call resonance_sweep(datafile, 'resonance_sweep')
      end if
      call vac_init(vac, mesh%nedge, mesh%ntri, mesh%m_res_min, mesh%m_res_max)
      call generate_vacfield(vac)
      call vac_write(vac, datafile, 'vac')
      ! pass effective toroidal mode number and runmode to FreeFem++
      call FEM_init(mesh%n, mesh%nedge, mesh%npoint, runmode_flags)
      call FEM_extend_mesh
    else
      ! initialize equilibrium field
      call read_field_input
      call geqdsk_import_hdf5(equil, datafile, 'equil')
      call init_field(equil)
      ! read in preprocessed data
      call load_symfluxcoord(datafile, 'symfluxcoord')
      call mesh_read(mesh, datafile, 'mesh')
      call read_cache
      call read_profiles_hdf5(datafile, 'equil/profiles')
      call vac_init(vac, mesh%nedge, mesh%ntri, mesh%m_res_min, mesh%m_res_max)
      call vac_read(vac, datafile, 'vac')
      ! reload config parameters here in case they changed since the meshing phase
      call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
      call conf_arr%export_hdf5(datafile, 'config')
      ! pass effective toroidal mode number and runmode to FreeFem++
      call FEM_init(mesh%n, mesh%nedge, mesh%npoint, runmode)
    end if
    if (preconditioner .or. iterations) then
      call perteq_init(perteq)
      if (preconditioner) then
        call FDM_compute_matrix(fdm)
        call FDM_write(fdm, datafile, 'iter/FDM')
        call compute_FLR2_coeff(flr2)
        call FLR2_write(flr2, datafile, 'iter/FLR2')
        call debug_initial_iteration(perteq, fdm, flr2)
        call precond_compute(precond, perteq, fdm, flr2)
        call precond_write(precond, datafile, 'iter')
      else
        call FDM_read(fdm, datafile, 'iter/FDM')
        call FLR2_read(flr2, datafile, 'iter/FLR2')
        call precond_read(precond, datafile, 'iter')
      end if
      if (iterations) then
        call perteq_iterate(perteq, precond, fdm, flr2)
      end if
      call FDM_deinit(fdm)
      call FLR2_deinit(flr2)
      call precond_deinit(precond)
      call perteq_deinit(perteq)
    end if
    call FEM_deinit
    call mephit_deinit
  end subroutine mephit_run

  subroutine mephit_deinit
    use magdata_in_symfluxcoor_mod, only: unload_magdata_in_symfluxcoord
    use hdf5_tools, only: h5_deinit
    use geqdsk_tools, only: geqdsk_deinit
    use mephit_conf, only: conf_arr, logger
    use mephit_util, only: deinit_field
    use mephit_mesh, only: equil, fs, fs_half, psi_fine, &
      mesh, cache, mesh_deinit, cache_deinit, deinit_profiles
    use mephit_pert, only: vac, vac_deinit

    if (allocated(psi_fine)) deallocate(psi_fine)  ! intermediate step, to be removed
    call cache_deinit(cache)
    call fs%deinit
    call fs_half%deinit
    call mesh_deinit(mesh)
    call unload_magdata_in_symfluxcoord
    call vac_deinit(vac)
    call deinit_field
    call deinit_profiles
    call geqdsk_deinit(equil)
    call conf_arr%deinit
    call logger%deinit
    call h5_deinit
  end subroutine mephit_deinit

  subroutine perteq_init(perteq)
    use mephit_conf, only: conf, currn_model_kilca
    use mephit_mesh, only: mesh
    use mephit_pert, only: RT0_init, L1_init
    type(perteq_t), intent(inout) :: perteq

    call L1_init(perteq%pn, mesh%npoint)
    call RT0_init(perteq%Bn, mesh%nedge, mesh%ntri)
    call RT0_init(perteq%Bnplas, mesh%nedge, mesh%ntri)
    call RT0_init(perteq%jn, mesh%nedge, mesh%ntri)
    call L1_init(perteq%jnpar_B0, mesh%npoint)
    call L1_init(perteq%AnR, mesh%npoint)
    call L1_init(perteq%AnZ, mesh%npoint)
    if (conf%currn_model == currn_model_kilca) then
      allocate(perteq%Phi_mn(0:mesh%nflux, mesh%m_res_min:mesh%m_res_max))
      allocate(perteq%Phi_aligned_mn(0:mesh%nflux, mesh%m_res_min:mesh%m_res_max))
    end if
  end subroutine perteq_init

  subroutine perteq_deinit(perteq)
    use mephit_pert, only: L1_deinit, RT0_deinit
    type(perteq_t), intent(inout) :: perteq

    call L1_deinit(perteq%pn)
    call RT0_deinit(perteq%Bn)
    call RT0_deinit(perteq%Bnplas)
    call RT0_deinit(perteq%jn)
    call L1_deinit(perteq%jnpar_B0)
    call L1_deinit(perteq%AnR)
    call L1_deinit(perteq%AnZ)
    if (allocated(perteq%Phi_mn)) deallocate(perteq%Phi_mn)
    if (allocated(perteq%Phi_aligned_mn)) deallocate(perteq%Phi_aligned_mn)
  end subroutine perteq_deinit

  subroutine perteq_read(perteq)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use mephit_conf, only: datafile, conf, currn_model_kilca
    use mephit_pert, only: L1_read, RT0_read
    type(perteq_t), intent(inout) :: perteq
    integer(HID_T) :: h5id_root

    call L1_read(perteq%pn, datafile, 'iter/pn')
    call RT0_read(perteq%Bn, datafile, 'iter/Bn')
    call RT0_read(perteq%Bnplas, datafile, 'iter/Bnplas')
    call RT0_read(perteq%jn, datafile, 'iter/jn')
    call L1_read(perteq%jnpar_B0, datafile, 'iter/jnpar_Bmod')
    call L1_read(perteq%AnR, datafile, 'iter/AnR')
    call L1_read(perteq%AnZ, datafile, 'iter/AnZ')
    if (conf%currn_model == currn_model_kilca) then
      call h5_open(datafile, h5id_root)
      call h5_get(h5id_root, 'iter/Phi_mn', perteq%Phi_mn)
      call h5_get(h5id_root, 'iter/Phi_aligned_mn', perteq%Phi_aligned_mn)
      call h5_close(h5id_root)
    end if
  end subroutine perteq_read

  subroutine precond_deinit(precond)
    type(precond_t), intent(inout) :: precond

    precond%nritz = 0
    if (allocated(precond%Lr)) deallocate(precond%Lr)
    if (allocated(precond%eigvecs)) deallocate(precond%eigvecs)
  end subroutine precond_deinit

  subroutine precond_write(precond, datafile, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_pert, only: RT0_init, RT0_deinit, RT0_write, RT0_tor_comp_from_zero_div
    use mephit_mesh, only: mesh
    use mephit_conf, only: conf
    type(precond_t), intent(in) :: precond
    character(len = *), intent(in) :: datafile
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"
    type(RT0_t) :: eigvec
    integer(HID_T) :: h5id_root
    integer :: i

    grp = trim(group)
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/nritz', precond%nritz, comment = 'number of Ritz values')
    if (precond%nritz > 0) then
      call h5_add(h5id_root, grp // '/Lr', precond%Lr, &
        lbound(precond%Lr), ubound(precond%Lr), &
        comment = 'matrix used in preconditioner', unit = '1')
      call h5_add(h5id_root, grp // '/eigvecs', precond%eigvecs, &
        lbound(precond%eigvecs), ubound(precond%eigvecs), &
        comment = 'iteration eigenvectors', unit = '1')
      call h5_close(h5id_root)
      call RT0_init(eigvec, mesh%nedge, mesh%ntri)
      do i = 1, min(precond%nritz, conf%max_eig_out)
        write (postfix, postfix_fmt) i
        eigvec%DOF(:) = precond%eigvecs(:, i)
        call RT0_tor_comp_from_zero_div(eigvec)
        call RT0_write(eigvec, datafile, grp // '/eigvec' // postfix, &
          'iteration eigenvector', 'G', 1)
      end do
      call RT0_deinit(eigvec)
    else
      call h5_close(h5id_root)
    end if
  end subroutine precond_write

  subroutine precond_read(precond, datafile, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use mephit_pert, only: RT0_init, RT0_deinit, RT0_read
    use mephit_mesh, only: mesh
    type(precond_t), intent(inout) :: precond
    character(len = *), intent(in) :: datafile
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call precond_deinit(precond)
    call h5_open(datafile, h5id_root)
    call h5_get(h5id_root, grp // '/nritz', precond%nritz)
    if (precond%nritz > 0) then
      allocate(precond%Lr(precond%nritz, precond%nritz))
      call h5_get(h5id_root, grp // '/Lr', precond%Lr)
      allocate(precond%eigvecs(mesh%nedge, precond%nritz))
      call h5_get(h5id_root, grp // '/eigvecs', precond%eigvecs)
    end if
    call h5_close(h5id_root)
  end subroutine precond_read

  pure function precond_apply(precond, vec)
    type(precond_t), intent(in) :: precond
    complex(dp), intent(in) :: vec(:)
    complex(dp) :: precond_apply(size(vec))

    precond_apply = matmul(precond%eigvecs(:, 1:precond%nritz), matmul(precond%Lr, &
      matmul(transpose(conjg(precond%eigvecs(:, 1:precond%nritz))), vec)))
  end function precond_apply

  subroutine precond_compute(precond, perteq, fdm, flr2)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, logger, datafile, cmplx_fmt, runmode_precon
    use mephit_util, only: arnoldi_break
    use mephit_mesh, only: mesh
    use mephit_pert, only: vac, RT0_init, RT0_deinit, RT0_write, RT0_tor_comp_from_zero_div, RT0_L2int
    use mephit_flr2, only: flr2_t
    type(precond_t), intent(inout) :: precond
    type(perteq_t), intent(inout) :: perteq
    type(fdm_t), intent(in) :: fdm
    type(flr2_t), intent(in) :: flr2
    integer :: i, j, info
    integer(HID_T) :: h5id_root
    complex(dp), allocatable :: start_vector(:), eigvals(:), Yr(:, :)
    integer, allocatable :: ipiv(:)
    type(RT0_t) :: Bn_prev

    call precond_deinit(precond)
    if (runmode_precon == conf%runmode) then
      ! construct start vector K * B_vac
      perteq%Bn%DOF(:) = vac%Bn%DOF
      perteq%Bn%comp_phi(:) = vac%Bn%comp_phi
      call compute_presn(perteq, fdm, conf%damp)
      call compute_currn(perteq, fdm, flr2, conf%damp, .false.)
      call compute_magfn(perteq)
      allocate(start_vector, source = perteq%Bn%DOF)
      ! calculate eigenvectors -- system dimension: number of non-redundant edges in core plasma
      call arnoldi_break(mesh%nedge, conf%nkrylov, 3, conf%ritz_threshold, conf%ritz_rel_err, &
        next_iteration_arnoldi, start_vector, info, precond%nritz, eigvals, precond%eigvecs)
      deallocate(start_vector)
      if (info < 0) then
        write (logger%msg, '("Error ", i0, " in routine arnoldi_break")') info
        if (logger%err) call logger%write_msg
        error stop
      end if
      call h5_open_rw(datafile, h5id_root)
      call h5_create_parent_groups(h5id_root, 'iter/')
      call h5_add(h5id_root, 'iter/eigvals', eigvals, lbound(eigvals), ubound(eigvals), &
        comment = 'iteration eigenvalues')
      call h5_close(h5id_root)
      if (logger%info) then
        write (logger%msg, '("Arnoldi method yields ", i0, " Ritz eigenvalues > ", es24.16e3)') &
          precond%nritz, conf%ritz_threshold
        call logger%write_msg
        do i = 1, precond%nritz
          write (logger%msg, '("lambda ", i0, ": ", ' // cmplx_fmt // ')') i, eigvals(i)
          call logger%write_msg
        end do
      end if
      if (info > 0) then
        logger%msg = 'Warning: not all eigenvalues converged, consider changing ' // &
          'your configuration (nkrylov, ritz_threshold, ritz_rel_err)'
        if (logger%warn) call logger%write_msg
      end if
      if (precond%nritz > 0) then
        allocate(precond%Lr(precond%nritz, precond%nritz), Yr(precond%nritz, precond%nritz))
        Yr = (0d0, 0d0)
        do i = 1, precond%nritz
          Yr(i, i) = (1d0, 0d0)
          do j = 1, precond%nritz
            precond%Lr(i, j) = sum(conjg(precond%eigvecs(:, i)) * &
              precond%eigvecs(:, j)) * (eigvals(j) - (1d0, 0d0))
          end do
        end do
        allocate(ipiv(precond%nritz))
        call zgesv(precond%nritz, precond%nritz, precond%Lr, &
          precond%nritz, ipiv, Yr, precond%nritz, info)
        deallocate(ipiv)
        if (info == 0) then
          logger%msg = 'Successfully inverted matrix for preconditioner'
          if (logger%info) call logger%write_msg
        else
          write (logger%msg, '("Matrix inversion for preconditioner failed: ' // &
            'zgesv returns error ", i0)') info
          if (logger%err) call logger%write_msg
          error stop
        end if
        do i = 1, precond%nritz
          precond%Lr(i, :) = eigvals(i) * Yr(i, :)
        end do
        deallocate(Yr)
        ! debug kernel of linear operator
        perteq%Bn%DOF(:) = 0d0
        do i = 1, precond%nritz
          if (abs(eigvals(i)) < 1d0) exit
          perteq%Bn%DOF(:) = perteq%Bn%DOF + precond%eigvecs(:, i)
        end do
        perteq%Bn%DOF(:) = RT0_L2int(perteq%Bn) / RT0_L2int(vac%Bn) * perteq%Bn%DOF
        call RT0_tor_comp_from_zero_div(perteq%Bn)
        call RT0_init(Bn_prev, mesh%nedge, mesh%ntri)
        Bn_prev%DOF(:) = perteq%Bn%DOF
        call compute_presn(perteq, fdm, conf%damp)
        call compute_currn(perteq, fdm, flr2, conf%damp, .false.)
        call compute_magfn(perteq)
        perteq%Bn%DOF(:) = perteq%Bn%DOF - precond_apply(precond, perteq%Bn%DOF - Bn_prev%DOF)
        call RT0_tor_comp_from_zero_div(perteq%Bn)
        call RT0_write(perteq%Bn, datafile, 'iter/debug_kernel', &
          'preconditioner applied to its kernel', 'G', 1)
        call RT0_deinit(Bn_prev)
      end if
    end if

  contains

    ! computes B_(n+1) = K * (B_n + B_vac) for Arnoldi iterations
    subroutine next_iteration_arnoldi(old_val, new_val)
      complex(dp), intent(in) :: old_val(:)
      complex(dp), intent(out) :: new_val(:)

      perteq%Bn%DOF(:) = old_val
      call RT0_tor_comp_from_zero_div(perteq%Bn)
      call compute_presn(perteq, fdm, conf%damp)
      call compute_currn(perteq, fdm, flr2, conf%damp, .false.)
      call compute_magfn(perteq)
      new_val(:) = perteq%Bn%DOF
    end subroutine next_iteration_arnoldi
  end subroutine precond_compute

  subroutine perteq_iterate(perteq, precond, fdm, flr2)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, logger, datafile, runmode_single, currn_model_kilca
    use mephit_mesh, only: mesh
    use mephit_pert, only: vac, L1_write, RT0_init, RT0_deinit, RT0_tor_comp_from_zero_div, RT0_L2int
    use mephit_flr2, only: flr2_t
    type(perteq_t), intent(inout) :: perteq
    type(precond_t), intent(inout) :: precond
    type(fdm_t), intent(in) :: fdm
    type(flr2_t), intent(in) :: flr2
    integer :: kiter, niter, maxiter
    integer(HID_T) :: h5id_root
    real(dp), allocatable :: L2int_Bn_diff(:)
    real(dp) :: L2int_Bnvac, rel_err
    type(RT0_t) :: Bn_prev, Bn_diff
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"

    L2int_Bnvac = RT0_L2int(vac%Bn)
    write (logger%msg, '("L2int_Bnvac = ", es24.16e3)') L2int_Bnvac
    if (logger%info) call logger%write_msg
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, 'iter/')
    call h5_add(h5id_root, 'iter/L2int_Bnvac', L2int_Bnvac, &
      comment = 'L2 integral of magnetic field (vacuum)', unit = 'Mx')
    call h5_close(h5id_root)
    if (runmode_single == conf%runmode) then
      maxiter = 0
    else
      maxiter = conf%niter
    end if
    allocate(L2int_Bn_diff(0:maxiter))
    L2int_Bn_diff = ieee_value(0d0, ieee_quiet_nan)
    call RT0_init(Bn_prev, mesh%nedge, mesh%ntri)
    call RT0_init(Bn_diff, mesh%nedge, mesh%ntri)
    perteq%Bn%DOF(:) = vac%Bn%DOF
    perteq%Bn%comp_phi(:) = vac%Bn%comp_phi
    call perteq_write('("iter/", a, "_vac")', ' (vacuum)', magfmn = perteq%Bn)
    if (precond%nritz > 0) then
      perteq%Bn%DOF(:) = perteq%Bn%DOF - precond_apply(precond, perteq%Bn%DOF)
      call RT0_tor_comp_from_zero_div(perteq%Bn)
    end if
    niter = maxiter
    do kiter = 0, maxiter
      write (logger%msg, '("Iteration ", i2, " of ", i2)') kiter, maxiter
      if (logger%info) call logger%write_msg
      write (postfix, postfix_fmt) kiter
      Bn_prev%DOF(:) = perteq%Bn%DOF
      Bn_prev%comp_phi(:) = perteq%Bn%comp_phi
#ifdef USE_MFEM
      if (kiter <= 1) then
        call MFEM_test(perteq%pn)
        call perteq_write('("iter/", a, "MFEM_' // postfix // '")', &
            ' (after MFEM iteration)', presn = perteq%pn, presmn = perteq%pn)
      end if
#endif
      ! compute B_(n+1) = K * B_n + B_vac ... different from next_iteration_arnoldi
      call compute_presn(perteq, fdm, conf%damp)
      if (kiter <= 1) then
        call perteq_write('("iter/", a, "' // postfix // '")', &
          ' (after iteration)', presn = perteq%pn, presmn = perteq%pn)
      end if
      call compute_currn(perteq, fdm, flr2, conf%damp, .false.)
      call compute_magfn(perteq)
      perteq%Bn%DOF(:) = perteq%Bn%DOF + vac%Bn%DOF
      if (precond%nritz > 0) then
        perteq%Bn%DOF(:) = perteq%Bn%DOF - precond_apply(precond, perteq%Bn%DOF - Bn_prev%DOF)
      end if
      call RT0_tor_comp_from_zero_div(perteq%Bn)
      Bn_diff%DOF(:) = perteq%Bn%DOF - Bn_prev%DOF
      Bn_diff%comp_phi(:) = perteq%Bn%comp_phi - Bn_prev%comp_phi
      L2int_Bn_diff(kiter) = RT0_L2int(Bn_diff)
      write (logger%msg, '("L2int_Bn_diff = ", es24.16e3)') L2int_Bn_diff(kiter)
      if (logger%info) call logger%write_msg
      if (kiter <= 1) then
        call perteq_write('("iter/", a, "' // postfix // '")', ' (after iteration)', &
          parcurrn = perteq%jnpar_B0, currn = perteq%jn, magfn = perteq%Bn)
        call perteq_write('("iter/", a, "_diff' // postfix // '")', &
          ' (difference between iterations)', magfn = Bn_diff)
      else
        if (L2int_Bn_diff(kiter) < conf%iter_rel_err * L2int_Bnvac) then
          niter = kiter
          exit
        end if
      end if
      call perteq_write('("iter/", a, "' // postfix // '")', ' (after iteration)', &
        presmn = perteq%pn, parcurrmn = perteq%jnpar_B0, Ires = perteq%jnpar_B0, &
        currmn = perteq%jn, magfmn = perteq%Bn)
    end do
    rel_err = L2int_Bn_diff(niter) / L2int_Bnvac
    write (logger%msg, '("Relative error after ", i0, " iterations: ", es24.16e3)') &
      niter, rel_err
    if (logger%info) call logger%write_msg
    if (maxiter > 0 .and. rel_err >= conf%iter_rel_err) then
      write (logger%msg, '("Requested relative error ", es24.16e3, ' // &
        '" could not be reached within the requested ", i0, " iterations")') &
        conf%iter_rel_err, conf%niter
      if (logger%warn) call logger%write_msg
    end if
    perteq%Bnplas%DOF(:) = perteq%Bn%DOF - vac%Bn%DOF
    perteq%Bnplas%comp_phi(:) = perteq%Bn%comp_phi - vac%Bn%comp_phi
    if (conf%kilca_scale_factor /= 0) then
      call check_furth(perteq%jn, perteq%Bnplas)
    end if
    ! save results
    call h5_open_rw(datafile, h5id_root)
    call h5_add(h5id_root, 'iter/niter', niter, comment = 'actual number of iterations')
    call h5_add(h5id_root, 'iter/rel_err', rel_err, comment = 'relative error of iterations')
    call h5_add(h5id_root, 'iter/L2int_Bn_diff', L2int_Bn_diff, &
      lbound(L2int_Bn_diff), ubound(L2int_Bn_diff), &
      comment = 'L2 integral of magnetic field (difference between iterations)', unit = 'Mx')
    if (conf%currn_model == currn_model_kilca) then
      call h5_add(h5id_root, 'iter/Phi_mn', perteq%Phi_mn, &
        lbound(perteq%Phi_mn), ubound(perteq%Phi_mn), &
        comment = 'electric potential perturbation', unit = 'statV')
      call h5_add(h5id_root, 'iter/Phi_aligned_mn', perteq%Phi_aligned_mn, &
        lbound(perteq%Phi_aligned_mn), ubound(perteq%Phi_aligned_mn), &
        comment = 'aligned electric potential perturbation', unit = 'statV')
    end if
    call h5_close(h5id_root)
    deallocate(L2int_Bn_diff)
    call L1_write(perteq%AnR, datafile, 'iter/AnR', &
      'R component of vector potential for GORILLA (full perturbation)', 'G cm')
    call L1_write(perteq%AnZ, datafile, 'iter/AnZ', &
      'Z component of vector potential for GORILLA (full perturbation)', 'G cm')
    call perteq_write('("iter/", a)', ' (full perturbation)', &
      perteq%pn, perteq%pn, perteq%jnpar_B0, perteq%jnpar_B0, perteq%jnpar_B0, &
      perteq%jn, perteq%jn, perteq%Bn, perteq%Bn)
    call perteq_write('("iter/", a, "_plas")', ' (plasma response)', &
      magfn = perteq%Bnplas, magfmn = perteq%Bnplas)

    call RT0_deinit(Bn_prev)
    call RT0_deinit(Bn_diff)
  end subroutine perteq_iterate

  subroutine debug_initial_iteration(perteq, fdm, flr2)
    use mephit_conf, only: conf
    use mephit_pert, only: vac
    use mephit_flr2, only: flr2_t
    type(perteq_t), intent(inout) :: perteq
    type(fdm_t), intent(in) :: fdm
    type(flr2_t), intent(in) :: flr2

    if (conf%debug_initial) then
      perteq%Bn%DOF(:) = vac%Bn%DOF
      perteq%Bn%comp_phi(:) = vac%Bn%comp_phi
#ifdef USE_MFEM
      call MFEM_test(perteq%pn)
      call perteq_write('("debug_MFEM_initial/MFEM_", a)', &
        ' (initial MFEM iteration)', presn = perteq%pn, presmn = perteq%pn)
#endif
      call compute_presn(perteq, fdm, .false.)
#ifdef USE_MFEM
      call perteq_write('("debug_MFEM_initial/", a)', &
        ' (initial iteration)', presn = perteq%pn, presmn = perteq%pn)
#endif
      call compute_currn(perteq, fdm, flr2, .false., .true.)
      perteq%Bn%DOF(:) = vac%Bn%DOF
      perteq%Bn%comp_phi(:) = vac%Bn%comp_phi
      call compute_presn(perteq, fdm, .true.)
      call compute_currn(perteq, fdm, flr2, .true., .true.)
    end if
  end subroutine debug_initial_iteration

  ! This subroutine calls a C function that pipes the data to/from FreeFem.
  subroutine compute_magfn(perteq)
    use mephit_mesh, only: mesh
    use mephit_pert, only: RT0_tor_comp_from_zero_div
    type(perteq_t), intent(inout) :: perteq

    call FEM_compute_magfn(mesh%nedge, mesh%npoint, perteq%jn%DOF, &
      perteq%Bn%DOF, perteq%AnR%DOF, perteq%AnZ%DOF)
    call RT0_tor_comp_from_zero_div(perteq%Bn)
  end subroutine compute_magfn

  subroutine unit_B0(R, Z, vector) bind(C, name = 'unit_B0')
    use iso_c_binding, only: c_double
    use field_sub, only : field
    real(c_double), intent(in), value :: R
    real(c_double), intent(in), value :: Z
    real(c_double), intent(out) :: vector(3)
    real(dp) :: dum

    call field(R, 0d0, Z, vector(1), vector(3), vector(2), &
      dum, dum, dum, dum, dum, dum, dum, dum, dum)
    dum = sqrt(sum(vector * vector))
    vector = vector / dum
  end subroutine unit_B0

  subroutine presn_inhom(R, Z, scalar) bind(C, name = 'presn_inhom')
    use iso_c_binding, only: c_double, c_double_complex
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use field_eq_mod, only: psib
    use field_sub, only: psif
    use mephit_util, only: interp1d
    use mephit_mesh, only: equil, mesh, point_location
    use mephit_pert, only: RT0_interp, vac
    use field_sub, only : field
    real(c_double), intent(in), value :: R
    real(c_double), intent(in), value :: Z
    complex(c_double_complex), intent(out) :: scalar
    real(dp) :: B_0(3), dum, psi, dp0_dpsi
    integer :: ktri
    complex(dp) :: B_n(3)

    call field(R, 0d0, Z, B_0(1), B_0(2), B_0(3), &
      dum, dum, dum, dum, dum, dum, dum, dum, dum)
    psi = psif - psib  ! see interp_psi_pol in mephit_util
    dp0_dpsi = interp1d(equil%psi_eqd, equil%pprime, psi, 3)
    ktri = point_location(R, Z)
    if (ktri <= 0 .or. ktri > mesh%ntri) then
      scalar = cmplx(ieee_value(0d0, ieee_quiet_nan), ieee_value(0d0, ieee_quiet_nan), dp)
      return
    end if
    call RT0_interp(vac%Bn, ktri, R, Z, B_n)
    scalar = -dp0_dpsi * (B_n(1) * B_0(3) - B_n(3) * B_0(1)) * R / sqrt(sum(B_0 * B_0))
  end subroutine presn_inhom

#ifdef USE_MFEM
  subroutine MFEM_test(pn)
    use iso_c_binding, only: c_int, c_null_char, c_loc, c_funloc
    use mephit_conf, only: conf, logger, basename_suffix, decorate_filename
    type(L1_t), intent(inout) :: pn
    character(len = 1024) :: mesh_file
    integer(c_int) :: status

    mesh_file = decorate_filename('core_plasma.mesh', '', basename_suffix)
    status = FEM_test(trim(mesh_file) // c_null_char, conf%n, size(pn%DOF), &
      pn%DOF, c_funloc(unit_B0), c_funloc(presn_inhom))
    if (logger%debug) then
      write (logger%msg, '("FEM_test return status: ", i0)') status
      call logger%write_msg
    end if
  end subroutine MFEM_test
#endif

  subroutine FDM_init(fdm, nnz)
    type(FDM_t), intent(inout) :: fdm
    integer, intent(in) :: nnz

    call FDM_deinit(fdm)
    fdm%nnz = nnz
    allocate(fdm%irow(nnz))
    allocate(fdm%icol(nnz))
    allocate(fdm%aval_MDE(nnz))
    allocate(fdm%aval_MDE_damped(nnz))
    allocate(fdm%aval_jnperp(nnz))
  end subroutine FDM_init

  subroutine FDM_deinit(fdm)
    type(FDM_t), intent(inout) :: fdm

    fdm%nnz = 0
    if (allocated(fdm%irow)) deallocate(fdm%irow)
    if (allocated(fdm%icol)) deallocate(fdm%icol)
    if (allocated(fdm%aval_MDE)) deallocate(fdm%aval_MDE)
    if (allocated(fdm%aval_MDE_damped)) deallocate(fdm%aval_MDE_damped)
    if (allocated(fdm%aval_jnperp)) deallocate(fdm%aval_jnperp)
  end subroutine FDM_deinit

  subroutine FDM_write(fdm, datafile, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(FDM_t), intent(in) :: fdm
    character(len = *), intent(in) :: datafile
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/nnz', fdm%nnz, &
      comment = 'number of nonzero entries in sparse matrix')
    call h5_add(h5id_root, grp // '/irow', fdm%irow, lbound(fdm%irow), ubound(fdm%irow), &
      comment = 'row indices of sparse matrix entries (COO format)')
    call h5_add(h5id_root, grp // '/icol', fdm%icol, lbound(fdm%icol), ubound(fdm%icol), &
      comment = 'column indices of sparse matrix entries (COO format)')
    call h5_add(h5id_root, grp // '/aval_MDE', fdm%aval_MDE, &
      lbound(fdm%aval_MDE), ubound(fdm%aval_MDE), &
      comment = 'sparse matrix entries for MDE (COO format)', unit = 'G cm^-1')
    call h5_add(h5id_root, grp // '/aval_MDE_damped', fdm%aval_MDE_damped, &
      lbound(fdm%aval_MDE_damped), ubound(fdm%aval_MDE_damped), &
      comment = 'sparse matrix entries for damped MDE (COO format)', unit = 'G cm^-1')
    call h5_add(h5id_root, grp // '/aval_jnperp', fdm%aval_jnperp, &
      lbound(fdm%aval_jnperp), ubound(fdm%aval_jnperp), &
      comment = 'sparse matrix entries for helical current ODE (COO format)', unit = '1')
    call h5_close(h5id_root)
  end subroutine FDM_write

  subroutine FDM_read(fdm, datafile, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(FDM_t), intent(inout) :: fdm
    character(len = *), intent(in) :: datafile
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer :: nnz

    grp = trim(group)
    call h5_open(datafile, h5id_root)
    call h5_get(h5id_root, grp // '/nnz', nnz)
    call FDM_init(fdm, nnz)
    call h5_get(h5id_root, grp // '/irow', fdm%irow)
    call h5_get(h5id_root, grp // '/icol', fdm%icol)
    call h5_get(h5id_root, grp // '/aval_MDE', fdm%aval_MDE)
    call h5_get(h5id_root, grp // '/aval_MDE_damped', fdm%aval_MDE_damped)
    call h5_get(h5id_root, grp // '/aval_jnperp', fdm%aval_jnperp)
    call h5_close(h5id_root)
  end subroutine FDM_read

  subroutine FDM_compute_matrix(fdm)
    use mephit_util, only: imun
    use mephit_mesh, only: mesh, cache, fs
    type(FDM_t), intent(inout) :: fdm
    complex(dp), dimension(:), allocatable :: a, b, b_damped, a_jnperp, b_jnperp
    integer :: kf, kp, kedge, ndim

    ! discretised ODE: $a_{k} (y_{k+1} - y_{k}) + b_{k} (y_{k+1} + y_{k}), k \mod N$
    ! resulting (upper) diagonal elements:
    ! d = -a + b
    ! du = a + b
    ndim = mesh%npoint - 1
    allocate(a(ndim), b(ndim), b_damped(ndim))
    a(:) = (0d0, 0d0)
    b(:) = (0d0, 0d0)
    b_damped(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        ! use midpoint of poloidal edge
        associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge))
          a(kedge) = (f%B0(1) * mesh%edge_R(kedge) + f%B0(3) * mesh%edge_Z(kedge)) / &
            (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
          b(kedge) = 0.5d0 * imun * mesh%n * f%B0(2) / R
          b_damped(kedge) = 0.5d0 * imun * (mesh%n + imun * mesh%damping(kf)) * f%B0(2) / R
        end associate
      end do
    end do
    allocate(a_jnperp(ndim), b_jnperp(ndim))
    a_jnperp(:) = (0d0, 0d0)
    b_jnperp(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        associate (s => cache%sample_jnperp(kedge))
          a_jnperp(kedge) = (s%dR_dtheta * mesh%edge_R(kedge) + s%dZ_dtheta * mesh%edge_Z(kedge)) / &
            (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
          b_jnperp(kedge) = -0.5d0 * imun * mesh%n / fs%F(kf) * &
            (s%dR_dtheta * s%B0_R + s%dZ_dtheta * s%B0_Z)
        end associate
      end do
    end do
    ! assemble sparse matrix (COO format) from blocks
    ! icol should be non-decreasing, and irow should be increasing for fixed icol
    call FDM_init(fdm, 2 * ndim)
    do kf = 1, mesh%nflux
      kedge = mesh%kp_low(kf)
      ! first column, diagonal
      fdm%irow(2 * kedge - 1) = kedge
      fdm%icol(2 * kedge - 1) = kedge
      fdm%aval_MDE(2 * kedge - 1) = -a(kedge) + b(kedge)
      fdm%aval_MDE_damped(2 * kedge - 1) = -a(kedge) + b_damped(kedge)
      fdm%aval_jnperp(2 * kedge - 1) = -a_jnperp(kedge) + b_jnperp(kedge)
      ! first column, off-diagonal (lower left corner)
      fdm%irow(2 * kedge) = kedge + mesh%kp_max(kf) - 1
      fdm%icol(2 * kedge) = kedge
      fdm%aval_MDE(2 * kedge) = &
        a(kedge + mesh%kp_max(kf) - 1) + b(kedge + mesh%kp_max(kf) - 1)
      fdm%aval_MDE_damped(2 * kedge) = &
        a(kedge + mesh%kp_max(kf) - 1) + b_damped(kedge + mesh%kp_max(kf) - 1)
      fdm%aval_jnperp(2 * kedge) = &
        a_jnperp(kedge + mesh%kp_max(kf) - 1) + b_jnperp(kedge + mesh%kp_max(kf) - 1)
      do kp = 2, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        ! off-diagonal
        fdm%irow(2 * kedge - 1) = kedge - 1
        fdm%icol(2 * kedge - 1) = kedge
        fdm%aval_MDE(2 * kedge - 1) = a(kedge - 1) + b(kedge - 1)
        fdm%aval_MDE_damped(2 * kedge - 1) = a(kedge - 1) + b_damped(kedge - 1)
        fdm%aval_jnperp(2 * kedge - 1) = a_jnperp(kedge - 1) + b_jnperp(kedge - 1)
        ! diagonal
        fdm%irow(2 * kedge) = kedge
        fdm%icol(2 * kedge) = kedge
        fdm%aval_MDE(2 * kedge) = -a(kedge) + b(kedge)
        fdm%aval_MDE_damped(2 * kedge) = -a(kedge) + b_damped(kedge)
        fdm%aval_jnperp(2 * kedge) = -a_jnperp(kedge) + b_jnperp(kedge)
      end do
    end do
    deallocate(a, b, b_damped, a_jnperp, b_jnperp)
  end subroutine FDM_compute_matrix

  subroutine FDM_solve(fdm, aval, inhom, solution)
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: logger
    use mephit_mesh, only: mesh
    type(FDM_t), intent(in) :: fdm
    complex(dp), dimension(:), intent(in) :: aval
    complex(dp), dimension(:), intent(in) :: inhom
    complex(dp), dimension(:), intent(out) :: solution
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(:), allocatable :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: ndim
    real(dp), parameter :: small = tiny(0d0)

    if (size(aval) /= fdm%nnz) then
      call logger%msg_arg_size('FDM_solve', &
        'size(aval)', 'fdm%nnz', size(aval), fdm%nnz)
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (size(inhom) /= mesh%npoint) then
      call logger%msg_arg_size('FDM_solve', &
        'size(inhom)', 'mesh%npoint', size(inhom), mesh%npoint)
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (size(solution) /= mesh%npoint) then
      call logger%msg_arg_size('FDM_solve', &
        'size(solution)', 'mesh%npoint', size(solution), mesh%npoint)
      if (logger%err) call logger%write_msg
      error stop
    end if
    ndim = mesh%npoint - 1
    solution(:) = inhom
    call sparse_solve(ndim, ndim, fdm%nnz, fdm%irow, fdm%icol, aval, solution(2:))
    ! first point on axis - average over enclosing flux surface
    solution(1) = sum(solution(2:mesh%kp_max(1) + 1)) / dble(mesh%kp_max(1))
    ! check relative error
    call sparse_matmul(ndim, ndim, fdm%irow, fdm%icol, aval, solution(2:), resid)
    resid(:) = resid - inhom(2:)
    allocate(rel_err(ndim))
    where (abs(inhom(2:)) >= small)
      rel_err(:) = abs(resid) / abs(inhom(2:))
    elsewhere
      rel_err(:) = 0d0
    end where
    max_rel_err = maxval(rel_err)
    avg_rel_err = sum(rel_err) / ndim
    write (logger%msg, '("FDM_solve: diagonalization max_rel_err = ", ' // &
      'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (logger%debug) call logger%write_msg
    deallocate(resid, rel_err)
  end subroutine FDM_solve

  subroutine compute_presn(perteq, fdm, apply_damping)
    use mephit_mesh, only: fs, mesh, cache
    type(perteq_t), intent(inout) :: perteq
    logical, intent(in) :: apply_damping
    type(FDM_t), intent(in) :: fdm
    complex(dp), dimension(mesh%npoint) :: inhom
    integer :: kf, kp, kedge

    inhom = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        ! use midpoint of poloidal edge
        associate (f => cache%mid_fields(kedge))
          inhom(kedge + 1) = -fs%dp_dpsi(kf) * perteq%Bn%DOF(kedge) * &
            (f%B0(1) * mesh%edge_R(kedge) + f%B0(3) * mesh%edge_Z(kedge)) / &
            (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
        end associate
      end do
    end do
    if (apply_damping) then
      call FDM_solve(fdm, fdm%aval_MDE_damped, inhom, perteq%pn%DOF)
    else
      call FDM_solve(fdm, fdm%aval_MDE, inhom, perteq%pn%DOF)
    end if
    perteq%pn%DOF(1) = (0d0, 0d0)
  end subroutine compute_presn

  subroutine debug_MDE(group, presn, magfn, currn_perp, currn_par)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: datafile
    use mephit_util, only: pi, clight, zd_cross
    use mephit_mesh, only: mesh, cache, equilibrium_field, curr0_geqdsk
    use mephit_pert, only: L1_t, L1_interp, RT0_t, RT0_interp
    character(len = *), intent(in) :: group
    type(L1_t), intent(in) :: presn
    type(RT0_t), intent(in) :: magfn
    type(RT0_t), intent(in) :: currn_perp
    type(L1_t), intent(in) :: currn_par
    integer(HID_T) :: h5id_root
    character(len = len_trim(group)) :: grp
    integer :: ndim, kf, kpol, k
    real(dp) :: psi, B0(3), dB0_dR(3), dB0_dZ(3), Bmod, dBmod_dR, dBmod_dZ, &
      j0(3), dj0_dR(3), dj0_dZ(3), grad_j0B0(3), B0_grad_B0(3)
    complex(dp) :: Bn(3), dBn_dR(3), dBn_dphi(3), dBn_dZ(3), grad_BnB0(3), jnperp(3, 0:3)
    complex(dp), allocatable :: lorentz(:, :), pn(:), grad_pn(:, :), Bn_psi_contravar(:), &
      jnpar_Bmod(:), grad_jnpar_Bmod(:, :), div_jnperp(:), div_jnperp_RT0(:)

    grp = trim(adjustl(group))
    ndim = mesh%npoint
    allocate(lorentz(3, ndim), pn(ndim), grad_pn(3, ndim), Bn_psi_contravar(ndim), &
      jnpar_Bmod(ndim), grad_jnpar_Bmod(3, ndim), div_jnperp(ndim), div_jnperp_RT0(ndim))
    lorentz(:, :) = (0d0, 0d0)
    pn(:) = (0d0, 0d0)
    grad_pn(:, :) = (0d0, 0d0)
    Bn_psi_contravar(:) = (0d0, 0d0)
    jnpar_Bmod(:) = (0d0, 0d0)
    grad_jnpar_Bmod(:, :) = (0d0, 0d0)
    div_jnperp(:) = (0d0, 0d0)
    div_jnperp_RT0(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kpol = 1, mesh%kp_max(kf)
        k = mesh%kp_low(kf) + kpol
        associate (s => cache%sample_polmodes(k))
          call equilibrium_field(s%R, s%Z, B0, dB0_dR, dB0_dZ, psi, Bmod, dBmod_dR, dBmod_dZ)
          call curr0_geqdsk(s%R, psi, B0, dB0_dR, dB0_dZ, j0, dj0_dR, dj0_dZ)
          call L1_interp(presn, s%ktri, s%R, s%Z, pn(k), grad_pn(:, k))
          call RT0_interp(magfn, s%ktri, s%R, s%Z, Bn, dBn_dR, dBn_dphi, dBn_dZ)
          call L1_interp(currn_par, s%ktri, s%R, s%Z, jnpar_Bmod(k), grad_jnpar_Bmod(:, k))
          call RT0_interp(currn_perp, s%ktri, s%R, s%Z, jnperp(:, 0), &
            jnperp(:, 1), jnperp(:, 2), jnperp(:, 3))
          Bn_psi_contravar(k) = s%R * (Bn(1) * B0(3) - Bn(3) * B0(1))
          lorentz(:, k) = (zd_cross(jnperp(:, 0), B0) - zd_cross(Bn, j0)) / clight
          div_jnperp_RT0(k) = jnperp(1, 1) + jnperp(2, 2) + jnperp(3, 3) + jnperp(1, 0) / s%R
          grad_j0B0 = [sum(dj0_dR * B0 + dB0_dR * j0), 0d0, sum(dj0_dZ * B0 + dB0_dZ * j0)]
          grad_BnB0 = [sum(dBn_dR * B0 + dB0_dR * Bn), sum(dBn_dphi * B0), sum(dBn_dZ * B0 + dB0_dZ * Bn)]
          B0_grad_B0 = [sum(dB0_dR * B0), 0d0, sum(dB0_dZ * B0)]
          div_jnperp(k) = (-2d0 / Bmod ** 2 * (clight * sum(zd_cross(grad_pn(:, k), B0) * B0_grad_B0) + &
            sum(Bn * B0) * sum(j0 * B0_grad_B0) - sum(Bn * B0_grad_B0) * sum(j0 * B0)) + &
            sum(grad_BnB0 * j0 - Bn * grad_j0B0) + 4d0 * pi * sum(grad_pn(:, k) * j0)) / Bmod ** 2
        end associate
      end do
    end do

    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/pn', pn, lbound(pn), ubound(pn), &
      comment = 'perturbation pressure', unit = 'dyn cm^-2')
    call h5_add(h5id_root, grp // '/grad_pn', grad_pn, lbound(grad_pn), ubound(grad_pn), &
      comment = 'perturbation pressure gradient', unit = 'dyn cm^-3')
    call h5_add(h5id_root, grp // '/Bn_psi_contravar', Bn_psi_contravar, &
      lbound(Bn_psi_contravar), ubound(Bn_psi_contravar), unit = 'G^2 cm', &
      comment = 'perturbation magnetic field (contravariant psi component)')
    call h5_add(h5id_root, grp // '/jnpar_Bmod', jnpar_Bmod, &
      lbound(jnpar_Bmod), ubound(jnpar_Bmod), unit = 's^-1', &
      comment = 'parallel perturbation current density')
    call h5_add(h5id_root, grp // '/grad_jnpar_Bmod', grad_jnpar_Bmod, &
      lbound(grad_jnpar_Bmod), ubound(grad_jnpar_Bmod), unit = 'cm^-1 s^-1', &
      comment = 'parallel perturbation current density gradient')
    call h5_add(h5id_root, grp // '/div_jnperp', div_jnperp, &
      lbound(div_jnperp), ubound(div_jnperp), unit = 'statA cm^-3', &
      comment = 'perpendicular perturbation current density divergence')
    call h5_add(h5id_root, grp // '/div_jnperp_RT0', div_jnperp_RT0, &
      lbound(div_jnperp_RT0), ubound(div_jnperp_RT0), unit = 'statA cm^-3', &
      comment = 'perpendicular perturbation current density divergence (RT0 interpolation)')
    call h5_add(h5id_root, grp // '/lorentz', lorentz, lbound(lorentz), ubound(lorentz), &
      comment = 'perturbation Lorentz force density', unit = 'dyn cm^-3')
    call h5_close(h5id_root)

    deallocate(lorentz, pn, grad_pn, Bn_psi_contravar, jnpar_Bmod, grad_jnpar_Bmod, div_jnperp, div_jnperp_RT0)
  end subroutine debug_MDE

  subroutine compute_currn(perteq, fdm, flr2, apply_damping, debug_initial)
    use mephit_conf, only: conf, currn_model_kilca, currn_model_mhd, logger, datafile
    use mephit_mesh, only: mesh
    use mephit_pert, only: polmodes_t, polmodes_init, polmodes_write, polmodes_deinit, &
      L1_t, L1_init, L1_deinit, RT0_t, RT0_init, RT0_deinit
    use mephit_flr2, only: flr2_t
    type(perteq_t), intent(inout) :: perteq
    type(fdm_t), intent(in) :: fdm
    type(flr2_t), intent(in) :: flr2
    logical, intent(in) :: apply_damping
    logical, intent(in) :: debug_initial
    type(polmodes_t) :: resonant_jmnpar_over_Bmod
    type(L1_t) :: resonant_jnpar_over_Bmod
    type(RT0_t) :: resonant_jn

    call polmodes_init(resonant_jmnpar_over_Bmod, conf%m_max, mesh%nflux)
    call L1_init(resonant_jnpar_over_Bmod, mesh%npoint)
    call RT0_init(resonant_jn, mesh%nedge, mesh%ntri)

    call current_from_pressure_balance(perteq%pn, perteq%Bn, fdm, apply_damping, perteq%jnpar_B0, perteq%jn)
    if (debug_initial) then
      call debug_MDE("debug_MDE_initial", perteq%pn, perteq%Bn, perteq%jn, perteq%jnpar_B0)
      if (apply_damping) then
        call perteq_write('("debug_KiLCA/", a, "_incl")', &
          ' from iMHD including damping', parcurrmn = perteq%jnpar_B0)
      else
        call perteq_write('("debug_KiLCA/", a, "_excl")', &
          ' from iMHD excluding damping', parcurrmn = perteq%jnpar_B0)
      end if
    end if
    select case (conf%currn_model)
    case (currn_model_mhd)
      call compute_shielding_current(perteq%pn, resonant_jmnpar_over_Bmod)
    case (currn_model_kilca)
      call compute_flr2_current(perteq%Bn, resonant_jmnpar_over_Bmod, &
        flr2, perteq%Phi_mn, perteq%Phi_aligned_mn)
    case default
      write (logger%msg, '("unknown response current model selection", i0)') conf%currn_model
      if (logger%err) call logger%write_msg
      error stop
    end select
    call helical_current_from_parallel_current(resonant_jmnpar_over_Bmod, fdm, &
      resonant_jnpar_over_Bmod, resonant_jn)
    perteq%jnpar_B0%DOF(:) = perteq%jnpar_B0%DOF + resonant_jnpar_over_Bmod%DOF
    perteq%jn%DOF(:) = perteq%jn%DOF + resonant_jn%DOF
    if (debug_initial) then
      call polmodes_write(resonant_jmnpar_over_Bmod, datafile, 'debug_KiLCA/jmnpar_Bmod_KiLCA', &
        'parallel current density from KiLCA', 's^-1')  ! SI: H^-1
      call perteq_write('("debug_KiLCA/", a, "_total")', &
        ' including KiLCA current', parcurrmn = perteq%jnpar_B0)
    end if
    ! hack: overwrite to save only last iteration step
    call debug_MDE("debug_MDE_final", perteq%pn, perteq%Bn, perteq%jn, perteq%jnpar_B0)

    call polmodes_deinit(resonant_jmnpar_over_Bmod)
    call L1_deinit(resonant_jnpar_over_Bmod)
    call RT0_deinit(resonant_jn)
  end subroutine compute_currn

  subroutine current_from_pressure_balance(pn, Bn, fdm, apply_damping, jnpar_B0, jn)
    use mephit_util, only: pi, clight, zd_cross
    use mephit_mesh, only: mesh, cache, fs
    use mephit_pert, only: L1_t, L1_interp, RT0_t, RT0_interp
    type(L1_t), intent(in) :: pn
    type(RT0_t), intent(in) :: Bn
    type(fdm_t), intent(in) :: fdm
    logical, intent(in) :: apply_damping
    type(L1_t), intent(inout) :: jnpar_B0
    type(RT0_t), intent(inout) :: jn
    integer :: kf, kp, ktri, kedge, k
    real(dp), dimension(3) :: grad_j0B0, B0_grad_B0, edge_perp
    complex(dp) :: zdum, B0_jnpar, avg_Bn_tor
    complex(dp), dimension(3) :: grad_pn, B_n, dBn_dR, dBn_dZ, dBn_dphi, grad_BnB0
    complex(dp), dimension(mesh%npoint) :: inhom

    jn%DOF(:) = (0d0, 0d0)
    jn%comp_phi(:) = (0d0, 0d0)
    jnpar_B0%DOF(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        ktri = mesh%edge_tri(1, kedge)
        ! use midpoint of poloidal edge
        associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge), Z => mesh%mid_Z(kedge))
          call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
          call RT0_interp(Bn, ktri, R, Z, B_n, dBn_dR, dBn_dphi, dBn_dZ)
          grad_j0B0 = [sum(f%dj0_dR * f%B0 + f%dB0_dR * f%j0), 0d0, sum(f%dj0_dZ * f%B0 + f%dB0_dZ * f%j0)]
          grad_BnB0 = [sum(dBn_dR * f%B0 + f%dB0_dR * B_n), sum(dBn_dphi * f%B0), sum(dBn_dZ * f%B0 + f%dB0_dZ * B_n)]
          B0_grad_B0 = [sum(f%dB0_dR * f%B0), 0d0, sum(f%dB0_dZ * f%B0)]
          inhom(kedge + 1) = (-2d0 / f%Bmod ** 2 * (clight * sum(zd_cross(grad_pn, f%B0) * B0_grad_B0) + &
            sum(B_n * f%B0) * sum(f%j0 * B0_grad_B0) - sum(B_n * B0_grad_B0) * sum(f%j0 * f%B0)) + &
            sum(grad_BnB0 * f%j0 - B_n * grad_j0B0) + 4d0 * pi * sum(grad_pn * f%j0)) / f%Bmod ** 2
        end associate
      end do
    end do
    if (apply_damping) then
      call FDM_solve(fdm, fdm%aval_MDE_damped, inhom, jnpar_B0%DOF)
    else
      call FDM_solve(fdm, fdm%aval_MDE, inhom, jnpar_B0%DOF)
    end if
    avg_Bn_tor = sum(Bn%comp_phi(1:mesh%kp_max(1)) * mesh%area(1:mesh%kp_max(1))) / sum(mesh%area(1:mesh%kp_max(1)))
    jnpar_B0%DOF(1) = clight * (fs%dp_dpsi(0) * mesh%R_O + fs%FdF_dpsi(0) / (4d0 * pi * mesh%R_O)) * &
      (mesh%R_O / fs%F(0)) ** 2 * (avg_Bn_tor + 4d0 * pi * pn%DOF(1) * mesh%R_O / fs%F(0))
    ! project to poloidal current density
    do kedge = 1, mesh%nedge
      ktri = mesh%edge_tri(1, kedge)
      edge_perp = [mesh%edge_Z(kedge), 0d0, -mesh%edge_R(kedge)]
      do k = 1, mesh%GL_order
        associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
          call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
          call RT0_interp(Bn, ktri, R, Z, B_n)
          call L1_interp(jnpar_B0, ktri, R, Z, B0_jnpar)
          B0_jnpar = B0_jnpar * f%Bmod ** 2
          jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * R * &
            (B0_jnpar * sum(f%B0 * edge_perp) - clight * sum(zd_cross(grad_pn, f%B0) * edge_perp) + &
            sum(f%j0 * f%B0) * sum(B_n * edge_perp) - sum(B_n * f%B0) * sum(f%j0 * edge_perp)) / f%Bmod ** 2
        end associate
      end do
    end do
    ! project to current density
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
          call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
          call RT0_interp(Bn, ktri, R, Z, B_n)
          call L1_interp(jnpar_B0, ktri, R, Z, B0_jnpar)
          B0_jnpar = B0_jnpar * f%Bmod ** 2
          jn%comp_phi(ktri) = jn%comp_phi(ktri) + mesh%GL2_weights(k) * &
            (B0_jnpar * f%B0(2) - clight * (grad_pn(3) * f%B0(1) - grad_pn(3) * f%B0(3)) + &
            sum(f%j0 * f%B0) * B_n(2) - sum(B_n * f%B0) * f%j0(2)) / f%Bmod ** 2
        end associate
      end do
    end do
  end subroutine current_from_pressure_balance

  subroutine helical_current_from_parallel_current(jmnpar_over_Bmod, fdm, jnpar_over_Bmod, jn)
    use mephit_util, only: imun
    use mephit_mesh, only: equil, mesh, cache, fs
    use mephit_pert, only: RT0_t, polmodes_t, L1_t, L1_init, L1_interp, L1_deinit
    type(polmodes_t), intent(in) :: jmnpar_over_Bmod
    type(fdm_t), intent(in) :: fdm
    type(L1_t), intent(inout) :: jnpar_over_Bmod
    type(RT0_t), intent(inout) :: jn
    integer :: m, m_res, kf, kedge, ktri, kp, k, kpoi
    real(dp) :: edge_perp(2)
    complex(dp) :: inhom_jnperp(mesh%npoint), coeff_jnperp_interp, jnpar_over_Bmod_interp
    type(L1_t) :: coeff_jnperp

    inhom_jnperp(:) = (0d0, 0d0)
    call L1_init(coeff_jnperp, mesh%npoint)
    jnpar_over_Bmod%DOF(:) = (0d0, 0d0)
    do m = mesh%m_res_min, mesh%m_res_max
      m_res = -equil%cocos%sgn_q * m
      do kf = 1, mesh%nflux
        do kp = 1, mesh%kp_max(kf)
          ! expand Fourier modes
          kpoi = mesh%kp_low(kf) + kp
          jnpar_over_Bmod%DOF(kpoi) = jnpar_over_Bmod%DOF(kpoi) + &
            jmnpar_over_Bmod%coeff(m_res, kf) * exp(imun * m_res * mesh%node_theta_flux(kpoi))
          ! compute inhomogeneity of ODE for jnperp
          kedge = kpoi - 1
          inhom_jnperp(kpoi) = inhom_jnperp(kpoi) + imun * (mesh%n * fs%q(kf) + m_res) / fs%F(kf) * &
            jmnpar_over_Bmod%coeff(m_res, kf) * exp(imun * m_res * cache%sample_jnperp(kedge)%theta)
        end do
      end do
    end do
    ! compute perpendicular current density coefficient
    call FDM_solve(fdm, fdm%aval_jnperp, inhom_jnperp, coeff_jnperp%DOF)
    ! add poloidal current density
    do kedge = 1, mesh%nedge
      ktri = mesh%edge_tri(1, kedge)
      edge_perp = [mesh%edge_Z(kedge), -mesh%edge_R(kedge)]
      do k = 1, mesh%GL_order
        associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
          call L1_interp(coeff_jnperp, ktri, R, Z, coeff_jnperp_interp)
          call L1_interp(jnpar_over_Bmod, ktri, R, Z, jnpar_over_Bmod_interp)
          jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * R * &
            (coeff_jnperp_interp * R * f%B0(2) + jnpar_over_Bmod_interp) * &
            sum([f%B0(1), f%B0(3)] * edge_perp)
        end associate
      end do
    end do
    ! add toroidal current density
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
          call L1_interp(coeff_jnperp, ktri, R, Z, coeff_jnperp_interp)
          call L1_interp(jnpar_over_Bmod, ktri, R, Z, jnpar_over_Bmod_interp)
          jn%comp_phi(ktri) = jn%comp_phi(ktri) + mesh%GL2_weights(k) * &
            (coeff_jnperp_interp * R * (f%B0(1) ** 2 + f%B0(3) ** 2) + jnpar_over_Bmod_interp * f%B0(2))
        end associate
      end do
    end do
    call L1_deinit(coeff_jnperp)
  end subroutine helical_current_from_parallel_current

  subroutine compute_flr2_coeff(flr2)
    use mephit_conf, only: conf
    use mephit_util, only: ev2erg
    use mephit_mesh, only: equil, mesh, fs, dens_e, temp_e, temp_i, Phi0, dPhi0_dpsi, nu_i, nu_e
    use mephit_flr2, only: flr2_t, flr2_init, flr2_coeff
    type(flr2_t), intent(inout) :: flr2
    integer :: m_min, m_max

    m_min = min(-equil%cocos%sgn_q * mesh%m_res_min, -equil%cocos%sgn_q * mesh%m_res_max)
    m_max = max(-equil%cocos%sgn_q * mesh%m_res_min, -equil%cocos%sgn_q * mesh%m_res_max)
    call flr2_init(flr2, mesh%nflux, m_min, m_max)
    call flr2_coeff(flr2, m_min, m_max, mesh%n, mesh%nflux, conf%m_i, conf%Z_i, &
      mesh%R_O, fs%psi(1:), fs%q(1:), fs%F(1:), Phi0%y(1:), mesh%avg_R2gradpsi2(1:), &
      dens_e%y(1:), temp_e%y(1:) * ev2erg, temp_i%y(1:) * ev2erg, nu_e%y(1:), nu_i%y(1:))
  end subroutine compute_flr2_coeff

  subroutine compute_flr2_current(Bn, jmnpar_over_Bmod, flr2, Phi_mn, Phi_aligned_mn)
    use mephit_conf, only: conf
    use mephit_util, only: imun, resample1d
    use mephit_mesh, only: equil, mesh, fs, fs_half, dPhi0_dpsi
    use mephit_pert, only: polmodes_t, &
      vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, RT0_poloidal_modes
    use mephit_flr2, only: flr2_t, flr2_response_current
    type(RT0_t), intent(in) :: Bn
    type(polmodes_t), intent(inout) :: jmnpar_over_Bmod
    type(flr2_t), intent(in) :: flr2
    complex(dp), intent(out) :: Phi_mn(0:, mesh%m_res_min:)
    complex(dp), intent(out) :: Phi_aligned_mn(0:, mesh%m_res_min:)
    integer :: m, m_res, kf_min
    real(dp) :: Bmnrho_over_B0theta(1:mesh%nflux), rtemp(1:mesh%nflux)
    complex(dp) :: Bmnpsi_over_B0phi(1:mesh%nflux)
    type(vec_polmodes_t) :: Bmn

    kf_min = 0 ! max(1, mesh%res_ind(mesh%m_res_min) / 4)
    Bmnpsi_over_B0phi(:) = (0d0, 0d0)
    jmnpar_over_Bmod%coeff(:, :) = (0d0, 0d0)
    Phi_mn(:, :) = (0d0, 0d0)
    Phi_aligned_mn(:, :) = (0d0, 0d0)
    call vec_polmodes_init(Bmn, conf%m_max, mesh%nflux)
    call RT0_poloidal_modes(Bn, Bmn)
    do m = mesh%m_res_min, mesh%m_res_max
      m_res = -equil%cocos%sgn_q * m
      Bmnrho_over_B0theta(:) = Bmn%coeff_rad(m_res, :)%Re
      call resample1d(fs_half%psi, Bmnrho_over_B0theta / fs_half%q * &
        equil%cocos%sgn_dpsi, fs%psi(1:), rtemp, 3)
      Bmnpsi_over_B0phi%Re = rtemp
      Bmnrho_over_B0theta(:) = Bmn%coeff_rad(m_res, :)%Im
      call resample1d(fs_half%psi, Bmnrho_over_B0theta / fs_half%q * &
        equil%cocos%sgn_dpsi, fs%psi(1:), rtemp, 3)
      Bmnpsi_over_B0phi%Im = rtemp
      call flr2_response_current(flr2, 1, m_res, mesh%nflux, &
        Bmnpsi_over_B0phi, jmnpar_over_Bmod%coeff(m_res, 1:), Phi_mn(1:, m))
      Phi_aligned_mn(1:, m) = imun * Bmnpsi_over_B0phi * &
        fs%q(1:) * dPhi0_dpsi%y(1:) / (m_res + mesh%n * fs%q(1:))
      jmnpar_over_Bmod%coeff(m_res, :kf_min) = (0d0, 0d0)  ! suppress spurious current near axis
    end do
    call vec_polmodes_deinit(Bmn)
  end subroutine compute_flr2_current

  subroutine compute_shielding_current(pn, jmnpar_over_Bmod)
    use mephit_conf, only: conf, conf_arr
    use mephit_mesh, only: equil, mesh, cache
    use mephit_pert, only: L1_t, polmodes_t, polmodes_init, polmodes_deinit, L1_poloidal_modes
    type(L1_t), intent(in) :: pn
    type(polmodes_t), intent(inout) :: jmnpar_over_Bmod
    type(polmodes_t) :: pmn
    integer :: m, m_res, kf

    jmnpar_over_Bmod%coeff(:, :) = (0d0, 0d0)
    call polmodes_init(pmn, conf%m_max, mesh%nflux)
    call L1_poloidal_modes(pn, pmn)
    do m = mesh%m_res_min, mesh%m_res_max
      m_res = -equil%cocos%sgn_q * m
      if (abs(conf_arr%refinement(m)) > 1d0) then
        kf = mesh%res_ind(m)
        jmnpar_over_Bmod%coeff(m_res, kf-1:kf) = &
          (pmn%coeff(m_res, kf) - pmn%coeff(m_res, kf-1)) * &
          cache%shielding(m)%coeff * conf_arr%sheet_current_factor(m)
      end if
    end do
    call polmodes_deinit(pmn)
  end subroutine compute_shielding_current

  subroutine perteq_write(name_fmt, comment, &
    presn, presmn, parcurrn, parcurrmn, Ires, currn, currmn, magfn, magfmn)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_mesh, only: mesh, cache
    use mephit_pert, only: vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, vec_polmodes_write, &
      polmodes_t, polmodes_init, polmodes_deinit, polmodes_write, &
      L1_poloidal_modes, RT0_poloidal_modes, L1_write, RT0_write, compute_Ires
    character(len = *), intent(in) :: name_fmt
    character(len = *), intent(in) :: comment
    type(L1_t), intent(in), optional :: presn
    type(L1_t), intent(in), optional :: presmn
    type(L1_t), intent(in), optional :: parcurrn
    type(L1_t), intent(in), optional :: parcurrmn
    type(L1_t), intent(in), optional :: Ires
    type(RT0_t), intent(in), optional :: currn
    type(RT0_t), intent(in), optional :: currmn
    type(RT0_t), intent(in), optional :: magfn
    type(RT0_t), intent(in), optional :: magfmn
    character(len = 1024) :: name
    integer(HID_T) :: h5id_root
    type(polmodes_t) :: polmodes
    type(vec_polmodes_t) :: vec_polmodes
    integer :: m
    complex(dp), allocatable :: Imn_res(:)

    if (present(presmn) .or. present(parcurrmn)) then
      call polmodes_init(polmodes, conf%m_max, mesh%nflux)
    end if
    if (present(currmn) .or. present(magfmn)) then
      call vec_polmodes_init(vec_polmodes, conf%m_max, mesh%nflux)
    end if

    if (present(presn)) then
      write (name, name_fmt) 'pn'
      call L1_write(presn, datafile, trim(name), &
        'pressure' // comment, 'dyn cm^-2')
    end if
    if (present(presmn)) then
      write (name, name_fmt) 'pmn'
      call L1_poloidal_modes(presmn, polmodes)
      call polmodes_write(polmodes, datafile, trim(name), &
        'poloidal modes of pressure' // comment, 'dyn cm^-2')
    end if
    if (present(parcurrn)) then
      write (name, name_fmt) 'jnpar_Bmod'
      call L1_write(parcurrn, datafile, trim(name), &
        'parallel current density' // comment, 's^-1')  ! SI: H^-1
    end if
    if (present(parcurrmn)) then
      write (name, name_fmt) 'jmnpar_Bmod'
      call L1_poloidal_modes(parcurrmn, polmodes)
      call polmodes_write(polmodes, datafile, trim(name), &
        'poloidal modes of parallel current density' // comment, 's^-1')  ! SI: H^-1
    end if
    if (present(currn)) then
      write (name, name_fmt) 'jn'
      call RT0_write(currn, datafile, trim(name), &
        'current density' // comment, 'statA cm^-2', 1)
    end if
    if (present(currmn)) then
      write (name, name_fmt) 'jmn'
      call RT0_poloidal_modes(currmn, vec_polmodes)
      call vec_polmodes_write(vec_polmodes, datafile, trim(name), &
        'poloidal modes of current density' // comment, 'statA cm^-2')
    end if
    if (present(magfn)) then
      write (name, name_fmt) 'Bn'
      call RT0_write(magfn, datafile, trim(name), &
        'magnetic field' // comment, 'G', 1)
    end if
    if (present(magfmn)) then
      write (name, name_fmt) 'Bmn'
      call RT0_poloidal_modes(magfmn, vec_polmodes)
      call vec_polmodes_write(vec_polmodes, datafile, trim(name), &
        'poloidal modes of magnetic field' // comment, 'G')
    end if

    if (present(presmn) .or. present(parcurrmn)) then
      call polmodes_deinit(polmodes)
    end if
    if (present(currmn) .or. present(magfmn)) then
      call vec_polmodes_deinit(vec_polmodes)
    end if

    ! resonant currents
    if (present(Ires)) then
      allocate(Imn_res(mesh%m_res_min:mesh%m_res_max))
      do m = mesh%m_res_min, mesh%m_res_max
        call compute_Ires(cache%shielding(m)%sample_Ires, cache%shielding(m)%GL_weights, &
          Ires, m, Imn_res(m))
      end do
      write (name, name_fmt) 'Ires'
      call h5_open_rw(datafile, h5id_root)
      call h5_add(h5id_root, trim(name), Imn_res, lbound(Imn_res), ubound(Imn_res), &
        comment = 'resonant currents' // comment, unit = 'statA')
      call h5_close(h5id_root)
    end if
  end subroutine perteq_write

  subroutine resonance_sweep(file, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use netlib_mod, only: zeroin
    use mephit_conf, only: conf, logger
    use mephit_util, only: ev2erg, linspace, interp1d
    use mephit_mesh, only: equil, mesh, cache, fs, dens_e, temp_e, temp_i, Phi0, nu_i, nu_e, E_r
    use mephit_pert, only: L1_t, L1_init, L1_deinit, polmodes_t, polmodes_init, polmodes_deinit, &
      L1_sum_poloidal_modes, compute_Ires
    use mephit_flr2, only: flr2_t, flr2_init, flr2_deinit, flr2_coeff, flr2_response_current
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    type(polmodes_t) :: jmnpar_over_Bmod
    type(L1_t) :: jnpar_over_Bmod
    type(flr2_t) :: flr2
    real(dp) :: Delta_E_r(conf%resonance_sweep), E_r_zero(conf%resonance_sweep)
    complex(dp) :: Bmnpsi_over_B0phi(0:mesh%nflux), &
      Phi_mn(0:mesh%nflux, mesh%m_res_min:mesh%m_res_max), &
      Imn_res(mesh%m_res_min:mesh%m_res_max, conf%resonance_sweep)
    integer :: k, m, m_res, m_min, m_max, kf_min, kf_max
    integer(HID_T) :: h5id_root

    if (logger%info) then
      write (logger%msg, '("Start resonance sweep with ", i0, " steps.")') &
        conf%resonance_sweep
      call logger%write_msg
    end if
    Bmnpsi_over_B0phi = (1d0, 0d0)
    call polmodes_init(jmnpar_over_Bmod, conf%m_max, mesh%nflux)
    jmnpar_over_Bmod%coeff(:, :) = (0d0, 0d0)
    call L1_init(jnpar_over_Bmod, mesh%npoint)
    m_min = min(-equil%cocos%sgn_q * mesh%m_res_min, -equil%cocos%sgn_q * mesh%m_res_max)
    m_max = max(-equil%cocos%sgn_q * mesh%m_res_min, -equil%cocos%sgn_q * mesh%m_res_max)
    call flr2_init(flr2, mesh%nflux, m_min, m_max)
    kf_min = mesh%res_ind(mesh%m_res_min) - 1
    kf_max = mesh%res_ind(mesh%m_res_max)
    Delta_E_r = linspace(-E_r%y(kf_min), -E_r%y(kf_max), conf%resonance_sweep, 0, 0)
    do k = 1, conf%resonance_sweep
      E_r_zero(k) = zeroin(fs%psi(kf_min), fs%psi(kf_max), E_r_shifted, 1d-9)
      call flr2_coeff(flr2, m_min, m_max, mesh%n, mesh%nflux, conf%m_i, conf%Z_i, &
        mesh%R_O, fs%psi(1:), fs%q(1:), fs%F(1:), Phi0%y(1:) - fs%rsmall(1:) * Delta_E_r(k), mesh%avg_R2gradpsi2(1:), &
        dens_e%y(1:), temp_e%y(1:) * ev2erg, temp_i%y(1:) * ev2erg, nu_e%y(1:), nu_i%y(1:))
      do m = mesh%m_res_min, mesh%m_res_max
        m_res = -equil%cocos%sgn_q * m
        call flr2_response_current(flr2, 1, m_res, mesh%nflux, &
          Bmnpsi_over_B0phi(1:), jmnpar_over_Bmod%coeff(m_res, 1:), Phi_mn(1:, m))
      end do
      call L1_sum_poloidal_modes(jmnpar_over_Bmod, jnpar_over_Bmod)
      do m = mesh%m_res_min, mesh%m_res_max
        call compute_Ires(cache%shielding(m)%sample_Ires, cache%shielding(m)%GL_weights, &
          jnpar_over_Bmod, m, Imn_res(m, k))
      end do
    end do
    call polmodes_deinit(jmnpar_over_Bmod)
    call L1_deinit(jnpar_over_Bmod)
    call flr2_deinit(flr2)
    grp = trim(group)
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/Delta_E_r', Delta_E_r, lbound(Delta_E_r), ubound(Delta_E_r), &
      comment = 'shift of electrical field', unit = 'statV cm^-1')
    call h5_add(h5id_root, grp // '/E_r_zero', E_r_zero, lbound(E_r_zero), ubound(E_r_zero), &
      comment = 'zero of electrical field in terms of poloidal flux', unit = 'Mx')
    call h5_add(h5id_root, grp // '/Imn_res', Imn_res, lbound(Imn_res), ubound(Imn_res), &
      comment = 'resonant current', unit = 'statA')
    call h5_close(h5id_root)
    if (logger%info) then
      write (logger%msg, '("Finished resonance sweep with ", i0, " steps.")') &
        conf%resonance_sweep
      call logger%write_msg
    end if

  contains
    function E_r_shifted(psi)
      real(dp), intent(in) :: psi
      real(dp) :: E_r_shifted

      E_r_shifted = interp1d(fs%psi, E_r%y, psi, 3) + Delta_E_r(k)
    end function E_r_shifted
  end subroutine resonance_sweep

  subroutine check_furth(currn, Bn_plas)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_util, only: imun, clight
    use mephit_mesh, only: equil, fs_half, mesh, cache
    use mephit_pert, only: RT0_t, vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, &
      RT0_poloidal_modes
    type(RT0_t), intent(in) :: currn
    type(RT0_t), intent(in) :: Bn_plas
    type(vec_polmodes_t) :: Bmn_plas
    character(len = *), parameter :: dataset = 'debug_furth'
    integer(HID_T) :: h5id_root
    integer :: kf, kt, ktri, kilca_m_res
    complex(dp) :: sheet_flux(mesh%nflux)
    real(dp) :: k_z, k_theta(mesh%nflux)

    call vec_polmodes_init(Bmn_plas, conf%m_max, mesh%nflux)
    call RT0_poloidal_modes(Bn_plas, Bmn_plas)
    kilca_m_res = -equil%cocos%sgn_q * abs(conf%kilca_pol_mode)
    k_z = mesh%n / mesh%R_O
    k_theta(:) = kilca_m_res / fs_half%rad
    sheet_flux(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        associate (s => cache%sample_polmodes_half(ktri))
          sheet_flux(kf) = sheet_flux(kf) + mesh%area(s%ktri) * currn%comp_phi(s%ktri) * &
            exp(-imun * kilca_m_res * s%theta)
        end associate
      end do
      sheet_flux(kf) = sheet_flux(kf) / dble(mesh%kt_max(kf))
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
    call vec_polmodes_deinit(Bmn_plas)
  end subroutine check_furth

end module mephit_iter
