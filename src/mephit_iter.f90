module mephit_iter
  use iso_fortran_env, only: dp => real64
  use mephit_pert, only: L1_t, RT0_t

  implicit none

  private

  public :: mephit_run, mephit_deinit

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
     subroutine FEM_init(tormode, nedge, runmode) bind(C, name = 'FEM_init')
       use iso_c_binding, only: c_int
       integer(c_int), intent(in), value :: tormode, nedge, runmode
     end subroutine FEM_init

     subroutine FEM_extend_mesh() bind(C, name = 'FEM_extend_mesh')
     end subroutine FEM_extend_mesh

     subroutine FEM_compute_Bn(nedge, Jn, Bn) bind(C, name = 'FEM_compute_Bn')
       use iso_c_binding, only: c_int, c_double_complex
       integer(c_int), intent(in), value :: nedge
       complex(c_double_complex), intent(in) :: Jn(1:nedge)
       complex(c_double_complex), intent(out) :: Bn(1:nedge)
     end subroutine FEM_compute_Bn

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

  subroutine mephit_run(runmode, config) bind(C, name = 'mephit_run')
    use iso_c_binding, only: c_int, c_ptr
    use magdata_in_symfluxcoor_mod, only: load_magdata_in_symfluxcoord
    use mephit_util, only: C_F_string, get_field_filenames, init_field
    use mephit_conf, only: conf, config_read, config_export_hdf5, conf_arr, logger, datafile
    use mephit_mesh, only: equil, mesh, generate_mesh, mesh_write, mesh_read, write_cache, read_cache, psi_interpolator
    use mephit_pert, only: generate_vacfield, vac, vac_init, vac_write, vac_read
    use hdf5_tools, only: h5_init, h5overwrite
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
    call config_read(conf, config_filename)
    call logger%init('-', conf%log_level, conf%quiet)
    call h5_init
    h5overwrite = .true.
    call config_export_hdf5(conf, datafile, 'config')
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
       call mesh_write(mesh, datafile, 'mesh')
       call write_cache
       call vac_init(vac, mesh%nedge, mesh%ntri, mesh%m_res_min, mesh%m_res_max)
       call generate_vacfield(vac)
       call vac_write(vac, datafile, 'vac')
       ! pass effective toroidal mode number and runmode to FreeFem++
       call FEM_init(mesh%n, mesh%nedge, runmode_flags)
       call FEM_extend_mesh
    else
       ! initialize equilibrium field
       call equil%import_hdf5(datafile, 'equil')
       call init_field(equil)
       call psi_interpolator%init(4, equil%psi_eqd)
       ! read in preprocessed data
       call mesh_read(mesh, datafile, 'mesh')
       call read_cache
       call load_magdata_in_symfluxcoord
       call vac_init(vac, mesh%nedge, mesh%ntri, mesh%m_res_min, mesh%m_res_max)
       call vac_read(vac, datafile, 'vac')
       ! reload config parameters here in case they changed since the meshing phase
       call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
       call conf_arr%export_hdf5(datafile, 'config')
       ! pass effective toroidal mode number and runmode to FreeFem++
       call FEM_init(mesh%n, mesh%nedge, runmode)
    end if
    if (iterations .or. analysis) then
       call iter_init
       if (iterations) then
          call MFEM_test
          call mephit_iterate
          call FEM_deinit
       else
          if (analysis) then
             call iter_read
          end if
          ! FEM_deinit is not needed because scripts/maxwell_daemon.edp exits
          ! if iterations are not requested
       end if
       if (analysis) then
          call mephit_postprocess
       end if
       call iter_deinit
    end if
    call mephit_deinit
  end subroutine mephit_run

  subroutine mephit_deinit
    use magdata_in_symfluxcoor_mod, only: unload_magdata_in_symfluxcoord
    use hdf5_tools, only: h5_deinit
    use mephit_conf, only: conf_arr, logger
    use mephit_util, only: deinit_field
    use mephit_mesh, only: equil, fs, fs_half, psi_interpolator, psi_fine_interpolator, &
         mesh, cache, mesh_deinit, cache_deinit
    use mephit_pert, only: vac, vac_deinit

    call psi_interpolator%deinit
    call psi_fine_interpolator%deinit
    call cache_deinit(cache)
    call fs%deinit
    call fs_half%deinit
    call mesh_deinit(mesh)
    call unload_magdata_in_symfluxcoord
    call vac_deinit(vac)
    call deinit_field
    call equil%deinit
    call conf_arr%deinit
    call logger%deinit
    call h5_deinit
  end subroutine mephit_deinit

  subroutine iter_init
    use mephit_mesh, only: mesh
    use mephit_pert, only: RT0_init, L1_init

    call L1_init(pn, mesh%npoint)
    call RT0_init(Bn, mesh%nedge, mesh%ntri)
    call RT0_init(Bnplas, mesh%nedge, mesh%ntri)
    call RT0_init(jn, mesh%nedge, mesh%ntri)
    call L1_init(jnpar_B0, mesh%npoint)
  end subroutine iter_init

  subroutine iter_deinit
    use mephit_pert, only: L1_deinit, RT0_deinit

    call L1_deinit(pn)
    call RT0_deinit(Bn)
    call RT0_deinit(Bnplas)
    call RT0_deinit(jn)
    call L1_deinit(jnpar_B0)
  end subroutine iter_deinit

  subroutine iter_read
    use mephit_conf, only: datafile
    use mephit_pert, only: L1_read, RT0_read

    call L1_read(pn, datafile, 'iter/pn')
    call RT0_read(Bn, datafile, 'iter/Bn')
    call RT0_read(Bnplas, datafile, 'iter/Bnplas')
    call RT0_read(jn, datafile, 'iter/jn')
    call L1_read(jnpar_B0, datafile, 'iter/jnpar_B0')
  end subroutine iter_read

  subroutine iter_step
    use mephit_conf, only: conf, currn_eps, currn_gl, currn_mde, logger
    ! compute pressure based on previous perturbation field
    call compute_presn
    ! compute currents based on previous perturbation field
    select case (conf%currn)
    case (currn_eps)
       call compute_currn
    case (currn_gl)
       call compute_currn_GL
    case (currn_mde)
       call compute_currn_MDE
    case default
       write (logger%msg, '("unknown current perturbation selection: ", i0)') conf%curr_prof
       if (logger%err) call logger%write_msg
       error stop
    end select
    ! use field code to generate new field from currents
    call compute_Bn
  end subroutine iter_step

  subroutine mephit_iterate
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_util, only: arnoldi_break
    use mephit_mesh, only: mesh
    use mephit_pert, only: L1_write, &
         RT0_init, RT0_deinit, RT0_write, RT0_compute_tor_comp, RT0_L2int, &
         vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, vec_polmodes_write, &
         polmodes_t, polmodes_init, polmodes_deinit, polmodes_write, &
         L1_poloidal_modes, RT0_poloidal_modes, vac
    use mephit_conf, only: conf, logger, runmode_precon, runmode_single, datafile, cmplx_fmt

    logical :: preconditioned
    integer :: kiter, niter, maxiter, ndim, nritz, i, j, info
    integer(HID_T) :: h5id_root
    real(dp), allocatable :: L2int_Bn_diff(:)
    real(dp) :: L2int_Bnvac, rel_err
    type(RT0_t) :: Bn_prev, Bn_diff
    complex(dp), allocatable :: eigvals(:), eigvecs(:, :), Lr(:, :), Yr(:, :)
    integer, allocatable :: ipiv(:)
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"
    integer, parameter :: m_max = 24
    type(polmodes_t) :: pmn
    type(vec_polmodes_t) :: jmn

    ! system dimension: number of non-redundant edges in core plasma
    ndim = mesh%nedge
    ! runmodes
    preconditioned = runmode_precon == conf%runmode
    if (runmode_single == conf%runmode) then
       maxiter = 0
    else
       maxiter = conf%niter
    end if
    if (preconditioned) then
       ! calculate eigenvectors
       call arnoldi_break(ndim, conf%nkrylov, 3, conf%ritz_threshold, conf%ritz_rel_err, &
            next_iteration_arnoldi, info, nritz, eigvals, eigvecs)
       if (info /= 0) then
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
               nritz, conf%ritz_threshold
          call logger%write_msg
          do i = 1, nritz
             write (logger%msg, '("lambda ", i0, ": ", ' // cmplx_fmt // ')') i, eigvals(i)
             call logger%write_msg
          end do
       end if
       if (nritz > 0) then
          do i = 1, min(nritz, conf%max_eig_out)
             write (postfix, postfix_fmt) i
             Bn%DOF(:) = eigvecs(:, i)
             call RT0_compute_tor_comp(Bn)
             call RT0_write(Bn, datafile, 'iter/eigvec' // postfix, &
                  'iteration eigenvector', 'G', 1)
          end do
          allocate(Lr(nritz, nritz), Yr(nritz, nritz))
          Yr = (0d0, 0d0)
          do i = 1, nritz
             Yr(i, i) = (1d0, 0d0)
             do j = 1, nritz
                Lr(i, j) = sum(conjg(eigvecs(:, i)) * eigvecs(:, j)) * (eigvals(j) - (1d0, 0d0))
             end do
          end do
          allocate(ipiv(nritz))
          call zgesv(nritz, nritz, Lr, nritz, ipiv, Yr, nritz, info)
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
          do i = 1, nritz
             Lr(i, :) = eigvals(i) * Yr(i, :)
          end do
          deallocate(Yr)
       else
          preconditioned = .false.
       end if
    end if

    L2int_Bnvac = RT0_L2int(vac%Bn)
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, 'iter/')
    call h5_add(h5id_root, 'iter/L2int_Bnvac', L2int_Bnvac, &
         comment = 'L2 integral of magnetic field (vacuum)', unit = 'Mx')
    call h5_close(h5id_root)
    allocate(L2int_Bn_diff(0:maxiter))
    L2int_Bn_diff = ieee_value(0d0, ieee_quiet_nan)
    call polmodes_init(pmn, m_max, mesh%nflux)
    call vec_polmodes_init(jmn, m_max, mesh%nflux)
    call RT0_init(Bn_prev, mesh%nedge, mesh%ntri)
    call RT0_init(Bn_diff, mesh%nedge, mesh%ntri)
    Bn%DOF(:) = vac%Bn%DOF
    Bn%comp_phi(:) = vac%Bn%comp_phi
    if (preconditioned) then
       Bn%DOF(:) = Bn%DOF - matmul(eigvecs(:, 1:nritz), matmul(Lr, &
            matmul(transpose(conjg(eigvecs(:, 1:nritz))), Bn%DOF)))
       call RT0_compute_tor_comp(Bn)
    end if
    niter = maxiter
    do kiter = 0, maxiter
       write (logger%msg, '("Iteration ", i2, " of ", i2)') kiter, maxiter
       if (logger%info) call logger%write_msg
       write (postfix, postfix_fmt) kiter
       Bn_prev%DOF(:) = Bn%DOF
       Bn_prev%comp_phi(:) = Bn%comp_phi
       ! compute B_(n+1) = K * B_n + B_vac ... different from next_iteration_arnoldi
       call iter_step
       Bn%DOF(:) = Bn%DOF + vac%Bn%DOF
       if (preconditioned) then
          Bn%DOF(:) = Bn%DOF - matmul(eigvecs(:, 1:nritz), matmul(Lr, &
               matmul(transpose(conjg(eigvecs(:, 1:nritz))), Bn%DOF - Bn_prev%DOF)))
       end if
       call RT0_compute_tor_comp(Bn)
       Bn_diff%DOF(:) = Bn%DOF - Bn_prev%DOF
       Bn_diff%comp_phi(:) = Bn%comp_phi - Bn_prev%comp_phi
       L2int_Bn_diff(kiter) = RT0_L2int(Bn_diff)
       if (kiter <= 1) then
          call L1_write(pn, datafile, 'iter/pn' // postfix, &
               'pressure (after iteration)', 'dyn cm^-2')
          call L1_write(pn, datafile, 'iter/jnpar_B0' // postfix, &
               'parallel current density (after iteration)', 's^-1')  ! SI: H^-1
          call RT0_write(jn, datafile, 'iter/jn' // postfix, &
               'current density (after iteration)', 'statA cm^-2', 1)
          call RT0_write(Bn, datafile, 'iter/Bn' // postfix, &
               'magnetic field (after iteration)', 'G', 1)
          call RT0_write(Bn_diff, datafile, 'iter/Bn_diff' // postfix, &
               'magnetic field (difference between iterations)', 'G', 1)
          call L1_poloidal_modes(pn, pmn)
          call polmodes_write(pmn, datafile, 'postprocess/pmn' // postfix, &
               'pressure (after iteration)', 'dyn cm^-2')
          call RT0_poloidal_modes(jn, jmn)
          call vec_polmodes_write(jmn, datafile, 'postprocess/jmn' // postfix, &
               'current density (after iteration)', 'statA cm^-2')
       else
          if (L2int_Bn_diff(kiter) > conf%ritz_threshold ** kiter * L2int_Bnvac .or. &
               L2int_Bn_diff(kiter) < conf%iter_rel_err * L2int_Bnvac) then
             niter = kiter
             exit
          end if
       end if
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
    if (preconditioned) then
       deallocate(Lr, eigvals, eigvecs)
    end if
    call debug_MDE("debug_MDE", pn, Bn, jn, jnpar_B0)
    call debug_currn("debug_currn", pn, Bn, jn)
    Bnplas%DOF(:) = Bn%DOF - vac%Bn%DOF
    Bnplas%comp_phi(:) = Bn%comp_phi - vac%Bn%comp_phi
    call h5_open_rw(datafile, h5id_root)
    call h5_add(h5id_root, 'iter/niter', niter, comment = 'actual number of iterations')
    call h5_add(h5id_root, 'iter/rel_err', rel_err, comment = 'relative error of iterations')
    call h5_add(h5id_root, 'iter/L2int_Bn_diff', L2int_Bn_diff, &
         lbound(L2int_Bn_diff), ubound(L2int_Bn_diff), &
         comment = 'L2 integral of magnetic field (difference between iterations)', unit = 'Mx')
    call h5_close(h5id_root)
    deallocate(L2int_Bn_diff)
    call L1_write(pn, datafile, 'iter/pn', &
         'pressure (full perturbation)', 'dyn cm^-2')
    call L1_write(pn, datafile, 'iter/jnpar_B0', &
         'parallel current density (after iteration)', 's^-1')  ! SI: H^-1
    call RT0_write(Bn, datafile, 'iter/Bn', &
         'magnetic field (full perturbation)', 'G', 2)
    call RT0_write(Bnplas, datafile, 'iter/Bnplas', &
         'magnetic field (plasma response)', 'G', 2)
    call RT0_write(jn, datafile, 'iter/jn', &
         'current density (full perturbation)', 'statA cm^-2', 1)
    call RT0_deinit(Bn_prev)
    call RT0_deinit(Bn_diff)
    call polmodes_deinit(pmn)
    call vec_polmodes_deinit(jmn)

  contains

    ! computes B_(n+1) = K * (B_n + B_vac) for Arnoldi iterations
    subroutine next_iteration_arnoldi(old_val, new_val)
      complex(dp), intent(in) :: old_val(:)
      complex(dp), intent(out) :: new_val(:)
      Bn%DOF(:) = old_val + vac%Bn%DOF
      call RT0_compute_tor_comp(Bn)
      call iter_step
      new_val(:) = Bn%DOF
    end subroutine next_iteration_arnoldi
  end subroutine mephit_iterate

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi.
  !>
  !> This subroutine calls a C function that pipes the data to/from FreeFem.
  subroutine compute_Bn
    use mephit_mesh, only: mesh
    use mephit_pert, only: RT0_compute_tor_comp

    call FEM_compute_Bn(mesh%nedge, jn%DOF, Bn%DOF)
    call RT0_compute_tor_comp(Bn)
  end subroutine compute_Bn

  subroutine unit_B0(R, Z, vector) bind(C, name = 'unit_B0')
    use iso_c_binding, only: c_double
    use mephit_conf, only: logger
    real(c_double), intent(in), value :: R
    real(c_double), intent(in), value :: Z
    real(c_double), intent(out) :: vector(3)
    real(dp) :: dum

    if (logger%debug) then
       logger%msg = 'called unit_B0'
       call logger%write_msg
    end if
    call field(R, 0d0, Z, vector(1), vector(3), vector(2), &
         dum, dum, dum, dum, dum, dum, dum, dum, dum)
    dum = sqrt(sum(vector * vector))
    vector = vector / dum
  end subroutine unit_B0

  subroutine presn_inhom(R, Z, scalar) bind(C, name = 'presn_inhom')
    use iso_c_binding, only: c_double, c_double_complex
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use field_eq_mod, only: psif, psib
    use mephit_conf, only: logger
    use mephit_mesh, only: equil, psi_interpolator, mesh, point_location
    use mephit_pert, only: RT0_interp, vac
    real(c_double), intent(in), value :: R
    real(c_double), intent(in), value :: Z
    complex(c_double_complex), intent(out) :: scalar
    real(dp) :: B_0(3), dum, psi, dp0_dpsi
    integer :: ktri
    complex(dp) :: B_n(3)

    if (logger%debug) then
       logger%msg = 'called presn_inhom'
       call logger%write_msg
    end if
    call field(R, 0d0, Z, B_0(1), B_0(2), B_0(3), &
         dum, dum, dum, dum, dum, dum, dum, dum, dum)
    psi = psif - psib  ! see intperp_psi_pol in mephit_util
    dp0_dpsi = psi_interpolator%eval(equil%pprime, psi)
    ktri = point_location(R, Z)
    if (ktri <= 0 .or. ktri > mesh%ntri) then
       scalar = cmplx(ieee_value(0d0, ieee_quiet_nan), ieee_value(0d0, ieee_quiet_nan), dp)
       return
    end if
    call RT0_interp(vac%Bn, ktri, R, Z, B_n)
    scalar = -dp0_dpsi * (B_n(1) * B_0(3) - B_n(3) * B_0(1)) * R / sqrt(sum(B_0 * B_0))
  end subroutine presn_inhom

  subroutine MFEM_test()
    use iso_c_binding, only: c_int, c_null_char, c_loc, c_funloc
    use mephit_conf, only: conf, logger
    character(len = *), parameter :: mesh_file = 'core_plasma.mesh'
    integer(c_int) :: status

    status = FEM_test(trim(mesh_file) // c_null_char, conf%n, size(pn%DOF), &
         pn%DOF, c_funloc(unit_B0), c_funloc(presn_inhom))
    if (logger%debug) then
       write (logger%msg, '("FEM_test return status: ", i0)') status
       call logger%write_msg
    end if
  end subroutine MFEM_test

  subroutine solve_MDE(inhom, solution)
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: conf, logger
    use mephit_util, only: imun
    use mephit_mesh, only: mesh, cache
    complex(dp), dimension(:), intent(in) :: inhom
    complex(dp), dimension(:), intent(out) :: solution
    complex(dp), dimension(maxval(mesh%kp_max)) :: a, b, x, d, du
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(maxval(mesh%kp_max)) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kp, kedge, k
    integer, dimension(2 * maxval(mesh%kp_max)) :: irow, icol
    complex(dp), dimension(2 * maxval(mesh%kp_max)) :: aval
    real(dp), parameter :: small = tiny(0d0)

    if (size(inhom) /= mesh%npoint) then
       call logger%msg_arg_size('solve_MDE', 'size(inhom)', 'mesh%npoint', size(inhom), mesh%npoint)
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (size(solution) /= mesh%npoint) then
       call logger%msg_arg_size('solve_MDE', 'size(solution)', 'mesh%npoint', size(solution), mesh%npoint)
       if (logger%err) call logger%write_msg
       error stop
    end if
    max_rel_err = 0d0
    avg_rel_err = 0d0
    solution(:) = (0d0, 0d0)
    a = (0d0, 0d0)
    b = (0d0, 0d0)
    x = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ! use midpoint of poloidal edge
          associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge))
            a(kp) = (f%B0(1) * mesh%edge_R(kedge) + f%B0(3) * mesh%edge_Z(kedge)) / &
                 (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
            b(kp) = imun * (mesh%n + imun * conf%damp) * f%B0(2) / R
            x(kp) = inhom(kedge + 1)
          end associate
       end do
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       associate (ndim => mesh%kp_max(kf), nz => 2 * mesh%kp_max(kf), &
            k_min => mesh%kp_low(kf) + 1, k_max => mesh%kp_low(kf) + mesh%kp_max(kf))
         ! assemble sparse matrix (COO format)
         ! first column, diagonal
         irow(1) = 1
         icol(1) = 1
         aval(1) = d(1)
         ! first column, off-diagonal
         irow(2) = ndim
         icol(2) = 1
         aval(2) = du(ndim)
         do k = 2, ndim
            ! off-diagonal
            irow(2*k-1) = k-1
            icol(2*k-1) = k
            aval(2*k-1) = du(k-1)
            ! diagonal
            irow(2*k) = k
            icol(2*k) = k
            aval(2*k) = d(k)
         end do
         call sparse_solve(ndim, ndim, nz, irow(:nz), icol(:nz), aval(:nz), x(:ndim))
         call sparse_matmul(ndim, ndim, irow(:nz), icol(:nz), aval(:nz), x(:ndim), resid)
         resid(:) = resid - inhom(k_min:k_max)
         where (abs(inhom(k_min:k_max)) >= small)
            rel_err(:ndim) = abs(resid(:ndim)) / abs(inhom(k_min:k_max))
         elsewhere
            rel_err(:ndim) = 0d0
         end where
         max_rel_err = max(max_rel_err, maxval(rel_err(:ndim)))
         avg_rel_err = avg_rel_err + sum(rel_err(:ndim))
       end associate
       if (kf == 1) then ! first point on axis - average over enclosing flux surface
          solution(1) = sum(x(:mesh%kp_max(1))) / dble(mesh%kp_max(1))
       end if
       do kp = 1, mesh%kp_max(kf)
          solution(mesh%kp_low(kf) + kp) = x(kp)
       end do
    end do
    avg_rel_err = avg_rel_err / dble(mesh%npoint - 1)
    write (logger%msg, '("solve_MDE: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (logger%debug) call logger%write_msg
    if (allocated(resid)) deallocate(resid)
  end subroutine solve_MDE

  subroutine compute_presn
    use mephit_mesh, only: fs, mesh, cache
    complex(dp), dimension(mesh%npoint) :: inhom
    integer :: kf, kp, kedge

    inhom = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ! use midpoint of poloidal edge
          associate (f => cache%mid_fields(kedge))
            inhom(kedge + 1) = -fs%dp_dpsi(kf) * Bn%DOF(kedge) * &
                 (f%B0(1) * mesh%edge_R(kedge) + f%B0(3) * mesh%edge_Z(kedge)) / &
                 (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
          end associate
       end do
    end do
    call solve_MDE(inhom, pn%DOF)
  end subroutine compute_presn

  subroutine debug_MDE(group, presn, magfn, currn_perp, currn_par)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use field_eq_mod, only: psif, psib
    use mephit_conf, only: conf, datafile
    use mephit_util, only: imun, pi, clight, zd_cross
    use mephit_mesh, only: equil, mesh, cache, psi_interpolator, point_location
    use mephit_pert, only: L1_t, L1_interp, RT0_t, RT0_interp
    character(len = *), intent(in) :: group
    type(L1_t), intent(in) :: presn
    type(RT0_t), intent(in) :: magfn
    type(RT0_t), intent(in) :: currn_perp
    type(L1_t), intent(in) :: currn_par
    integer(HID_T) :: h5id_root
    character(len = len_trim(group)) :: grp
    integer :: ndim, kf, kp, kedge, ktri, k
    real(dp) :: dum, psi, B_0(3), j_0(3), stencil_R(5), stencil_Z(5), &
         dp0_dpsi, dF_dpsi, FdF_dpsi, fprime(equil%nw)
    complex(dp) :: zdum, B_n(3), jnperp(3, 5), grad_p_n(3)
    complex(dp), allocatable :: grad_pn(:, :), Bn_psi_contravar(:), grad_jnpar(:, :), &
         div_jnperp(:), div_jnperp_RT0(:)

    grp = trim(adjustl(group))
    fprime = equil%ffprim / equil%fpol
    ndim = mesh%npoint - 1
    allocate(grad_pn(3, ndim), Bn_psi_contravar(ndim), grad_jnpar(3, ndim), &
         div_jnperp(ndim), div_jnperp_RT0(ndim))
    div_jnperp(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ktri = mesh%edge_tri(1, kedge)
          associate (R => mesh%mid_R(kedge), Z => mesh%mid_Z(kedge), c => cache%mid_fields(kedge))
            call L1_interp(presn, ktri, R, Z, zdum, grad_pn(:, kedge))
            call RT0_interp(magfn, ktri, R, Z, B_n)
            call L1_interp(currn_par, ktri, R, Z, zdum, grad_jnpar(:, kedge))
            call RT0_interp(currn_perp, ktri, R, Z, jnperp(:, 4), &
                 jnperp(:, 1), jnperp(:, 2), jnperp(:, 3))
            Bn_psi_contravar(kedge) = R * (B_n(1) * c%B0(3) - B_n(3) * c%B0(1))
            div_jnperp_RT0(kedge) = jnperp(1, 1) + jnperp(2, 2) + jnperp(3, 3) + &
                 jnperp(1, 4) / R
          end associate
          if (kf < mesh%nflux) then
             stencil_R = mesh%mid_R(kedge) + [0, 1, -1, 0, 0] * 0.2d0 * conf%max_Delta_rad
             stencil_Z = mesh%mid_Z(kedge) + [0, 0, 0, 1, -1] * 0.2d0 * conf%max_Delta_rad
             do k = 1, 5
                call field(stencil_R(k), 0d0, stencil_Z(k), B_0(1), B_0(2), B_0(3), &
                     dum, dum, dum, dum, dum, dum, dum, dum, dum)
                psi = psif - psib  ! see intperp_psi_pol in mephit_util
                dp0_dpsi = psi_interpolator%eval(equil%pprime, psi)
                FdF_dpsi = psi_interpolator%eval(equil%ffprim, psi)
                dF_dpsi = psi_interpolator%eval(fprime, psi)
                j_0(1) = 0.25d0 / pi * clight * dF_dpsi * B_0(1)
                j_0(2) = clight * (dp0_dpsi * stencil_R(k) + 0.25d0 / (pi * stencil_R(k)) * FdF_dpsi)
                j_0(3) = 0.25d0 / pi * clight * dF_dpsi * B_0(3)
                ktri = point_location(stencil_R(k), stencil_Z(k), psi)
                call L1_interp(presn, ktri, stencil_R(k), stencil_Z(k), zdum, grad_p_n)
                call RT0_interp(magfn, ktri, stencil_R(k), stencil_Z(k), B_n)
                jnperp(:, k) = (sum(B_0 * J_0) * B_n - sum(B_0 * B_n) * J_0 + clight * &
                     zd_cross(grad_p_n, B_0)) / sum(B_0 ** 2)
             end do
             div_jnperp(kedge) = imun * mesh%n / stencil_R(1) * jnperp(2, 1) + &
                  (jnperp(1, 2) - jnperp(1, 3) + jnperp(3, 4) - jnperp(3, 5)) / &
                  (0.4d0 * conf%max_Delta_rad) + jnperp(1, 1) / stencil_R(1)
          end if
       end do
    end do

    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/grad_pn', grad_pn, lbound(grad_pn), ubound(grad_pn), &
         comment = 'perturbation pressure gradient', unit = 'dyn cm^-2')
    call h5_add(h5id_root, grp // '/Bn_psi_contravar', Bn_psi_contravar, &
         lbound(Bn_psi_contravar), ubound(Bn_psi_contravar), &
         comment = 'perturbation pressure gradient', unit = 'G^2 cm')
    call h5_add(h5id_root, grp // '/grad_jnpar', grad_jnpar, &
         lbound(grad_jnpar), ubound(grad_jnpar), unit = 'cm^-1 s^-1', &
         comment = 'parallel perturbation current density gradient')
    call h5_add(h5id_root, grp // '/div_jnperp', div_jnperp, &
         lbound(div_jnperp), ubound(div_jnperp), unit = 'statA cm^-3', &
         comment = 'perpendicular perturbation current density divergence (five-point stencil)')
    call h5_add(h5id_root, grp // '/div_jnperp_RT0', div_jnperp_RT0, &
         lbound(div_jnperp_RT0), ubound(div_jnperp_RT0), unit = 'statA cm^-3', &
         comment = 'perpendicular perturbation current density divergence (RT0 interpolation)')
    call h5_close(h5id_root)

    deallocate(grad_pn, Bn_psi_contravar, grad_jnpar, div_jnperp, div_jnperp_RT0)
  end subroutine debug_MDE

  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #mephit_conf::currn_file.
  subroutine compute_currn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: conf, logger
    use mephit_mesh, only: fs, mesh, cache
    use mephit_pert, only: RT0_compute_tor_comp
    use mephit_util, only: imun, clight
    complex(dp), dimension(maxval(mesh%kt_max)) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(maxval(mesh%kt_max)) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kp, kt, ktri, ke, kedge, k
    integer, dimension(2 * maxval(mesh%kt_max)) :: irow, icol
    complex(dp), dimension(2 * maxval(mesh%kt_max)) :: aval
    complex(dp) :: Bnphi_edge, Delta_pn
    real(dp) :: Delta_p0
    real(dp), parameter :: small = tiny(0d0)
    ! hack: first call in mephit_iter() is always without plasma response
    ! and without additional shielding currents
    logical, save :: first_call = .true.

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          associate (f => cache%mid_fields(kedge), nodes => mesh%edge_node(:, kedge))
            jn%DOF(kedge) = (f%j0(2) * Bn%DOF(kedge) + clight * mesh%mid_R(kedge) * &
                 (pn%DOF(nodes(2)) - pn%DOF(nodes(1)))) / f%B0(2)
          end associate
       end do
    end do
    do kf = 1, mesh%nflux
       Delta_p0 = fs%p(kf) - fs%p(kf-1)
       x = (0d0, 0d0)
       do kt = 1, mesh%kt_max(kf)
          ! radial edge
          kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
          Bnphi_edge = 0.5d0 * sum(Bn%comp_phi(mesh%edge_tri(:, kedge)))
          Delta_pn = pn%DOF(mesh%edge_node(2, kedge)) - pn%DOF(mesh%edge_node(1, kedge))
          ! first term on source side: flux through edge f
          ktri = mesh%edge_tri(2, kedge)
          ke = kt
          if (mesh%orient(ktri)) then
             x(ke) = x(ke) - jn%DOF(mesh%tri_edge(1, ktri))
          else
             x(ke) = x(ke) + jn%DOF(mesh%tri_edge(1, ktri))
          end if
          associate (f => cache%mid_fields(kedge), B0_flux => cache%B0_flux(kedge))
            ! diagonal matrix element - edge i
            d(ke) = -1d0 - imun * (mesh%n + imun * conf%damp) * &
                 mesh%area(ktri) * 0.5d0 * f%B0(2) / B0_flux
            ! additional term from edge i on source side
            x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * mesh%mid_R(kedge) / B0_flux * &
                 (Bnphi_edge / f%B0(2) * Delta_p0 - Delta_pn) + f%j0(2) * &
                 (Bnphi_edge / f%B0(2) - Bn%DOF(kedge) / B0_flux))
            ! superdiagonal matrix element - edge o
            ktri = mesh%edge_tri(1, kedge)
            ke = mod(kt + mesh%kt_max(kf) - 2, mesh%kt_max(kf)) + 1
            du(ke) = 1d0 + imun * (mesh%n + imun * conf%damp) * &
                 mesh%area(ktri) * 0.5d0 * f%B0(2) / (-B0_flux)
            ! additional term from edge o on source side
            x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * mesh%mid_R(kedge) / (-B0_flux) * &
                 (Bnphi_edge / f%B0(2) * (-Delta_p0) - (-Delta_pn)) + f%j0(2) * &
                 (Bnphi_edge / f%B0(2) - (-Bn%DOF(kedge)) / (-B0_flux)))
          end associate
       end do
       associate (ndim => mesh%kt_max(kf), nz => 2 * mesh%kt_max(kf))
         ! assemble sparse matrix (COO format)
         ! first column, diagonal
         irow(1) = 1
         icol(1) = 1
         aval(1) = d(1)
         ! first column, off-diagonal
         irow(2) = ndim
         icol(2) = 1
         aval(2) = du(ndim)
         do k = 2, ndim
            ! off-diagonal
            irow(2*k-1) = k-1
            icol(2*k-1) = k
            aval(2*k-1) = du(k-1)
            ! diagonal
            irow(2*k) = k
            icol(2*k) = k
            aval(2*k) = d(k)
         end do
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
          kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
          jn%DOF(kedge) = -x(kt)
       end do
    end do
    avg_rel_err = avg_rel_err / dble(mesh%ntri)
    write (logger%msg, '("compute_currn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (logger%debug) call logger%write_msg
    if (allocated(resid)) deallocate(resid)

    if (first_call) then
       first_call = .false.
       call RT0_compute_tor_comp(jn)
       call debug_currn("debug_currn_000", pn, Bn, jn)
    end if
    call add_sheet_current
    call RT0_compute_tor_comp(jn)
  end subroutine compute_currn

  subroutine compute_currn_GL
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: logger
    use mephit_mesh, only: mesh, cache
    use mephit_pert, only: L1_interp, RT0_interp, RT0_compute_tor_comp
    use mephit_util, only: imun, clight
    real(dp), parameter :: small = tiny(0d0)
    ! hack: first call in mephit_iter() is always without plasma response
    ! and without additional shielding currents
    logical, save :: first_call = .true.
    integer :: kf, kp, kt, ktri, kedge, k, nodes(3)
    complex(dp) :: series, zdum, grad_pn(3), B_n(3)
    complex(dp), dimension(maxval(mesh%kt_max)) :: x, d, du, inhom
    integer, dimension(2 * maxval(mesh%kt_max)) :: irow, icol
    complex(dp), dimension(2 * maxval(mesh%kt_max)) :: aval
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(maxval(mesh%kt_max)) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err

    max_rel_err = 0d0
    avg_rel_err = 0d0
    jn%DOF(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ktri = mesh%edge_tri(1, kedge)
          do k = 1, mesh%GL_order
             call L1_interp(pn, ktri, mesh%GL_R(k, kedge), mesh%GL_Z(k, kedge), zdum, grad_pn)
             call RT0_interp(Bn, ktri, mesh%GL_R(k, kedge), mesh%GL_Z(k, kedge), B_n)
             associate (f => cache%edge_fields(k, kedge))
               jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * mesh%GL_R(k, kedge) * &
                    ((clight * grad_pn(3) + f%j0(2) * B_n(1)) * mesh%edge_Z(kedge) - &
                    (-clight * grad_pn(1) + f%j0(2) * B_n(3)) * mesh%edge_R(kedge)) / f%B0(2)
             end associate
          end do
       end do
    end do
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ! indices F, O, I of nodes as defined in mephit_pert::RT0_interp
          if (mesh%orient(ktri)) then
             nodes = mesh%tri_node([3, 1, 2], ktri)
          else
             nodes = mesh%tri_node([1, 3, 2], ktri)
          end if
          ! area contribution to inhomogeneity
          series = (0d0, 0d0)
          do k = 1, mesh%GL2_order
             call L1_interp(pn, ktri, mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), zdum, grad_pn)
             call RT0_interp(Bn, ktri, mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), B_n)
             associate (f => cache%area_fields(k, ktri))
               series = series + mesh%GL2_weights(k) / (f%B0(1) ** 2 + f%B0(3) ** 2) * &
                    ((B_n(2) * f%j0(1) - f%j0(2) * B_n(1) - clight * grad_pn(3)) * f%B0(1) + &
                    (B_n(2) * f%j0(3) - f%j0(2) * B_n(3) + clight * grad_pn(1)) * f%B0(3))
             end associate
          end do
          x(kt) = -imun * mesh%n * mesh%area(ktri) * series
          ! edge f contribution to inhomogeneity
          kedge = mesh%tri_edge(1, ktri)
          series = (0d0, 0d0)
          do k = 1, mesh%GL2_order
             associate (f => cache%area_fields(k, ktri))
               series = series + mesh%GL2_weights(k) * f%B0(2) / mesh%GL2_R(k, ktri) * &
                    ((mesh%GL2_R(k, ktri) - mesh%node_R(nodes(1))) * f%B0(1) + &
                    (mesh%GL2_Z(k, ktri) - mesh%node_Z(nodes(1))) * f%B0(3)) / &
                    (f%B0(1) ** 2 + f%B0(3) ** 2)
             end associate
          end do
          if (mesh%orient(ktri)) then
             x(kt) = x(kt) - jn%DOF(ktri) * (1d0 + imun * mesh%n * 0.5d0 * series)
          else
             x(kt) = x(kt) + jn%DOF(ktri) * (1d0 + imun * mesh%n * 0.5d0 * series)
          end if
          ! edge i contribution to diagonal element
          kedge = mesh%tri_edge(3, ktri)
          series = (0d0, 0d0)
          do k = 1, mesh%GL2_order
             associate (f => cache%area_fields(k, ktri))
               series = series + mesh%GL2_weights(k) * f%B0(2) / mesh%GL2_R(k, ktri) * &
                    ((mesh%GL2_R(k, ktri) - mesh%node_R(nodes(3))) * f%B0(1) + &
                    (mesh%GL2_Z(k, ktri) - mesh%node_Z(nodes(3))) * f%B0(3)) / &
                    (f%B0(1) ** 2 + f%B0(3) ** 2)
             end associate
          end do
          d(kt) = -1d0 - imun * mesh%n * 0.5d0 * series
          ! edge o contribution to upper diagonal element
          kedge = mesh%tri_edge(2, ktri)
          series = (0d0, 0d0)
          do k = 1, mesh%GL2_order
             associate (f => cache%area_fields(k, ktri))
               series = series + mesh%GL2_weights(k) * f%B0(2) / mesh%GL2_R(k, ktri) * &
                    ((mesh%GL2_R(k, ktri) - mesh%node_R(nodes(2))) * f%B0(1) + &
                    (mesh%GL2_Z(k, ktri) - mesh%node_Z(nodes(2))) * f%B0(3)) / &
                    (f%B0(1) ** 2 + f%B0(3) ** 2)
             end associate
          end do
          du(kt) = 1d0 + imun * mesh%n * 0.5d0 * series
       end do
       associate (ndim => mesh%kt_max(kf), nz => 2 * mesh%kt_max(kf))
         ! assemble sparse matrix (COO format)
         ! first column, diagonal
         irow(1) = 1
         icol(1) = 1
         aval(1) = d(1)
         ! first column, off-diagonal
         irow(2) = ndim
         icol(2) = 1
         aval(2) = du(ndim)
         do k = 2, ndim
            ! off-diagonal
            irow(2*k-1) = k-1
            icol(2*k-1) = k
            aval(2*k-1) = du(k-1)
            ! diagonal
            irow(2*k) = k
            icol(2*k) = k
            aval(2*k) = d(k)
         end do
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
          kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
          jn%DOF(kedge) = -x(kt)
       end do
    end do
    avg_rel_err = avg_rel_err / dble(mesh%ntri)
    write (logger%msg, '("compute_currn_GL: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (logger%debug) call logger%write_msg
    if (allocated(resid)) deallocate(resid)

    if (first_call) then
       first_call = .false.
       call RT0_compute_tor_comp(jn)
       call debug_currn("debug_currn_000", pn, Bn, jn)
    end if
    call add_sheet_current
    call RT0_compute_tor_comp(jn)
  end subroutine compute_currn_GL

  subroutine compute_currn_MDE
    use mephit_mesh, only: mesh, cache
    use mephit_pert, only: L1_interp, RT0_interp
    use mephit_util, only: pi, clight, zd_cross
    ! hack: first call in mephit_iter() is always without plasma response
    ! and without additional shielding currents
    logical, save :: first_call = .true.
    integer :: kf, kp, ktri, kedge, k
    real(dp), dimension(3) :: n_f, grad_j0B0, B0_grad_B0
    complex(dp) :: zdum, B0_jnpar
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
    call solve_MDE(inhom, jnpar_B0%DOF)
    call add_shielding_current
    do kedge = 1, mesh%nedge
       do k = 1, mesh%GL_order
          ktri = mesh%edge_tri(1, kedge)
          n_f = [mesh%edge_Z(kedge), 0d0, -mesh%edge_R(kedge)]
          associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
            call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
            call RT0_interp(Bn, ktri, R, Z, B_n)
            call L1_interp(jnpar_B0, ktri, R, Z, B0_jnpar)
            B0_jnpar = B0_jnpar * f%Bmod ** 2
            jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * R * &
                 (B0_jnpar * sum(f%B0 * n_f) - clight * sum(zd_cross(grad_pn, f%B0) * n_f) + &
                 sum(f%j0 * f%B0) * sum(B_n * n_f) - sum(B_n * f%B0) * sum(f%j0 * n_f)) / f%Bmod ** 2
          end associate
       end do
    end do
    n_f = [0d0, 1d0, 0d0]
    do ktri = 1, mesh%ntri
       do k = 1, mesh%GL2_order
          associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
            call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
            call RT0_interp(Bn, ktri, R, Z, B_n)
            call L1_interp(jnpar_B0, ktri, R, Z, B0_jnpar)
            B0_jnpar = B0_jnpar * f%Bmod ** 2
            jn%comp_phi(ktri) = jn%comp_phi(ktri) + mesh%GL2_weights(k) * &
                 (B0_jnpar * sum(f%B0 * n_f) - clight * sum(zd_cross(grad_pn, f%B0) * n_f) + &
                 sum(f%j0 * f%B0) * sum(B_n * n_f) - sum(B_n * f%B0) * sum(f%j0 * n_f)) / f%Bmod ** 2
          end associate
       end do
    end do

    if (first_call) then
       first_call = .false.
       call debug_currn("debug_currn_000", pn, Bn, jn, inhom)
       call debug_MDE("debug_MDE_000", pn, Bn, jn, jnpar_B0)
    end if
  end subroutine compute_currn_MDE

  subroutine add_shielding_current
    use mephit_conf, only: conf_arr
    use mephit_util, only: imun
    use mephit_mesh, only: equil, mesh
    integer :: m, m_res, kf, kpoi_min, kpoi_max
    complex(dp) :: pn_outer, pn_inner

    do m = mesh%m_res_min, mesh%m_res_max
       m_res = -equil%cocos%sgn_q * m
       if (abs(conf_arr%sheet_current_factor(m)) > 0d0) then
          kf = mesh%res_ind(m)
          kpoi_min = mesh%kp_low(kf) + 1
          kpoi_max = mesh%kp_low(kf) + mesh%kp_max(kf)
          pn_outer = sum(exp(-imun * m_res * mesh%node_theta_flux(kpoi_min:kpoi_max)) * &
               pn%DOF(kpoi_min:kpoi_max)) / dble(mesh%kp_max(kf))
          kf = mesh%res_ind(m) - 1
          kpoi_min = mesh%kp_low(kf) + 1
          kpoi_max = mesh%kp_low(kf) + mesh%kp_max(kf)
          pn_inner = sum(exp(-imun * m_res * mesh%node_theta_flux(kpoi_min:kpoi_max)) * &
               pn%DOF(kpoi_min:kpoi_max)) / dble(mesh%kp_max(kf))
          kf = mesh%res_ind(m)
          kpoi_min = mesh%kp_low(kf) + 1
          kpoi_max = mesh%kp_low(kf) + mesh%kp_max(kf)
          jnpar_B0%DOF(kpoi_min:kpoi_max) = jnpar_B0%DOF(kpoi_min:kpoi_max) + &
               mesh%shielding_coeff(m) * conf_arr%sheet_current_factor(m) * &
               (pn_outer - pn_inner) * exp(imun * m_res * mesh%node_theta_flux(kpoi_min:kpoi_max))
       end if
    end do
  end subroutine add_shielding_current

  subroutine debug_currn(group, presn, magfn, currn, inhom)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: datafile
    use mephit_util, only: clight, zd_cross
    use mephit_mesh, only: mesh, cache
    use mephit_pert, only: RT0_interp, L1_interp
    character(len = *), intent(in) :: group
    type(L1_t), intent(in) :: presn
    type(RT0_t), intent(in) :: magfn, currn
    complex(dp), intent(in), optional :: inhom(:)
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer :: kf, kt, ktri
    complex(dp) :: zdum, Bn(3), jn(3)
    complex(dp), dimension(3, mesh%ntri) :: grad_pn, lorentz

    grp = trim(adjustl(group))
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          associate (R => mesh%cntr_R(ktri), Z => mesh%cntr_Z(ktri), c => cache%cntr_fields(ktri))
            call L1_interp(presn, ktri, R, Z, zdum, grad_pn(:, ktri))
            call RT0_interp(magfn, ktri, R, Z, Bn)
            call RT0_interp(currn, ktri, R, Z, jn)
            lorentz(:, ktri) = (zd_cross(jn, c%B0) - zd_cross(Bn, c%j0)) / clight
          end associate
       end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/grad_pn', grad_pn, lbound(grad_pn), ubound(grad_pn), &
         comment = 'perturbation pressure gradient', unit = 'dyn cm^-2')
    call h5_add(h5id_root, grp // '/lorentz', lorentz, lbound(lorentz), ubound(lorentz), &
         comment = 'perturbation Lorentz force density', unit = 'dyn cm^-2')
    call h5_add(h5id_root, grp // '/I', currn%DOF, lbound(currn%DOF), ubound(currn%DOF), &
         comment = 'perturbation current density degrees of freedom', unit = 'statA')
    if (present(inhom)) then
       call h5_add(h5id_root, grp // '/inhom', inhom, lbound(inhom), ubound(inhom), &
            comment = 'inhomogeneities of linear systems of equations', unit = 'statA cm^-3')
    end if
    call h5_close(h5id_root)
  end subroutine debug_currn

  subroutine add_sheet_current
    use mephit_conf, only: conf_arr
    use mephit_mesh, only: mesh, cache
    integer :: m, kf, kt, kedge, ke, k
    complex(dp) :: pn_outer, pn_inner

    do m = mesh%m_res_min, mesh%m_res_max
       kf = mesh%res_ind(m)
       if (abs(conf_arr%sheet_current_factor(m)) > 0d0) then
          do kt = 1, mesh%kt_max(kf)
             kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
             ke = mesh%shielding_kt_low(m) + kt
             do k = 1, mesh%GL_order
                pn_outer = sum(pn%DOF(mesh%kp_low(kf) + [mesh%shielding_L1_kp(0, k, ke), &
                     mod(mesh%shielding_L1_kp(0, k, ke), mesh%kp_max(kf)) + 1]) * &
                     [mesh%shielding_L1_weight(0, k, ke), &
                     1d0 - mesh%shielding_L1_weight(0, k, ke)])
                pn_inner = sum(pn%DOF(mesh%kp_low(kf-1) + [mesh%shielding_L1_kp(-1, k, ke), &
                     mod(mesh%shielding_L1_kp(-1, k, ke), mesh%kp_max(kf-1)) + 1]) * &
                     [mesh%shielding_L1_weight(-1, k, ke), &
                     1d0 - mesh%shielding_L1_weight(-1, k, ke)])
                jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * &
                     mesh%shielding_coeff(m) * conf_arr%sheet_current_factor(m) * &
                     (pn_outer - pn_inner) * mesh%GL_R(k, kedge) * &
                     (cache%edge_fields(k, kedge)%B0(1) * mesh%edge_Z(kedge) - &
                     cache%edge_fields(k, kedge)%B0(3) * mesh%edge_R(kedge))
             end do
          end do
       end if
    end do
  end subroutine add_sheet_current

  subroutine mephit_postprocess
    use magdata_in_symfluxcoor_mod, only: ntheta
    use mephit_conf, only: conf, datafile
    use mephit_mesh, only: mesh, coord_cache_ext_t, compute_sample_Ipar
    use mephit_pert, only: vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, vec_polmodes_write, &
         polmodes_t, polmodes_init, polmodes_deinit, polmodes_write, &
         L1_poloidal_modes, RT0_poloidal_modes, vac
    integer, parameter :: m_max = 24
    integer :: k
    character(len = 24) :: dataset
    type(polmodes_t) :: polmodes
    type(vec_polmodes_t) :: vec_polmodes
    type(coord_cache_ext_t) :: sample_Ipar(ntheta, conf%nrad_Ipar)

    ! poloidal modes
    call polmodes_init(polmodes, m_max, mesh%nflux)
    call L1_poloidal_modes(pn, polmodes)
    call polmodes_write(polmodes, datafile, 'postprocess/pmn', &
         'poloidal modes of pressure perturbation', 'dyn cm^-2')
    call polmodes_deinit(polmodes)
    call vec_polmodes_init(vec_polmodes, m_max, mesh%nflux)
    call RT0_poloidal_modes(Bn, vec_polmodes)
    call vec_polmodes_write(vec_polmodes, datafile, 'postprocess/Bmn', &
         'poloidal modes of magnetic field (full perturbation)', 'G')
    call RT0_poloidal_modes(vac%Bn, vec_polmodes)
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
    do k = lbound(mesh%res_modes, 1), ubound(mesh%res_modes, 1)
       call compute_sample_Ipar(sample_Ipar, mesh%res_modes(k))
       write (dataset, '("postprocess/Imn_par_", i0)') mesh%res_modes(k)
       call write_Ipar_symfluxcoord(mesh%res_modes(k), sample_Ipar, datafile, dataset)
       if (conf%kilca_scale_factor /= 0) then
          call write_Ipar(mesh%res_modes(k), sample_Ipar, datafile, 'postprocess/Imn_par_KiLCA')
       end if
    end do
  end subroutine mephit_postprocess

  subroutine check_furth(currn, Bmn_plas)
    use mephit_conf, only: conf, datafile
    use mephit_util, only: imun, clight
    use mephit_mesh, only: equil, fs_half, mesh, cache
    use mephit_pert, only: RT0_t, vec_polmodes_t
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(RT0_t), intent(in) :: currn
    type(vec_polmodes_t), intent(in) :: Bmn_plas
    character(len = *), parameter :: dataset = 'debug_furth'
    integer(HID_T) :: h5id_root
    integer :: kf, kpol, kilca_m_res
    complex(dp) :: sheet_flux(mesh%nflux)
    real(dp) :: k_z, k_theta(mesh%nflux)

    kilca_m_res = -equil%cocos%sgn_q * abs(conf%kilca_pol_mode)
    k_z = mesh%n / mesh%R_O
    k_theta(:) = kilca_m_res / fs_half%rad
    sheet_flux(:) = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kpol = 1, cache%npol
          associate (s => cache%sample_polmodes_half(kpol, kf))
            sheet_flux(kf) = sheet_flux(kf) + mesh%area(s%ktri) * currn%comp_phi(s%ktri) * &
                 exp(-imun * kilca_m_res * s%theta)
          end associate
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
  subroutine write_Ipar_symfluxcoord(m, sample_Ipar, file, dataset)
    use mephit_conf, only: conf
    use mephit_util, only: imun, pi, clight
    use mephit_mesh, only: coord_cache_ext_t
    use mephit_pert, only: RT0_interp
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    integer, intent(in) :: m
    type(coord_cache_ext_t), dimension(:, :), intent(in) :: sample_Ipar
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root
    integer :: nrad, npol, krad, kpol
    real(dp) :: B0_psi, dB0R_dpsi, dB0phi_dpsi, dB0Z_dpsi, &
         B0_dB0_dpsi, dhphi2_dpsi, B0_theta, dB0R_dtheta, dB0phi_dtheta, dB0Z_dtheta, &
         B0_dB0_dtheta, dhphi2_dtheta, common_term, dB0theta_dpsi, dB0psi_dtheta, &
         dhphihtheta_dpsi, dhphihpsi_dtheta
    complex(dp) :: vec(3), dvec_dR(3), dvec_dZ(3), &
         jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
         dBnR_dR, dBnR_dZ, dBnZ_dR, dBnZ_dZ, &
         Bn_psi, dBnpsi_dpsi, Bn_theta, Delta_mn, part_int, bndry
    real(dp), dimension(size(sample_Ipar, 2)) :: rad, psi, I_char
    complex(dp), dimension(size(sample_Ipar, 2)) :: jmn_par_neg, jmn_par_pos, &
         part_int_neg, part_int_pos, bndry_neg, bndry_pos, Delta_mn_neg, Delta_mn_pos

    npol = size(sample_Ipar, 1)
    nrad = size(sample_Ipar, 2)
    psi(:) = sample_Ipar(1, :)%psi
    rad(:) = sample_Ipar(1, :)%rad
    jmn_par_neg = (0d0, 0d0)
    jmn_par_pos = (0d0, 0d0)
    part_int_neg = (0d0, 0d0)
    part_int_pos = (0d0, 0d0)
    bndry_neg = (0d0, 0d0)
    bndry_pos = (0d0, 0d0)
    I_char = 0d0
    Delta_mn_neg = (0d0, 0d0)
    Delta_mn_pos = (0d0, 0d0)
    do krad = 1, nrad
       do kpol = 1, npol
          associate (s => sample_Ipar(kpol, krad))
            call RT0_interp(jn, s%ktri, s%R, s%Z, vec)
            jn_R = vec(1)
            jn_phi = vec(2)
            jn_Z = vec(3)
            ! include h^phi in current density
            jn_par = (jn_R * s%B0_R + jn_Z * s%B0_Z + jn_phi * s%B0_phi) * &
                 s%B0_phi / s%B0_2 * s%sqrt_g / s%R
            jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * exp(imun * abs(m) * s%theta)
            jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * exp(-imun * abs(m) * s%theta)
            ! comparison with indirect calculation (Boozer and Nuehrenberg 2006)
            call RT0_interp(Bn, s%ktri, s%R, s%Z, vec, dvec_dR = dvec_dR, dvec_dZ = dvec_dZ)
            Bn_R = vec(1)
            Bn_phi = vec(2)
            Bn_Z = vec(3)
            dBnR_dR = dvec_dR(1)
            dBnZ_dR = dvec_dR(3)
            dBnR_dZ = dvec_dZ(1)
            dBnZ_dZ = dvec_dZ(3)
            dB0phi_dpsi = (s%dB0phi_dR * s%dR_dpsi + s%dB0phi_dZ * s%dZ_dpsi) / s%R - &
                 s%B0_phi * s%dR_dpsi / s%R ** 2
            Bn_psi = s%R * (Bn_R * s%B0_Z - Bn_Z * s%B0_R)
            dBnpsi_dpsi = s%dR_dpsi * (Bn_R * s%B0_Z - Bn_Z * s%B0_R) + s%R * ( &
                 Bn_R * (s%dB0Z_dR * s%dR_dpsi + s%dB0Z_dZ * s%dZ_dpsi) - &
                 Bn_Z * (s%dB0R_dR * s%dR_dpsi + s%dB0R_dZ * s%dZ_dpsi) + &
                 (dBnR_dR * s%dR_dpsi + dBnR_dZ * s%dZ_dpsi) * s%B0_Z - &
                 (dBnZ_dR * s%dR_dpsi + dBnZ_dZ * s%dZ_dpsi) * s%B0_R)
            Delta_mn = s%dq_dpsi / s%q * ((Bn_psi + dBnpsi_dpsi) / s%B0_phi * s%R - &
                 Bn_psi * dB0phi_dpsi / s%B0_phi ** 2 * s%R ** 2)
            Delta_mn_neg(krad) = Delta_mn_neg(krad) + Delta_mn * exp(imun * abs(m) * s%theta)
            Delta_mn_pos(krad) = Delta_mn_pos(krad) + Delta_mn * exp(-imun * abs(m) * s%theta)
            I_char(krad) = I_char(krad) + s%B0_2 / &
                 ((s%B0_R ** 2 + s%B0_Z ** 2) * s%q * s%q * s%R * s%B0_phi)
            ! comparison with indirect calculation (Ampere's law and integration by parts)
            B0_psi = s%B0_R * s%dR_dpsi + s%B0_Z * s%dZ_dpsi  ! covariant component
            dB0R_dpsi = s%dB0R_dR * s%dR_dpsi + s%dB0R_dZ * s%dZ_dpsi
            dB0phi_dpsi = s%dB0phi_dR * s%dR_dpsi + s%dB0phi_dZ * s%dZ_dpsi
            dB0Z_dpsi = s%dB0Z_dR * s%dR_dpsi + s%dB0Z_dZ * s%dZ_dpsi
            B0_dB0_dpsi = s%B0_R * dB0R_dpsi + s%B0_phi * dB0phi_dpsi + s%B0_Z * dB0Z_dpsi
            dhphi2_dpsi = 2d0 * s%B0_phi * (dB0phi_dpsi / s%B0_2 - &
                 s%B0_phi * B0_dB0_dpsi / s%B0_2 ** 2)
            B0_theta = s%B0_R * s%dR_dtheta + s%B0_Z * s%dZ_dtheta  ! covariant component
            dB0R_dtheta = s%dB0R_dR * s%dR_dtheta + s%dB0R_dZ * s%dZ_dtheta
            dB0phi_dtheta = s%dB0phi_dR * s%dR_dtheta + s%dB0phi_dZ * s%dZ_dtheta
            dB0Z_dtheta = s%dB0Z_dR * s%dR_dtheta + s%dB0Z_dZ * s%dZ_dtheta
            B0_dB0_dtheta = s%B0_R * dB0R_dtheta + s%B0_phi * dB0phi_dtheta + s%B0_Z * dB0Z_dtheta
            dhphi2_dtheta = 2d0 * s%B0_phi * (dB0phi_dtheta / s%B0_2 - &
                 s%B0_phi * B0_dB0_dtheta / s%B0_2 ** 2)
            common_term = s%B0_R * s%d2R_dpsi_dtheta + s%B0_Z * s%d2Z_dpsi_dtheta + &
                 s%dB0R_dR * s%dR_dpsi * s%dR_dtheta + s%dB0Z_dZ * s%dZ_dpsi * s%dZ_dtheta
            dB0theta_dpsi = common_term + &
                 s%dB0R_dZ * s%dZ_dpsi * s%dR_dtheta + s%dB0Z_dR * s%dR_dpsi * s%dZ_dtheta
            dB0psi_dtheta = common_term + &
                 s%dB0R_dZ * s%dZ_dtheta * s%dR_dpsi + s%dB0Z_dR * s%dR_dtheta * s%dZ_dpsi
            dhphihtheta_dpsi = ((dB0phi_dpsi / s%R - s%B0_phi * s%dR_dpsi / s%R ** 2) * B0_theta + &
                 s%B0_phi / s%R * dB0theta_dpsi) / s%B0_2 - 2d0 * s%B0_phi / s%R * B0_theta * &
                 B0_dB0_dpsi / s%B0_2 ** 2
            dhphihpsi_dtheta = ((dB0phi_dtheta / s%R - s%B0_phi * s%dR_dtheta / s%R ** 2) * B0_psi + &
                 s%B0_phi / s%R * dB0psi_dtheta) / s%B0_2 - 2d0 * s%B0_phi / s%R * B0_psi * &
                 B0_dB0_dtheta / s%B0_2 ** 2
            Bn_psi = Bn_R * s%dR_dpsi + Bn_Z * s%dZ_dpsi
            Bn_theta = Bn_R * s%dR_dtheta + Bn_Z * s%dZ_dtheta
            part_int = dhphi2_dpsi * Bn_theta - dhphi2_dtheta * Bn_psi + Bn_phi / s%R * &
                 (dhphihtheta_dpsi - dhphihpsi_dtheta) + imun * conf%n * s%B0_phi / s%R * &
                 (B0_psi * Bn_theta - B0_theta * Bn_psi) / s%B0_2
            part_int_neg(krad) = part_int_neg(krad) + (part_int - imun * abs(m) * &
                 s%B0_phi * (s%B0_phi * Bn_psi - B0_psi * Bn_phi)) * &
                 exp(imun * abs(m) * s%theta)
            part_int_pos(krad) = part_int_pos(krad) + (part_int + imun * abs(m) * &
                 s%B0_phi * (s%B0_phi * Bn_psi - B0_psi * Bn_phi)) * &
                 exp(-imun * abs(m) * s%theta)
            bndry = s%B0_phi * (Bn_phi * B0_theta - Bn_theta * s%B0_phi) / s%B0_2
            bndry_neg(krad) = bndry_neg(krad) + bndry * exp(imun * abs(m) * s%theta)
            bndry_pos(krad) = bndry_pos(krad) + bndry * exp(-imun * abs(m) * s%theta)
          end associate
       end do
       jmn_par_neg(krad) = jmn_par_neg(krad) / dble(npol)
       jmn_par_pos(krad) = jmn_par_pos(krad) / dble(npol)
       part_int_neg(krad) = part_int_neg(krad) / dble(npol) * 0.25d0 * clight / pi
       part_int_pos(krad) = part_int_pos(krad) / dble(npol) * 0.25d0 * clight / pi
       bndry_neg(krad) = bndry_neg(krad) / dble(npol) * 0.25d0 * clight / pi
       bndry_pos(krad) = bndry_pos(krad) / dble(npol) * 0.25d0 * clight / pi
       Delta_mn_neg(krad) = Delta_mn_neg(krad) / dble(npol)
       Delta_mn_pos(krad) = Delta_mn_pos(krad) / dble(npol)
       I_char(krad) = dble(npol) * 0.5d0 * clight / I_char(krad)
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
  subroutine write_Ipar(m, sample_Ipar, file, dataset)
    use mephit_util, only: imun, pi, clight, bent_cyl2straight_cyl
    use mephit_mesh, only: mesh, coord_cache_ext_t
    use mephit_pert, only: RT0_interp
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    integer, intent(in) :: m
    type(coord_cache_ext_t), dimension(:, :), intent(in) :: sample_Ipar
    character(len = *), intent(in) :: file, dataset
    integer(HID_T) :: h5id_root
    integer :: nrad, npol, krad, kpol
    real(dp) :: dB0R_drad, dB0phi_drad, dB0Z_drad, B0_dB0_drad, &
         B0_theta, dB0theta_drad, dhz2_drad, dradhthetahz_drad
    complex(dp) :: vec(3), jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
         Bn_rad, Bn_pol, Bn_tor, part_int, bndry
    real(dp), dimension(size(sample_Ipar, 2)) :: rad, psi
    complex(dp), dimension(size(sample_Ipar, 2)) :: jmn_par_neg, jmn_par_pos, &
         part_int_neg, part_int_pos, bndry_neg, bndry_pos

    npol = size(sample_Ipar, 1)
    nrad = size(sample_Ipar, 2)
    psi(:) = sample_Ipar(1, :)%psi
    rad(:) = sample_Ipar(1, :)%rad
    jmn_par_neg = (0d0, 0d0)
    jmn_par_pos = (0d0, 0d0)
    part_int_neg = (0d0, 0d0)
    part_int_pos = (0d0, 0d0)
    bndry_neg = (0d0, 0d0)
    bndry_pos = (0d0, 0d0)
    do krad = 1, nrad
       do kpol = 1, npol
          associate (s => sample_Ipar(kpol, krad))
            call RT0_interp(jn, s%ktri, s%R, s%Z, vec)
            jn_R = vec(1)
            jn_phi = vec(2)
            jn_Z = vec(3)
            ! include h^z in current density
            jn_par = (jn_R * s%B0_R + jn_Z * s%B0_Z + jn_phi * s%B0_phi) * s%B0_phi / s%B0_2
            jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * exp(imun * abs(m) * s%theta)
            jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * exp(-imun * abs(m) * s%theta)
            ! comparison with indirect calculation (Ampere's law and integration by parts)
            dB0R_drad = s%dB0R_dR * cos(s%theta) + s%dB0R_dZ * sin(s%theta)
            dB0phi_drad = s%dB0phi_dR * cos(s%theta) + s%dB0phi_dZ * sin(s%theta)
            dB0Z_drad = s%dB0Z_dR * cos(s%theta) + s%dB0Z_dZ * sin(s%theta)
            B0_dB0_drad = s%B0_R * dB0R_drad + s%B0_phi * dB0phi_drad + s%B0_Z * dB0Z_drad
            B0_theta = s%B0_Z * cos(s%theta) - s%B0_R * sin(s%theta)
            dB0theta_drad = s%dB0Z_dR * cos(s%theta) ** 2 - s%dB0R_dZ * sin(s%theta) ** 2 + &
                 (s%dB0Z_dZ - s%dB0R_dR) * sin(s%theta) * cos(s%theta)
            dhz2_drad = 2d0 * s%B0_phi * (s%B0_2 * dB0phi_drad - s%B0_phi * B0_dB0_drad) / &
                 s%B0_2 ** 2
            dradhthetahz_drad = B0_theta * s%B0_phi / s%B0_2 + rad(krad) * &
                 ((B0_theta * dB0phi_drad + dB0theta_drad * s%B0_phi) / s%B0_2 - &
                 2d0 * B0_theta * s%B0_phi * B0_dB0_drad / s%B0_2 ** 2)
            call RT0_interp(Bn, s%ktri, s%R, s%Z, vec)
            Bn_R = vec(1)
            Bn_phi = vec(2)
            Bn_Z = vec(3)
            call bent_cyl2straight_cyl(Bn_R, Bn_phi, Bn_Z, s%theta, &
                 Bn_rad, Bn_pol, Bn_tor)
            part_int = -rad(krad) * Bn_pol * dhz2_drad + Bn_tor * dradhthetahz_drad
            part_int_neg(krad) = part_int_neg(krad) + (part_int + imun * s%B0_phi * &
                 (dble(mesh%n) / mesh%R_O * rad(krad) * B0_theta + &
                 abs(m) * s%B0_phi) * Bn_rad / s%B0_2) * exp(imun * abs(m) * s%theta)
            part_int_pos(krad) = part_int_pos(krad) + (part_int + imun * s%B0_phi * &
                 (dble(mesh%n) / mesh%R_O * rad(krad) * B0_theta - &
                 abs(m) * s%B0_phi) * Bn_rad / s%B0_2) * exp(-imun * abs(m) * s%theta)
            bndry = s%B0_phi * rad(krad) * (s%B0_phi * Bn_pol - B0_theta * Bn_tor) / s%B0_2
            bndry_neg(krad) = bndry_neg(krad) + bndry * exp(imun * abs(m) * s%theta)
            bndry_pos(krad) = bndry_pos(krad) + bndry * exp(-imun * abs(m) * s%theta)
          end associate
       end do
       jmn_par_neg(krad) = jmn_par_neg(krad) / dble(npol)
       jmn_par_pos(krad) = jmn_par_pos(krad) / dble(npol)
       part_int_neg(krad) = part_int_neg(krad) / dble(npol) * 0.25d0 * clight / pi
       part_int_pos(krad) = part_int_pos(krad) / dble(npol) * 0.25d0 * clight / pi
       bndry_neg(krad) = bndry_neg(krad) / dble(npol) * 0.25d0 * clight / pi
       bndry_pos(krad) = bndry_pos(krad) / dble(npol) * 0.25d0 * clight / pi
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
end module mephit_iter
