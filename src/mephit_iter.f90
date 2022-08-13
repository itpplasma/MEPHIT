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
    use mephit_mesh, only: equil, mesh, cache, fft, psi_interpolator, &
         generate_mesh, mesh_write, mesh_read, write_cache, read_cache
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
       call fft%init(cache%npol)
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
    call fft%deinit
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
    call L1_read(jnpar_B0, datafile, 'iter/jnpar_Bmod')
  end subroutine iter_read

  subroutine iter_step
    ! compute pressure based on previous perturbation field
    call compute_presn
    ! compute currents based on previous perturbation field
    call compute_currn
    ! use field code to generate new field from currents
    call compute_Bn
  end subroutine iter_step

  subroutine mephit_iterate
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_util, only: arnoldi_break
    use mephit_mesh, only: mesh, cache
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
    type(polmodes_t) :: pmn, jmnpar_Bmod
    type(vec_polmodes_t) :: jmn, Bmn
    complex(dp) :: Ires(mesh%m_res_min:mesh%m_res_max)

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
    call polmodes_init(jmnpar_Bmod, m_max, mesh%nflux)
    call vec_polmodes_init(jmn, m_max, mesh%nflux)
    call vec_polmodes_init(Bmn, m_max, mesh%nflux)
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
          call L1_write(jnpar_B0, datafile, 'iter/jnpar_Bmod' // postfix, &
               'parallel current density (after iteration)', 's^-1')  ! SI: H^-1
          call RT0_write(jn, datafile, 'iter/jn' // postfix, &
               'current density (after iteration)', 'statA cm^-2', 1)
          call RT0_write(Bn, datafile, 'iter/Bn' // postfix, &
               'magnetic field (after iteration)', 'G', 1)
          call RT0_write(Bn_diff, datafile, 'iter/Bn_diff' // postfix, &
               'magnetic field (difference between iterations)', 'G', 1)
       else
          if (L2int_Bn_diff(kiter) > conf%ritz_threshold ** kiter * L2int_Bnvac .or. &
               L2int_Bn_diff(kiter) < conf%iter_rel_err * L2int_Bnvac) then
             niter = kiter
             exit
          end if
       end if
       call L1_poloidal_modes(pn, pmn)
       call polmodes_write(pmn, datafile, 'iter/pmn' // postfix, &
            'pressure (after iteration)', 'dyn cm^-2')
       call L1_poloidal_modes(jnpar_B0, jmnpar_Bmod)
       call polmodes_write(jmnpar_Bmod, datafile, 'iter/jmnpar_Bmod' // postfix, &
            'parallel current density (after iteration)', 's^-1')  ! SI: H^-1
       call compute_Ires(cache%sample_Ires, cache%GL_weights, jnpar_B0, Ires)
       call h5_open_rw(datafile, h5id_root)
       call h5_add(h5id_root, 'iter/Ires' // postfix, Ires, lbound(Ires), ubound(Ires), &
            comment = 'resonant currents (after iteration)', unit = 'statA')
       call h5_close(h5id_root)
       call RT0_poloidal_modes(jn, jmn)
       call vec_polmodes_write(jmn, datafile, 'iter/jmn' // postfix, &
            'current density (after iteration)', 'statA cm^-2')
       call RT0_poloidal_modes(Bn, Bmn)
       call vec_polmodes_write(Bmn, datafile, 'iter/Bmn' // postfix, &
            'magnetic field (after iteration)', 'G')
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
    call L1_write(jnpar_B0, datafile, 'iter/jnpar_Bmod', &
         'parallel current density (full perturbation)', 's^-1')  ! SI: H^-1
    call RT0_write(Bn, datafile, 'iter/Bn', &
         'magnetic field (full perturbation)', 'G', 2)
    call RT0_write(Bnplas, datafile, 'iter/Bnplas', &
         'magnetic field (plasma response)', 'G', 2)
    call RT0_write(jn, datafile, 'iter/jn', &
         'current density (full perturbation)', 'statA cm^-2', 1)
    call RT0_deinit(Bn_prev)
    call RT0_deinit(Bn_diff)
    call polmodes_deinit(pmn)
    call polmodes_deinit(jmnpar_Bmod)
    call vec_polmodes_deinit(jmn)
    call vec_polmodes_deinit(Bmn)

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

  subroutine compute_currn
    use mephit_mesh, only: mesh, cache, field_cache_t
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
    if (first_call) then
       do kedge = 1, mesh%nedge
          do k = 1, mesh%GL_order
             ktri = mesh%edge_tri(1, kedge)
             n_f = [mesh%edge_Z(kedge), 0d0, -mesh%edge_R(kedge)]
             call project_combined(jn%DOF(kedge), ktri, cache%edge_fields(k, kedge), &
                  mesh%GL_weights(k) * mesh%GL_R(k, kedge), mesh%GL_R(k, kedge), mesh%GL_Z(k, kedge), n_f)
          end do
       end do
       n_f = [0d0, 1d0, 0d0]
       do ktri = 1, mesh%ntri
          do k = 1, mesh%GL2_order
             call project_combined(jn%comp_phi(ktri), ktri, cache%area_fields(k, ktri), &
                  mesh%GL2_weights(k), mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), n_f)
          end do
       end do
       call debug_currn("debug_currn_initial", pn, Bn, jn, inhom)
       call debug_MDE("debug_MDE_initial", pn, Bn, jn, jnpar_B0)
       jn%DOF(:) = (0d0, 0d0)
       jn%comp_phi(:) = (0d0, 0d0)
    end if
    call add_shielding_current
    do kedge = 1, mesh%nedge
       do k = 1, mesh%GL_order
          ktri = mesh%edge_tri(1, kedge)
          n_f = [mesh%edge_Z(kedge), 0d0, -mesh%edge_R(kedge)]
          call project_combined(jn%DOF(kedge), ktri, cache%edge_fields(k, kedge), &
               mesh%GL_weights(k) * mesh%GL_R(k, kedge), mesh%GL_R(k, kedge), mesh%GL_Z(k, kedge), n_f)
       end do
    end do
    n_f = [0d0, 1d0, 0d0]
    do ktri = 1, mesh%ntri
       do k = 1, mesh%GL2_order
          call project_combined(jn%comp_phi(ktri), ktri, cache%area_fields(k, ktri), &
               mesh%GL2_weights(k), mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), n_f)
       end do
    end do

    if (first_call) then
       first_call = .false.
    else
       ! hack: overwrite to save only last iteration step
       call debug_currn("debug_currn_final", pn, Bn, jn, inhom)
       call debug_MDE("debug_MDE_final", pn, Bn, jn, jnpar_B0)
    end if

  contains
    subroutine project_combined(comp, ktri, f, weight, R, Z, n_f)
      complex(dp), intent(inout) :: comp
      integer, intent(in) :: ktri
      type(field_cache_t), intent(in) :: f
      real(dp), intent(in) :: weight, R, Z, n_f(3)

      call L1_interp(pn, ktri, R, Z, zdum, grad_pn)
      call RT0_interp(Bn, ktri, R, Z, B_n)
      call L1_interp(jnpar_B0, ktri, R, Z, B0_jnpar)
      B0_jnpar = B0_jnpar * f%Bmod ** 2
      comp = comp + weight * &
           (B0_jnpar * sum(f%B0 * n_f) - clight * sum(zd_cross(grad_pn, f%B0) * n_f) + &
           sum(f%j0 * f%B0) * sum(B_n * n_f) - sum(B_n * f%B0) * sum(f%j0 * n_f)) / f%Bmod ** 2
    end subroutine project_combined
  end subroutine compute_currn

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

  subroutine mephit_postprocess
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_mesh, only: mesh, cache
    use mephit_pert, only: vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, vec_polmodes_write, &
         polmodes_t, polmodes_init, polmodes_deinit, polmodes_write, &
         L1_poloidal_modes, RT0_poloidal_modes, vac
    integer, parameter :: m_max = 24
    integer(HID_T) :: h5id_root
    type(polmodes_t) :: polmodes
    type(vec_polmodes_t) :: vec_polmodes
    complex(dp), allocatable :: Ires(:)

    ! poloidal modes
    call polmodes_init(polmodes, m_max, mesh%nflux)
    call L1_poloidal_modes(pn, polmodes)
    call polmodes_write(polmodes, datafile, 'postprocess/pmn', &
         'poloidal modes of pressure perturbation', 'dyn cm^-2')
    call L1_poloidal_modes(jnpar_B0, polmodes)
    call polmodes_write(polmodes, datafile, 'iter/jmnpar_Bmod', &
         'parallel current density', 's^-1')  ! SI: H^-1
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
    ! resonant currents
    allocate(Ires(mesh%m_res_min:mesh%m_res_max))
    call compute_Ires(cache%sample_Ires, cache%GL_weights, jnpar_B0, Ires)
    call h5_open_rw(datafile, h5id_root)
    call h5_add(h5id_root, 'postprocess/Ires', Ires, lbound(Ires), ubound(Ires), &
         comment = 'resonant currents', unit = 'statA')
    call h5_close(h5id_root)
  end subroutine mephit_postprocess

  subroutine check_furth(currn, Bmn_plas)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_util, only: imun, clight
    use mephit_mesh, only: equil, fs_half, mesh, cache
    use mephit_pert, only: RT0_t, vec_polmodes_t
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

  !> compute resonant currents
  subroutine compute_Ires(sample_Ires, GL_weights, jnpar_Bmod, Ires)
    use mephit_conf, only: logger
    use mephit_util, only: imun
    use mephit_mesh, only: coord_cache_t, mesh, equil
    use mephit_pert, only: L1_t, L1_interp
    type(coord_cache_t), intent(in) :: sample_Ires(:, :, mesh%m_res_min:)
    real(dp), intent(in) :: GL_weights(:, mesh%m_res_min:)
    type(L1_t), intent(in) :: jnpar_Bmod
    complex(dp), intent(out) :: Ires(mesh%m_res_min:)
    integer :: nrad, npol, krad, kpol, m
    complex(dp) :: jn_par, jmn_par

    if (size(GL_weights, 1) /= size(sample_Ires, 2)) then
       call logger%msg_arg_size('compute_Ires', &
            'size(GL_weights, 1)', 'size(sample_Ires, 2)', &
            size(GL_weights, 1), size(sample_Ires, 2))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (ubound(GL_weights, 2) /= mesh%m_res_max) then
       call logger%msg_arg_size('compute_Ires', &
            'ubound(GL_weights, 2)', 'mesh%m_res_max', &
            ubound(GL_weights, 2), mesh%m_res_max)
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (ubound(Ires, 1) /= mesh%m_res_max) then
       call logger%msg_arg_size('compute_Ires', &
            'ubound(Ires, 1)', 'mesh%m_res_max', &
            ubound(Ires, 1), mesh%m_res_max)
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (ubound(sample_Ires, 3) /= mesh%m_res_max) then
       call logger%msg_arg_size('compute_Ires', &
            'ubound(sample_Ires, 3)', 'mesh%m_res_max', &
            ubound(sample_Ires, 3), mesh%m_res_max)
       if (logger%err) call logger%write_msg
       error stop
    end if
    npol = size(sample_Ires, 1)
    nrad = size(sample_Ires, 2)
    do m = mesh%m_res_min, mesh%m_res_max
       Ires(m) = (0d0, 0d0)
       do krad = 1, nrad
          jmn_par = (0d0, 0d0)
          do kpol = 1, npol
             associate (s => sample_Ires(kpol, krad, m))
               call L1_interp(jnpar_Bmod, s%ktri, s%R, s%Z, jn_par)
               jn_par = jn_par * sqrt(s%B0_R ** 2 + s%B0_phi ** 2 + s%B0_Z ** 2)
               ! jmn_par includes the area differential
               jmn_par = jmn_par + s%sqrt_g / s%R * jn_par * &
                    exp(imun * equil%cocos%sgn_q * m * s%theta)
             end associate
          end do
          Ires(m) = Ires(m) + jmn_par / dble(npol) * GL_weights(krad, m)
       end do
    end do
  end subroutine compute_Ires

end module mephit_iter
