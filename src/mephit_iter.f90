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
  end interface

contains

  subroutine mephit_run(runmode, config) bind(C, name = 'mephit_run')
    use iso_c_binding, only: c_int, c_ptr
    use magdata_in_symfluxcoor_mod, only: load_magdata_in_symfluxcoord
    use mephit_util, only: C_F_string, get_field_filenames, init_field
    use mephit_conf, only: conf, config_read, config_export_hdf5, conf_arr, logger, datafile
    use mephit_mesh, only: equil, mesh, generate_mesh, mesh_write, mesh_read, write_cache, read_cache
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
         mesh, cache, mesh_deinit, cache_deinit, &
         B0R_edge, B0phi_edge, B0Z_edge, B0R_Omega, B0phi_Omega, B0Z_Omega, B0_flux, j0phi_edge
    use mephit_pert, only: vac, vac_deinit

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
  end subroutine iter_init

  subroutine iter_deinit
    use mephit_pert, only: L1_deinit, RT0_deinit

    call L1_deinit(pn)
    call RT0_deinit(Bn)
    call RT0_deinit(Bnplas)
    call RT0_deinit(jn)
  end subroutine iter_deinit

  subroutine iter_read
    use mephit_conf, only: datafile
    use mephit_pert, only: L1_read, RT0_read

    call L1_read(pn, datafile, 'iter/pn')
    call RT0_read(Bn, datafile, 'iter/Bn')
    call RT0_read(Bnplas, datafile, 'iter/Bnplas')
    call RT0_read(jn, datafile, 'iter/jn')
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
    use mephit_mesh, only: mesh
    use mephit_pert, only: L1_write, RT0_init, RT0_deinit, RT0_write, RT0_compute_tor_comp, &
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
       call arnoldi_break(ndim, conf%nkrylov, conf%ritz_threshold, conf%ritz_rel_err, &
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

    call FEM_compute_L2int(mesh%nedge, vac%Bn%DOF, L2int_Bnvac)
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
       call FEM_compute_L2int(mesh%nedge, Bn_diff%DOF, L2int_Bn_diff(kiter))
       if (kiter <= 1) then
          call L1_write(pn, datafile, 'iter/pn' // postfix, &
               'pressure (after iteration)', 'dyn cm^-2')
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
    if (rel_err >= conf%iter_rel_err) then
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

  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  subroutine compute_presn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: conf, logger
    use mephit_util, only: imun
    use mephit_mesh, only: fs, mesh, B0R_edge, B0phi_edge, B0Z_edge
    complex(dp), dimension(maxval(mesh%kp_max)) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(maxval(mesh%kp_max)) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kp, kedge, k
    integer, dimension(2 * maxval(mesh%kp_max)) :: irow, icol
    complex(dp), dimension(2 * maxval(mesh%kp_max)) :: aval
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    a = (0d0, 0d0)
    b = (0d0, 0d0)
    x = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          ! use midpoint of poloidal edge
          a(kp) = (B0R_edge(kedge) * mesh%edge_R(kedge) + B0Z_edge(kedge) * mesh%edge_Z(kedge)) / &
               (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2)
          x(kp) = -fs%dp_dpsi(kf) * a(kp) * Bn%DOF(kedge)
          b(kp) = imun * (mesh%n + imun * conf%damp) * B0phi_edge(kedge) / mesh%mid_R(kedge)
       end do
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       associate (ndim => mesh%kp_max(kf), nz => 2 * mesh%kp_max(kf))
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
       if (kf == 1) then ! first point on axis - average over enclosing flux surface
          pn%DOF(1) = sum(x(:mesh%kp_max(1))) / dble(mesh%kp_max(1))
       end if
       do kp = 1, mesh%kp_max(kf)
          pn%DOF(mesh%kp_low(kf) + kp) = x(kp)
       end do
    end do

    avg_rel_err = avg_rel_err / dble(mesh%npoint - 1)
    write (logger%msg, '("compute_presn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (logger%debug) call logger%write_msg
    if (allocated(resid)) deallocate(resid)
  end subroutine compute_presn

  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #mephit_conf::currn_file.
  subroutine compute_currn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mephit_conf, only: conf, logger
    use mephit_mesh, only: fs, mesh, B0phi_edge, B0_flux, j0phi_edge
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

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          kedge = mesh%kp_low(kf) + kp - 1
          jn%DOF(kedge) = j0phi_edge(kedge) / B0phi_edge(kedge) * Bn%DOF(kedge) + &
               clight * mesh%mid_R(kedge) / B0phi_edge(kedge) * &
               (pn%DOF(mesh%edge_node(2, kedge)) - pn%DOF(mesh%edge_node(1, kedge)))
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
          ! diagonal matrix element - edge i
          d(ke) = -1d0 - imun * (mesh%n + imun * conf%damp) * &
               mesh%area(ktri) * 0.5d0 * B0phi_edge(kedge) / B0_flux(kedge)
          ! additional term from edge i on source side
          x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * mesh%mid_R(kedge) / B0_flux(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) * Delta_p0 - Delta_pn) + j0phi_edge(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) - Bn%DOF(kedge) / B0_flux(kedge)))
          ! superdiagonal matrix element - edge o
          ktri = mesh%edge_tri(1, kedge)
          ke = mod(kt + mesh%kt_max(kf) - 2, mesh%kt_max(kf)) + 1
          du(ke) = 1d0 + imun * (mesh%n + imun * conf%damp) * &
               mesh%area(ktri) * 0.5d0 * B0phi_edge(kedge) / (-B0_flux(kedge))
          ! additional term from edge o on source side
          x(ke) = x(ke) - imun * mesh%n * mesh%area(ktri) * 0.5d0 * (clight * mesh%mid_R(kedge) / (-B0_flux(kedge)) * &
               (Bnphi_edge / B0phi_edge(kedge) * (-Delta_p0) - (-Delta_pn)) + j0phi_edge(kedge) * &
               (Bnphi_edge / B0phi_edge(kedge) - (-Bn%DOF(kedge)) / (-B0_flux(kedge))))
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

    call add_sheet_current
    call RT0_compute_tor_comp(jn)
  end subroutine compute_currn

  subroutine add_sheet_current
    use mephit_conf, only: conf_arr
    use mephit_mesh, only: mesh
    integer :: m, kf, kt, kedge, ke, k
    real(dp) :: B0_R, B0_Z, dum
    complex(dp) :: pn_outer, pn_inner

    do m = mesh%m_res_min, mesh%m_res_max
       kf = mesh%res_ind(m)
       if (abs(conf_arr%sheet_current_factor(m)) > 0d0) then
          do kt = 1, mesh%kt_max(kf)
             kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
             ke = mesh%shielding_kt_low(m) + kt
             do k = 1, mesh%GL_order
                call field(mesh%GL_R(k, kedge), 0d0, mesh%GL_Z(k, kedge), B0_R, dum, B0_Z, &
                     dum, dum, dum, dum, dum, dum, dum, dum, dum)
                pn_outer = sum(pn%DOF(mesh%kp_low(kf) + [mesh%shielding_L1_kp(0, k, ke), &
                     mod(mesh%shielding_L1_kp(0, k, ke), mesh%kp_max(kf)) + 1]) * &
                     [mesh%shielding_L1_weight(0, k, ke), &
                     1d0 - mesh%shielding_L1_weight(0, k, ke)])
                pn_inner = sum(pn%DOF(mesh%kp_low(kf-1) + [mesh%shielding_L1_kp(-1, k, ke), &
                     mod(mesh%shielding_L1_kp(-1, k, ke), mesh%kp_max(kf-1)) + 1]) * &
                     [mesh%shielding_L1_weight(-1, k, ke), &
                     1d0 - mesh%shielding_L1_weight(-1, k, ke)])
                jn%DOF(kedge) = jn%DOF(kedge) + mesh%GL_weights(k) * &
                     mesh%shielding_coeff(m) * conf_arr%sheet_current_factor(m) * (pn_outer - pn_inner) * &
                     (B0_R * mesh%edge_Z(kedge) - B0_Z * mesh%edge_R(kedge)) * mesh%GL_R(k, kedge)
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

  subroutine check_furth(jn, Bmn_plas)
    use mephit_conf, only: conf, datafile
    use mephit_util, only: imun, clight
    use mephit_mesh, only: equil, fs_half, mesh, cache
    use mephit_pert, only: RT0_t, vec_polmodes_t
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(RT0_t), intent(in) :: jn
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
            sheet_flux(kf) = sheet_flux(kf) + mesh%area(s%ktri) * jn%comp_phi(s%ktri) * &
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
    complex(dp) :: jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
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
            call RT0_interp(s%ktri, jn, s%R, s%Z, jn_R, jn_Z, jn_phi)
            ! include h^phi in current density
            jn_par = (jn_R * s%B0_R + jn_Z * s%B0_Z + jn_phi * s%B0_phi) * &
                 s%B0_phi / s%B0_2 * s%sqrt_g / s%R
            jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * exp(imun * abs(m) * s%theta)
            jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * exp(-imun * abs(m) * s%theta)
            ! comparison with indirect calculation (Boozer and Nuehrenberg 2006)
            call RT0_interp(s%ktri, Bn, s%R, s%Z, Bn_R, Bn_Z, Bn_phi, &
                 dBnR_dR, dBnR_dZ, dBnZ_dR, dBnZ_dZ)
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
    complex(dp) :: jn_R, jn_Z, jn_phi, jn_par, Bn_R, Bn_Z, Bn_phi, &
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
            call RT0_interp(s%ktri, jn, s%R, s%Z, jn_R, jn_Z, jn_phi)
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
            call RT0_interp(s%ktri, Bn, s%R, s%Z, Bn_R, Bn_Z, Bn_phi)
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
