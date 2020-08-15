module magdif
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: Bnflux, Bnphi, magdif_init, magdif_cleanup, magdif_single, magdif_iterated

  character(len = 1024) :: magdif_bin_dir = '.'

  !> Pressure perturbation \f$ p_{n} \f$ in dyn cm^-1.
  !>
  !> Values are taken at each mesh point and the indexing scheme is the same as for
  !> #mesh_mod::mesh_point.
  complex(dp), allocatable :: presn(:)

  !> Edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$ in G cm^2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: Bnflux(:,:)

  !> Physical toroidal component of magnetic perturbation \f$ B_{n (\phi)} \f$ in G.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: Bnphi(:)

  !> Vacuum perturbation edge fluxes \f$ R \vec{B}_{n} \cdot \vec{n} \f$ in G cm^2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: Bnflux_vac(:,:)

  !> Physical toroidal component of vacuum magnetic perturbation \f$ B_{n(\phi)} \f$ in G.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: Bnphi_vac(:)

  !> Edge currents \f$ R \vec{j}_{n} \cdot \vec{n} \f$ in statampere.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  complex(dp), allocatable :: jnflux(:,:)

  !> Physical toroidal component of current perturbation \f$ j_{n (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are taken at each triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element.
  complex(dp), allocatable :: jnphi(:)

contains

  !> Initialize magdif module
  subroutine magdif_init(bin_dir)
    use input_files, only: gfile
    use magdif_conf, only: conf, log, magdif_log, decorate_filename
    use magdif_mesh, only: equil, cache_equilibrium_field, flux_func_cache_check, &
         init_flux_variables, compute_j0phi, check_curr0, check_safety_factor
    use magdif_pert, only: compute_Bn_nonres, write_vector_dof, write_vector_plot
    character(len = *), intent(in) :: bin_dir

    magdif_bin_dir = bin_dir
    log = magdif_log('-', conf%log_level, conf%quiet)

    ! only depends on config variables
    call read_mesh
    call load_magdata_in_symfluxcoord

    ! depends on mesh data
    call cache_equilibrium_field

    ! needs initialized field_eq
    call equil%read(gfile)
    call equil%classify
    if (equil%cocos%index /= 3) then
       write (log%msg, '("GEQDSK file ", a, " is not conforming to COCOS 3")') trim(gfile)
       if (log%err) call log%write
       error stop
    end if

    ! depends on mesh data, equilibrium field and G EQDSK profiles
    call flux_func_cache_check
    call init_flux_variables

    ! depends on flux variables
    call compute_j0phi
    call check_curr0
    call check_safety_factor

    if (conf%nonres) then
       call compute_Bn_nonres(Bnflux, Bnphi)
    else
       call read_Bn(conf%Bn_vac_file)
    end if
    Bnflux_vac = Bnflux
    Bnphi_vac = Bnphi
    call write_vector_dof(Bnflux_vac, Bnphi_vac, conf%Bn_vacout_file)
    call write_vector_plot(Bnflux_vac, Bnphi_vac, &
         decorate_filename(conf%Bn_vacout_file, 'plot_', ''))
    log%msg = 'magdif initialized'
    if (log%info) call log%write
  end subroutine magdif_init

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp
    use magdif_mesh, only: B0R, B0phi, B0Z, B0R_Omega, B0phi_Omega, B0Z_Omega, B0flux, &
         j0phi
    use magdif_conf, only: log
    if (allocated(B0r)) deallocate(B0r)
    if (allocated(B0phi)) deallocate(B0phi)
    if (allocated(B0z)) deallocate(B0z)
    if (allocated(B0r_Omega)) deallocate(B0r_Omega)
    if (allocated(B0phi_Omega)) deallocate(B0phi_Omega)
    if (allocated(B0z_Omega)) deallocate(B0z_Omega)
    if (allocated(B0flux)) deallocate(B0flux)
    if (allocated(presn)) deallocate(presn)
    if (allocated(jnflux)) deallocate(jnflux)
    if (allocated(Bnflux)) deallocate(Bnflux)
    if (allocated(Bnphi)) deallocate(Bnphi)
    if (allocated(Bnflux_vac)) deallocate(Bnflux_vac)
    if (allocated(Bnphi_vac)) deallocate(Bnphi_vac)
    if (allocated(jnphi)) deallocate(jnphi)
    if (allocated(j0phi)) deallocate(j0phi)
    if (allocated(mesh_point)) deallocate(mesh_point)
    if (allocated(mesh_element)) deallocate(mesh_element)
    if (allocated(mesh_element_rmp)) deallocate(mesh_element_rmp)
    log%msg = 'magdif cleanup finished'
    if (log%info) call log%write
  end subroutine magdif_cleanup

  subroutine magdif_single
    use magdif_conf, only: conf, decorate_filename
    use magdif_pert, only: write_vector_plot
    ! compute pressure based on previous perturbation field
    call compute_presn
    ! compute currents based on previous perturbation field
    call compute_currn
    ! use field code to generate new field from currents
    call compute_Bn
    ! read new bnflux from field code
    call read_Bn(conf%Bn_file)
    call write_vector_plot(Bnflux, Bnphi, decorate_filename(conf%Bn_file, 'plot_', ''))
  end subroutine magdif_single

  subroutine magdif_iterated
    use arnoldi_mod, only: ieigen, ngrow, tol, eigvecs  ! arnoldi.f90
    use mesh_mod, only: ntri
    use magdif_mesh, only: mesh
    use magdif_pert, only: check_div_free, check_redundant_edges, &
         write_vector_dof, write_vector_plot, write_vector_plot_rect, write_scalar_dof
    use magdif_conf, only: conf, log, runmode_precon, decorate_filename, cmplx_fmt

    logical :: preconditioned
    integer :: kiter, ndim, i, j, info
    complex(dp) :: Bnflux_diff(ntri, 3)
    complex(dp) :: Bnphi_diff(ntri)
    complex(dp) :: Bn(mesh%nedge), Bn_prev(mesh%nedge)
    complex(dp) :: eigvals(conf%nritz)
    complex(dp), allocatable :: Lr(:,:), Yr(:,:)
    integer, allocatable :: ipiv(:)
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"

    ! system dimension: number of non-redundant edges in core plasma
    ndim = mesh%nedge
    preconditioned = runmode_precon == conf%runmode
    if (preconditioned) then
       tol = conf%ritz_threshold
       ! calculate eigenvectors
       ieigen = 1
       call arnoldi(ndim, conf%nritz, eigvals, next_iteration_arnoldi)
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
             call unpack_dof(Bnflux, Bnphi, eigvecs(:, i))
             call write_vector_dof(Bnflux, Bnphi, &
                  decorate_filename(conf%eigvec_file, '', postfix))
             call write_vector_plot(Bnflux, Bnphi, &
                  decorate_filename(conf%eigvec_file, 'plot_', postfix))
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

    call write_vector_dof(Bnflux_vac, Bnphi_vac, conf%Bn_diff_file)
    call compute_L2int
    call pack_dof(Bnflux_vac, Bn_prev)
    if (preconditioned) then
       Bn_prev = Bn_prev - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
            matmul(transpose(conjg(eigvecs(:, 1:ngrow))), Bn_prev)))
    end if
    do kiter = 0, conf%niter - 1
       write (log%msg, '("Iteration ", i2, " of ", i2)') kiter, conf%niter - 1
       if (log%info) call log%write
       write (postfix, postfix_fmt) kiter

       call next_iteration(ndim, Bn_prev, Bn)
       if (preconditioned) then
          Bn = Bn - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
               matmul(transpose(conjg(eigvecs(:, 1:ngrow))), Bn - Bn_prev)))
          call unpack_dof(Bnflux, Bnphi, Bn)
          call check_redundant_edges(Bnflux, .false., 'B_n')
          call check_div_free(Bnflux, Bnphi, conf%n, conf%rel_err_Bn, 'B_n')
       end if

       call unpack_dof(Bnflux_diff, Bnphi_diff, Bn - Bn_prev)
       call write_vector_dof(Bnflux_diff, Bnphi_diff, conf%Bn_diff_file)
       call compute_L2int
       call write_vector_dof(Bnflux, Bnphi, &
            decorate_filename(conf%Bn_file, '', postfix))
       if (kiter <= 1) then
          call write_vector_plot(Bnflux_diff, Bnphi_diff, &
               decorate_filename(conf%Bn_diff_file, 'plot_', postfix))
          call write_vector_plot(Bnflux, Bnphi, &
               decorate_filename(conf%Bn_file, 'plot_', postfix))
          call write_vector_plot(jnflux, jnphi, &
               decorate_filename(conf%currn_file, 'plot_', postfix))
          call write_poloidal_modes(jnflux, jnphi, &
               decorate_filename('currmn.dat', '', postfix))
          call write_scalar_dof(presn, decorate_filename(conf%presn_file, '', postfix))
       end if

       call pack_dof(Bnflux, Bn_prev)
    end do
    call write_vector_dof(Bnflux, Bnphi, conf%Bn_file)
    call write_vector_plot(Bnflux, Bnphi, decorate_filename(conf%Bn_file, 'plot_', ''))
    call write_vector_plot_rect(Bnflux, Bnphi, &
         decorate_filename(conf%Bn_file, 'rect_', ''))
    call write_poloidal_modes(Bnflux, Bnphi, 'Bmn.dat')
    call write_poloidal_modes(Bnflux_vac, Bnphi_vac, 'Bmn_vac.dat')
    call write_poloidal_modes(Bnflux - Bnflux_vac, Bnphi - Bnphi_vac, 'Bmn_plas.dat')
    call write_poloidal_modes(jnflux, jnphi, 'currmn.dat')
    call write_Ipar(2048)

    if (allocated(Lr)) deallocate(Lr)

  contains

    pure subroutine pack_dof(pol_flux, packed)
      complex(dp), intent(in) :: pol_flux(ntri, 3)
      complex(dp), intent(out) :: packed(mesh%nedge)
      integer :: kedge
      do kedge = 1, mesh%nedge
         packed(kedge) = pol_flux(mesh%edge_map2ktri(kedge, 1), mesh%edge_map2ke(kedge, 1))
      end do
    end subroutine pack_dof

    pure subroutine unpack_dof(pol_flux, tor_comp, packed)
      use magdif_util, only: imun
      use mesh_mod, only: mesh_element_rmp
      complex(dp), intent(out) :: pol_flux(ntri, 3), tor_comp(ntri)
      complex(dp), intent(in) :: packed(mesh%nedge)
      integer :: kedge, ktri
      do kedge = 1, mesh%nedge
         pol_flux(mesh%edge_map2ktri(kedge, 1), mesh%edge_map2ke(kedge, 1)) = &
              packed(kedge)
         if (mesh%edge_map2ktri(kedge, 2) > 0) then
            pol_flux(mesh%edge_map2ktri(kedge, 2), mesh%edge_map2ke(kedge, 2)) = &
                 -packed(kedge)
         end if
      end do
      do ktri = 1, mesh%kt_low(conf%nflux + 1)
         tor_comp(ktri) = sum(pol_flux(ktri, :)) * &
              imun / conf%n / mesh_element_rmp(ktri)%area
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
      call unpack_dof(Bnflux, Bnphi, xold)
      call magdif_single
      Bnflux = Bnflux + Bnflux_vac
      Bnphi = Bnphi + Bnphi_vac
      call pack_dof(Bnflux, xnew)
    end subroutine next_iteration

    ! computes B_(n+1) = K*(B_n + B_vac) ... as in kin2d.f90
    ! next_iteration in arnoldi_mod is still declared external and has no interface,
    ! so we use explicit-shape arrays here
    subroutine next_iteration_arnoldi(n, xold, xnew)
      integer, intent(in) :: n
      complex(dp), intent(in) :: xold(n)
      complex(dp), intent(out) :: xnew(n)
      call unpack_dof(Bnflux, Bnphi, xold)
      Bnflux = Bnflux + Bnflux_vac
      Bnphi = Bnphi + Bnphi_vac
      call magdif_single
      call pack_dof(Bnflux, xnew)
    end subroutine next_iteration_arnoldi
  end subroutine magdif_iterated

  !> Reads mesh points and triangles.
  !>
  !> #mesh_mod::npoint and #mesh_mod::mesh_point are read directly from an unformatted
  !> #magdif_config::point_file, while #mesh_mod::ntri and #mesh_mod::mesh_element are
  !> read directly from an unformatted #magdif_config::tri_file. #presn, #bnflux, #bnphi,
  !> #bnflux_vac, #bnphi_vac, #j0phi, #jnphi and #jnflux are allocated and initialized to
  !> zero. Deallocation is done in magdif_cleanup().
  subroutine read_mesh
    use magdif_conf, only: conf, log
    use magdif_mesh, only: fs, fs_half, mesh, &
         B0R, B0phi, B0Z, B0R_Omega, B0phi_Omega, B0Z_Omega, B0flux, j0phi
    use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, mesh_element_rmp
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    integer(HID_T) :: h5id_magdif
    integer :: ktri
    integer, dimension(:, :), allocatable :: tri_node, adj_tri, adj_edge, tri_li, tri_lo, tri_lf
    integer, dimension(:), allocatable :: tri_orient

    call h5_open('magdif.h5', h5id_magdif)
    call h5_get(h5id_magdif, 'mesh/R_O', mesh%R_O)
    call h5_get(h5id_magdif, 'mesh/Z_O', mesh%Z_O)
    call h5_get(h5id_magdif, 'mesh/R_min', mesh%R_min)
    call h5_get(h5id_magdif, 'mesh/Z_min', mesh%Z_min)
    call h5_get(h5id_magdif, 'mesh/R_max', mesh%R_max)
    call h5_get(h5id_magdif, 'mesh/Z_max', mesh%Z_max)
    call h5_get(h5id_magdif, 'mesh/nflux', conf%nflux)
    call h5_get(h5id_magdif, 'mesh/npoint', npoint)
    call h5_get(h5id_magdif, 'mesh/ntri', ntri)
    call h5_get(h5id_magdif, 'mesh/nedge', mesh%nedge)
    call h5_get(h5id_magdif, 'mesh/m_res_min', mesh%m_res_min)
    call h5_get(h5id_magdif, 'mesh/m_res_max', mesh%m_res_max)
    write (log%msg, '("nflux = ", i0, ", npoint = ", i0, ", ntri = ", i0)') &
         conf%nflux, npoint, ntri
    if (log%info) call log%write
    call fs%init(conf%nflux, .false.)
    call fs_half%init(conf%nflux, .true.)
    ! TODO: allocate deferred-shape arrays in hdf5_tools and skip allocation here
    allocate(mesh%kp_max(conf%nflux))
    allocate(mesh%kt_max(conf%nflux))
    allocate(mesh%kp_low(conf%nflux + 1))
    allocate(mesh%kt_low(conf%nflux + 1))
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(mesh_element_rmp(ntri))
    allocate(tri_node(3, ntri))
    allocate(adj_tri(3, ntri))
    allocate(adj_edge(3, ntri))
    allocate(tri_li(2, ntri))
    allocate(tri_lo(2, ntri))
    allocate(tri_lf(2, ntri))
    allocate(tri_orient(ntri))
    allocate(mesh%edge_map2global(ntri, 3))
    allocate(mesh%edge_map2ktri(mesh%nedge, 2))
    allocate(mesh%edge_map2ke(mesh%nedge, 2))
    call h5_get(h5id_magdif, 'mesh/kp_max', mesh%kp_max)
    call h5_get(h5id_magdif, 'mesh/kp_low', mesh%kp_low)
    call h5_get(h5id_magdif, 'mesh/kt_max', mesh%kt_max)
    call h5_get(h5id_magdif, 'mesh/kt_low', mesh%kt_low)
    call h5_get(h5id_magdif, 'mesh/node_R', mesh_point%rcoord)
    call h5_get(h5id_magdif, 'mesh/node_Z', mesh_point%zcoord)
    call h5_get(h5id_magdif, 'mesh/tri_node', tri_node)
    call h5_get(h5id_magdif, 'mesh/tri_node_F', mesh_element%knot_h)
    call h5_get(h5id_magdif, 'mesh/tri_li', tri_li)
    call h5_get(h5id_magdif, 'mesh/tri_lo', tri_lo)
    call h5_get(h5id_magdif, 'mesh/tri_lf', tri_lf)
    call h5_get(h5id_magdif, 'mesh/tri_ei', mesh_element_rmp%ei)
    call h5_get(h5id_magdif, 'mesh/tri_eo', mesh_element_rmp%eo)
    call h5_get(h5id_magdif, 'mesh/tri_ef', mesh_element_rmp%ef)
    call h5_get(h5id_magdif, 'mesh/tri_orient', tri_orient)
    call h5_get(h5id_magdif, 'mesh/adj_tri', adj_tri)
    call h5_get(h5id_magdif, 'mesh/adj_edge', adj_edge)
    call h5_get(h5id_magdif, 'mesh/edge_map2kedge', mesh%edge_map2global)
    call h5_get(h5id_magdif, 'mesh/edge_map2ktri', mesh%edge_map2ktri)
    call h5_get(h5id_magdif, 'mesh/edge_map2ke', mesh%edge_map2ke)
    call h5_get(h5id_magdif, 'mesh/tri_centr_R', mesh_element_rmp%R_Omega)
    call h5_get(h5id_magdif, 'mesh/tri_centr_Z', mesh_element_rmp%Z_Omega)
    call h5_get(h5id_magdif, 'mesh/det', mesh_element%det_3)
    call h5_get(h5id_magdif, 'cache/fs/psi', fs%psi)
    call h5_get(h5id_magdif, 'cache/fs/rad', fs%rad)
    call h5_get(h5id_magdif, 'cache/fs_half/psi', fs_half%psi)
    call h5_get(h5id_magdif, 'cache/fs_half/rad', fs_half%rad)
    call h5_close(h5id_magdif)
    ! TODO: remove intermediate when mesh_mod is refactored into magdif_mesh
    mesh_element_rmp%area = 0.5d0 * mesh_element%det_3
    do ktri = 1, ntri
       mesh_element(ktri)%i_knot = tri_node(:, ktri)
       mesh_element(ktri)%neighbour = adj_tri(:, ktri)
       mesh_element(ktri)%neighbour_edge = adj_edge(:, ktri)
       mesh_element_rmp(ktri)%li = tri_li(:, ktri)
       mesh_element_rmp(ktri)%lo = tri_lo(:, ktri)
       mesh_element_rmp(ktri)%lf = tri_lf(:, ktri)
    end do
    where (tri_orient == 0)
       mesh_element_rmp%orient = .false.
    elsewhere
       mesh_element_rmp%orient = .true.
    end where
    deallocate(tri_node)
    deallocate(adj_tri)
    deallocate(adj_edge)
    deallocate(tri_li)
    deallocate(tri_lo)
    deallocate(tri_lf)
    deallocate(tri_orient)

    allocate(B0r(ntri, 3))
    allocate(B0phi(ntri, 3))
    allocate(B0z(ntri, 3))
    allocate(B0r_Omega(ntri))
    allocate(B0phi_Omega(ntri))
    allocate(B0z_Omega(ntri))
    allocate(B0flux(ntri, 3))
    allocate(presn(npoint))
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    allocate(Bnflux_vac(ntri, 3))
    allocate(Bnphi_vac(ntri))
    allocate(jnphi(ntri))
    allocate(j0phi(ntri, 3))
    allocate(jnflux(ntri, 3))

    B0r = 0d0
    B0phi = 0d0
    B0z = 0d0
    B0r_Omega = 0d0
    B0phi_Omega = 0d0
    B0z_Omega = 0d0
    B0flux = 0d0
    presn = 0d0
    Bnflux = 0d0
    Bnphi = 0d0
    Bnflux_vac = 0d0
    Bnphi_vac = 0d0
    jnphi = 0d0
    j0phi = 0d0
    jnflux = 0d0
  end subroutine read_mesh

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi via an external program. No data
  !> is read yet; this is done by read_bn().
  !>
  !> Currently, this subroutine Calls FreeFem++ via shell script maxwell.sh in the current
  !> directory and handles the exit code. For further information see maxwell.sh and the
  !> script called therein.
  subroutine compute_Bn
    use magdif_conf, only: conf, log
    use magdif_mesh, only: mesh
    integer :: stat = 0, dummy = 0
    character(len = 1024) :: fem_cmd
    write (fem_cmd, '(a, "/maxwell.sh -n ", i0, " -N ", i0, " -J ", a, " -B ", a)') &
         trim(magdif_bin_dir), conf%n, mesh%kt_low(conf%nflux + 1), trim(conf%currn_file), &
         trim(conf%Bn_file)
    call execute_command_line(fem_cmd, exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       write (log%msg, '("FreeFem++ failed with exit code ", i0)') stat
       if (log%err) call log%write
       error stop
    end if
  end subroutine compute_Bn

  subroutine compute_L2int
    use magdif_conf, only: conf, log
    use magdif_mesh, only: mesh
    integer :: stat = 0, dummy = 0
    character(len = 1024) :: L2int_cmd
    write (L2int_cmd, '(a, "/L2int.sh -N ", i0, " -B ", a, " -C ", a)') &
         trim(magdif_bin_dir), mesh%kt_low(conf%nflux + 1), trim(conf%Bn_diff_file), &
         trim(conf%conv_file)
    call execute_command_line(L2int_cmd, exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       write (log%msg, '("FreeFem++ failed with exit code ", i0)') stat
       if (log%err) call log%write
       error stop
    end if
  end subroutine compute_L2int

  !> Reads fluxes of perturbation field and checks divergence-freeness.
  !>
  !> @param filename name of the formatted file containing the data
  !>
  !> Line numbers in the given file correspond to the global triangle index of
  !> #mesh_mod::mesh_element. The content of each line is read into #bnflux and #bnphi
  !> with numbering of edges in #bnflux as in #mesh_mod::mesh_element and the imaginary
  !> part of each value immediately following the real part.
  subroutine read_Bn(filename)
    use magdif_conf, only: conf, longlines
    use mesh_mod, only: ntri, mesh_element_rmp
    use magdif_pert, only: check_div_free, check_redundant_edges
    character(len = 1024) :: filename
    integer :: ktri, fid
    real(dp) :: dummy8(8)

    open(newunit = fid, file = filename, recl = longlines, status = 'old')
    do ktri = 1, ntri
       read (fid, *) dummy8
       Bnflux(ktri,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(ktri,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(ktri,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(ktri) = cmplx(dummy8(7), dummy8(8), dp) / mesh_element_rmp(ktri)%area
    end do
    close(fid)

    call check_redundant_edges(Bnflux, .false., 'B_n')
    call check_div_free(Bnflux, Bnphi, conf%n, conf%rel_err_Bn, 'B_n')
  end subroutine read_Bn

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
    use magdif_pert, only: write_scalar_dof
    use mesh_mod, only: knot, triangle_rmp, mesh_point, mesh_element_rmp
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
    type(triangle_rmp) :: tri
    type(knot) :: base, tip
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, conf%nflux
       inner: do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          if (.not. tri%orient) cycle inner
          ! use midpoint of edge f
          base = mesh_point(tri%lf(1))
          tip = mesh_point(tri%lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord

          kp = tri%lf(1) - mesh%kp_low(kf)

          a(kp) = (B0r(ktri, tri%ef) * lr + B0z(ktri, tri%ef) * lz) / (lr ** 2 + lz ** 2)

          x(kp) = -fs%dp_dpsi(kf) * a(kp) * Bnflux(ktri, tri%ef)

          b(kp) = imun * (conf%n + imun * conf%damp) * B0phi(ktri, tri%ef) / r
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
          presn(1) = sum(x) / size(x)
       end if
       do kp = 1, mesh%kp_max(kf)
          presn(mesh%kp_low(kf) + kp) = x(kp)
       end do
    end do

    avg_rel_err = avg_rel_err / sum(mesh%kp_max(1:conf%nflux))
    write (log%msg, '("compute_presn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (log%debug) call log%write
    if (allocated(resid)) deallocate(resid)

    call write_scalar_dof(presn, conf%presn_file)
  end subroutine compute_presn

  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #magdif_conf::currn_file.
  subroutine compute_currn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mesh_mod, only: triangle, triangle_rmp, mesh_element, mesh_element_rmp, mesh_point
    use magdif_conf, only: conf, log, decorate_filename
    use magdif_mesh, only: fs, mesh, B0phi, B0flux, j0phi
    use magdif_pert, only: check_div_free, check_redundant_edges, write_vector_dof, &
         write_vector_plot
    use magdif_util, only: imun, clight
    complex(dp), dimension(2 * conf%nkpol) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(2 * conf%nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, ktri
    integer :: nz
    integer, dimension(4 * conf%nkpol) :: irow, icol
    complex(dp), dimension(4 * conf%nkpol) :: aval
    type(triangle) :: elem
    type(triangle_rmp) :: tri
    complex(dp) :: Bnphi_Gamma
    real(dp) :: r
    real(dp), parameter :: small = tiny(0d0)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, conf%nflux ! loop through flux surfaces
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          elem = mesh_element(ktri)
          tri = mesh_element_rmp(ktri)
          ! first term on source side: flux through edge f
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          jnflux(ktri, tri%ef) = j0phi(ktri, tri%ef) / B0phi(ktri, tri%ef) * &
               Bnflux(ktri, tri%ef) + &
               clight * r / B0phi(ktri, tri%ef) * (presn(tri%lf(2)) - presn(tri%lf(1)))
          x(kt) = -jnflux(ktri, tri%ef)
          ! diagonal matrix element - edge i
          d(kt) = -1d0 - imun * (conf%n + imun * conf%damp) * tri%area * 0.5d0 * &
               B0phi(ktri, tri%ei) / B0flux(ktri, tri%ei)
          ! additional term from edge i on source side
          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(tri%ei)))
          x(kt) = x(kt) - imun * conf%n * tri%area * 0.5d0 * (clight * r / B0flux(ktri, tri%ei) * &
               (Bnphi_Gamma / B0phi(ktri, tri%ei) * (fs%p(kf) - fs%p(kf-1)) - &
               (presn(tri%li(2)) - presn(tri%li(1)))) + j0phi(ktri, tri%ei) * &
               (Bnphi_Gamma / B0phi(ktri, tri%ei) - Bnflux(ktri, tri%ei) / B0flux(ktri, tri%ei)))
          ! superdiagonal matrix element - edge o
          du(kt) = 1d0 + imun * (conf%n + imun * conf%damp) * tri%area * 0.5d0 * &
               B0phi(ktri, tri%eo) / B0flux(ktri, tri%eo)
          ! additional term from edge o on source side
          r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(tri%eo)))
          x(kt) = x(kt) - imun * conf%n * tri%area * 0.5d0 * (clight * r / B0flux(ktri, tri%eo) * &
               (Bnphi_Gamma / B0phi(ktri, tri%eo) * (fs%p(kf-1) - fs%p(kf)) - &
               (presn(tri%lo(2)) - presn(tri%lo(1)))) + j0phi(ktri, tri%eo) * &
               (Bnphi_Gamma / B0phi(ktri, tri%eo) - Bnflux(ktri, tri%eo) / B0flux(ktri, tri%eo)))
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
          tri = mesh_element_rmp(ktri)
          jnflux(ktri, tri%ei) = -x(kt)
          jnflux(ktri, tri%eo) = x(mod(kt, mesh%kt_max(kf)) + 1)
          jnphi(ktri) = sum(jnflux(ktri, :)) * imun / conf%n / tri%area
       end do
    end do
    avg_rel_err = avg_rel_err / sum(mesh%kt_max(1:conf%nflux))
    write (log%msg, '("compute_currn: diagonalization max_rel_err = ", ' // &
         'es24.16e3, ", avg_rel_err = ", es24.16e3)') max_rel_err, avg_rel_err
    if (log%debug) call log%write
    if (allocated(resid)) deallocate(resid)

    call add_sheet_current
    call check_redundant_edges(jnflux, .false., 'jnflux')
    call check_div_free(jnflux, jnphi, conf%n, conf%rel_err_currn, 'jnflux')
    call write_vector_dof(jnflux, jnphi, conf%currn_file)
    call write_vector_plot(jnflux, jnphi, decorate_filename(conf%currn_file, 'plot_', ''))
  end subroutine compute_currn

  subroutine add_sheet_current
    use mesh_mod, only: triangle_rmp, mesh_element_rmp, mesh_element
    use magdif_conf, only: conf, conf_arr
    use magdif_util, only: imun
    use magdif_mesh, only: mesh, B0flux
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri, tri_adj
    complex(dp) :: presn_half

    do kf = 1, conf%nflux
       if (mesh%m_res(kf) > 0) then
          if (abs(conf_arr%sheet_current_factor(mesh%m_res(kf))) > 0d0) then
             do kt = 1, mesh%kt_max(kf)
                ktri = mesh%kt_low(kf) + kt
                tri = mesh_element_rmp(ktri)
                ! add sheet current on edge i
                if (mod(kt, 2) == 0) then
                   ! edge i is diagonal
                   if (tri%orient) then
                      presn_half = 0.5d0 * sum(presn(tri%lf(:)))
                   else
                      tri_adj = mesh_element_rmp(mesh_element(ktri)%neighbour(tri%ei))
                      presn_half = 0.5d0 * sum(presn(tri_adj%lf(:)))
                   end if
                else
                   presn_half = presn(tri%li(2))
                end if
                jnflux(ktri, tri%ei) = jnflux(ktri, tri%ei) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * &
                     B0flux(ktri, tri%ei) * presn_half
                ! add sheet current on edge o
                if (mod(kt, 2) == 1) then
                   ! edge o is diagonal
                   if (tri%orient) then
                      presn_half = 0.5d0 * sum(presn(tri%lf(:)))
                   else
                      tri_adj = mesh_element_rmp(mesh_element(ktri)%neighbour(tri%eo))
                      presn_half = 0.5d0 * sum(presn(tri_adj%lf(:)))
                   end if
                else
                   presn_half = presn(tri%lo(1))
                end if
                jnflux(ktri, tri%eo) = jnflux(ktri, tri%eo) + &
                     conf_arr%sheet_current_factor(mesh%m_res(kf)) * &
                     B0flux(ktri, tri%eo) * presn_half
                ! adjust toroidal current density
                jnphi(ktri) = sum(jnflux(ktri, :)) * imun / conf%n / tri%area
             end do
          end if
       end if
    end do
  end subroutine add_sheet_current

  subroutine write_poloidal_modes(pol_flux, tor_comp, outfile)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, longlines, decorate_filename
    use magdif_util, only: imun, clight, bent_cyl2straight_cyl
    use magdif_mesh, only: equil, fs, fs_half, mesh, point_location
    use magdif_pert, only: interp_RT0
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer, parameter :: mmax = 24
    character(len = 20) :: fmt
    complex(dp), dimension(-mmax:mmax) :: coeff_rad, coeff_pol, coeff_tor, fourier_basis
    integer :: kf, kt, ktri, m, fid_rad, fid_pol, fid_tor
    integer :: kilca_m_res, fid_furth  ! only KiLCA geometry
    type(triangle_rmp) :: tri
    complex(dp) :: comp_R, comp_Z, comp_rad, comp_pol, comp_tor
    complex(dp) :: sheet_flux  ! only KiLCA geometry
    real(dp) :: R, Z, theta, B0_R, B0_Z, dum
    real(dp) :: k_z, k_theta  ! only KiLCA geometry
    real(dp) :: sqrt_g, dR_dtheta, dZ_dtheta, q, q_sum  ! only ASDEX geometry

    write (fmt, '(a, i3, a)') '(', 4 * mmax + 2 + 2, '(es23.15e3, 1x))'
    if (conf%kilca_scale_factor /= 0) then
       kilca_m_res = -equil%cocos%sgn_q * abs(conf%kilca_pol_mode)
       k_z = conf%n / equil%rcentr
       open(newunit = fid_rad, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_r'))
       open(newunit = fid_pol, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_theta'))
       open(newunit = fid_tor, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_z'))
       open(newunit = fid_furth, recl = longlines, status = 'replace', &
            file = decorate_filename(outfile, 'furth_', ''))
    else
       kilca_m_res = 0  ! to suppress compiler warning
       open(newunit = fid_rad, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_psi'))
       open(newunit = fid_pol, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_theta'))
       open(newunit = fid_tor, recl = 3 * longlines, status = 'replace', &
            file = decorate_filename(outfile, '', '_phi'))
    end if
    do kf = 1, conf%nflux
       coeff_rad = 0d0
       coeff_pol = 0d0
       coeff_tor = 0d0
       q_sum = 0d0
       sheet_flux = 0d0
       do kt = 1, mesh%kt_max(kf)
          theta = (dble(kt) - conf%debug_pol_offset) / dble(mesh%kt_max(kf)) * 2d0 * pi
          fourier_basis = [(exp(-imun * m * theta), m = -mmax, mmax)]
          if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
             R = equil%rmaxis + fs_half%rad(kf) * cos(theta)
             Z = equil%zmaxis + fs_half%rad(kf) * sin(theta)
          else
             ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
             call magdata_in_symfluxcoord_ext(2, dum, fs_half%psi(kf) - fs%psi(0), &
                  theta, q, dum, sqrt_g, dum, dum, R, dum, dR_dtheta, Z, dum, dZ_dtheta)
          end if
          call field(R, 0d0, Z, B0_R, dum, B0_Z, dum, dum, dum, dum, dum, dum, dum, &
               dum, dum)
          ktri = point_location(R, Z)
          tri = mesh_element_rmp(ktri)
          call interp_RT0(ktri, pol_flux, R, Z, comp_R, comp_Z)
          if (conf%kilca_scale_factor /= 0) then
             call bent_cyl2straight_cyl(comp_R, tor_comp(ktri), comp_Z, theta, &
                  comp_rad, comp_pol, comp_tor)
             sheet_flux = sheet_flux + tri%area * comp_tor * fourier_basis(kilca_m_res)
          else
             comp_rad = (comp_R * B0_Z - comp_Z * B0_R) * R * sqrt_g * q
             comp_pol = comp_R * dR_dtheta + comp_Z * dZ_dtheta
             comp_tor = tor_comp(ktri)
          end if
          coeff_rad = coeff_rad + comp_rad * fourier_basis
          coeff_pol = coeff_pol + comp_pol * fourier_basis
          coeff_tor = coeff_tor + comp_tor * fourier_basis
          q_sum = q_sum + q
       end do
       coeff_rad = coeff_rad / mesh%kt_max(kf)
       coeff_pol = coeff_pol / mesh%kt_max(kf)
       coeff_tor = coeff_tor / mesh%kt_max(kf)
       q = q_sum / mesh%kt_max(kf)
       if (conf%kilca_scale_factor /= 0) then
          write (fid_rad, fmt) fs_half%rad(kf), q, real(coeff_rad), aimag(coeff_rad)
          write (fid_pol, fmt) fs_half%rad(kf), q, real(coeff_pol), aimag(coeff_pol)
          write (fid_tor, fmt) fs_half%rad(kf), q, real(coeff_tor), aimag(coeff_tor)
          k_theta = kilca_m_res / fs_half%rad(kf)
          sheet_flux = -2d0 * imun / clight / k_theta * sheet_flux
          write (fid_furth, '(7(1x, es24.16e3))') fs_half%rad(kf), k_z, k_theta, &
               coeff_rad(-kilca_m_res), sheet_flux
       else
          write (fid_rad, fmt) fs_half%psi(kf), q, real(coeff_rad), aimag(coeff_rad)
          write (fid_pol, fmt) fs_half%psi(kf), q, real(coeff_pol), aimag(coeff_pol)
          write (fid_tor, fmt) fs_half%psi(kf), q, real(coeff_tor), aimag(coeff_tor)
       end if
    end do
    close(fid_rad)
    close(fid_pol)
    close(fid_tor)
    if (conf%kilca_pol_mode /= 0) then
       close(fid_furth)
    end if
  end subroutine write_poloidal_modes


  !> calculate parallel current (density) on a finer grid
  subroutine write_Ipar(rad_resolution)
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: ntheta
    use magdif_conf, only: conf, conf_arr, longlines, decorate_filename
    use magdif_util, only: imun, linspace, bent_cyl2straight_cyl
    use magdif_mesh, only: equil, fs, mesh, point_location
    use magdif_pert, only: interp_RT0
    integer, intent(in) :: rad_resolution
    integer :: kf_min, kf_max, krad, kt, ktri, fid_jpar
    real(dp) :: theta, R, Z, dum, B0_R, B0_phi, B0_Z, B0_2, dB0R_dR, dB0R_dZ, dB0phi_dR, &
         dB0phi_dZ, dB0Z_dR, dB0Z_dZ, dB0R_drad, dB0phi_drad, dB0Z_drad, B0_dB0_drad, &
         B0_theta, dB0theta_drad, dhz2_drad, dradhthetahz_drad
    complex(dp) :: jn_R, jn_Z, jn_par, Bn_R, Bn_Z, Bn_rad, Bn_pol, Bn_tor, part_int, bndry
    real(dp), dimension(rad_resolution) :: rad
    complex(dp), dimension(rad_resolution) :: jmn_par_neg, jmn_par_pos, &
         part_int_neg, part_int_pos, bndry_neg, bndry_pos

    if (conf%kilca_scale_factor /= 0) return
    kf_min = mesh%res_ind(abs(conf%kilca_pol_mode)) &
         - conf_arr%additions(abs(conf%kilca_pol_mode)) &
         - conf_arr%deletions(abs(conf%kilca_pol_mode))
    if (kf_min < 1) kf_min = 1
    kf_max = mesh%res_ind(abs(conf%kilca_pol_mode)) &
         + conf_arr%additions(abs(conf%kilca_pol_mode)) &
         + conf_arr%deletions(abs(conf%kilca_pol_mode))
    if (kf_max > conf%nflux) kf_max = conf%nflux
    rad = linspace(fs%rad(kf_min), fs%rad(kf_max), rad_resolution, 0, 0)
    jmn_par_neg = (0d0, 0d0)
    jmn_par_pos = (0d0, 0d0)
    part_int_neg = (0d0, 0d0)
    part_int_pos = (0d0, 0d0)
    bndry_neg = (0d0, 0d0)
    bndry_pos = (0d0, 0d0)
    open(newunit = fid_jpar, status = 'replace', recl = longlines, &
         file = 'currn_par.dat')
    do krad = 1, rad_resolution
       do kt = 1, ntheta
          theta = 2d0 * pi * (dble(kt) - 0.5d0) / dble(ntheta)
          R = equil%rmaxis + rad(krad) * cos(theta)
          Z = equil%zmaxis + rad(krad) * sin(theta)
          ktri = point_location(R, Z)
          call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dB0R_dR, dum, dB0R_dZ, dB0phi_dR, &
               dum, dB0phi_dZ, dB0Z_dR, dum, dB0Z_dZ)
          B0_2 = B0_R * B0_R + B0_Z * B0_Z + B0_phi * B0_phi
          call interp_RT0(ktri, jnflux, R, Z, jn_R, jn_Z)
          ! include h^z in current density
          jn_par = (jn_R * B0_R + jn_Z * B0_Z + jnphi(ktri) * B0_phi) * B0_phi / B0_2
          jmn_par_neg(krad) = jmn_par_neg(krad) + jn_par * &
               exp(imun * abs(conf%kilca_pol_mode) * theta)
          jmn_par_pos(krad) = jmn_par_pos(krad) + jn_par * &
               exp(-imun * abs(conf%kilca_pol_mode) * theta)
          ! comparison with indirect calculation (Ampere's law and integration by parts)
          dB0R_drad = dB0R_dR * cos(theta) + dB0R_dZ + sin(theta)
          dB0phi_drad = dB0phi_dR * cos(theta) + dB0phi_dZ + sin(theta)
          dB0Z_drad = dB0Z_dR * cos(theta) + dB0Z_dZ + sin(theta)
          B0_dB0_drad = B0_R * dB0R_drad + B0_phi * dB0phi_drad + B0_Z * dB0Z_drad
          B0_theta = B0_Z * cos(theta) - B0_R * sin(theta)
          dB0theta_drad = dB0Z_dR * cos(theta) ** 2 - dB0R_dZ * sin(theta) ** 2 + &
               (dB0Z_dZ - dB0R_dR) * sin(theta) * cos(theta)
          dhz2_drad = 2d0 * B0_phi * (B0_2 * dB0phi_drad - B0_phi * B0_dB0_drad) / &
               B0_2 ** 2
          dradhthetahz_drad = B0_theta * B0_phi / B0_2 + rad(krad) * &
               ((B0_theta * dB0phi_drad + dB0theta_drad * B0_phi) / B0_2 - &
               2d0 * B0_theta * B0_phi * B0_dB0_drad / B0_2 ** 2)
          call interp_RT0(ktri, Bnflux, R, Z, Bn_R, Bn_Z)
          call bent_cyl2straight_cyl(Bn_R, Bnphi(ktri), Bn_Z, theta, &
               Bn_rad, Bn_pol, Bn_tor)
          part_int = -rad(krad) * Bn_pol * dhz2_drad + Bn_tor * dradhthetahz_drad
          part_int_neg(krad) = part_int_neg(krad) + (part_int + imun * B0_phi * &
               (dble(conf%n) / equil%rmaxis * rad(krad) * B0_theta + &
               abs(conf%kilca_pol_mode) * B0_phi) * Bn_rad / B0_2) * &
               exp(imun * abs(conf%kilca_pol_mode) * theta)
          part_int_pos(krad) = part_int_pos(krad) + (part_int + imun * B0_phi * &
               (dble(conf%n) / equil%rmaxis * rad(krad) * B0_theta - &
               abs(conf%kilca_pol_mode) * B0_phi) * Bn_rad / B0_2) * &
               exp(-imun * abs(conf%kilca_pol_mode) * theta)
          bndry = B0_phi * rad(krad) * (B0_phi * Bn_pol - B0_theta * Bn_tor) / B0_2
          bndry_neg(krad) = bndry_neg(krad) + bndry * &
               exp(imun * abs(conf%kilca_pol_mode) * theta)
          bndry_pos(krad) = bndry_pos(krad) + bndry * &
               exp(-imun * abs(conf%kilca_pol_mode) * theta)
       end do
       jmn_par_neg(krad) = jmn_par_neg(krad) / dble(ntheta)
       jmn_par_pos(krad) = jmn_par_pos(krad) / dble(ntheta)
       part_int_neg(krad) = part_int_neg(krad) / dble(ntheta)
       part_int_pos(krad) = part_int_pos(krad) / dble(ntheta)
       bndry_neg(krad) = bndry_neg(krad) / dble(ntheta)
       bndry_pos(krad) = bndry_pos(krad) / dble(ntheta)
       write (fid_jpar, '(13(1x, es24.16e3))') rad(krad), &
            jmn_par_neg(krad), jmn_par_pos(krad), &
            part_int_neg(krad), part_int_pos(krad), bndry_neg(krad), bndry_pos(krad)
    end do
    close(fid_jpar)
  end subroutine write_Ipar
end module magdif
