module magdif
  use from_nrtype, only: dp                            ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi                              ! PRELOAD/SRC/orbit_mod.f90
  use mesh_mod, only: npoint, ntri, knot, triangle, &  ! PRELOAD/SRC/mesh_mod.f90
       mesh_point, mesh_element, triangle_rmp, mesh_element_rmp
  use sparse_mod, only: sparse_solve, sparse_matmul    ! MHD/SRC/sparse_mod.f90
  use magdif_config                                    ! MHD/SRC/magdif_config.f90
  use magdif_util                                      ! MHD/SRC/magdif_util.f90

  implicit none

  type(g_eqdsk), allocatable :: efit
  type(flux_func), allocatable :: fluxvar
  type(flux_func_cache), allocatable :: fs
  type(flux_func_cache), allocatable :: fs_half

  !> Poloidal mode number \f$ m \f$ (dimensionless) in resonance at given flux surface.
  !>
  !> Indexing is the same as for #q, on which the values depend. If no resonances are
  !> expected at a given index, #m_res is 0.
  integer, allocatable :: m_res(:)

  !> \f$ R \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0r(:,:)

  !> \f$ \phi \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0phi(:,:)

  !> \f$ Z \f$ component of equilibrium magnetic field \f$ B_{0} \f$.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: B0z(:,:)

  !> \f$ \phi \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0r_Omega(:)

  !> \f$ Z \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0phi_Omega(:)

  !> \f$ R \f$ component of equilibrium magnetic field \f$ B_{0} (\Omega) \f$.
  !>
  !> Values are stored seprately for each triangle and the indexing scheme is the same as
  !> for #mesh_mod::mesh_element.
  real(dp), allocatable :: B0z_Omega(:)

  real(dp), allocatable :: B0flux(:,:)

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

  !> Physical toroidal component of equilibrium current \f$ j_{0 (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: j0phi(:,:)

  !> Number of knots on the flux surface given by the array index.
  !>
  !> The array index ranges from 1 for the innermost flux surface to
  !> #magdif_config::nflux+2 for consistency with #kp_low.
  integer, allocatable :: kp_max(:)

  !> Number of triangles inside the flux surface given by the array index.
  !>
  !> The array index ranges from 1 for the innermost flux surface to
  !> #magdif_config::nflux+1 for consistency with #kt_low.
  integer, allocatable :: kt_max(:)

  !> Global index of the last knot of the previous flux surface given by the array index.
  !>
  !> The global index of knots in #mesh_mod::mesh_point on the flux surface kf runs from
  !> #kp_low (kf)+1 to #kp_low (kf)+#kp_max (kf), so #kp_low is determined by cumulatively
  !> adding consecutive values of #kp_max. The array index ranges from 1, giving the
  !> global index of the knot on the magnetic axis (which has to be 1), to
  !> #magdif_config::nflux+1, giving the last knot on the flux surface just outside the
  !> last closed flux surface. The latter is needed for some interpolations.
  integer, allocatable :: kp_low(:)

  !> Global index of the last triangle of the previous flux surface given by the array
  !> index.
  !>
  !> The global index of triangles in #mesh_mod::mesh_element inside the flux surface kf
  !> runs from #kt_low (kf)+1 to #kt_low (kf)+#kt_max (kf), so #kt_low is determined by
  !> cumulatively adding consecutive values of #kt_max. The array index ranges from 1,
  !> giving the global index of the non-existent triangle on the magnetic axis (which is
  !> therefore 0), to #magdif_config::nflux+1, giving the last triangle just outside the
  !> last closed flux surface. The latter is needed for some interpolations.
  integer, allocatable :: kt_low(:)

contains

  !> Initialize magdif module
  subroutine magdif_init
    use input_files, only: gfile
    integer :: m

    call log_open

    ! only depends on config variables
    call init_indices
    if (kilca_scale_factor /= 0) then
       R0 = R0 * kilca_scale_factor
       n = n * kilca_scale_factor
    end if

    call read_mesh

    ! depends on mesh data
    call cache_mesh_data
    call cache_equilibrium_field

    ! needs initialized field_eq
    allocate(efit)
    call efit%read(gfile)

    ! depends on mesh data, equilibrium field and EFIT profiles
    call init_flux_variables

    ! depends on safety factor
    call read_delayed_config
    if (log_info) then
       log_msg = 'poloidal mode number, sheet current factor'
       call log_write
       do m = lbound(sheet_current_factor, 1), ubound(sheet_current_factor, 1)
          write (log_msg, '(i2, 1x, ' // cmplx_fmt // ')') m, sheet_current_factor(m)
          call log_write
       end do
    end if

    ! depends on flux variables
    call compute_j0phi
    call check_curr0

    if (nonres) then
       call compute_Bn_nonres
    else
       call read_Bn(Bn_vac_file)
    end if
    Bnflux_vac = Bnflux
    Bnphi_vac = Bnphi
    call write_vector_dof(Bnflux_vac, Bnphi_vac, Bn_vacout_file)
    call write_vector_plot(Bnflux_vac, Bnphi_vac, &
         decorate_filename(Bn_vacout_file, 'plot_', ''))
    log_msg = 'magdif initialized'
    if (log_info) call log_write
  end subroutine magdif_init

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    if (allocated(efit)) deallocate(efit)
    if (allocated(fluxvar)) deallocate(fluxvar)
    if (allocated(fs)) deallocate(fs)
    if (allocated(fs_half)) deallocate(fs_half)
    if (allocated(m_res)) deallocate(m_res)
    if (allocated(sheet_current_factor)) deallocate(sheet_current_factor)
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
    if (allocated(kp_max)) deallocate(kp_max)
    if (allocated(kt_max)) deallocate(kt_max)
    if (allocated(kp_low)) deallocate(kp_low)
    if (allocated(kt_low)) deallocate(kt_low)
    log_msg = 'magdif cleanup finished'
    if (log_info) call log_write
    call log_close
  end subroutine magdif_cleanup

  subroutine magdif_single
    call compute_presn     ! compute pressure based on previous perturbation field
    call compute_currn     ! compute currents based on previous perturbation field
    call compute_Bn        ! use field code to generate new field from currents
    call read_Bn(Bn_file)  ! read new bnflux from field code
    call write_vector_plot(Bnflux, Bnphi, decorate_filename(Bn_file, 'plot_', ''))
  end subroutine magdif_single

  subroutine magdif_iterated
    use arnoldi_mod

    logical :: preconditioned
    integer :: kiter, ndim, i, j, info
    complex(dp) :: Bnflux_diff(ntri, 3)
    complex(dp) :: Bnphi_diff(ntri)
    complex(dp) :: Bn(ntri * 4), Bn_prev(ntri * 4)
    complex(dp) :: eigvals(nritz)
    complex(dp), allocatable :: Lr(:,:), Yr(:,:)
    integer, allocatable :: ipiv(:)
    character(len = 4) :: postfix
    character(len = *), parameter :: postfix_fmt = "('_', i0.3)"

    ! system dimension N ! - no need to include extended mesh
    ndim = ntri * 4 ! kt_low(nflux+1) * 4
    preconditioned = runmode_precon == runmode
    if (preconditioned) then
       ! calculate eigenvectors
       ieigen = 1
       call arnoldi(ndim, nritz, eigvals, next_iteration_arnoldi)
       if (log_info) then
          write (log_msg, '("Arnoldi method yields ", i0, " Ritz eigenvalues > ", f0.2)') &
               ngrow, tol
          call log_write
          do i = 1, ngrow
             write (log_msg, '("lambda ", i0, ": ", ' // cmplx_fmt // ')') i, eigvals(i)
             call log_write
          end do
       end if
       if (ngrow > 0) then
          do i = 1, min(ngrow, max_eig_out)
             write (postfix, postfix_fmt) i
             call unpack2(Bnflux, Bnphi, eigvecs(:, i))
             call write_vector_dof(Bnflux, Bnphi, &
                  decorate_filename(eigvec_file, '', postfix))
             call write_vector_plot(Bnflux, Bnphi, &
                  decorate_filename(eigvec_file, 'plot_', postfix))
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
             log_msg = 'Successfully inverted matrix for preconditioner'
             if (log_info) call log_write
          else
             write (log_msg, '("Matrix inversion for preconditioner failed: ' // &
                  'zgesv returns error ", i0)') info
             if (log_err) call log_write
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

    call write_vector_dof(Bnflux_vac, Bnphi_vac, Bn_diff_file)
    call compute_L2int
    call pack2(Bnflux_vac, Bnphi_vac, Bn_prev)
    if (preconditioned) then
       Bn_prev = Bn_prev - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
            matmul(transpose(conjg(eigvecs(:, 1:ngrow))), Bn_prev)))
    end if
    do kiter = 0, niter-1
       write (log_msg, '("Iteration ", i2, " of ", i2)') kiter, niter - 1
       if (log_info) call log_write
       write (postfix, postfix_fmt) kiter

       call next_iteration(ndim, Bn_prev, Bn)
       if (preconditioned) then
          Bn = Bn - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
               matmul(transpose(conjg(eigvecs(:, 1:ngrow))), Bn - Bn_prev)))
          call unpack2(Bnflux, Bnphi, Bn)
          call check_redundant_edges(Bnflux, .false., 'B_n')
          call check_div_free(Bnflux, Bnphi, n, rel_err_Bn, 'B_n')
       end if

       call unpack2(Bnflux_diff, Bnphi_diff, Bn - Bn_prev)
       call write_vector_dof(Bnflux_diff, Bnphi_diff, Bn_diff_file)
       call compute_L2int
       call write_vector_dof(Bnflux, Bnphi, &
            decorate_filename(Bn_file, '', postfix))
       if (kiter <= 1) then
          call write_vector_plot(Bnflux_diff, Bnphi_diff, &
               decorate_filename(Bn_diff_file, 'plot_', postfix))
          call write_vector_plot(Bnflux, Bnphi, &
               decorate_filename(Bn_file, 'plot_', postfix))
          call write_vector_plot(jnflux, jnphi, &
               decorate_filename(currn_file, 'plot_', postfix))
          call write_scalar_dof(presn, decorate_filename(presn_file, '', postfix))
       end if

       call pack2(Bnflux, Bnphi, Bn_prev)
    end do
    call write_vector_dof(Bnflux, Bnphi, Bn_file)
    call write_vector_plot(Bnflux, Bnphi, decorate_filename(Bn_file, 'plot_', ''))
    if (kilca_scale_factor /= 0) then
       call write_kilca_modes(Bnflux, Bnphi, kilca_pol_mode_file)
       call write_kilca_modes(Bnflux_vac, Bnphi_vac, &
            decorate_filename(kilca_pol_mode_file, '', '_vac'))
       call write_kilca_modes(Bnflux_vac - Bnflux, Bnphi_vac - Bnphi, &
            decorate_filename(kilca_pol_mode_file, '', '_plas'))
    end if

    if (allocated(Lr)) deallocate(Lr)

  contains

    pure subroutine pack2(pol_flux, tor_comp, packed)
      complex(dp), intent(in) :: pol_flux(ntri, 3), tor_comp(ntri)
      complex(dp), intent(out) :: packed(ndim)
      packed(1:ntri*3) = reshape(pol_flux, (/ ntri * 3 /))
      packed(ntri*3+1:ndim) = tor_comp
    end subroutine pack2

    pure subroutine unpack2(pol_flux, tor_comp, packed)
      complex(dp), intent(out) :: pol_flux(ntri, 3), tor_comp(ntri)
      complex(dp), intent(in) :: packed(ndim)
      pol_flux = reshape(packed(1:ntri*3), (/ ntri, 3 /))
      tor_comp = packed(ntri*3+1:ndim)
    end subroutine unpack2

    subroutine next_iteration(n, xold, xnew)
      ! computes B_(n+1) = K*B_n + B_vac ... different from kin2d.f90
      integer, intent(in) :: n
      complex(dp), intent(in) :: xold(n)
      complex(dp), intent(out) :: xnew(n)
      call unpack2(Bnflux, Bnphi, xold)
      call magdif_single
      Bnflux = Bnflux + Bnflux_vac
      Bnphi = Bnphi + Bnphi_vac
      call pack2(Bnflux, Bnphi, xnew)
    end subroutine next_iteration

    subroutine next_iteration_arnoldi(n, xold, xnew)
      ! computes B_(n+1) = K*(B_n + B_vac) ... as in kin2d.f90
      integer, intent(in) :: n
      complex(dp), intent(in) :: xold(n)
      complex(dp), intent(out) :: xnew(n)
      call unpack2(Bnflux, Bnphi, xold)
      Bnflux = Bnflux + Bnflux_vac
      Bnphi = Bnphi + Bnphi_vac
      call magdif_single
      call pack2(Bnflux, Bnphi, xnew)
    end subroutine next_iteration_arnoldi
  end subroutine magdif_iterated

  !> Allocates and initializes #kp_low, #kp_max, #kt_low and #kt_max based on the values
  !> of #magdif_config::nflux and #magdif_config::nkpol. Deallocation is done in
  !> magdif_cleanup().
  subroutine init_indices
    integer :: kf

    allocate (kp_max(nflux+1))
    allocate (kt_max(nflux+1))
    allocate (kp_low(nflux+1))
    allocate (kt_low(nflux+1))

    kp_max = nkpol
    kt_max = 2 * nkpol
    kt_max(1) = nkpol

    kp_low(1) = 1
    do kf = 2, nflux+1
       kp_low(kf) = kp_low(kf-1) + kp_max(kf-1)
    end do
    kt_low(1) = 0
    do kf = 2, nflux+1
       kt_low(kf) = kt_low(kf-1) + kt_max(kf-1)
    end do

    write (log_msg, '("Number of points up to LCFS: ", i0, ' // nl_fmt // &
         ', "Number of triangles up to LCFS: ", i0)') kp_low(nflux+1), kt_low(nflux+1)
    if (log_info) call log_write
  end subroutine init_indices

  !> Reads mesh points and triangles.
  !>
  !> #mesh_mod::npoint and #mesh_mod::mesh_point are read directly from an unformatted
  !> #magdif_config::point_file, while #mesh_mod::ntri and #mesh_mod::mesh_element are
  !> read directly from an unformatted #magdif_config::tri_file. #presn, #bnflux, #bnphi,
  !> #bnflux_vac, #bnphi_vac, #j0phi, #jnphi and #jnflux are allocated and initialized to
  !> zero. Deallocation is done in magdif_cleanup().
  subroutine read_mesh
    use mesh_mod, only: bphicovar
    integer :: fid

    open(newunit = fid, file = point_file, form = 'unformatted', status = 'old')
    read (fid) npoint
    write (log_msg, '("npoint = ", i0)') npoint
    if (log_info) call log_write
    allocate(mesh_point(npoint))
    read (fid) mesh_point
    close(fid)

    open(newunit = fid, file = tri_file, form = 'unformatted', status = 'old')
    read(fid) ntri
    write (log_msg, '("ntri   = ", i0)') ntri
    if (log_info) call log_write
    allocate(mesh_element(ntri))
    read (fid) mesh_element
    read (fid) bphicovar
    close(fid)

    allocate(mesh_element_rmp(ntri))

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

  subroutine cache_mesh_data
    integer :: kf, kt, ktri
    type(triangle) :: elem
    type(triangle_rmp) :: tri

    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          elem = mesh_element(ktri)
          tri%area = 0.5d0 * elem%det_3
          call get_labeled_edges(elem, tri%li, tri%lo, tri%lf, tri%ei, tri%eo, tri%ef, &
               tri%orient)
          call ring_centered_avg_coord(elem, tri%R_Omega, tri%Z_Omega)
          mesh_element_rmp(ktri) = tri
       end do
    end do
  end subroutine cache_mesh_data

  !> Computes #bnflux and #bnphi from #jnflux and #jnphi via an external program. No data
  !> is read yet; this is done by read_bn().
  !>
  !> Currently, this subroutine Calls FreeFem++ via shell script maxwell.sh in the current
  !> directory and handles the exit code. For further information see maxwell.sh and the
  !> script called therein.
  subroutine compute_Bn
    integer :: stat = 0, dummy = 0
    character(len = 1024) :: fem_cmd
    write (fem_cmd, '("./maxwell.sh -n ", i0, " -N ", i0, " -J ", a, " -B ", a)') &
         n, kt_low(nflux + 1), trim(currn_file), trim(Bn_file)
    call execute_command_line(fem_cmd, exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       write (log_msg, '("FreeFem++ failed with exit code ", i0)') stat
       if (log_err) call log_write
       error stop
    end if
  end subroutine compute_Bn

  subroutine compute_L2int
    integer :: stat = 0, dummy = 0
    character(len = 1024) :: L2int_cmd
    write (L2int_cmd, '("./L2int.sh -N ", i0, " -B ", a, " -C ", a)') &
         kt_low(nflux + 1), trim(Bn_diff_file), trim(conv_file)
    call execute_command_line(L2int_cmd, exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       write (log_msg, '("FreeFem++ failed with exit code ", i0)') stat
       if (log_err) call log_write
       error stop
    end if
  end subroutine compute_L2int

  !> Checks if divergence-freeness of the given vector field is fulfilled on each
  !> triangle, otherwise halts the program.
  !>
  !> @param pol_flux poloidal flux components, e.g. #jnflux - first index refers to the
  !> triangle, second index refers to the edge on that triangle, as per
  !> #mesh_mod::mesh_element
  !> @param tor_comp toroidal field components, e.g. #jnphi - index refers to the triangle
  !> @param rel_err relative error threshold, i.e. maximum acceptable value for divergence
  !> after division by absolute flux through the triangle
  !> @param field_name name given to the vector field in the error message
  !>
  !> This subroutine calculates the divergence via the divergence theorem, i.e. by adding
  !> up the fluxes of the vector field through triangle edges. If this sum, divided by the
  !> absolute flux, is higher than \p rel_err on any triangle, it halts the program.
  subroutine check_div_free(pol_flux, tor_comp, n, rel_err, field_name)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: rel_err
    character(len = *), intent(in) :: field_name

    integer :: kt
    real(dp) :: div, abs_flux

    do kt = 1, ntri
       abs_flux = sum(abs(pol_flux(kt,:))) + abs(imun * n * tor_comp(kt) * &
            mesh_element(kt)%det_3 * 0.5d0)
       if (abs_flux > 0d0) then
          div = abs((sum(pol_flux(kt,:)) + imun * n * tor_comp(kt) * &
               mesh_element(kt)%det_3 * 0.5d0)) / abs_flux
          if (div > rel_err) then
             write (log_msg, '("divergence of ", a, ' // &
                  '" above threshold in triangle ", i0, ": ", es23.16)') &
                  trim(field_name), kt, div
              if (log_err) call log_write
              error stop
          end if
       end if
    end do
  end subroutine check_div_free

  subroutine check_redundant_edges(pol_quant, same_sign, name)
    complex(dp), intent(in) :: pol_quant(:,:)
    logical, intent(in) :: same_sign
    character(len = *), intent(in) :: name
    integer :: ktri, ktri_adj, ke, ke_adj
    logical :: checked(ntri, 3), inconsistent
    real(dp), parameter :: eps = epsilon(1d0)

    checked = .false.
    do ktri = 1, kt_low(nflux+1)
       do ke = 1, 3
          if (ktri > kt_low(nflux) .and. ke == mesh_element_rmp(ktri)%ef) cycle
          if (.not. checked(ktri, ke)) then
             checked(ktri, ke) = .true.
             ktri_adj = mesh_element(ktri)%neighbour(ke)
             ke_adj = mesh_element(ktri)%neighbour_edge(ke)
             checked(ktri_adj, ke_adj) = .true.
             inconsistent = .false.
             if (real(pol_quant(ktri, ke)) == 0d0) then
                inconsistent = inconsistent .or. real(pol_quant(ktri_adj, ke_adj)) /= 0d0
             else
                if (same_sign) then
                   inconsistent = inconsistent .or. real(pol_quant(ktri_adj, ke_adj)) / &
                        real(pol_quant(ktri, ke)) - 1d0 > eps
                else
                   inconsistent = inconsistent .or. real(pol_quant(ktri_adj, ke_adj)) / &
                        (-real(pol_quant(ktri, ke))) - 1d0 > eps
                end if
             end if
             if (aimag(pol_quant(ktri, ke)) == 0d0) then
                inconsistent = inconsistent .or. aimag(pol_quant(ktri_adj, ke_adj)) /= 0d0
             else
                if (same_sign) then
                   inconsistent = inconsistent .or. aimag(pol_quant(ktri_adj, ke_adj)) / &
                        aimag(pol_quant(ktri, ke)) - 1d0 > eps
                else
                   inconsistent = inconsistent .or. aimag(pol_quant(ktri_adj, ke_adj)) / &
                        (-aimag(pol_quant(ktri, ke))) - 1d0 > eps
                end if
             end if
             if (inconsistent) then
                write (log_msg, '("inconsistent redundant edges: ", ' // &
                     'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ', ", ", ' // &
                     'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ')') &
                     trim(name), ktri, ke, pol_quant(ktri, ke), &
                     trim(name), ktri_adj, ke_adj, pol_quant(ktri_adj, ke_adj)
                if (log_err) call log_write
                error stop
             end if
          end if
       end do
    end do
  end subroutine check_redundant_edges

  !> Reads fluxes of perturbation field and checks divergence-freeness.
  !>
  !> @param filename name of the formatted file containing the data
  !>
  !> Line numbers in the given file correspond to the global triangle index of
  !> #mesh_mod::mesh_element. The content of each line is read into #bnflux and #bnphi
  !> with numbering of edges in #bnflux as in #mesh_mod::mesh_element and the imaginary
  !> part of each value immediately following the real part.
  subroutine read_Bn(filename)
    character(len = 1024) :: filename
    integer :: kt, fid
    real(dp) :: dummy8(8)

    open(newunit = fid, file = filename, recl = longlines, status = 'old')
    do kt = 1, ntri
       read (fid, *) dummy8
       Bnflux(kt,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(kt,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(kt,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(kt) = cmplx(dummy8(7), dummy8(8), dp) / mesh_element(kt)%det_3 * 2d0
    end do
    close(fid)

    call check_div_free(Bnflux, Bnphi, n, rel_err_Bn, 'B_n')
    call check_redundant_edges(Bnflux, .false., 'B_n')
  end subroutine read_Bn

  !> Allocates and computes the safety factor #q and #m_res.
  !>
  !> Also allocates #magdif_config::sheet_current_factor, to be read in via
  !> magdif_config::read_delayed_config() in magdif_init(). All deallocation is done in
  !> magdif_cleanup().
  subroutine compute_safety_factor
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri
    integer :: m_res_min, m_res_max, m
    real(dp), dimension(nflux) :: abs_err

    select case (q_prof)
    case (q_prof_flux)
       fs_half%q = 0d0
       do kf = 1, nflux
          do kt = 1, kt_max(kf)
             ktri = kt_low(kf) + kt
             tri = mesh_element_rmp(ktri)
             fs_half%q(kf) = fs_half%q(kf) + B0phi_Omega(ktri) * tri%area
          end do
          fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
       end do
       ! use linear interpolation for full-grid steps for now
       fs%q(1:nflux-1) = 0.5d0 * (fs_half%q(1:nflux-1) + fs_half%q(2:nflux))
       ! linear extrapolation for values at separatrix and magnetic axis
       fs%q(nflux) = fs%q(nflux-1) + (fs%q(nflux-1) - fs%q(nflux-2))
       fs%q(0) = fs%q(1) - (fs%q(2) - fs%q(1))
    case (q_prof_efit)
       do kf = 0, nflux
          ! use positive q
          fs%q(kf) = abs(fluxvar%interp(efit%qpsi, kf, .false.))
       end do
       do kf = 1, nflux
          ! use positive q
          fs_half%q(kf) = abs(fluxvar%interp(efit%qpsi, kf, .true.))
       end do
    case default
       write (log_msg, '("unknown q profile selection: ", i0)') q_prof
       if (log_err) call log_write
       error stop
    end select

    allocate(m_res(nflux))
    m_res = 0
    if (kilca_scale_factor /= 0) then
       m_res_min = max(ceiling(minval(fs_half%q) * n), n / kilca_scale_factor + 1)
    else
       m_res_min = max(ceiling(minval(fs_half%q) * n), n + 1)
    end if
    m_res_max = floor(maxval(fs_half%q) * n)
    do m = m_res_max, m_res_min, -1
       abs_err = [(abs(fs_half%q(kf) - dble(m) / dble(n)), kf = 1, nflux)]
       m_res(minloc(abs_err, 1)) = m
    end do
    write (log_msg, '("resonant m: ", i0, " .. ", i0)') m_res_min, m_res_max
    if (log_debug) call log_write
    allocate(sheet_current_factor(m_res_min:m_res_max))
    sheet_current_factor = (0d0, 0d0)
  end subroutine compute_safety_factor

  !> Initializes quantities that are constant on each flux surface. The safety factor #q
  !> is initialized separately in compute_safety_factor().
  !>
  !> #psimin, #psimax and #psi are initialized from #mesh_mod::mesh_point::psi_pol.
  !> #pres0, #dpres0_dpsi, #dens and #temp are given an assumed profile based on #psimin,
  !> #psimax, #psi, #magdif_conf::di0, #magdif_conf::d_min, #magdif_conf::ti0 and
  !> #magdif_conf::t_min.
  subroutine init_flux_variables
    integer :: kf

    allocate(fs)
    call fs%init(nflux, .false.)

    ! magnetic axis at k == 0 is not counted as flux surface
    fs%psi(0) = mesh_point(1)%psi_pol
    do kf = 1, nflux
       ! average over the loop to smooth out numerical errors
       fs%psi(kf) = sum(mesh_point((kp_low(kf) + 1):kp_low(kf+1))%psi_pol) / kp_max(kf)
    end do

    allocate(fs_half)
    call fs_half%init(nflux, .true.)

    ! use linear interpolation for half-grid steps for now
    fs_half%psi = 0.5d0 * (fs%psi(0:nflux-1) + fs%psi(1:nflux))

    allocate(fluxvar)
    call fluxvar%init(4, efit%nw, nflux, fs%psi, fs_half%psi)

    call compute_pres_prof
    call compute_safety_factor
    do kf = 0, nflux
       fs%F(kf) = fluxvar%interp(efit%fpol, kf, .false.)
       fs%FdF_dpsi(kf) = fluxvar%interp(efit%ffprim, kf, .false.)
    end do
    do kf = 1, nflux
       fs_half%F(kf) = fluxvar%interp(efit%fpol, kf, .true.)
       fs_half%FdF_dpsi(kf) = fluxvar%interp(efit%ffprim, kf, .true.)
    end do
    call write_fluxvar
  end subroutine init_flux_variables

  subroutine compute_pres_prof
    use for_macrostep, only: t_min, d_min  ! PRELOAD/SRC/orbit_mod.f90
    use constants, only: ev2erg            ! PRELOAD/SRC/orbit_mod.f90
    integer :: kf
    real(dp) :: ddens_dpsi, dtemp_dpsi, psimin, psimax

    !> Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
    real(dp) :: dens(0:nflux)

    !> Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
    real(dp) :: temp(0:nflux)

    dens = 0d0
    temp = 0d0
    psimin = minval(mesh_point%psi_pol)
    psimax = maxval(mesh_point%psi_pol)
    select case (pres_prof)
    case (pres_prof_eps)
       ddens_dpsi = di0 / psimax
       dtemp_dpsi = ti0 / psimax
       dens = (fs%psi - psimin) / psimax * di0 + d_min
       temp = (fs%psi - psimin) / psimax * ti0 + t_min
       write (log_msg, '("temp@axis: ", es23.16, ", dens@axis: ", es23.16)') temp(0), dens(0)
       if (log_info) call log_write
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psimin) / psimax * di0 + d_min
       temp(1:) = (fs_half%psi - psimin) / psimax * ti0 + t_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_par)
       ddens_dpsi = (di0 - d_min) / (psimax - psimin)
       dtemp_dpsi = (ti0 - t_min) / (psimax - psimin)
       dens = (fs%psi - psimin) / (psimax - psimin) * (di0 - d_min) + d_min
       temp = (fs%psi - psimin) / (psimax - psimin) * (ti0 - t_min) + t_min
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psimin) / (psimax - psimin) * (di0 - d_min) + d_min
       temp(1:) = (fs_half%psi - psimin) / (psimax - psimin) * (ti0 - t_min) + t_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_efit)
       do kf = 0, nflux
          fs%p(kf) = fluxvar%interp(efit%pres, kf, .false.)
          fs%dp_dpsi(kf) = fluxvar%interp(efit%pprime, kf, .false.)
       end do
       do kf = 1, nflux
          fs_half%p(kf) = fluxvar%interp(efit%pres, kf, .true.)
          fs_half%dp_dpsi(kf) = fluxvar%interp(efit%pprime, kf, .true.)
       end do
    case default
       write (log_msg, '("unknown pressure profile selection", i0)') pres_prof
       if (log_err) call log_write
       error stop
    end select
  end subroutine compute_pres_prof


  subroutine cache_equilibrium_field
    real(dp) :: r, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kf, kt, ktri, ke, fid
    type(triangle) :: elem
    type(knot) :: base, tip
    real(dp) :: n_r, n_z

    open(newunit = fid, file = 'plot_B0.dat', recl = longlines, status = 'replace')
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          elem = mesh_element(ktri)
          do ke = 1, 3
             base = mesh_point(elem%i_knot(ke))
             tip = mesh_point(elem%i_knot(mod(ke, 3) + 1))
             r = (base%rcoord + tip%rcoord) * 0.5d0
             z = (base%zcoord + tip%zcoord) * 0.5d0
             call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
                  dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
             B0r(ktri, ke) = Br
             B0phi(ktri, ke) = Bp
             B0z(ktri, ke) = Bz
             n_r = tip%zcoord - base%zcoord
             n_z = base%rcoord - tip%rcoord
             B0flux(ktri, ke) = r * (Br * n_r + Bz * n_z)
          end do
          call ring_centered_avg_coord(elem, r, z)
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          B0r_Omega(ktri) = Br
          B0phi_Omega(ktri) = Bp
          B0z_Omega(ktri) = Bz
          write (fid, '(5(1x, es23.16))') r, z, Br, Bz, Bp
       end do
    end do
    do kt = kt_low(nflux+1) + 1, ntri
       write (fid, '(5(1x, es23.16))') R0, 0d0, 0d0, 0d0, 0d0
    end do
    close(fid)
  end subroutine cache_equilibrium_field

  !> Computes equilibrium current density #j0phi from given equilibrium magnetic field and
  !> assumed equilibrium pressure #pres0.
  !>
  !> This step is necessary because equilibrium pressure is not given experimentally as is
  !> \f$ \vec{B}_{0} \f$; arbitrary values are assumed. Consistency of MHD equilibrium is
  !> necessary in the derivation, while Ampere's equation is not used.
  subroutine compute_j0phi
    integer :: kf, kt, ktri, fid
    real(dp) :: r, z
    real(dp) :: Btor2
    real(dp), dimension(nflux) :: B2avg, B2avg_half
    real(dp) :: plot_j0phi
    type(triangle_rmp) :: tri

    open(newunit = fid, file = j0phi_file, recl = longlines, status = 'replace')
    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0r(ktri, tri%ef) ** 2 + &
               B0phi(ktri, tri%ef) ** 2 + B0z(ktri, tri%ef) ** 2
          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0r(ktri, tri%ei) ** 2 + &
               B0phi(ktri, tri%ei) ** 2 + B0z(ktri, tri%ei) ** 2
       end do
       B2avg(kf) = B2avg(kf) / kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / kt_max(kf)

       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)

          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             if (tri%orient) then
                j0phi(ktri, tri%ef) = clight * (fs%dp_dpsi(kf) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf) / r)
             else
                j0phi(ktri, tri%ef) = clight * (fs%dp_dpsi(kf-1) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf-1) / r)
             end if
          case (curr_prof_rot)
             z = sum(mesh_point(tri%lf(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%ef) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%ef) ** 2
             if (.not. tri%orient) then
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf-1) * (1d0 - &
                     Btor2 / B2avg(kf-1))
             else
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf) * (1d0 - &
                     Btor2 / B2avg(kf))
             end if
          end select

          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             j0phi(ktri, tri%ei) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(tri%li(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%ei) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%ei) ** 2
             j0phi(ktri, tri%ei) = clight * r * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             j0phi(ktri, tri%eo) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(tri%lo(:))%zcoord) * 0.5d0
             j0phi(ktri, tri%eo) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(ktri, tri%eo) ** 2
             j0phi(ktri, tri%eo) = clight * r * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          select case (curr_prof)
          case (curr_prof_efit)
             plot_j0phi = clight * (fs_half%dp_dpsi(kf) * tri%r_Omega + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / tri%r_Omega)
          case (curr_prof_rot)
             plot_j0phi = j0phi_ampere(tri%r_Omega, tri%z_Omega)
          case (curr_prof_ps)
             Btor2 = B0phi_Omega(ktri) ** 2
             plot_j0phi = clight * tri%r_Omega * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          write (fid, '(4(1x, es23.16))') &
               j0phi(ktri, 1), j0phi(ktri, 2), j0phi(ktri, 3), plot_j0phi
       end do
    end do
    do kt = kt_low(nflux+1) + 1, ntri
       write (fid, '(4(1x, es23.16))') j0phi(kt, 1), j0phi(kt, 2), j0phi(kt, 3), 0d0
    end do
    close(fid)

    call check_redundant_edges(cmplx(j0phi, 0d0, dp), .true., 'j0phi')

  contains
    function j0phi_ampere(r, z) result (rotB_phi)
      real(dp), intent(in) :: r, z
      real(dp) :: rotB_phi
      real(dp) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
           dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
      call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
           dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
      rotB_phi = 0.25d0 / pi * clight * (dBrdZ - dBzdR)
    end function j0phi_ampere

  end subroutine compute_j0phi


  subroutine check_curr0
    integer :: kf, kt, ktri, fid_amp, fid_gs, fid_prof
    real(dp) :: r, z, n_r, n_z, cmp_amp, cmp_gs, Jr, Jp, Jz
    real(dp) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    type(triangle_rmp) :: tri

    open(newunit = fid_prof, file = 'cmp_prof.dat', recl = longlines, status = 'replace')
    open(newunit = fid_amp, file = 'j0_amp.dat', recl = longlines, status = 'replace')
    open(newunit = fid_gs, file = 'j0_gs.dat', recl = longlines, status = 'replace')
    do kf = 1, nflux
       cmp_amp = 0d0
       cmp_gs = 0d0
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          call radially_outward_normal(ktri, n_r, n_z)
          r = tri%r_Omega
          z = tri%z_Omega
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          dBpdR = dBpdR + fs_half%FdF_dpsi(kf) / fs_half%F(kf) * Bz
          dBpdZ = -fs_half%FdF_dpsi(kf) / fs_half%F(kf) * Br
          Jr = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * &
               B0r_Omega(ktri)
          Jz = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * &
               B0z_Omega(ktri)
          Jp = clight * (fs_half%dp_dpsi(kf) * r + 0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          write (fid_gs, '(3(1x, es23.16))') Jr, Jp, Jz
          cmp_gs = cmp_gs + ((Jp * Bz - Jz * Bp) * n_r + (Jr * Bp - Jp * Br) * n_z) / &
               (n_r**2 + n_z**2) / (fs%psi(kf) - fs%psi(kf-1)) * 2d0 * tri%area / clight
          Jr = 0.25d0 / pi * clight * (-dBpdZ)
          Jp = 0.25d0 / pi * clight * (dBrdZ - dBzdR)
          Jz = 0.25d0 / pi * clight * (dBpdR + Bp / r)
          write (fid_amp, '(3(1x, es23.16))') Jr, Jp, Jz
          cmp_amp = cmp_amp + ((Jp * Bz - Jz * Bp) * n_r + (Jr * Bp - Jp * Br) * n_z) / &
               (n_r**2 + n_z**2) / (fs%psi(kf) - fs%psi(kf-1)) * 2d0 * tri%area / clight
       end do
       cmp_amp = cmp_amp / kt_max(kf)
       cmp_gs = cmp_gs / kt_max(kf)
       write (fid_prof, '(2(1x, es23.16))') cmp_amp, cmp_gs
    end do
    close(fid_prof)
    close(fid_amp)
    close(fid_gs)
  end subroutine check_curr0

  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  subroutine compute_presn
    real(dp) :: r, Deltapsi
    real(dp) :: lr, lz  ! edge vector components
    complex(dp) :: Bnpsi
    complex(dp), dimension(nkpol) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, ktri, kp
    integer :: nz
    integer, dimension(2*nkpol) :: irow, icol
    complex(dp), dimension(2*nkpol) :: aval
    type(triangle_rmp) :: tri
    type(knot) :: base, tip
    integer :: common_tri(2)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux
       if (kf < nflux) then
          Deltapsi = fs%psi(kf+1) - fs%psi(kf-1)
       else
          ! use linear interpolation for values outside LCFS
          Deltapsi = 2d0 * (fs%psi(kf) - fs%psi(kf-1))
       end if
       inner: do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          if (.not. tri%orient) cycle inner
          if (kf < nflux) then
             common_tri = [ktri, mesh_element(ktri)%neighbour(tri%ef)]
          else
             ! triangles outside LCFS are comparatively distorted
             ! use linear extrapolation instead, i.e. a triangle of the same size
             common_tri = [ktri, ktri]
          end if
          ! use midpoint of edge f
          base = mesh_point(tri%lf(1))
          tip = mesh_point(tri%lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord

          Bnpsi = Bnflux(ktri, tri%ef) / r * Deltapsi / &
               sum(mesh_element_rmp(common_tri(:))%area) * 0.5d0

          kp = tri%lf(1) - kp_low(kf)

          x(kp) = -fs%dp_dpsi(kf) * Bnpsi

          a(kp) = (B0r(ktri, tri%ef) * lr + B0z(ktri, tri%ef) * lz) / (lr ** 2 + lz ** 2)

          b(kp) = imun * (n + imun * damp) * B0phi(ktri, tri%ef) / r
       end do inner

       ! solve linear system
       d = -a + b * 0.5d0
       du = a + b * 0.5d0
       call assemble_sparse(nkpol, d, du, nz, irow, icol, aval)
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(nkpol, nkpol, nz, irow, icol, aval, x)
       call sparse_matmul(nkpol, nkpol, irow, icol, aval, x, resid)
       resid = resid - inhom
       where (inhom /= 0d0)
          rel_err = abs(resid) / abs(inhom)
       elsewhere
          rel_err = 0d0
       end where
       max_rel_err = max(max_rel_err, maxval(rel_err))
       avg_rel_err = avg_rel_err + sum(rel_err)

       if (kf == 1) then ! first point on axis - average over enclosing flux surface
          presn(1) = sum(x) / size(x)
       end if
       do kp = 1, kp_max(kf)
          presn(kp_low(kf) + kp) = x(kp)
       end do
    end do

    avg_rel_err = avg_rel_err / sum(kp_max(1:nflux))
    write (log_msg, '("compute_presn: diagonalization max_rel_err = ", ' // &
         'es23.16, ", avg_rel_err = ", es23.16)') max_rel_err, avg_rel_err
    if (log_debug) call log_write
    if (allocated(resid)) deallocate(resid)

    call write_scalar_dof(presn, presn_file)
  end subroutine compute_presn

  !> Map edge symbols to integer indices.
  !>
  !> @param elem the triangle for which indices are to be obtained
  !> @param li knot indices for base (1) and tip (2) of edge i
  !> @param lo knot indices for base (1) and tip (2) of edge o
  !> @param lf knot indices for base (1) and tip (2) of edge f
  !> @param ei index of edge i, e.g. for #bnflux
  !> @param eo index of edge o, e.g. for #bnflux
  !> @param ef index of edge f, e.g. for #bnflux
  !> @param orient true if edge f lies on the outer flux surface, false otherwise
  !>
  !> It is assumed that the knots are globally numbered in ascending order starting from
  !> the magnetic axis and going counter-clockwise around each flux surface. Furthermore
  !> it is assumed that edge 1 goes from knot 1 to knot 2, edge 2 from knot 2 to knot 3
  !> and edge 3 from knot 3 to knot 1. It is not assumed that knots are orderd
  !> counter-clockwise locally.

  subroutine get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    logical, intent(out) :: orient
    integer, dimension(3) :: i_knot_diff
    integer :: knot_i, knot_o, knot_f
    integer :: i1, i2
    logical :: closing_loop

    log_msg = 'cannot find correct label for triangle edges'

    ! initialize to suppress compiler warnings
    i1 = 0
    i2 = 0
    knot_f = elem%knot_h
    select case (knot_f)
    case (1)
       i1 = 2
       i2 = 3
    case (2)
       i1 = 3
       i2 = 1
    case (3)
       i1 = 1
       i2 = 2
    end select
    if (elem%i_knot(i1) == elem%i_knot(i2)) then
       if (log_err) call log_write
       error stop
    end if
    ! last triangle in strip if indices not next to each other
    closing_loop = abs(elem%i_knot(i1) - elem%i_knot(i2)) /= 1
    i_knot_diff = elem%i_knot - elem%i_knot(knot_f)
    if (all(i_knot_diff >= 0)) then
       ! knot_f lies on inner surface
       orient = .true.
       if ((elem%i_knot(i1) < elem%i_knot(i2)) .neqv. closing_loop) then
          ! i1 is next after knot_f counter-clockwise
          knot_o = i1
          knot_i = i2
       else
          ! i2 is next after knot_f counter-clockwise
          knot_o = i2
          knot_i = i1
       end if
       ei = knot_f
       eo = knot_i
       ef = knot_o
       li = (/ elem%i_knot(knot_f), elem%i_knot(knot_o) /)
       lo = (/ elem%i_knot(knot_i), elem%i_knot(knot_f) /)
       lf = (/ elem%i_knot(knot_o), elem%i_knot(knot_i) /)
    else if (all(i_knot_diff <= 0)) then
       ! knot_f lies on outer surface
       orient = .false.
       if ((elem%i_knot(i1) > elem%i_knot(i2)) .neqv. closing_loop) then
          ! i1 is next after knot_f counter-clockwise
          knot_i = i1
          knot_o = i2
       else
          ! i2 is next after knot_f counter-clockwise
          knot_i = i2
          knot_o = i1
       end if
       ei = knot_o
       eo = knot_f
       ef = knot_i
       li = (/ elem%i_knot(knot_o), elem%i_knot(knot_f) /)
       lo = (/ elem%i_knot(knot_f), elem%i_knot(knot_i) /)
       lf = (/ elem%i_knot(knot_i), elem%i_knot(knot_o) /)
    else
       if (log_err) call log_write
       error stop
    end if
  end subroutine get_labeled_edges


  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #magdif_conf::currn_file.
  subroutine compute_currn
    complex(dp), dimension(2 * nkpol) :: x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(2 * nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, ktri
    integer :: nz
    integer, dimension(4 * nkpol) :: irow, icol
    complex(dp), dimension(4 * nkpol) :: aval
    type(triangle) :: elem
    type(triangle_rmp) :: tri
    complex(dp) :: Bnphi_Gamma
    real(dp) :: r

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux ! loop through flux surfaces
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          elem = mesh_element(ktri)
          tri = mesh_element_rmp(ktri)
          ! first term on source side: flux through edge f
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          jnflux(ktri, tri%ef) = j0phi(ktri, tri%ef) / B0phi(ktri, tri%ef) * &
               Bnflux(ktri, tri%ef) + &
               clight * r / B0phi(ktri, tri%ef) * (presn(tri%lf(2)) - presn(tri%lf(1)))
          x(kt) = -jnflux(ktri, tri%ef)
          ! diagonal matrix element - edge i
          d(kt) = -1d0 - imun * (n + imun * damp) * tri%area * 0.5d0 * &
               B0phi(ktri, tri%ei) / B0flux(ktri, tri%ei)
          ! additional term from edge i on source side
          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(tri%ei)))
          x(kt) = x(kt) - imun * n * tri%area * 0.5d0 * (clight * r / B0flux(ktri, tri%ei) * &
               (Bnphi_Gamma / B0phi(ktri, tri%ei) * (fs%p(kf) - fs%p(kf-1)) - &
               (presn(tri%li(2)) - presn(tri%li(1)))) + j0phi(ktri, tri%ei) * &
               (Bnphi_Gamma / B0phi(ktri, tri%ei) - Bnflux(ktri, tri%ei) / B0flux(ktri, tri%ei)))
          ! superdiagonal matrix element - edge o
          du(kt) = 1d0 + imun * (n + imun * damp) * tri%area * 0.5d0 * &
               B0phi(ktri, tri%eo) / B0flux(ktri, tri%eo)
          ! additional term from edge o on source side
          r = sum(mesh_point(tri%lo(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(tri%eo)))
          x(kt) = x(kt) - imun * n * tri%area * 0.5d0 * (clight * r / B0flux(ktri, tri%eo) * &
               (Bnphi_Gamma / B0phi(ktri, tri%eo) * (fs%p(kf-1) - fs%p(kf)) - &
               (presn(tri%lo(2)) - presn(tri%lo(1)))) + j0phi(ktri, tri%eo) * &
               (Bnphi_Gamma / B0phi(ktri, tri%eo) - Bnflux(ktri, tri%eo) / B0flux(ktri, tri%eo)))
       end do
       call assemble_sparse(kt_max(kf), d(:kt_max(kf)), du(:kt_max(kf)), nz, &
            irow(:(2*kt_max(kf))), icol(:(2*kt_max(kf))), aval(:(2*kt_max(kf))))
       inhom = x  ! remember inhomogeneity before x is overwritten with the solution
       call sparse_solve(kt_max(kf), kt_max(kf), nz, irow(:nz), icol(:nz), aval(:nz), &
            x(:kt_max(kf)))
       call sparse_matmul(kt_max(kf), kt_max(kf), irow(:nz), icol(:nz), aval(:nz), &
            x(:kt_max(kf)), resid)
       resid = resid - inhom(:kt_max(kf))
       where (inhom(:kt_max(kf)) /= 0d0)
          rel_err(:kt_max(kf)) = abs(resid(:kt_max(kf))) / abs(inhom(:kt_max(kf)))
       elsewhere
          rel_err(:kt_max(kf)) = 0d0
       end where
       max_rel_err = max(max_rel_err, maxval(rel_err(:kt_max(kf))))
       avg_rel_err = avg_rel_err + sum(rel_err(:kt_max(kf)))
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          jnflux(ktri, tri%ei) = -x(kt)
          jnflux(ktri, tri%eo) = x(mod(kt, kt_max(kf)) + 1)
          jnphi(ktri) = sum(jnflux(ktri, :)) * imun / n / tri%area
       end do
    end do
    avg_rel_err = avg_rel_err / sum(kt_max(1:nflux))
    write (log_msg, '("compute_currn: diagonalization max_rel_err = ", ' // &
         'es23.16, ", avg_rel_err = ", es23.16)') max_rel_err, avg_rel_err
    if (log_debug) call log_write
    if (allocated(resid)) deallocate(resid)

    call add_sheet_current
    call check_redundant_edges(jnflux, .false., 'jnflux')
    call check_div_free(jnflux, jnphi, n, rel_err_currn, 'jnflux')
    call write_vector_dof(jnflux, jnphi, currn_file)
    call write_vector_plot(jnflux, jnphi, decorate_filename(currn_file, 'plot_', ''))
  end subroutine compute_currn

  subroutine add_sheet_current
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri, tri_adj
    complex(dp) :: presn_half

    do kf = 1, nflux
       if (m_res(kf) > 0) then
          if (sheet_current_factor(m_res(kf)) /= (0d0, 0d0)) then
             do kt = 1, kt_max(kf)
                ktri = kt_low(kf) + kt
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
                     sheet_current_factor(m_res(kf)) * B0flux(ktri, tri%ei) * presn_half
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
                     sheet_current_factor(m_res(kf)) * B0flux(ktri, tri%eo) * presn_half
                ! adjust toroidal current density
                jnphi(ktri) = sum(jnflux(ktri, :)) * imun / n / tri%area
             end do
          end if
       end if
    end do
  end subroutine add_sheet_current

  subroutine avg_flux_on_quad(pol_flux, tor_comp)
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)

    integer :: kf, kt, ktri1, ktri2
    complex(dp) :: tor_flux_avg, tor_flux_diff
    type(triangle_rmp) :: tri1, tri2

    do kf = 2, nflux
       do kt = 1, kt_max(kf), 2
          ktri1 = kt_low(kf) + kt
          ktri2 = kt_low(kf) + kt + 1
          tri1 = mesh_element_rmp(ktri1)
          tri2 = mesh_element_rmp(ktri2)
          tor_flux_avg = 0.5d0 * (tor_comp(ktri2) * tri2%area + &
               tor_comp(ktri1) * tri1%area)
          tor_flux_diff = 0.5d0 * (tor_comp(ktri2) * tri2%area - &
               tor_comp(ktri1) * tri1%area)
          tor_comp(ktri1) = tor_flux_avg / tri1%area
          pol_flux(ktri1, tri1%eo) = pol_flux(ktri1, tri1%eo) - imun * n * tor_flux_diff
          tor_comp(ktri2) = tor_flux_avg / tri2%area
          pol_flux(ktri2, tri2%ei) = pol_flux(ktri2, tri2%ei) + imun * n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  subroutine compute_Bn_nonres
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri
    integer :: common_tri(2)
    real(dp) :: Deltapsi
    complex(dp) :: Bnpsi
    real(dp) :: r

    do kf = 1, nflux ! loop through flux surfaces
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          if (tri%orient) then
             if (kf < nflux) then
                Deltapsi = fs%psi(kf+1) - fs%psi(kf-1)
             else
                ! use linear interpolation for values outside LCFS
                Deltapsi = 2d0 * (fs%psi(kf) - fs%psi(kf-1))
             end if
          else
             Deltapsi = fs%psi(kf-2) - fs%psi(kf)
          end if
          if (kf == nflux .and. tri%orient) then
             ! triangles outside LCFS are comparatively distorted
             ! use linear extrapolation instead, i.e. a triangle of the same size
             common_tri = [ktri, ktri]
          else
             common_tri = [ktri, mesh_element(ktri)%neighbour(tri%ef)]
          end if
          ! use midpoint of edge f
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          Bnpsi = -R0 * B0phi(ktri, tri%ef) / r
          Bnflux(ktri, tri%ef) = Bnpsi * r / Deltapsi * 2d0 * &
               sum(mesh_element_rmp(common_tri(:))%area)
          Bnphi(ktri) = imun / n * Bnflux(ktri, tri%ef) / tri%area
          Bnflux(ktri, tri%ei) = (0d0, 0d0)
          Bnflux(ktri, tri%eo) = (0d0, 0d0)
       end do
    end do
    if (quad_avg) call avg_flux_on_quad(Bnflux, Bnphi)
    call check_div_free(Bnflux, Bnphi, n, rel_err_Bn, 'non-resonant B_n')
    call check_redundant_edges(Bnflux, .false., 'non-resonant B_n')
  end subroutine compute_Bn_nonres

  subroutine interp_RT0(ktri, pol_flux, r, z, pol_comp_r, pol_comp_z)
    integer, intent(in) :: ktri
    complex(dp), intent(in) :: pol_flux(:,:)
    real(dp), intent(in) :: r, z
    complex(dp), intent(out) :: pol_comp_r, pol_comp_z
    type(triangle) :: elem
    type(knot) :: node(3)

    elem = mesh_element(ktri)
    node = mesh_point(elem%i_knot(:))
    ! edge 1 lies opposite to knot 3, etc.
    pol_comp_r = 1d0 / elem%det_3 * ( &
         pol_flux(ktri, 1) * (r - node(3)%rcoord) + &
         pol_flux(ktri, 2) * (r - node(1)%rcoord) + &
         pol_flux(ktri, 3) * (r - node(2)%rcoord))
    pol_comp_z = 1d0 / elem%det_3 * ( &
         pol_flux(ktri, 1) * (z - node(3)%zcoord) + &
         pol_flux(ktri, 2) * (z - node(1)%zcoord) + &
         pol_flux(ktri, 3) * (z - node(2)%zcoord))
  end subroutine interp_RT0

  pure function sqrt_g(kf, kt, r) result(metric_det)
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: r
    real(dp) :: metric_det
    metric_det = -fs_half%q(kf) * r / B0phi_Omega(kt_low(kf) + kt)
  end function sqrt_g

  subroutine radially_outward_normal(ktri, n_r, n_z)
    integer, intent(in) :: ktri
    real(dp), intent(out) :: n_r
    real(dp), intent(out) :: n_z
    type(triangle_rmp) :: tri

    tri = mesh_element_rmp(ktri)
    if (tri%orient) then
       n_r = mesh_point(tri%lf(2))%zcoord - mesh_point(tri%lf(1))%zcoord
       n_z = mesh_point(tri%lf(1))%rcoord - mesh_point(tri%lf(2))%rcoord
    else
       n_r = mesh_point(tri%lf(1))%zcoord - mesh_point(tri%lf(2))%zcoord
       n_z = mesh_point(tri%lf(2))%rcoord - mesh_point(tri%lf(1))%rcoord
    end if
  end subroutine radially_outward_normal

  subroutine write_vector_plot(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: kf, kt, ktri, n_cutoff, fid
    type(triangle_rmp) :: tri
    complex(dp) :: pol_comp_r, pol_comp_z, dens_psi_contravar, proj_theta_covar
    real(dp) :: r, z, n_r, n_z

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    if (nonres) then
       n_cutoff = nflux - 1
    else
       n_cutoff = nflux
    end if
    do kf = 1, n_cutoff
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          r = tri%R_Omega
          z = tri%Z_Omega
          call interp_RT0(ktri, pol_flux, r, z, pol_comp_r, pol_comp_z)
          ! projection to contravariant psi component
          call radially_outward_normal(ktri, n_r, n_z)
          dens_psi_contravar = (pol_comp_r * n_r + pol_comp_z * n_z) * &
               (fs%psi(kf) - fs%psi(kf-1)) / tri%area * 0.5d0 * sqrt_g(kf, kt, r)
          ! projection to covariant theta component
          proj_theta_covar = -(pol_comp_r * B0r_Omega(ktri) + &
               pol_comp_z * B0z_Omega(ktri)) * sqrt_g(kf, kt, r)
          write (fid, '(12(1x, es23.16))') r, z, real(pol_comp_r), aimag(pol_comp_r), &
               real(pol_comp_z), aimag(pol_comp_z), &
               real(tor_comp(ktri)), aimag(tor_comp(ktri)), &
               real(dens_psi_contravar), aimag(dens_psi_contravar), &
               real(proj_theta_covar), aimag(proj_theta_covar)
       end do
    end do
    r = mesh_point(1)%rcoord
    z = mesh_point(1)%zcoord
    do kt = kt_low(n_cutoff+1) + 1, ntri
       write (fid, '(12(1x, es23.16))') r, z, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(fid)
  end subroutine write_vector_plot

  subroutine write_vector_dof(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: ktri, fid

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do ktri = 1, kt_low(nflux+1)
       write (fid, '(8(1x, es23.16))') &
            real(pol_flux(ktri, 1)), aimag(pol_flux(ktri, 1)), &
            real(pol_flux(ktri, 2)), aimag(pol_flux(ktri, 2)), &
            real(pol_flux(ktri, 3)), aimag(pol_flux(ktri, 3)), &
            real(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0), &
            aimag(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0)
    end do
    do ktri = kt_low(nflux+1) + 1, ntri
       write (fid, '(8(1x, es23.16))') 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(fid)
  end subroutine write_vector_dof

  !> Real and imaginary part of \p scalar_dof (e.g. #presn) are written, in that order,
  !> to \p outfile (e.g. #magdif_conf::presn_file), where line number corresponds to the
  !> knot index in #mesh_mod::mesh_point.
  subroutine write_scalar_dof(scalar_dof, outfile)
    complex(dp), intent(in) :: scalar_dof(:)
    character(len = *), intent(in) :: outfile

    integer :: kpoint, fid

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do kpoint = 1, kp_low(nflux+1)
       write (fid, '(2(1x, es23.16))') real(scalar_dof(kpoint)), aimag(scalar_dof(kpoint))
    end do
    do kpoint = kp_low(nflux+1) + 1, npoint ! write zeroes in remaining points until end
       write (fid, '(2(1x, es23.16))') 0d0, 0d0
    end do
    close(fid)
  end subroutine write_scalar_dof

  subroutine write_fluxvar
    integer :: kf, fid
    real(dp) :: rho_r, rho_z

    open(newunit = fid, file = fluxvar_file, recl = longlines, status = 'replace')
    write (fid, '(5(1x, es23.16))') 0d0, fs%psi(0), fs%q(0), fs%p(0), fs%dp_dpsi(0)
    do kf = 1, nflux
       rho_r = mesh_point(kp_low(kf) + 1)%rcoord - mesh_point(1)%rcoord
       rho_z = mesh_point(kp_low(kf) + 1)%zcoord - mesh_point(1)%zcoord
       write (fid, '(5(1x, es23.16))') &
            hypot(rho_r, rho_z), fs%psi(kf), fs%q(kf), fs%p(kf), fs%dp_dpsi(kf)
    end do
    close(fid)
  end subroutine write_fluxvar

  subroutine write_kilca_modes(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer, parameter :: mmax = 24
    character(len = 19) :: fmt
    complex(dp), dimension(-mmax:mmax) :: coeff_rho, coeff_theta, coeff_zeta, fourier_basis
    integer :: kf, kt, ktri, m, fid_rho, fid_theta, fid_zeta, fid_furth
    type(triangle_rmp) :: tri
    complex(dp) :: pol_comp_r, pol_comp_z, pol_comp_rho, pol_comp_theta, sheet_current
    real(dp) :: r, z, rmaxis, zmaxis, rho_max, rho, theta, k_zeta, k_theta

    write (fmt, '(a, i3, a)') '(', 4 * mmax + 2 + 2, '(1es22.15, 1x))'
    rmaxis = mesh_point(1)%rcoord
    zmaxis = mesh_point(1)%zcoord
    rho_max = hypot(mesh_point(kp_low(nflux)+1)%rcoord - rmaxis, &
         mesh_point(kp_low(nflux)+1)%zcoord - zmaxis)
    k_zeta = n / R0
    open(newunit = fid_rho, recl = 3 * longlines, status = 'replace', &
         file = decorate_filename(outfile, '', '_r'))
    open(newunit = fid_theta, recl = 3 * longlines, status = 'replace', &
         file = decorate_filename(outfile, '', '_theta'))
    open(newunit = fid_zeta, recl = 3 * longlines, status = 'replace', &
         file = decorate_filename(outfile, '', '_z'))
    open(newunit = fid_furth, recl = longlines, status = 'replace', &
         file = decorate_filename(outfile, 'furth_', ''))
    do kf = 1, nflux
       coeff_rho = 0d0
       coeff_theta = 0d0
       coeff_zeta = 0d0
       sheet_current = 0d0
       rho = rho_max * (dble(kf) - 0.5d0) / dble(nflux)
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          theta = (dble(kt) - 0.5d0) / dble(kt_max(kf)) * 2d0 * pi
          ! result_spectrum.f90 uses negative q, so poloidal modes are switched
          ! we invert the sign here to keep post-processing consistent
          fourier_basis = [(exp(imun * m * theta), m = -mmax, mmax)]
          r = rmaxis + rho * cos(theta)
          z = zmaxis + rho * sin(theta)
          tri = mesh_element_rmp(ktri)
          call interp_RT0(ktri, pol_flux, r, z, pol_comp_r, pol_comp_z)
          pol_comp_rho = pol_comp_r * cos(theta) + pol_comp_z * sin(theta)
          pol_comp_theta = pol_comp_z * cos(theta) - pol_comp_r * sin(theta)
          coeff_rho = coeff_rho + pol_comp_rho * fourier_basis
          coeff_theta = coeff_theta + pol_comp_theta * fourier_basis
          coeff_zeta = coeff_zeta + tor_comp(ktri) * fourier_basis
          sheet_current = sheet_current + tri%area * jnphi(ktri) * &
               fourier_basis(kilca_pol_mode)
       end do
       coeff_rho = coeff_rho / kt_max(kf)
       coeff_theta = coeff_theta / kt_max(kf)
       coeff_zeta = coeff_zeta / kt_max(kf)
       write (fid_rho, fmt) rho, fs_half%q(kf), real(coeff_rho), aimag(coeff_rho)
       write (fid_theta, fmt) rho, fs_half%q(kf), real(coeff_theta), aimag(coeff_theta)
       write (fid_zeta, fmt) rho, fs_half%q(kf), real(coeff_zeta), aimag(coeff_zeta)
       k_theta = kilca_pol_mode / rho
       sheet_current = -2d0 * imun / clight / k_theta * sheet_current
       write (fid_furth, *) rho, k_zeta, k_theta, real(coeff_rho(-kilca_pol_mode)), &
            aimag(coeff_rho(-kilca_pol_mode)), real(sheet_current), aimag(sheet_current)
    end do
    close(fid_furth)
    do kt = kt_low(nflux+1) + 1, ntri
       write (fid_rho, fmt) rho_max, fs_half%q(nflux), [(0d0, m = 1, 4 * mmax + 2)]
       write (fid_theta, fmt) rho_max, fs_half%q(nflux), [(0d0, m = 1, 4 * mmax + 2)]
       write (fid_zeta, fmt) rho_max, fs_half%q(nflux), [(0d0, m = 1, 4 * mmax + 2)]
    end do
    close(fid_rho)
    close(fid_theta)
    close(fid_zeta)
  end subroutine write_kilca_modes
end module magdif
