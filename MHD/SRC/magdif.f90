module magdif
  use from_nrtype, only: dp                            ! PRELOAD/SRC/from_nrtype.f90
  use constants, only: pi                              ! PRELOAD/SRC/orbit_mod.f90
  use mesh_mod, only: npoint, ntri, knot, triangle, &  ! PRELOAD/SRC/mesh_mod.f90
       mesh_point, mesh_element
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

  real(dp), parameter :: clight = 2.99792458d10      !< Speed of light in cm sec^-1.
  complex(dp), parameter :: imun = (0.0_dp, 1.0_dp)  !< Imaginary unit in double precision.

contains

  !> Initialize magdif module
  subroutine magdif_init
    use input_files, only: gfile
    integer :: m

    ! only depends on config variables
    call init_indices
    if (kilca_scale_factor /= 0) then
       R0 = R0 * kilca_scale_factor
       n = n * kilca_scale_factor
    end if

    call read_mesh

    ! depends on mesh data
    call cache_equilibrium_field

    ! needs initialized field_eq
    allocate(efit)
    call efit%read(gfile)

    ! depends on mesh data, equilibrium field and EFIT profiles
    call init_flux_variables

    ! depends on safety factor
    call read_delayed_config
    if (log_info) then
       write (logfile, *) 'poloidal mode number, sheet current factor'
       do m = lbound(sheet_current_factor, 1), ubound(sheet_current_factor, 1)
          write (logfile, *) m, sheet_current_factor(m)
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
    if (log_info) write(logfile, *) 'magdif initialized'
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
    if (allocated(kp_max)) deallocate(kp_max)
    if (allocated(kt_max)) deallocate(kt_max)
    if (allocated(kp_low)) deallocate(kp_low)
    if (allocated(kt_low)) deallocate(kt_low)
    if (log_info) write(logfile, *) 'magdif cleanup finished'
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
          write (logfile, *) 'Arnoldi method yields ', ngrow, ' Ritz eigenvalues > ', tol
          do i = 1, ngrow
             write (logfile, *) 'lambda ', i, ': ', eigvals(i)
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
          if (log_info) then
             if (info == 0) then
                write (logfile, *) 'Successfully inverted matrix for preconditioner'
             else
                write (logfile, *) 'Matrix inversion for preconditioner failed: zgesv returns error', info
                stop
             end if
          end if
          do i = 1, ngrow
             Lr(i, :) = eigvals(i) * Yr(i, :)
          end do
          if (allocated(Yr)) deallocate(Yr)
       else
          preconditioned = .false.
       end if
    end if

    call pack2(Bnflux_vac, Bnphi_vac, Bn_prev)
    if (preconditioned) then
       Bn_prev = Bn_prev - matmul(eigvecs(:, 1:ngrow), matmul(Lr, &
            matmul(transpose(conjg(eigvecs(:, 1:ngrow))), Bn_prev)))
    end if
    do kiter = 0, niter-1
       if (log_info) write(logfile, *) 'Iteration ', kiter, ' of ', niter-1
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
       call write_vector_dof(Bnflux_diff, Bnphi_diff, &
            decorate_filename(Bn_diff_file, '', postfix))
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

    if (log_info) then
       write (logfile, *) "Number of points up to LCFS: ", kp_low(nflux+1), &
            "Number of triangles up to LCFS: ", kt_low(nflux+1)
    end if
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

    open(1, file = point_file, form = 'unformatted')
    read(1) npoint
    if (log_info) write(logfile, *) 'npoint = ', npoint
    allocate(mesh_point(npoint))
    read(1) mesh_point
    close(1)

    open(1, file = tri_file, form = 'unformatted')
    read(1) ntri
    if (log_info) write(logfile, *) 'ntri   = ', ntri
    allocate(mesh_element(ntri))
    read(1) mesh_element
    read(1) bphicovar
    close(1)

    allocate(B0r(ntri, 3))
    allocate(B0phi(ntri, 3))
    allocate(B0z(ntri, 3))
    allocate(B0r_Omega(ntri))
    allocate(B0phi_Omega(ntri))
    allocate(B0z_Omega(ntri))
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
    integer :: stat = 0, dummy = 0
    call execute_command_line("./maxwell.sh", exitstat = stat, cmdstat = dummy)
    if (stat /= 0) then
       if (log_err) write(logfile, *) 'FreeFem++ failed with exit code ', stat
       stop 'FreeFem++ failed'
    end if
  end subroutine compute_Bn

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
    character(len = len_trim(field_name) + 20) :: err_msg
    err_msg = trim(field_name) // ' not divergence-free'

    do kt = 1, ntri
       abs_flux = sum(abs(pol_flux(kt,:))) + abs(imun * n * tor_comp(kt) * &
            mesh_element(kt)%det_3 * 0.5d0)
       if (abs_flux > 0d0) then
          div = abs((sum(pol_flux(kt,:)) + imun * n * tor_comp(kt) * &
               mesh_element(kt)%det_3 * 0.5d0)) / abs_flux
          if (div > rel_err) then
             if (log_err) write(logfile, *) err_msg, ' in triangle ', kt, ': ', div
             stop err_msg
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
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    real(dp), parameter :: eps = epsilon(1d0)
    character(len = len_trim(name) + 30) :: err_msg
    err_msg = trim(name) // ': inconsistent redundant edges'

    checked = .false.
    do ktri = 1, kt_low(nflux+1)
       elem = mesh_element(ktri)
       call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
       do ke = 1, 3
          if (ktri > kt_low(nflux) .and. ke == ef) cycle
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
                if (log_err) write (logfile, *) err_msg, ' - ', &
                     name, '(', ktri, ',', ke, ') = ', pol_quant(ktri, ke), &
                     name, '(', ktri_adj, ',', ke_adj, ') = ', pol_quant(ktri_adj, ke_adj)
                stop err_msg
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
    integer :: kt
    real(dp) :: dummy8(8)

    open(1, file = filename, recl = longlines)
    do kt = 1, ntri
       read(1, *) dummy8
       Bnflux(kt,1) = cmplx(dummy8(1), dummy8(2), dp)
       Bnflux(kt,2) = cmplx(dummy8(3), dummy8(4), dp)
       Bnflux(kt,3) = cmplx(dummy8(5), dummy8(6), dp)
       Bnphi(kt) = cmplx(dummy8(7), dummy8(8), dp) / mesh_element(kt)%det_3 * 2d0
    end do
    close(1)

    call check_div_free(Bnflux, Bnphi, n, rel_err_Bn, 'B_n')
    call check_redundant_edges(Bnflux, .false., 'B_n')
  end subroutine read_Bn

  !> Allocates and computes the safety factor #q and #m_res.
  !>
  !> Also allocates #magdif_config::sheet_current_factor, to be read in via
  !> magdif_config::read_delayed_config() in magdif_init(). All deallocation is done in
  !> magdif_cleanup().
  subroutine compute_safety_factor
    integer :: kf, kt
    type(triangle) :: elem
    integer :: m_res_min, m_res_max, m
    real(dp), dimension(nflux) :: abs_err

    fs_half%q = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          fs_half%q(kf) = fs_half%q(kf) + B0phi_Omega(kt_low(kf) + kt) * elem%det_3
       end do
       fs_half%q(kf) = fs_half%q(kf) * 0.25d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    ! use linear interpolation for full-grid steps for now
    fs%q(1:nflux-1) = 0.5d0 * (fs_half%q(1:nflux-1) + fs_half%q(2:nflux))
    ! linear extrapolation for values at separatrix and magnetic axis
    fs%q(nflux) = fs%q(nflux-1) + (fs%q(nflux-1) - fs%q(nflux-2))
    fs%q(0) = fs%q(1) - (fs%q(2) - fs%q(1))

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
    if (log_debug) write (logfile, *) 'resonant m: ', m_res_min, '..', m_res_max
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
    ! linear extrapolation for value just outside LCFS
    fs%psi(nflux+1) = fs%psi(nflux) + (fs%psi(nflux) - fs%psi(nflux-1))

    allocate(fs_half)
    call fs_half%init(nflux, .true.)

    ! use linear interpolation for half-grid steps for now
    fs_half%psi = 0.5d0 * (fs%psi(0:nflux-1) + fs%psi(1:nflux))

    allocate(fluxvar)
    call fluxvar%init(4, efit%nw, nflux, fs%psi(:nflux), fs_half%psi)

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
       dens = (fs%psi(:nflux) - psimin) / psimax * di0 + d_min
       temp = (fs%psi(:nflux) - psimin) / psimax * ti0 + t_min
       if (log_info)  write (logfile, *) 'temp@axis: ', temp(0), ' dens@axis: ', dens(0)
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psimin) / psimax * di0 + d_min
       temp(1:) = (fs_half%psi - psimin) / psimax * ti0 + t_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_par)
       ddens_dpsi = (di0 - d_min) / (psimax - psimin)
       dtemp_dpsi = (ti0 - t_min) / (psimax - psimin)
       dens = (fs%psi(:nflux) - psimin) / (psimax - psimin) * (di0 - d_min) + d_min
       temp = (fs%psi(:nflux) - psimin) / (psimax - psimin) * (ti0 - t_min) + t_min
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
          fs_half%p(kf) = fluxvar%interp(efit%pres, kf, .false.)
          fs_half%dp_dpsi(kf) = fluxvar%interp(efit%pprime, kf, .true.)
       end do
    case default
       stop 'Error: unknown pressure profile selection'
    end select
  end subroutine compute_pres_prof


  subroutine cache_equilibrium_field
    real(dp) :: r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kf, kt, ke
    type(triangle) :: elem
    type(knot) :: base, tip

    open(1, file = 'plot_B0.dat', recl = longlines)
    p = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          do ke = 1, 3
             base = mesh_point(elem%i_knot(ke))
             tip = mesh_point(elem%i_knot(mod(ke, 3) + 1))
             r = (base%rcoord + tip%rcoord) * 0.5d0
             z = (base%zcoord + tip%zcoord) * 0.5d0
             call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
                  dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
             B0r(kt_low(kf) + kt, ke) = Br
             B0phi(kt_low(kf) + kt, ke) = Bp
             B0z(kt_low(kf) + kt, ke) = Bz
          end do
          call ring_centered_avg_coord(elem, r, z)
          call field(r, p, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          B0r_Omega(kt_low(kf) + kt) = Br
          B0phi_Omega(kt_low(kf) + kt) = Bp
          B0z_Omega(kt_low(kf) + kt) = Bz
          write (1, *) r, z, Br, Bz, Bp
       end do
    end do
    do kt = kt_low(nflux+1) + 1, ntri
       write (1, *) R0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine cache_equilibrium_field

  !> Computes equilibrium current density #j0phi from given equilibrium magnetic field and
  !> assumed equilibrium pressure #pres0.
  !>
  !> This step is necessary because equilibrium pressure is not given experimentally as is
  !> \f$ \vec{B}_{0} \f$; arbitrary values are assumed. Consistency of MHD equilibrium is
  !> necessary in the derivation, while Ampere's equation is not used.
  subroutine compute_j0phi
    use input_files, only: gfile

    integer :: kf, kt
    real(dp) :: r, z
    real(dp) :: Btor2
    real(dp), dimension(nflux) :: B2avg, B2avg_half
    real(dp) :: plot_j0phi
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient

    open(1, file = j0phi_file, recl = longlines)
    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, nflux
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0r(kt_low(kf) + kt, ef) ** 2 + &
               B0phi(kt_low(kf) + kt, ef) ** 2 + B0z(kt_low(kf) + kt, ef) ** 2
          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0r(kt_low(kf) + kt, ei) ** 2 + &
               B0phi(kt_low(kf) + kt, ei) ** 2 + B0z(kt_low(kf) + kt, ei) ** 2
       end do
       B2avg(kf) = B2avg(kf) / kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / kt_max(kf)

       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)

          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             if (orient) then
                j0phi(kt_low(kf) + kt, ef) = clight * (fs%dp_dpsi(kf) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf) / r)
             else
                j0phi(kt_low(kf) + kt, ef) = clight * (fs%dp_dpsi(kf-1) * r + &
                     0.25d0 / pi * fs%FdF_dpsi(kf-1) / r)
             end if
          case (curr_prof_rot)
             z = sum(mesh_point(lf(:))%zcoord) * 0.5d0
             j0phi(kt_low(kf) + kt, ef) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(kt_low(kf) + kt, ef) ** 2
             if (.not. orient) then
                j0phi(kt_low(kf) + kt, ef) = clight * fs%dp_dpsi(kf-1) * (1d0 - Btor2 / &
                     B2avg(kf-1)) * r
             else
                j0phi(kt_low(kf) + kt, ef) = clight * fs%dp_dpsi(kf) * (1d0 - Btor2 / &
                     B2avg(kf)) * r
             end if
          end select

          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             j0phi(kt_low(kf) + kt, ei) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(li(:))%zcoord) * 0.5d0
             j0phi(kt_low(kf) + kt, ei) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(kt_low(kf) + kt, ei) ** 2
             j0phi(kt_low(kf) + kt, ei) = clight * (fs%p(kf) - fs%p(kf-1)) / &
                  (fs%psi(kf) - fs%psi(kf-1)) * (1d0 - Btor2 / B2avg_half(kf)) * r
          end select

          r = sum(mesh_point(lo(:))%rcoord) * 0.5d0
          select case (curr_prof)
          case (curr_prof_efit)
             j0phi(kt_low(kf) + kt, eo) = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             z = sum(mesh_point(lo(:))%zcoord) * 0.5d0
             j0phi(kt_low(kf) + kt, eo) = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi(kt_low(kf) + kt, eo) ** 2
             j0phi(kt_low(kf) + kt, eo) = clight * (fs%p(kf) - fs%p(kf-1)) / &
                  (fs%psi(kf) - fs%psi(kf-1)) * (1d0 - Btor2 / B2avg_half(kf)) * r
          end select

          call ring_centered_avg_coord(elem, r, z)
          select case (curr_prof)
          case (curr_prof_efit)
             plot_j0phi = clight * (fs_half%dp_dpsi(kf) * r + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          case (curr_prof_rot)
             plot_j0phi = j0phi_ampere(r, z)
          case (curr_prof_ps)
             Btor2 = B0phi_Omega(kt_low(kf) + kt) ** 2
             plot_j0phi = clight * (fs%p(kf) - fs%p(kf-1)) / &
                  (fs%psi(kf) - fs%psi(kf-1)) * (1d0 - Btor2 / B2avg_half(kf)) * r
          end select

          write (1, *) j0phi(kt_low(kf) + kt, 1), j0phi(kt_low(kf) + kt, 2), &
               j0phi(kt_low(kf) + kt, 3), plot_j0phi
       end do
    end do
    do kt = kt_low(nflux+1) + 1, ntri
       write (1, *) j0phi(kt, 1), j0phi(kt, 2), j0phi(kt, 3), 0d0
    end do
    close(1)

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
    integer :: kf, kt, fid_amp, fid_gs
    real(dp) :: r, z, n_r, n_z, cmp_amp, cmp_gs, Jr, Jp, Jz
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    real(dp) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    open(1, file = 'cmp_prof.dat', recl = longlines)
    open(newunit = fid_amp, file = 'j0_amp.dat', recl = longlines)
    open(newunit = fid_gs, file = 'j0_gs.dat', recl = longlines)
    do kf = 1, nflux
       cmp_amp = 0d0
       cmp_gs = 0d0
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient, n_r, n_z)
          call ring_centered_avg_coord(elem, r, z)
          call field(r, 0d0, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
               dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
          dBpdR = dBpdR + fs_half%FdF_dpsi(kf) / fs_half%F(kf) * Bz
          dBpdZ = -fs_half%FdF_dpsi(kf) / fs_half%F(kf) * Br
          Jr = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * &
               B0r_Omega(kt_low(kf) + kt)
          Jz = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * &
               B0z_Omega(kt_low(kf) + kt)
          Jp = clight * (fs_half%dp_dpsi(kf) * r + 0.25d0 / pi * fs_half%FdF_dpsi(kf) / r)
          write (fid_gs, *) Jr, Jp, Jz
          cmp_gs = cmp_gs + ((Jp * Bz - Jz * Bp) * n_r + (Jr * Bp - Jp * Br) * n_z) / &
               (n_r * n_r + n_z * n_z) / (fs%psi(kf) - fs%psi(kf-1)) * elem%det_3 / clight
          Jr = 0.25d0 / pi * clight * (-dBpdZ)
          Jp = 0.25d0 / pi * clight * (dBrdZ - dBzdR)
          Jz = 0.25d0 / pi * clight * (dBpdR + Bp / r)
          write (fid_amp, *) Jr, Jp, Jz
          cmp_amp = cmp_amp + ((Jp * Bz - Jz * Bp) * n_r + (Jr * Bp - Jp * Br) * n_z) / &
               (n_r * n_r + n_z * n_z) / (fs%psi(kf) - fs%psi(kf-1)) * elem%det_3 / clight
       end do
       cmp_amp = cmp_amp / kt_max(kf)
       cmp_gs = cmp_gs / kt_max(kf)
       write (1, *) cmp_amp, cmp_gs
    end do
    close(1)
    close(fid_amp)
    close(fid_gs)
  end subroutine check_curr0

  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  subroutine compute_presn
    real(dp) :: r
    real(dp) :: lr, lz  ! edge vector components
    complex(dp) :: Bnpsi
    complex(dp), dimension(nkpol) :: a, b, x, d, du, inhom
    complex(dp), dimension(:), allocatable :: resid
    real(dp), dimension(nkpol) :: rel_err
    real(dp) :: max_rel_err, avg_rel_err
    integer :: kf, kt, kp
    integer :: nz
    integer, dimension(2*nkpol) :: irow, icol
    complex(dp), dimension(2*nkpol) :: aval
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    type(knot) :: base, tip
    integer :: common_tri(2)

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux
inner: do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          if (.not. orient) cycle inner
          if (kf < nflux) then
             common_tri = (/ kt_low(kf) + kt, elem%neighbour(ef) /)
          else
             ! triangles outside LCFS are comparatively distorted
             ! use linear extrapolation instead, i.e. a triangle of the same size
             common_tri = (/ kt_low(kf) + kt, kt_low(kf) + kt /)
          end if
          ! use midpoint of edge f
          base = mesh_point(lf(1))
          tip = mesh_point(lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord

          Bnpsi = Bnflux(kt_low(kf) + kt, ef) / r * (fs%psi(kf+1) - fs%psi(kf-1)) / &
               sum(mesh_element(common_tri(:))%det_3)

          kp = lf(1) - kp_low(kf)

          x(kp) = -fs%dp_dpsi(kf) * Bnpsi

          a(kp) = (B0r(kt_low(kf) + kt, ef) * lr + B0z(kt_low(kf) + kt, ef) * lz) / &
               (lr ** 2 + lz ** 2)

          b(kp) = imun * (n + imun * damp) * B0phi(kt_low(kf) + kt, ef) / r
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
    if (log_debug) write (logfile, *) 'compute_presn: diagonalization' // &
         ' max_rel_err = ', max_rel_err, ' avg_rel_err = ', avg_rel_err
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

  subroutine get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient, nfr, nfz)
    type(triangle), intent(in) :: elem
    integer, dimension(2), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    logical, intent(out) :: orient
    real(dp), intent(out), optional :: nfr, nfz
    integer, dimension(3) :: i_knot_diff
    integer :: knot_i, knot_o, knot_f
    integer :: i1, i2
    logical :: closing_loop
    character(len = *), parameter :: errmsg = 'cannot find correct label for triangle edges'

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
       if (log_debug) write(logfile, *) errmsg
       stop errmsg
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
       if (log_debug) write(logfile, *) errmsg
       stop errmsg
    end if
    if (present(nfr)) then
       if (orient) then
          nfr = mesh_point(lf(2))%zcoord - mesh_point(lf(1))%zcoord
       else
          nfr = mesh_point(lf(1))%zcoord - mesh_point(lf(2))%zcoord
       end if
    end if
    if (present(nfz)) then
       if (orient) then
          nfz = mesh_point(lf(1))%rcoord - mesh_point(lf(2))%rcoord
       else
          nfz = mesh_point(lf(2))%rcoord - mesh_point(lf(1))%rcoord
       end if
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
    real(dp) :: area
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    real(dp) :: B0flux
    complex(dp) :: Bnphi_Gamma
    real(dp) :: r
    real(dp) :: n_r, n_z

    max_rel_err = 0d0
    avg_rel_err = 0d0
    do kf = 1, nflux ! loop through flux surfaces
       do kt = 1, kt_max(kf)
          ktri = kt_low(kf) + kt
          elem = mesh_element(ktri)
          area = elem%det_3 * 0.5d0
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          ! first term on source side: flux through edge f
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          jnflux(ktri, ef) = j0phi(ktri, ef) / B0phi(ktri, ef) * Bnflux(ktri, ef) + &
               clight * r / B0phi(ktri, ef) * (presn(lf(2)) - presn(lf(1)))
          x(kt) = -jnflux(ktri, ef)
          ! diagonal matrix element - edge i
          n_r = mesh_point(li(2))%zcoord - mesh_point(li(1))%zcoord
          n_z = mesh_point(li(1))%rcoord - mesh_point(li(2))%rcoord
          B0flux = r * (B0r(ktri, ei) * n_r + B0z(ktri, ei) * n_z)
          d(kt) = -1d0 - imun * (n + imun * damp) * area * 0.5d0 * &
               B0phi(ktri, ei) / B0flux
          ! additional term from edge i on source side
          r = sum(mesh_point(li(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(ei)))
          x(kt) = x(kt) - imun * n * area * 0.5d0 * (clight * r / B0flux * &
               (Bnphi_Gamma / B0phi(ktri, ei) * (fs%p(kf) - fs%p(kf-1)) - &
               (presn(li(2)) - presn(li(1)))) + j0phi(ktri, ei) * &
               (Bnphi_Gamma / B0phi(ktri, ei) - Bnflux(ktri, ei) / B0flux))
          ! superdiagonal matrix element - edge o
          n_r = mesh_point(lo(2))%zcoord - mesh_point(lo(1))%zcoord
          n_z = mesh_point(lo(1))%rcoord - mesh_point(lo(2))%rcoord
          B0flux = r * (B0r(ktri, eo) * n_r + B0z(ktri, eo) * n_z)
          du(kt) = 1d0 + imun * (n + imun * damp) * area * 0.5d0 * &
               B0phi(ktri, eo) / B0flux
          ! additional term from edge o on source side
          r = sum(mesh_point(lo(:))%rcoord) * 0.5d0
          Bnphi_Gamma = 0.5d0 * (Bnphi(ktri) + Bnphi(elem%neighbour(eo)))
          x(kt) = x(kt) - imun * n * area * 0.5d0 * (clight * r / B0flux * &
               (Bnphi_Gamma / B0phi(ktri, eo) * (fs%p(kf-1) - fs%p(kf)) - &
               (presn(lo(2)) - presn(lo(1)))) + j0phi(ktri, eo) * &
               (Bnphi_Gamma / B0phi(ktri, eo) - Bnflux(ktri, eo) / B0flux))
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
          elem = mesh_element(ktri)
          area = elem%det_3 * 0.5d0
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          jnflux(ktri, ei) = -x(kt)
          jnflux(ktri, eo) = x(mod(kt, kt_max(kf)) + 1)
          if (m_res(kf) > 0) then
             if (sheet_current_factor(m_res(kf)) /= (0d0, 0d0)) then
                ! add sheet current on edge i
                r = sum(mesh_point(li(:))%rcoord) * 0.5d0
                n_r = mesh_point(li(2))%zcoord - mesh_point(li(1))%zcoord
                n_z = mesh_point(li(1))%rcoord - mesh_point(li(2))%rcoord
                B0flux = r * (B0r(ktri, ei) * n_r + B0z(ktri, ei) * n_z)
                jnflux(ktri, ei) = jnflux(ktri, ei) + sheet_current_factor(m_res(kf)) * &
                     B0flux * sum(presn(li(:))) * 0.5d0
                ! add sheet current on edge o
                r = sum(mesh_point(lo(:))%rcoord) * 0.5d0
                n_r = mesh_point(lo(2))%zcoord - mesh_point(lo(1))%zcoord
                n_z = mesh_point(lo(1))%rcoord - mesh_point(lo(2))%rcoord
                B0flux = r * (B0r(ktri, eo) * n_r + B0z(ktri, eo) * n_z)
                jnflux(ktri, eo) = jnflux(ktri, eo) + sheet_current_factor(m_res(kf)) * &
                     B0flux * sum(presn(lo(:))) * 0.5d0
             end if
          end if
          jnphi(ktri) = sum(jnflux(ktri, :)) * imun / n / area
       end do
    end do
    avg_rel_err = avg_rel_err / sum(kt_max(1:nflux))
    if (log_debug) write (logfile, *) 'compute_currn: diagonalization' // &
         ' max_rel_err = ', max_rel_err, ' avg_rel_err = ', avg_rel_err
    if (allocated(resid)) deallocate(resid)

    call check_redundant_edges(jnflux, .false., 'jnflux')
    call check_div_free(jnflux, jnphi, n, rel_err_currn, 'jnflux')
    call write_vector_dof(jnflux, jnphi, currn_file)
    call write_vector_plot(jnflux, jnphi, decorate_filename(currn_file, 'plot_', ''))
  end subroutine compute_currn

  subroutine avg_flux_on_quad(pol_flux, tor_comp)
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)

    integer :: kf, kt
    complex(dp) :: tor_flux_avg, tor_flux_diff
    type(triangle) :: elem1, elem2
    integer, dimension(2) :: li, lo, lf
    integer :: ei1, eo1, ef1, ei2, eo2, ef2
    logical :: orient

    do kf = 2, nflux
       do kt = 1, kt_max(kf), 2
          elem1 = mesh_element(kt_low(kf) + kt)
          elem2 = mesh_element(kt_low(kf) + kt + 1)
          call get_labeled_edges(elem1, li, lo, lf, ei1, eo1, ef1, orient)
          call get_labeled_edges(elem2, li, lo, lf, ei2, eo2, ef2, orient)
          tor_flux_avg = 0.25d0 * (tor_comp(kt_low(kf) + kt) * elem1%det_3 + &
               tor_comp(kt_low(kf) + kt + 1) * elem2%det_3)
          tor_flux_diff = 0.25d0 * (-tor_comp(kt_low(kf) + kt) * elem1%det_3 + &
               tor_comp(kt_low(kf) + kt + 1) * elem2%det_3)
          tor_comp(kt_low(kf) + kt) = tor_flux_avg / elem1%det_3 * 2d0
          pol_flux(kt_low(kf) + kt, eo1) = pol_flux(kt_low(kf) + kt, eo1) - imun * n * &
                  tor_flux_diff
          tor_comp(kt_low(kf) + kt + 1) = tor_flux_avg / elem2%det_3 * 2d0
          pol_flux(kt_low(kf) + kt + 1, ei2) = pol_flux(kt_low(kf) + kt + 1, ei2) + &
               imun * n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  subroutine compute_Bn_nonres
    integer :: kf, kt
    type(triangle) :: elem
    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    integer :: common_tri(2)
    real(dp) :: Deltapsi
    complex(dp) :: Bnpsi
    real(dp) :: r

    do kf = 1, nflux ! loop through flux surfaces
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient)
          if (orient) then
             Deltapsi = fs%psi(kf+1) - fs%psi(kf-1)
          else
             Deltapsi = fs%psi(kf-2) - fs%psi(kf)
          end if
          if (kf == nflux .and. orient) then
             ! triangles outside LCFS are comparatively distorted
             ! use linear extrapolation instead, i.e. a triangle of the same size
             common_tri = (/ kt_low(kf) + kt, kt_low(kf) + kt /)
          else
             common_tri = (/ kt_low(kf) + kt, elem%neighbour(ef) /)
          end if
          ! use midpoint of edge f
          r = sum(mesh_point(lf(:))%rcoord) * 0.5d0
          Bnpsi = -R0 * B0phi(kt_low(kf) + kt, ef) / r
          Bnflux(kt_low(kf) + kt, ef) = Bnpsi * r / Deltapsi * &
               sum(mesh_element(common_tri(:))%det_3)
          Bnphi(kt_low(kf) + kt) = imun / n * Bnflux(kt_low(kf) + kt, ef) &
               / elem%det_3 * 2d0
          Bnflux(kt_low(kf) + kt, ei) = (0d0, 0d0)
          Bnflux(kt_low(kf) + kt, eo) = (0d0, 0d0)
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
    type(knot) :: tri(3)

    elem = mesh_element(ktri)
    tri = mesh_point(elem%i_knot(:))
    ! edge 1 lies opposite to knot 3, etc.
    pol_comp_r = 1d0 / elem%det_3 * ( &
         pol_flux(ktri, 1) * (r - tri(3)%rcoord) + &
         pol_flux(ktri, 2) * (r - tri(1)%rcoord) + &
         pol_flux(ktri, 3) * (r - tri(2)%rcoord))
    pol_comp_z = 1d0 / elem%det_3 * ( &
         pol_flux(ktri, 1) * (z - tri(3)%zcoord) + &
         pol_flux(ktri, 2) * (z - tri(1)%zcoord) + &
         pol_flux(ktri, 3) * (z - tri(2)%zcoord))
  end subroutine interp_RT0

  pure function sqrt_g(kf, kt, r) result(metric_det)
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: r
    real(dp) :: metric_det
    metric_det = -fs_half%q(kf) * r / B0phi_Omega(kt_low(kf) + kt)
  end function sqrt_g

  subroutine write_vector_plot(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: kf, kt, n_cutoff
    type(triangle) :: elem
    complex(dp) :: pol_comp_r, pol_comp_z, dens_psi_contravar, proj_theta_covar
    real(dp) :: r, z

    integer, dimension(2) :: li, lo, lf
    integer :: ei, eo, ef
    logical :: orient
    real(dp) :: n_r, n_z

    open(1, file = outfile, recl = longlines)
    if (nonres) then
       n_cutoff = nflux - 1
    else
       n_cutoff = nflux
    end if
    do kf = 1, n_cutoff
       do kt = 1, kt_max(kf)
          elem = mesh_element(kt_low(kf) + kt)
          call ring_centered_avg_coord(elem, r, z)
          call interp_RT0(kt_low(kf) + kt, pol_flux, r, z, pol_comp_r, pol_comp_z)
          call get_labeled_edges(elem, li, lo, lf, ei, eo, ef, orient, n_r, n_z)
          ! projection to contravariant psi component
          dens_psi_contravar = (pol_comp_r * n_r + pol_comp_z * n_z) * &
               (fs%psi(kf) - fs%psi(kf-1)) / elem%det_3 * sqrt_g(kf, kt, r)
          ! projection to covariant theta component
          proj_theta_covar = -(pol_comp_r * B0r_Omega(kt_low(kf) + kt) + &
               pol_comp_z * B0z_Omega(kt_low(kf) + kt)) * sqrt_g(kf, kt, r)
          write (1, *) r, z, real(pol_comp_r), aimag(pol_comp_r), real(pol_comp_z), &
               aimag(pol_comp_z), real(tor_comp(kt_low(kf) + kt)), &
               aimag(tor_comp(kt_low(kf) + kt)), real(dens_psi_contravar), &
               aimag(dens_psi_contravar), real(proj_theta_covar), aimag(proj_theta_covar)
       end do
    end do
    r = mesh_point(1)%rcoord
    z = mesh_point(1)%zcoord
    do kt = kt_low(n_cutoff+1) + 1, ntri
       write (1, *) r, z, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine write_vector_plot

  subroutine write_vector_dof(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: ktri

    open(1, file = outfile, recl = longlines)
    do ktri = 1, kt_low(nflux+1)
       write(1, *) &
            real(pol_flux(ktri, 1)), aimag(pol_flux(ktri, 1)), &
            real(pol_flux(ktri, 2)), aimag(pol_flux(ktri, 2)), &
            real(pol_flux(ktri, 3)), aimag(pol_flux(ktri, 3)), &
            real(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0), &
            aimag(tor_comp(ktri) * mesh_element(ktri)%det_3 * 0.5d0)
    end do
    do ktri = kt_low(nflux+1) + 1, ntri
       write (1, *) 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
    end do
    close(1)
  end subroutine write_vector_dof

  !> Real and imaginary part of \p scalar_dof (e.g. #presn) are written, in that order,
  !> to \p outfile (e.g. #magdif_conf::presn_file), where line number corresponds to the
  !> knot index in #mesh_mod::mesh_point.
  subroutine write_scalar_dof(scalar_dof, outfile)
    complex(dp), intent(in) :: scalar_dof(:)
    character(len = *), intent(in) :: outfile

    integer :: kpoint

    open(1, file = outfile, recl = longlines)
    do kpoint = 1, kp_low(nflux+1)
       write (1, *) real(scalar_dof(kpoint)), aimag(scalar_dof(kpoint))
    end do
    do kpoint = kp_low(nflux+1) + 1, npoint ! write zeroes in remaining points until end
       write (1, *) 0d0, 0d0
    end do
    close(1)
  end subroutine write_scalar_dof

  subroutine write_fluxvar
    integer :: kf
    real(dp) :: rho_r, rho_z

    open(1, file = fluxvar_file, recl = longlines)
    write (1, *) 0d0, fs%psi(0), fs_half%q(1), fs%p(0), fs%dp_dpsi(0)
    do kf = 1, nflux
       rho_r = mesh_point(kp_low(kf) + 1)%rcoord - mesh_point(1)%rcoord
       rho_z = mesh_point(kp_low(kf) + 1)%zcoord - mesh_point(1)%zcoord
       write (1, *) hypot(rho_r, rho_z), fs%psi(kf), fs_half%q(kf), fs%p(kf), &
            fs%dp_dpsi(kf)
    end do
    close(1)
  end subroutine write_fluxvar

  subroutine write_kilca_modes(pol_flux, tor_comp, outfile)
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer, parameter :: mmax = 24
    character(len = 19) :: fmt
    complex(dp), dimension(-mmax:mmax) :: coeff_rho, coeff_theta, coeff_zeta, fourier_basis
    integer :: kf, kt, m, fid_rho, fid_theta, fid_zeta, fid_furth
    type(triangle) :: elem
    complex(dp) :: pol_comp_r, pol_comp_z, pol_comp_rho, pol_comp_theta, sheet_current
    real(dp) :: r, z, rmaxis, zmaxis, rho_max, rho, theta, k_zeta, k_theta

    write (fmt, '(a, i3, a)') '(', 4 * mmax + 2 + 2, '(1es22.15, 1x))'
    rmaxis = mesh_point(1)%rcoord
    zmaxis = mesh_point(1)%zcoord
    rho_max = hypot(mesh_point(kp_low(nflux)+1)%rcoord - rmaxis, &
         mesh_point(kp_low(nflux)+1)%zcoord - zmaxis)
    k_zeta = n / R0
    open(newunit = fid_rho, recl = 3 * longlines, &
         file = decorate_filename(outfile, '', '_r'))
    open(newunit = fid_theta, recl = 3 * longlines, &
         file = decorate_filename(outfile, '', '_theta'))
    open(newunit = fid_zeta, recl = 3 * longlines, &
         file = decorate_filename(outfile, '', '_z'))
    open(newunit = fid_furth, recl = longlines, &
         file = decorate_filename(outfile, 'furth_', ''))
    do kf = 1, nflux
       coeff_rho = 0d0
       coeff_theta = 0d0
       coeff_zeta = 0d0
       sheet_current = 0d0
       rho = rho_max * (dble(kf) - 0.5d0) / dble(nflux)
       do kt = 1, kt_max(kf)
          theta = (dble(kt) - 0.5d0) / dble(kt_max(kf)) * 2d0 * pi
          ! result_spectrum.f90 uses negative q, so poloidal modes are switched
          ! we invert the sign here to keep post-processing consistent
          fourier_basis = [(exp(imun * m * theta), m = -mmax, mmax)]
          r = rmaxis + rho * cos(theta)
          z = zmaxis + rho * sin(theta)
          elem = mesh_element(kt_low(kf) + kt)
          call interp_RT0(kt_low(kf) + kt, pol_flux, r, z, pol_comp_r, pol_comp_z)
          pol_comp_rho = pol_comp_r * cos(theta) + pol_comp_z * sin(theta)
          pol_comp_theta = pol_comp_z * cos(theta) - pol_comp_r * sin(theta)
          coeff_rho = coeff_rho + pol_comp_rho * fourier_basis
          coeff_theta = coeff_theta + pol_comp_theta * fourier_basis
          coeff_zeta = coeff_zeta + tor_comp(kt_low(kf) + kt) * fourier_basis
          sheet_current = sheet_current + jnphi(kt_low(kf) + kt) * elem%det_3 * 0.5d0
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
