module magdif
  use iso_fortran_env, only: dp => real64
  use magdif_util, only: g_eqdsk, flux_func, flux_func_cache

  implicit none

  private

  public :: equil, fluxvar, fs, fs_half, magdif_mesh, mesh, &
       B0flux, Bnflux, Bnphi, init_indices, cache_mesh_data, check_redundant_edges, &
       check_div_free, point_location, interp_RT0, write_vector_dof, magdif_init, &
       magdif_cleanup, magdif_single, magdif_iterated, flux_func_cache_check

  character(len = 1024) :: magdif_bin_dir = '.'

  type(g_eqdsk) :: equil
  type(flux_func) :: fluxvar
  type(flux_func_cache) :: fs
  type(flux_func_cache) :: fs_half

  type :: magdif_mesh

     !> R coordinate of the O point in cm.
     real(dp) :: R_O

     !> Z coordinate of the O point in cm.
     real(dp) :: Z_O

     !> R coordinate of the X point in cm.
     real(dp) :: R_X

     !> Z coordinate of the X point in cm.
     real(dp) :: Z_X

     !> Minimal R value on computational grid in cm.
     real(dp) :: R_min

     !> Maximal Z value on computational grid in cm.
     real(dp) :: R_max

     !> Minimal R value on computational grid in cm.
     real(dp) :: Z_min

     !> Maximal Z value on computational grid in cm.
     real(dp) :: Z_max

     !> Number of flux surfaces. May differ from #magdif_conf::magdif_config::nflux due
     !> to refinement of flux surfaces.
     integer :: nflux

     !> Number of triangle edges. TODO: Add npoint and ntri along with mesh_element_rmp.
     integer :: nedge

     !> Toroidal mode number. May differ from #magdif_conf::magdif_config#n due to large
     !> aspect ratio scaling.
     integer :: n

     !> Minimal poloidal mode number in resonance.
     integer :: m_res_min

     !> Maximal poloidal mode number in resonance
     integer ::  m_res_max

     !> Number of unrefined flux surfaces to be replaced by refined ones.
     integer, allocatable :: deletions(:)

     !> Number of refined flux surfaces.
     integer, allocatable :: additions(:)

     !> Relative size of most refined flux surface.
     real(dp), allocatable :: refinement(:)

     !> Poloidal mode number \f$ m \f$ (dimensionless) in resonance at given flux surface.
     !>
     !> Indexing is the same as for #q, on which the values depend. If no resonances are
     !> expected at a given index, #m_res is 0.
     integer, allocatable :: m_res(:)

     !> Indices of flux surfaces where resonance corresponding to a poloidal mode (given as
     !> array index) occurs.
     integer, allocatable :: res_ind(:)

     !> Number of knots on the flux surface given by the array index.
     !>
     !> The array index ranges from 1 for the innermost flux surface to
     !> #magdif_config::nflux for the last closed flux surface.
     integer, allocatable :: kp_max(:)

     !> Number of triangles inside the flux surface given by the array index.
     !>
     !> The array index ranges from 1 for the innermost flux surface to
     !> #magdif_config::nflux for the last closed flux surface.
     integer, allocatable :: kt_max(:)

     !> Global index of the last knot of the previous flux surface given by the array index.
     !>
     !> The global index of knots in #mesh_mod::mesh_point on the flux surface kf runs from
     !> #kp_low (kf)+1 to #kp_low (kf)+#kp_max (kf), so #kp_low is determined by cumulatively
     !> adding consecutive values of #kp_max. The array index ranges from 1, giving the
     !> global index of the knot on the magnetic axis (which has to be 1), to
     !> #magdif_config::nflux+1, effectively giving the last knot on the last closed
     !> flux surface.
     integer, allocatable :: kp_low(:)

     !> Global index of the last triangle of the previous flux surface given by the array
     !> index.
     !>
     !> The global index of triangles in #mesh_mod::mesh_element inside the flux surface kf
     !> runs from #kt_low (kf)+1 to #kt_low (kf)+#kt_max (kf), so #kt_low is determined by
     !> cumulatively adding consecutive values of #kt_max. The array index ranges from 1,
     !> giving the global index of the non-existent triangle on the magnetic axis (which is
     !> therefore 0), to #magdif_config::nflux+1, giving the last triangle inside the last
     !> closed flux surface.
     integer, allocatable :: kt_low(:)

     integer, allocatable :: edge_map2global(:, :)
     integer, allocatable :: edge_map2ktri(:, :)
     integer, allocatable :: edge_map2ke(:, :)

   contains
     final :: magdif_mesh_destructor
  end type magdif_mesh

  type(magdif_mesh) :: mesh

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

contains

  !> Initialize magdif module
  subroutine magdif_init(bin_dir)
    use input_files, only: gfile
    use magdif_conf, only: conf, log, magdif_log, decorate_filename
    character(len = *), intent(in) :: bin_dir

    magdif_bin_dir = bin_dir
    log = magdif_log('-', conf%log_level, conf%quiet)

    ! only depends on config variables
    call read_mesh
    call load_magdata_in_symfluxcoord

    ! depends on mesh data
    call init_indices
    call cache_mesh_data
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
       call compute_Bn_nonres
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

  subroutine magdif_mesh_destructor(this)
    type(magdif_mesh), intent(inout) :: this
    if (allocated(this%deletions)) deallocate(this%deletions)
    if (allocated(this%additions)) deallocate(this%additions)
    if (allocated(this%refinement)) deallocate(this%refinement)
    if (allocated(this%m_res)) deallocate(this%m_res)
    if (allocated(this%res_ind)) deallocate(this%res_ind)
    if (allocated(this%kp_max)) deallocate(this%kp_max)
    if (allocated(this%kt_max)) deallocate(this%kt_max)
    if (allocated(this%kp_low)) deallocate(this%kp_low)
    if (allocated(this%kt_low)) deallocate(this%kt_low)
  end subroutine magdif_mesh_destructor

  !> Deallocates all previously allocated variables.
  subroutine magdif_cleanup
    use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp
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

  !> Allocates and initializes #kp_low, #kp_max, #kt_low and #kt_max based on the values
  !> of #magdif_config::nflux and #magdif_config::nkpol. Deallocation is done in
  !> magdif_cleanup().
  subroutine init_indices
    use magdif_conf, only: conf, log
    integer :: kf

    allocate(mesh%kp_max(conf%nflux))
    allocate(mesh%kt_max(conf%nflux))
    allocate(mesh%kp_low(conf%nflux + 1))
    allocate(mesh%kt_low(conf%nflux + 1))

    mesh%kp_max = conf%nkpol
    mesh%kt_max = 2 * conf%nkpol
    mesh%kt_max(1) = conf%nkpol

    mesh%kp_low(1) = 1
    do kf = 2, conf%nflux + 1
       mesh%kp_low(kf) = mesh%kp_low(kf-1) + mesh%kp_max(kf-1)
    end do
    mesh%kt_low(1) = 0
    do kf = 2, conf%nflux + 1
       mesh%kt_low(kf) = mesh%kt_low(kf-1) + mesh%kt_max(kf-1)
    end do

    write (log%msg, '("Number of points up to LCFS: ", i0)') mesh%kp_low(conf%nflux + 1)
    if (log%info) call log%write
    write (log%msg, '("Number of triangles up to LCFS: ", i0)') mesh%kt_low(conf%nflux + 1)
    if (log%info) call log%write
  end subroutine init_indices

  !> Reads mesh points and triangles.
  !>
  !> #mesh_mod::npoint and #mesh_mod::mesh_point are read directly from an unformatted
  !> #magdif_config::point_file, while #mesh_mod::ntri and #mesh_mod::mesh_element are
  !> read directly from an unformatted #magdif_config::tri_file. #presn, #bnflux, #bnphi,
  !> #bnflux_vac, #bnphi_vac, #j0phi, #jnphi and #jnflux are allocated and initialized to
  !> zero. Deallocation is done in magdif_cleanup().
  subroutine read_mesh
    use magdif_conf, only: conf, log
    use mesh_mod, only: npoint, ntri, mesh_point, mesh_element, mesh_element_rmp
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    integer(HID_T) :: h5id_magdif
    integer :: ktri
    integer, dimension(:, :), allocatable :: tri_node, adj_tri, adj_edge

    call h5_open('magdif.h5', h5id_magdif)
    call h5_get(h5id_magdif, 'mesh/nflux', conf%nflux)
    call h5_get(h5id_magdif, 'mesh/npoint', npoint)
    call h5_get(h5id_magdif, 'mesh/ntri', ntri)
    write (log%msg, '("nflux = ", i0, ", npoint = ", i0, ", ntri = ", i0)') &
         conf%nflux, npoint, ntri
    if (log%info) call log%write
    call fs%init(conf%nflux, .false.)
    call fs_half%init(conf%nflux, .true.)
    allocate(mesh_point(npoint))
    allocate(mesh_element(ntri))
    allocate(tri_node(3, ntri))
    allocate(adj_tri(3, ntri))
    allocate(adj_edge(3, ntri))
    call h5_get(h5id_magdif, 'mesh/m_res_min', mesh%m_res_min)
    call h5_get(h5id_magdif, 'mesh/m_res_max', mesh%m_res_max)
    call h5_get(h5id_magdif, 'mesh/node_R', mesh_point%rcoord)
    call h5_get(h5id_magdif, 'mesh/node_Z', mesh_point%zcoord)
    call h5_get(h5id_magdif, 'mesh/tri_node', tri_node)
    call h5_get(h5id_magdif, 'mesh/tri_node_F', mesh_element%knot_h)
    call h5_get(h5id_magdif, 'mesh/adj_tri', adj_tri)
    call h5_get(h5id_magdif, 'mesh/adj_edge', adj_edge)
    call h5_get(h5id_magdif, 'mesh/det', mesh_element%det_3)
    call h5_get(h5id_magdif, 'cache/fs/psi', fs%psi)
    call h5_get(h5id_magdif, 'cache/fs/rad', fs%rad)
    call h5_get(h5id_magdif, 'cache/fs_half/psi', fs_half%psi)
    call h5_get(h5id_magdif, 'cache/fs_half/rad', fs_half%rad)
    call h5_close(h5id_magdif)
    ! intermediate until mesh_mod is refactored into magdif_mesh
    do ktri = 1, ntri
       mesh_element(ktri)%i_knot = tri_node(:, ktri)
       mesh_element(ktri)%neighbour = adj_tri(:, ktri)
       mesh_element(ktri)%neighbour_edge = adj_edge(:, ktri)
    end do
    deallocate(tri_node)
    deallocate(adj_tri)
    deallocate(adj_edge)

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
    use mesh_mod, only: triangle, triangle_rmp, mesh_element, mesh_element_rmp, mesh_point
    use magdif_conf, only: conf
    use magdif_util, only: ring_centered_avg_coord
    integer :: ktri, kedge, ke, ke_adj, ktri_adj
    type(triangle) :: elem
    type(triangle_rmp) :: tri

    mesh%R_O = mesh_point(1)%rcoord
    mesh%Z_O = mesh_point(1)%zcoord
    mesh%R_min = minval(mesh_point(:mesh%kp_low(conf%nflux + 1))%rcoord)
    mesh%R_max = maxval(mesh_point(:mesh%kp_low(conf%nflux + 1))%rcoord)
    mesh%Z_min = minval(mesh_point(:mesh%kp_low(conf%nflux + 1))%zcoord)
    mesh%Z_max = maxval(mesh_point(:mesh%kp_low(conf%nflux + 1))%zcoord)

    mesh%nedge = (3 * mesh%kt_low(conf%nflux + 1) + mesh%kp_max(conf%nflux)) / 2
    allocate(mesh%edge_map2global(mesh%kt_low(conf%nflux + 1), 3))
    allocate(mesh%edge_map2ktri(mesh%nedge, 2))
    allocate(mesh%edge_map2ke(mesh%nedge, 2))
    mesh%edge_map2global = 0
    mesh%edge_map2ktri = 0
    mesh%edge_map2ke = 0
    kedge = 1

    do ktri = 1, mesh%kt_low(conf%nflux + 1)
       elem = mesh_element(ktri)
       tri%area = 0.5d0 * elem%det_3
       call get_labeled_edges(elem, tri%li, tri%lo, tri%lf, tri%ei, tri%eo, tri%ef, &
            tri%orient)
       call ring_centered_avg_coord(elem, tri%R_Omega, tri%Z_Omega)
       mesh_element_rmp(ktri) = tri

       do ke = 1, 3
          if (mesh%edge_map2global(ktri, ke) == 0) then
             ktri_adj = elem%neighbour(ke)
             ke_adj = elem%neighbour_edge(ke)
             if (ktri_adj > mesh%kt_low(conf%nflux + 1)) then
                mesh%edge_map2global(ktri, ke) = kedge
                mesh%edge_map2ktri(kedge, :) = [ktri, -1]
                mesh%edge_map2ke(kedge, :) = [ke, -1]
                kedge = kedge + 1
             else
                mesh%edge_map2global(ktri, ke) = kedge
                mesh%edge_map2global(ktri_adj, ke_adj) = kedge
                mesh%edge_map2ktri(kedge, :) = [ktri, ktri_adj]
                mesh%edge_map2ke(kedge, :) = [ke, ke_adj]
                kedge = kedge + 1
             end if
          end if
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
    use magdif_conf, only: conf, log
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
    use magdif_conf, only: log
    use magdif_util, only: imun
    use mesh_mod, only: ntri, mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    integer, intent(in) :: n
    real(dp), intent(in) :: rel_err
    character(len = *), intent(in) :: field_name

    integer :: ktri
    real(dp) :: div, abs_flux

    do ktri = 1, ntri
       abs_flux = sum(abs(pol_flux(ktri,:))) + abs(imun * n * tor_comp(ktri) * &
            mesh_element_rmp(ktri)%area)
       if (abs_flux > 0d0) then
          div = abs((sum(pol_flux(ktri,:)) + imun * n * tor_comp(ktri) * &
               mesh_element_rmp(ktri)%area)) / abs_flux
          if (div > rel_err) then
             write (log%msg, '("divergence of ", a, ' // &
                  '" above threshold in triangle ", i0, ": ", es24.16e3)') &
                  trim(field_name), ktri, div
              if (log%err) call log%write
              error stop
          end if
       end if
    end do
  end subroutine check_div_free

  subroutine check_redundant_edges(pol_quant, same_sign, name)
    use magdif_conf, only: log, cmplx_fmt
    complex(dp), intent(in) :: pol_quant(:,:)
    logical, intent(in) :: same_sign
    character(len = *), intent(in) :: name
    integer :: kedge, ktri, ktri_adj, ke, ke_adj
    logical :: inconsistent
    real(dp), parameter :: eps = epsilon(1d0), small = tiny(0d0)

    do kedge = 1, mesh%nedge
       ktri = mesh%edge_map2ktri(kedge, 1)
       ktri_adj = mesh%edge_map2ktri(kedge, 2)
       ke = mesh%edge_map2ke(kedge, 1)
       ke_adj = mesh%edge_map2ke(kedge, 2)
       if (ktri_adj <= 0) cycle
       inconsistent = .false.
       if (abs(real(pol_quant(ktri, ke))) < small) then
          inconsistent = inconsistent .or. abs(real(pol_quant(ktri_adj, ke_adj))) >= small
       else
          if (same_sign) then
             inconsistent = inconsistent .or. eps < abs(1d0 - &
                  real(pol_quant(ktri_adj, ke_adj)) / real(pol_quant(ktri, ke)))
          else
             inconsistent = inconsistent .or. eps < abs(1d0 + &
                  real(pol_quant(ktri_adj, ke_adj)) / real(pol_quant(ktri, ke)))
          end if
       end if
       if (abs(aimag(pol_quant(ktri, ke))) < small) then
          inconsistent = inconsistent .or. abs(aimag(pol_quant(ktri_adj, ke_adj))) >= small
       else
          if (same_sign) then
             inconsistent = inconsistent .or. eps < abs(1d0 - &
                  aimag(pol_quant(ktri_adj, ke_adj)) / aimag(pol_quant(ktri, ke)))
          else
             inconsistent = inconsistent .or. eps < abs(1d0 + &
                  aimag(pol_quant(ktri_adj, ke_adj)) / aimag(pol_quant(ktri, ke)))
          end if
       end if
       if (inconsistent) then
          write (log%msg, '("inconsistent redundant edges: ", ' // &
               'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ', ", ", ' // &
               'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ')') &
               trim(name), ktri, ke, pol_quant(ktri, ke), &
               trim(name), ktri_adj, ke_adj, pol_quant(ktri_adj, ke_adj)
          if (log%err) call log%write
          error stop
       end if
    end do
  end subroutine check_redundant_edges

  subroutine flux_func_cache_check
    use magdif_conf, only: log
    log%msg = 'checking flux_func_cache...'
    if (log%debug) call log%write
    write (log%msg, '("array bounds: fs%psi(", i0, ":", i0, "), ' // &
         ' fs%rad(", i0, ":", i0, "), fs_half%psi(", i0, ":", i0, "), ' // &
         'fs_half%rad(", i0, ":", i0, ")")') lbound(fs%psi, 1), ubound(fs%psi, 1), &
         lbound(fs%rad, 1), ubound(fs%rad, 1), lbound(fs_half%psi, 1), &
         ubound(fs_half%psi, 1), lbound(fs_half%rad, 1), ubound(fs_half%rad, 1)
    if (log%debug) call log%write
    write (log%msg, '("expected sign of psi''(r): ", sp, i0, ss)') equil%cocos%sgn_dpsi
    if (log%debug) call log%write
    write (log%msg, '(i0, " ordering violations for psi")') &
         count((fs%psi(1:) - fs_half%psi) * equil%cocos%sgn_dpsi <= 0d0) + &
         count([(fs_half%psi(1) - fs%psi(0)) * equil%cocos%sgn_dpsi] <= 0d0)
    if (log%debug) call log%write
    write (log%msg, '(i0, " ordering violations for rad")') &
         count(fs%rad(1:) <= fs_half%rad) + count([fs_half%rad(1)] <= [fs%rad(0)])
    if (log%debug) call log%write
  end subroutine flux_func_cache_check

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

  !> Allocates and computes the safety factor #q and #m_res.
  !>
  !> Also allocates #magdif_config::sheet_current_factor, to be read in via
  !> magdif_config::read_delayed_config() in magdif_init(). All deallocation is done in
  !> magdif_cleanup().
  subroutine compute_safety_factor
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, conf_arr, q_prof_flux, q_prof_rot, q_prof_geqdsk, &
         log, cmplx_fmt
    use magdif_util, only: flux_func
    use magdata_in_symfluxcoor_mod, only: psipol_max, psisurf, qsaf
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    integer :: kf, kt, ktri, m, kf_res
    type(triangle_rmp) :: tri
    type(flux_func) :: psi_interpolator
    real(dp), dimension(conf%nflux) :: abs_err

    select case (conf%q_prof)
    case (q_prof_flux)
       fs_half%q = 0d0
       do kf = 1, conf%nflux
          do kt = 1, mesh%kt_max(kf)
             ktri = mesh%kt_low(kf) + kt
             tri = mesh_element_rmp(ktri)
             fs_half%q(kf) = fs_half%q(kf) + B0phi_Omega(ktri) * tri%area
          end do
          fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
       end do
       call psi_interpolator%init(4, fs_half%psi)
       ! Lagrange polynomial extrapolation for values at separatrix and magnetic axis
       fs%q = [(psi_interpolator%interp(fs_half%q, fs%psi(kf)), kf = 0, conf%nflux)]
    case (q_prof_rot)
       ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
       ! load_magdata_in_symfluxcoord_ext divides by psipol_max
       call psi_interpolator%init(4, psisurf(1:) * psipol_max + fs%psi(0))
       ! Lagrange polynomial extrapolation for value at magnetic axis
       fs%q = [(psi_interpolator%interp(qsaf, fs%psi(kf)), kf = 0, conf%nflux)]
       fs_half%q = [(psi_interpolator%interp(qsaf, fs_half%psi(kf)), kf = 1, conf%nflux)]
    case (q_prof_geqdsk)
       fs%q = [(fluxvar%interp(equil%qpsi, fs%psi(kf)), kf = 0, conf%nflux)]
       fs_half%q = [(fluxvar%interp(equil%qpsi, fs_half%psi(kf)), kf = 1, conf%nflux)]
    case default
       write (log%msg, '("unknown q profile selection: ", i0)') conf%q_prof
       if (log%err) call log%write
       error stop
    end select

    allocate(mesh%m_res(conf%nflux))
    mesh%m_res = 0
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    mesh%res_ind = 0
    log%msg = 'resonance positions:'
    if (log%debug) call log%write
    do m = mesh%m_res_max, mesh%m_res_min, -1
       abs_err = [(abs(abs(fs_half%q(kf)) - dble(m) / dble(conf%n)), kf = 1, conf%nflux)]
       kf_res = minloc(abs_err, 1)
       mesh%res_ind(m) = kf_res
       mesh%m_res(kf_res) = m
       write (log%msg, '("m = ", i0, ", kf_res = ", i0, ' // &
            '", rho: ", f19.16, 2(" < ", f19.16))') m, kf_res, &
            [fs%rad(kf_res - 1), fs_half%rad(kf_res), fs%rad(kf_res)] / fs%rad(conf%nflux)
       if (log%debug) call log%write
    end do

    if (log%info) then
       log%msg = 'absolute poloidal mode number, sheet current factor'
       call log%write
       do m = conf%m_min, conf%m_max
          write (log%msg, '(i2, 1x, ' // cmplx_fmt // ')') m, &
               conf_arr%sheet_current_factor(m)
          call log%write
       end do
    end if
  end subroutine compute_safety_factor

  subroutine init_flux_variables
    use magdif_conf, only: conf
    integer :: kf

    ! initialize fluxvar with equidistant psi values
    call fluxvar%init(4, equil%psi_eqd)

    call compute_pres_prof
    call compute_safety_factor
    fs%F = [(fluxvar%interp(equil%fpol, fs%psi(kf)), kf = 0, conf%nflux)]
    fs%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs%psi(kf)), kf = 0, conf%nflux)]
    fs_half%F = [(fluxvar%interp(equil%fpol, fs_half%psi(kf)), kf = 1, conf%nflux)]
    fs_half%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs_half%psi(kf)), &
         kf = 1, conf%nflux)]
    call write_fluxvar
  end subroutine init_flux_variables

  subroutine compute_pres_prof
    use constants, only: ev2erg  ! orbit_mod.f90
    use magdif_conf, only: conf, pres_prof_eps, pres_prof_par, pres_prof_geqdsk, log
    integer :: kf
    real(dp) :: ddens_dpsi, dtemp_dpsi, psi_int, psi_ext

    !> Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
    real(dp) :: dens(0:conf%nflux)

    !> Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
    real(dp) :: temp(0:conf%nflux)

    dens = 0d0
    temp = 0d0
    if (equil%cocos%sgn_dpsi == -1) then
       psi_ext = minval(equil%psirz)
       psi_int = maxval(equil%psirz)
    else
       psi_ext = maxval(equil%psirz)
       psi_int = minval(equil%psirz)
    end if
    select case (conf%pres_prof)
    case (pres_prof_eps)
       ddens_dpsi = conf%dens_max / psi_int
       dtemp_dpsi = conf%temp_max / psi_int
       dens = (fs%psi - psi_ext) / psi_int * conf%dens_max + conf%dens_min
       temp = (fs%psi - psi_ext) / psi_int * conf%temp_max + conf%temp_min
       write (log%msg, '("temp@axis: ", es24.16e3, ", dens@axis: ", es24.16e3)') &
            temp(0), dens(0)
       if (log%info) call log%write
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psi_ext) / psi_int * conf%dens_max + conf%dens_min
       temp(1:) = (fs_half%psi - psi_ext) / psi_int * conf%temp_max + conf%temp_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_par)
       ddens_dpsi = (conf%dens_max - conf%dens_min) / (psi_int - psi_ext)
       dtemp_dpsi = (conf%temp_max - conf%temp_min) / (psi_int - psi_ext)
       dens = (fs%psi - psi_ext) / (psi_int - psi_ext) * (conf%dens_max - conf%dens_min) &
            + conf%dens_min
       temp = (fs%psi - psi_ext) / (psi_int - psi_ext) * (conf%temp_max - conf%temp_min) &
            + conf%temp_min
       fs%p = dens * temp * ev2erg
       fs%dp_dpsi = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
       dens(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
            (conf%dens_max - conf%dens_min) + conf%dens_min
       temp(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
            (conf%temp_max - conf%temp_min) + conf%temp_min
       fs_half%p = dens(1:) * temp(1:) * ev2erg
       fs_half%dp_dpsi = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
    case (pres_prof_geqdsk)
       fs%p = [(fluxvar%interp(equil%pres, fs%psi(kf)), kf = 0, conf%nflux)]
       fs%dp_dpsi = [(fluxvar%interp(equil%pprime, fs%psi(kf)), kf = 0, conf%nflux)]
       fs_half%p = [(fluxvar%interp(equil%pres, fs_half%psi(kf)), kf = 1, conf%nflux)]
       fs_half%dp_dpsi = [(fluxvar%interp(equil%pprime, fs_half%psi(kf)), &
            kf = 1, conf%nflux)]
    case default
       write (log%msg, '("unknown pressure profile selection", i0)') conf%pres_prof
       if (log%err) call log%write
       error stop
    end select
  end subroutine compute_pres_prof


  subroutine cache_equilibrium_field
    use mesh_mod, only: knot, triangle, mesh_point, mesh_element
    use magdif_conf, only: conf, longlines
    use magdif_util, only: ring_centered_avg_coord
    real(dp) :: r, z, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, &
         dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
    integer :: kf, kt, ktri, ke, fid
    type(triangle) :: elem
    type(knot) :: base, tip
    real(dp) :: n_r, n_z

    open(newunit = fid, file = 'plot_B0.dat', recl = longlines, status = 'replace')
    do kf = 1, conf%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
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
          write (fid, '(5(1x, es24.16e3))') r, z, Br, Bz, Bp
       end do
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
    use constants, only: pi  ! orbit_mod.f90
    use mesh_mod, only: triangle_rmp, mesh_element_rmp, mesh_point
    use magdif_conf, only: conf, curr_prof_ps, curr_prof_rot, curr_prof_geqdsk, longlines
    use magdif_util, only: clight
    integer :: kf, kt, ktri, fid
    real(dp) :: r, z
    real(dp) :: Btor2
    real(dp), dimension(conf%nflux) :: B2avg, B2avg_half
    real(dp) :: plot_j0phi
    type(triangle_rmp) :: tri

    open(newunit = fid, file = conf%j0phi_file, recl = longlines, status = 'replace')
    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, conf%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0r(ktri, tri%ef) ** 2 + &
               B0phi(ktri, tri%ef) ** 2 + B0z(ktri, tri%ef) ** 2
          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0r(ktri, tri%ei) ** 2 + &
               B0phi(ktri, tri%ei) ** 2 + B0z(ktri, tri%ei) ** 2
       end do
       B2avg(kf) = B2avg(kf) / mesh%kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / mesh%kt_max(kf)

       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)

          r = sum(mesh_point(tri%lf(:))%rcoord) * 0.5d0
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
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
             if (kf > 1 .and. .not. tri%orient) then
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf-1) * (1d0 - &
                     Btor2 / B2avg(kf-1))
             else
                j0phi(ktri, tri%ef) = clight * r * fs%dp_dpsi(kf) * (1d0 - &
                     Btor2 / B2avg(kf))
             end if
          end select

          r = sum(mesh_point(tri%li(:))%rcoord) * 0.5d0
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
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
          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
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

          select case (conf%curr_prof)
          case (curr_prof_geqdsk)
             plot_j0phi = clight * (fs_half%dp_dpsi(kf) * tri%r_Omega + &
                  0.25d0 / pi * fs_half%FdF_dpsi(kf) / tri%r_Omega)
          case (curr_prof_rot)
             plot_j0phi = j0phi_ampere(tri%r_Omega, tri%z_Omega)
          case (curr_prof_ps)
             Btor2 = B0phi_Omega(ktri) ** 2
             plot_j0phi = clight * tri%r_Omega * fs_half%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg_half(kf))
          end select

          write (fid, '(4(1x, es24.16e3))') &
               j0phi(ktri, 1), j0phi(ktri, 2), j0phi(ktri, 3), plot_j0phi
       end do
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
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, longlines
    use magdif_util, only: clight
    integer :: kf, kt, fid_amp, fid_gs, fid_prof
    real(dp) :: cmp_gradp, cmp_amp, cmp_gs, theta, R, Z, dum, B0_R, B0_phi, B0_Z, &
         dB0R_dZ, dB0phi_dR, dB0phi_dZ, dB0Z_dR, J0_R, J0_phi, J0_Z, grad_psi(3)

    open(newunit = fid_prof, file = 'cmp_prof.dat', recl = longlines, status = 'replace')
    open(newunit = fid_amp, file = 'j0_amp.dat', recl = longlines, status = 'replace')
    open(newunit = fid_gs, file = 'j0_gs.dat', recl = longlines, status = 'replace')
    do kf = 1, conf%nflux
       do kt = 1, mesh%kt_max(kf)
          theta = (dble(kt) - 0.5d0) / dble(mesh%kt_max(kf)) * 2d0 * pi
          ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
          call magdata_in_symfluxcoord_ext(2, dum, fs_half%psi(kf) - fs%psi(0), &
               theta, dum, dum, dum, dum, dum, R, dum, dum, Z, dum, dum)
          call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dum, dum, dB0R_dZ, dB0phi_dR, &
               dum, dB0phi_dZ, dB0Z_dR, dum, dum)
          ! left-hand side of iMHD force balance
          grad_psi = [R * B0_Z, 0d0, -R * B0_R]
          cmp_gradp = fs_half%dp_dpsi(kf) * dot_product(grad_psi, grad_psi) / &
               norm2(grad_psi)
          ! current density via Grad-Shafranov equation
          J0_R = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * B0_R
          J0_Z = 0.25d0 / pi * clight * fs_half%FdF_dpsi(kf) / fs_half%F(kf) * B0_Z
          J0_phi = clight * (fs_half%dp_dpsi(kf) * R + &
               0.25d0 / pi * fs_half%FdF_dpsi(kf) / R)
          write (fid_gs, '(3(1x, es24.16e3))') J0_R, J0_phi, J0_Z
          cmp_gs = dot_product([J0_phi * B0_Z - J0_Z * B0_phi, J0_Z * B0_R - J0_R * B0_Z, &
               J0_R * B0_phi - J0_phi * B0_R], grad_psi) / norm2(grad_psi) / clight
          ! current density via Ampere's equation
          J0_R = 0.25d0 / pi * clight * (-dB0phi_dZ)
          J0_phi = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          J0_Z = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
          write (fid_amp, '(3(1x, es24.16e3))') J0_R, J0_phi, J0_Z
          cmp_amp = dot_product([J0_phi * B0_Z - J0_Z * B0_phi, J0_Z * B0_R - J0_R * B0_Z, &
               J0_R * B0_phi - J0_phi * B0_R], grad_psi) / norm2(grad_psi) / clight
          write (fid_prof, '(3(1x, es24.16e3))') cmp_gradp, cmp_amp, cmp_gs
       end do
    end do
    close(fid_prof)
    close(fid_amp)
    close(fid_gs)
  end subroutine check_curr0


  subroutine check_safety_factor
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf
    use magdata_in_symfluxcoor_mod, only: nlabel, psipol_max, psisurf, rbeg, qsaf
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    integer :: kf, kt, ktri, fid
    type(triangle_rmp) :: tri
    real(dp), allocatable :: q(:)

    allocate(q(conf%nflux))
    q = 0d0
    do kf = 1, conf%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          q(kf) = q(kf) + B0phi_Omega(ktri) * tri%area
       end do
       q(kf) = q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    open(newunit = fid, file = 'check_q_step.dat', status = 'replace')
    do kf = 1, conf%nflux
       write (fid, '(3(1x, es24.16e3))') (fs_half%psi(kf) - fs%psi(0)) / &
            (fs%psi(conf%nflux) - fs%psi(0)), fs_half%rad(kf), q(kf)
    end do
    close(fid)
    deallocate(q)
    allocate(q(nlabel))
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    q = [(fluxvar%interp(equil%qpsi, psisurf(kf) * psipol_max + fs%psi(0)), &
         kf = 1, nlabel)]
    open(newunit = fid, file = 'check_q_cont.dat', status = 'replace')
    do kf = 1, nlabel
       write (fid, '(4(1x, es24.16e3))') psisurf(kf), rbeg(kf), qsaf(kf), q(kf)
    end do
    close(fid)
    deallocate(q)
  end subroutine check_safety_factor


  !> Computes pressure perturbation #presn from equilibrium quantities and #bnflux.
  subroutine compute_presn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use magdif_conf, only: conf, log
    use magdif_util, only: assemble_sparse, imun
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
    use mesh_mod, only: triangle
    use magdif_conf, only: log
    type(triangle), intent(in) :: elem
    integer, dimension(:), intent(out) :: li, lo, lf
    integer, intent(out) :: ei, eo, ef
    logical, intent(out) :: orient
    integer, dimension(3) :: i_knot_diff
    integer :: knot_i, knot_o, knot_f
    integer :: i1, i2
    logical :: closing_loop

    if (2 /= size(li)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(li)', 'actual size(li)', &
            2, size(li))
       if (log%err) call log%write
       error stop
    end if
    if (2 /= size(lo)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(lo)', 'actual size(lo)', &
            2, size(lo))
       if (log%err) call log%write
       error stop
    end if
    if (2 /= size(lf)) then
       call log%msg_arg_size('get_labeled_edges', 'expected size(lf)', 'actual size(lf)', &
            2, size(lf))
       if (log%err) call log%write
       error stop
    end if
    log%msg = 'cannot find correct label for triangle edges'

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
       if (log%err) call log%write
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
       li = [elem%i_knot(knot_f), elem%i_knot(knot_o)]
       lo = [elem%i_knot(knot_i), elem%i_knot(knot_f)]
       lf = [elem%i_knot(knot_o), elem%i_knot(knot_i)]
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
       li = [elem%i_knot(knot_o), elem%i_knot(knot_f)]
       lo = [elem%i_knot(knot_f), elem%i_knot(knot_i)]
       lf = [elem%i_knot(knot_i), elem%i_knot(knot_o)]
    else
       if (log%err) call log%write
       error stop
    end if
  end subroutine get_labeled_edges


  !> Computes current perturbation #jnflux and #jnphi from equilibrium quantities,
  !> #presn, #bnflux and #bnphi.
  !>
  !> This subroutine computes the fluxes through each triangle, separately for each flux
  !> surface. The result is written to #magdif_conf::currn_file.
  subroutine compute_currn
    use sparse_mod, only: sparse_solve, sparse_matmul
    use mesh_mod, only: triangle, triangle_rmp, mesh_element, mesh_element_rmp, mesh_point
    use magdif_conf, only: conf, log, decorate_filename
    use magdif_util, only: imun, clight, assemble_sparse
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

  subroutine avg_flux_on_quad(pol_flux, tor_comp)
    use magdif_conf, only: conf
    use magdif_util, only: imun
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)

    integer :: kf, kt, ktri1, ktri2
    complex(dp) :: tor_flux_avg, tor_flux_diff
    type(triangle_rmp) :: tri1, tri2

    do kf = 2, conf%nflux
       do kt = 1, mesh%kt_max(kf), 2
          ktri1 = mesh%kt_low(kf) + kt
          ktri2 = mesh%kt_low(kf) + kt + 1
          tri1 = mesh_element_rmp(ktri1)
          tri2 = mesh_element_rmp(ktri2)
          tor_flux_avg = 0.5d0 * (tor_comp(ktri2) * tri2%area + &
               tor_comp(ktri1) * tri1%area)
          tor_flux_diff = 0.5d0 * (tor_comp(ktri2) * tri2%area - &
               tor_comp(ktri1) * tri1%area)
          tor_comp(ktri1) = tor_flux_avg / tri1%area
          pol_flux(ktri1, tri1%eo) = pol_flux(ktri1, tri1%eo) - imun * conf%n * tor_flux_diff
          tor_comp(ktri2) = tor_flux_avg / tri2%area
          pol_flux(ktri2, tri2%ei) = pol_flux(ktri2, tri2%ei) + imun * conf%n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  subroutine compute_Bn_nonres
    use mesh_mod, only: knot, triangle_rmp, mesh_point, mesh_element_rmp
    use magdif_conf, only: conf
    use magdif_util, only: imun
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri
    type(knot) :: base, tip
    real(dp) :: r, lr, lz
    complex(dp) :: Bnpsi

    do kf = 1, conf%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          ! use midpoint of edge f
          base = mesh_point(tri%lf(1))
          tip = mesh_point(tri%lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord
          Bnpsi = -equil%rcentr * B0phi(ktri, tri%ef) / r
          Bnflux(ktri, tri%ef) = Bnpsi * (lr ** 2 + lz ** 2) / &
               (B0r(ktri, tri%ef) * lr + B0z(ktri, tri%ef) * lz)
          Bnphi(ktri) = imun / conf%n * Bnflux(ktri, tri%ef) / tri%area
          Bnflux(ktri, tri%ei) = (0d0, 0d0)
          Bnflux(ktri, tri%eo) = (0d0, 0d0)
       end do
    end do
    if (conf%quad_avg) call avg_flux_on_quad(Bnflux, Bnphi)
    call check_div_free(Bnflux, Bnphi, conf%n, conf%rel_err_Bn, 'non-resonant B_n')
    call check_redundant_edges(Bnflux, .false., 'non-resonant B_n')
  end subroutine compute_Bn_nonres

  subroutine interp_RT0(ktri, pol_flux, r, z, pol_comp_r, pol_comp_z)
    use mesh_mod, only: knot, triangle, mesh_point, mesh_element
    integer, intent(in) :: ktri
    complex(dp), intent(in) :: pol_flux(:,:)
    real(dp), intent(in) :: r, z
    complex(dp), intent(out) :: pol_comp_r, pol_comp_z
    type(triangle) :: elem
    type(knot) :: node(3)

    elem = mesh_element(ktri)
    node = mesh_point(elem%i_knot(:))
    ! edge 1 lies opposite to knot 3, etc.
    pol_comp_r = 1d0 / elem%det_3 / r * ( &
         pol_flux(ktri, 1) * (r - node(3)%rcoord) + &
         pol_flux(ktri, 2) * (r - node(1)%rcoord) + &
         pol_flux(ktri, 3) * (r - node(2)%rcoord))
    pol_comp_z = 1d0 / elem%det_3 / r * ( &
         pol_flux(ktri, 1) * (z - node(3)%zcoord) + &
         pol_flux(ktri, 2) * (z - node(1)%zcoord) + &
         pol_flux(ktri, 3) * (z - node(2)%zcoord))
  end subroutine interp_RT0

  pure function jacobian(kf, kt, r) result(metric_det)
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: r
    real(dp) :: metric_det
    metric_det = equil%cocos%sgn_dpsi * fs_half%q(kf) * r / B0phi_Omega(mesh%kt_low(kf) + kt)
  end function jacobian

  subroutine write_vector_plot(pol_flux, tor_comp, outfile)
    use magdif_conf, only: conf, longlines
    use mesh_mod, only: ntri, triangle_rmp, mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: k, kf, kt, ktri, n_cutoff, fid
    type(triangle_rmp) :: tri
    complex(dp) :: pol_comp_r, pol_comp_z, dens_psi_contravar, proj_theta_covar
    real(dp) :: r, z

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    if (conf%nonres) then
       n_cutoff = conf%nflux - 1
    else
       n_cutoff = conf%nflux
    end if
    do kf = 1, n_cutoff
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          r = tri%R_Omega
          z = tri%Z_Omega
          call interp_RT0(ktri, pol_flux, r, z, pol_comp_r, pol_comp_z)
          ! projection to contravariant psi component
          dens_psi_contravar = (pol_comp_r * B0z_Omega(ktri) - &
               pol_comp_z * B0r_Omega(ktri)) * r * jacobian(kf, kt, r)
          ! projection to covariant theta component
          proj_theta_covar = equil%cocos%sgn_dpsi * (pol_comp_r * B0r_Omega(ktri) + &
               pol_comp_z * B0z_Omega(ktri)) * jacobian(kf, kt, r)
          write (fid, '(12(1x, es24.16e3))') r, z, pol_comp_r, pol_comp_z, &
               tor_comp(ktri), dens_psi_contravar, proj_theta_covar
       end do
    end do
    r = equil%rmaxis
    z = equil%zmaxis
    do ktri = mesh%kt_low(n_cutoff+1) + 1, ntri
       write (fid, '(12(1x, es24.16e3))') r, z, (0d0, k = 1, 10)
    end do
    close(fid)
  end subroutine write_vector_plot

  ! based on http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
  function point_in_triangle(ktri, R, Z, thickness) result(probably)
    use mesh_mod, only: triangle, mesh_element, mesh_point
    integer, intent(in) :: ktri
    real(dp), intent(in) :: R, Z
    real(dp), intent(in), optional :: thickness
    logical :: probably
    real(dp), dimension(1:4) :: node_R, node_Z, dist_R, dist_Z
    real(dp), dimension(1:3) :: edge_R, edge_Z, edge_2, dist_2, dotprod
    type(triangle) :: elem

    probably = .false.
    if (ktri <= 0) return
    elem = mesh_element(ktri)
    node_R = mesh_point([elem%i_knot, elem%i_knot(1)])%rcoord
    node_Z = mesh_point([elem%i_knot, elem%i_knot(1)])%zcoord
    dist_R = R - node_R
    dist_Z = Z - node_Z
    edge_R = node_R(2:4) - node_R(1:3)
    edge_Z = node_Z(2:4) - node_Z(1:3)
    edge_2 = edge_R ** 2 + edge_Z ** 2
    ! perp_R = edge_Z, perp_Z = -edge_R
    dotprod = edge_Z * dist_R(1:3) - edge_R * dist_Z(1:3)
    probably = all(dotprod <= 0d0)
    if (probably .or. .not. present(thickness)) return
    ! reuse dotprod as parameter of edge vector in linear equation
    dotprod = edge_R * dist_R(1:3) + edge_Z * dist_Z(1:3)
    where (dotprod < 0d0)
       dist_2 = dist_R(1:3) ** 2 + dist_Z(1:3) ** 2
    elsewhere (dotprod > edge_2)
       dist_2 = dist_R(2:4) ** 2 + dist_Z(2:4) ** 2
    elsewhere
       dist_2 = dist_R(1:3) ** 2 + dist_Z(1:3) ** 2 - dotprod ** 2 / edge_2
    end where
    probably = any(dist_2 < thickness ** 2)
  end function point_in_triangle

  function point_location(r, z) result(location)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf
    use magdif_util, only: interp_psi_pol, binsearch
    use mesh_mod, only: mesh_point
    real(dp), intent(in) :: r, z
    integer :: location

    integer :: kf, kq, k, k_max, ktri, candidates(6)
    real(dp) :: psi, pol_frac, pol_offset, pol_interval(0:conf%nkpol), thickness

    location = -1
    if (R < mesh%R_min .or. R > mesh%R_max .or. Z < mesh%Z_min .or. Z > mesh%Z_max) return
    location = -2
    psi = interp_psi_pol(r, z)
    if (equil%cocos%sgn_dpsi == +1) then
       if (psi > fs%psi(conf%nflux)) return
    else
       if (psi < fs%psi(conf%nflux)) return
    end if
    call binsearch(fs%psi, lbound(fs%psi, 1), psi, kf)
    location = -3

    pol_interval = 0d0
    pol_interval(0:mesh%kp_max(kf)-1) = 0.5d0 / pi * atan2( &
         mesh_point((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))%zcoord - mesh%Z_O, &
         mesh_point((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))%rcoord - mesh%R_O)
    pol_offset = pol_interval(0)
    pol_interval = pol_interval - pol_offset
    pol_interval(mesh%kp_max(kf)) = 1d0
    where (pol_interval < 0d0) pol_interval = pol_interval + 1d0
    where (pol_interval > 1d0) pol_interval = pol_interval - 1d0
    pol_frac = 0.5d0 / pi * atan2(Z - mesh%Z_O, R - mesh%R_O) - pol_offset
    if (pol_frac < 0d0) pol_frac = pol_frac + 1d0
    if (pol_frac > 1d0) pol_frac = pol_frac - 1d0
    call binsearch(pol_interval, lbound(pol_interval, 1), pol_frac, kq)

    ! Triangle edges do not lie exactly on flux surfaces, so we include the two adjacent
    ! triangle strips in the search. The candidates are ordered by decreasing likelihood
    ! of being the correct guess, i.e., the current loop, the outer loop and the inner
    ! loop, or filler if any of these is not applicable.
    if (kf == 1) then
       k_max = 3
       candidates = [mesh%kt_low(kf) + kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            -1, -1, -1]
    elseif (kf == 2) then
       k_max = 5
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            mesh%kt_low(kf - 1) + kq, -1]
    elseif (kf == conf%nflux) then
       k_max = 4
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf - 1) + 2 * kq - 1, mesh%kt_low(kf - 1) + 2 * kq, &
            -1, -1]
    else
       k_max = 6
       candidates = [mesh%kt_low(kf) + 2 * kq - 1, mesh%kt_low(kf) + 2 * kq, &
            mesh%kt_low(kf + 1) + 2 * kq - 1, mesh%kt_low(kf + 1) + 2 * kq, &
            mesh%kt_low(kf - 1) + 2 * kq - 1, mesh%kt_low(kf - 1) + 2 * kq]
    end if

    thickness = fs%rad(conf%nflux) * sqrt(epsilon(1d0)) * 8d0
    do k = 1, k_max
       ktri = candidates(k)
       if (point_in_triangle(ktri, R, Z, thickness)) then
          location = ktri
          exit
       end if
    end do
  end function point_location

  subroutine write_vector_plot_rect(pol_flux, tor_comp, outfile)
    use magdif_conf, only: conf, longlines
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: kw, kh, ktri, fid
    type(triangle_rmp) :: tri
    complex(dp) :: comp_R, comp_Z, comp_phi
    real(dp) :: R, Z

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do kw = 1, equil%nw
       do kh = 1, equil%nh
          R = equil%R_eqd(kw)
          Z = equil%Z_eqd(kh)
          ktri = point_location(R, Z)
          if (ktri > mesh%kt_low(1) .and. ktri <= mesh%kt_low(conf%nflux + 1)) then
             tri = mesh_element_rmp(ktri)
             call interp_RT0(ktri, pol_flux, R, Z, comp_R, comp_Z)
             comp_phi = tor_comp(ktri)
          else
             comp_R = 0d0
             comp_Z = 0d0
             comp_phi = 0d0
          end if
          write (fid, '(i6, 12(1x, es24.16e3))') ktri, R, Z, comp_R, comp_Z, comp_phi
       end do
    end do
    close(fid)
  end subroutine write_vector_plot_rect

  subroutine write_vector_dof(pol_flux, tor_comp, outfile)
    use magdif_conf, only: conf, longlines
    use mesh_mod, only: mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile

    integer :: ktri, fid

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do ktri = 1, mesh%kt_low(conf%nflux + 1)
       write (fid, '(8(1x, es23.15e3))') pol_flux(ktri, :), &
            tor_comp(ktri) * mesh_element_rmp(ktri)%area
    end do
    close(fid)
  end subroutine write_vector_dof

  !> Real and imaginary part of \p scalar_dof (e.g. #presn) are written, in that order,
  !> to \p outfile (e.g. #magdif_conf::presn_file), where line number corresponds to the
  !> knot index in #mesh_mod::mesh_point.
  subroutine write_scalar_dof(scalar_dof, outfile)
    use magdif_conf, only: conf, longlines
    complex(dp), intent(in) :: scalar_dof(:)
    character(len = *), intent(in) :: outfile

    integer :: kpoint, fid

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do kpoint = 1, mesh%kp_low(conf%nflux + 1)
       write (fid, '(2(1x, es24.16e3))') scalar_dof(kpoint)
    end do
    close(fid)
  end subroutine write_scalar_dof

  subroutine write_fluxvar
    use magdif_conf, only: conf, longlines
    use mesh_mod, only: mesh_point
    integer :: kf, fid
    real(dp) :: rho_r, rho_z

    open(newunit = fid, file = conf%fluxvar_file, recl = longlines, status = 'replace')
    write (fid, '(5(1x, es24.16e3))') 0d0, fs%psi(0), fs%q(0), fs%p(0), fs%dp_dpsi(0)
    do kf = 1, conf%nflux
       rho_r = mesh_point(mesh%kp_low(kf) + 1)%rcoord - mesh_point(1)%rcoord
       rho_z = mesh_point(mesh%kp_low(kf) + 1)%zcoord - mesh_point(1)%zcoord
       write (fid, '(5(1x, es24.16e3))') &
            hypot(rho_r, rho_z), fs%psi(kf), fs%q(kf), fs%p(kf), fs%dp_dpsi(kf)
    end do
    close(fid)
  end subroutine write_fluxvar

  subroutine write_poloidal_modes(pol_flux, tor_comp, outfile)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, longlines, decorate_filename
    use magdif_util, only: imun, clight, bent_cyl2straight_cyl
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
