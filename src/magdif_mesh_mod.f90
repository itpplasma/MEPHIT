module magdif_mesh

  use iso_fortran_env, only: dp => real64
  use magdif_util, only: g_eqdsk, flux_func

  implicit none

  private

  public :: equil, fluxvar, flux_func_cache, fs, fs_half, mesh_t, mesh, B0r, B0phi, B0z, &
       B0r_Omega, B0phi_Omega, B0z_Omega, B0flux, j0phi, coord_cache, sample_polmodes, &
       coord_cache_ext, coord_cache_ext_init, compute_sample_Ipar, coord_cache_ext_deinit, &
       coord_cache_ext_write, coord_cache_ext_read, &
       flux_func_cache_init, flux_func_cache_check, flux_func_cache_destructor, generate_mesh, &
       compute_resonance_positions, refine_eqd_partition, refine_resonant_surfaces, write_kilca_convexfile, &
       create_mesh_points, init_indices, common_triangles, &
       connect_mesh_points, get_labeled_edges, write_mesh_cache, read_mesh_cache, &
       magdif_mesh_destructor, init_flux_variables, compute_pres_prof_eps, compute_pres_prof_par, &
       compute_pres_prof_geqdsk, compute_safety_factor_flux, compute_safety_factor_rot, &
       compute_safety_factor_geqdsk, check_safety_factor, cache_equilibrium_field, &
       compute_j0phi, check_curr0, point_location, point_in_triangle

  type(g_eqdsk) :: equil

  !> Structure containing flux functions evaluated at a specific flux surface, indicated
  !> by a common array index. For details see flux_func_cache_init().
  type :: flux_func_cache
     private
     integer :: nflux

     !> Magnetic flux surface label \f$ \psi \f$ in maxwell.
     !>
     !> \f$ \psi \f$ is the disc poloidal flux divided by \f$ 2 \pi \f$. Its sign is
     !> positive and its magnitude is growing in the radially inward direction.
     real(dp), dimension(:), allocatable, public :: psi

     !> Minor radius \f$ r \f$ in centimeter.
     real(dp), dimension(:), allocatable, public :: rad

     real(dp), dimension(:), allocatable, public :: F

     !> Equilibrium pressure \f$ p_{0} \f$ in barye.
     real(dp), dimension(:), allocatable, public :: p

     real(dp), dimension(:), allocatable, public :: FdF_dpsi

     !> Derivative of equilibrium pressure w.r.t. flux surface label,
     !> \f$ p_{0}'(\psi) \f$, in barye per maxwell.
     real(dp), dimension(:), allocatable, public :: dp_dpsi

     !> Safety factor \f$ q \f$ (dimensionless).
     real(dp), dimension(:), allocatable, public :: q
   contains
     procedure :: init => flux_func_cache_init
     final :: flux_func_cache_destructor
  end type flux_func_cache

  type(flux_func) :: fluxvar
  type(flux_func_cache) :: fs
  type(flux_func_cache) :: fs_half

  type :: mesh_t

     !> R coordinate of the O point in cm.
     real(dp) :: R_O

     !> Z coordinate of the O point in cm.
     real(dp) :: Z_O

     !> R coordinate of the X point in cm.
     real(dp) :: R_X

     !> Z coordinate of the X point in cm.
     real(dp) :: Z_X

     !> Geometric poloidal angle where symmetry flux poloidal angle is zero.
     !>
     !> When theta0_at_xpoint is true, this is atan2(Z_X - Z_O, R_X - R_O);
     !> otherwise, this is zero.
     real(dp) :: theta0

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

     !> Number of triangle nodes.
     integer :: npoint

     !> Number of triangles.
     integer :: ntri

     !> Number of triangle edges.
     integer :: nedge

     !> Toroidal mode number. May differ from #magdif_conf::magdif_config#n due to large
     !> aspect ratio scaling.
     integer :: n

     !> Minimal poloidal mode number in resonance.
     integer :: m_res_min

     !> Maximal poloidal mode number in resonance
     integer :: m_res_max

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

     !> Poloidal flux at resonance position corresponding to a poloidal mode (given as
     !> array index).
     real(dp), allocatable :: psi_res(:)

     !> Normalized minor radius (along line connecting X point and O point) at resonance
     !> position corresponding to a poloidal mode (given as array index).
     real(dp), allocatable :: rad_norm_res(:)

     !> Poloidal modes that are expected to be in resonance. This might be different from
     !> m_res_min:m_res_max for specially constructed vacuum perturbation fields.
     integer, allocatable :: res_modes(:)

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

     !> Symmetry flux poloidal angle where geometric poloidal angle is zero.
     real(dp), allocatable :: theta_offset(:)

     real(dp), allocatable :: node_R(:)
     real(dp), allocatable :: node_Z(:)
     real(dp), allocatable :: node_theta_flux(:)
     real(dp), allocatable :: node_theta_geom(:)

     integer, allocatable :: tri_node(:, :)
     integer, allocatable :: tri_node_F(:)

     integer, allocatable :: li(:, :)
     integer, allocatable :: lo(:, :)
     integer, allocatable :: lf(:, :)
     integer, allocatable :: ei(:)
     integer, allocatable :: eo(:)
     integer, allocatable :: ef(:)
     logical, allocatable :: orient(:)

     integer, allocatable :: adj_tri(:, :)
     integer, allocatable :: adj_edge(:, :)

     integer, allocatable :: edge_map2global(:, :)
     integer, allocatable :: edge_map2ktri(:, :)
     integer, allocatable :: edge_map2ke(:, :)

     !> Surface integral Jacobian used for normalization of poloidal modes in GPEC
     complex(dp), allocatable :: gpec_jacfac(:, :)

     real(dp), allocatable :: area(:)
     real(dp), allocatable :: R_Omega(:)
     real(dp), allocatable :: Z_Omega(:)

   contains
     final :: magdif_mesh_destructor
  end type mesh_t

  type(mesh_t) :: mesh

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

  !> Physical toroidal component of equilibrium current \f$ j_{0 (\phi)} \f$ in
  !> statampere cm^-2.
  !>
  !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
  !> refers to the triangle and the indexing scheme is the same as for
  !> #mesh_mod::mesh_element. The second index refers to the edge and can be interpreted
  !> by get_labeled_edges().
  real(dp), allocatable :: j0phi(:,:)

  type :: coord_cache
     integer :: n
     integer, allocatable, dimension(:) :: ktri
     real(dp), allocatable, dimension(:) :: R, Z, psi, theta, sqrt_g, B0_R, B0_Z, &
          dR_dtheta, dZ_dtheta
  end type coord_cache

  type(coord_cache) :: sample_polmodes

  type, extends(coord_cache) :: coord_cache_ext
     integer :: nrad, npol, m
     real(dp), allocatable, dimension(:) :: rad, q, dq_dpsi, dR_dpsi, dZ_dpsi, &
          d2R_dpsi_dtheta, d2Z_dpsi_dtheta, B0_phi, dB0R_dR, dB0R_dZ, dB0phi_dR, &
          dB0phi_dZ, dB0Z_dR, dB0Z_dZ, B0_2
  end type coord_cache_ext

contains

  !> Set up arrays of cached values of flux functions.
  !>
  !> @nflux number of flux surfaces
  !> @half_step values are taken at flux surfaces (false) or between flux surfaces (true)
  !>
  !> For full-grid quantities, values are taken on flux surfaces with indices running
  !> from 0 to \p nflux, i.e. from the magnetic axis to the separatrix. An exception is
  !> made for \psi, where the index runs up to \p nflux +1. This value is extrapolated for
  !> finite differences in magdif::compute_presn() and magdif::compute_bn_nonres().
  !> For half-grid quantities, values are taken between two flux surfaces with indices
  !> running from 1 to \p nflux, i.e. from the triangle strip surrounding the magnetic
  !> axis to the triangle strip just inside the separatrix.
  subroutine flux_func_cache_init(this, nflux, half_step)
    class(flux_func_cache), intent(inout) :: this
    integer, intent(in) :: nflux
    logical, intent(in) :: half_step

    call flux_func_cache_destructor(this)
    if (half_step) then
       allocate(this%psi(nflux))
       allocate(this%rad(nflux))
       allocate(this%F(nflux))
       allocate(this%p(nflux))
       allocate(this%FdF_dpsi(nflux))
       allocate(this%dp_dpsi(nflux))
       allocate(this%q(nflux))
    else
       allocate(this%psi(0:nflux))
       allocate(this%rad(0:nflux))
       allocate(this%F(0:nflux))
       allocate(this%p(0:nflux))
       allocate(this%FdF_dpsi(0:nflux))
       allocate(this%dp_dpsi(0:nflux))
       allocate(this%q(0:nflux))
    end if
    this%psi = 0d0
    this%rad = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
  end subroutine flux_func_cache_init

  subroutine flux_func_cache_destructor(this)
    type(flux_func_cache), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%F)) deallocate(this%F)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%FdF_dpsi)) deallocate(this%FdF_dpsi)
    if (allocated(this%dp_dpsi)) deallocate(this%dp_dpsi)
    if (allocated(this%q)) deallocate(this%q)
  end subroutine flux_func_cache_destructor

  subroutine flux_func_cache_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(flux_func_cache), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', &
         cache%psi, lbound(cache%psi), ubound(cache%psi), unit = 'Mx', &
         comment = 'poloidal flux ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad', &
         cache%rad, lbound(cache%rad), ubound(cache%rad), unit = 'cm', &
         comment = 'radial position on OX line ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/F', &
         cache%F, lbound(cache%F), ubound(cache%F), unit = 'G cm', &
         comment = 'covariant toroidal equilibrium field ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/p', &
         cache%p, lbound(cache%p), ubound(cache%p), unit = 'dyn cm^-2', &
         comment = 'equilibrium pressure ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/FdF_dpsi', &
         cache%FdF_dpsi, lbound(cache%FdF_dpsi), ubound(cache%FdF_dpsi), unit = 'G', &
         comment = 'FF''(psi) ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dp_dpsi', &
         cache%dp_dpsi, lbound(cache%dp_dpsi), ubound(cache%dp_dpsi), &
         comment = 'p''(psi) ' // trim(adjustl(comment)), unit = 'dyn cm^-2 Mx^-1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/q', &
         cache%q, lbound(cache%q), ubound(cache%q), unit = '1', &
         comment = 'safety factor ' // trim(adjustl(comment)))
    call h5_close(h5id_root)
  end subroutine flux_func_cache_write

  subroutine flux_func_cache_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(flux_func_cache), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi', cache%psi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rad', cache%rad)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/F', cache%F)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/p', cache%p)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/FdF_dpsi', cache%FdF_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dp_dpsi', cache%dp_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/q', cache%q)
    call h5_close(h5id_root)
  end subroutine flux_func_cache_read

  subroutine coord_cache_init(this, n)
    class(coord_cache), intent(inout) :: this
    integer, intent(in) :: n

    call coord_cache_deinit(this)
    this%n = n
    allocate(this%ktri(this%n))
    allocate(this%R(this%n))
    allocate(this%Z(this%n))
    allocate(this%psi(this%n))
    allocate(this%theta(this%n))
    allocate(this%sqrt_g(this%n))
    allocate(this%B0_R(this%n))
    allocate(this%B0_Z(this%n))
    allocate(this%dR_dtheta(this%n))
    allocate(this%dZ_dtheta(this%n))
  end subroutine coord_cache_init

  subroutine coord_cache_deinit(this)
    class(coord_cache), intent(inout) :: this

    this%n = 0
    if (allocated(this%ktri)) deallocate(this%ktri)
    if (allocated(this%R)) deallocate(this%R)
    if (allocated(this%Z)) deallocate(this%Z)
    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%theta)) deallocate(this%theta)
    if (allocated(this%sqrt_g)) deallocate(this%sqrt_g)
    if (allocated(this%B0_R)) deallocate(this%B0_R)
    if (allocated(this%B0_Z)) deallocate(this%B0_Z)
    if (allocated(this%dR_dtheta)) deallocate(this%dR_dtheta)
    if (allocated(this%dZ_dtheta)) deallocate(this%dZ_dtheta)
  end subroutine coord_cache_deinit

  subroutine coord_cache_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(coord_cache), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ktri', cache%ktri, &
         lbound(cache%ktri), ubound(cache%ktri), &
         comment = 'triangle index of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R', cache%R, &
         lbound(cache%R), ubound(cache%R), unit = 'cm', &
         comment = 'R coordinate of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z', cache%Z, &
         lbound(cache%Z), ubound(cache%Z), unit = 'cm', &
         comment = 'Z coordinate of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', cache%theta, &
         lbound(cache%psi), ubound(cache%psi), unit = 'Mx', &
         comment = 'poloidal flux at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/theta', cache%theta, &
         lbound(cache%theta), ubound(cache%theta), unit = 'rad', &
         comment = 'flux poloidal angle at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/sqrt_g', cache%sqrt_g, &
         lbound(cache%sqrt_g), ubound(cache%sqrt_g), unit = 'cm G^-1', &
         comment = 'Jacobian at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache%B0_R, &
         lbound(cache%B0_R), ubound(cache%B0_R), unit = 'G', &
         comment = 'R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache%B0_Z, &
         lbound(cache%B0_Z), ubound(cache%B0_Z), unit = 'G', &
         comment = 'Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dR_dtheta', cache%dR_dtheta, &
         lbound(cache%dR_dtheta), ubound(cache%dR_dtheta), unit = 'cm rad^-1', &
         comment = 'Jacobian element (R, theta) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dZ_dtheta', cache%dZ_dtheta, &
         lbound(cache%dZ_dtheta), ubound(cache%dZ_dtheta), unit = 'cm rad^-1', &
         comment = 'Jacobian element (Z, theta) at ' // trim(adjustl(comment)))
    call h5_close(h5id_root)
  end subroutine coord_cache_write

  subroutine coord_cache_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(coord_cache), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/ktri', cache%ktri)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R', cache%R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z', cache%Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi', cache%theta)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/theta', cache%theta)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/sqrt_g', cache%sqrt_g)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache%B0_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache%B0_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dR_dtheta', cache%dR_dtheta)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dZ_dtheta', cache%dZ_dtheta)
    call h5_close(h5id_root)
  end subroutine coord_cache_read

  subroutine coord_cache_ext_init(this, nrad, npol)
    type(coord_cache_ext), intent(inout) :: this
    integer, intent(in) :: nrad, npol

    call coord_cache_ext_deinit(this)
    call coord_cache_init(this, nrad * npol)
    this%nrad = nrad
    this%npol = npol
    this%m = 0
    allocate(this%rad(this%n))
    allocate(this%q(this%n))
    allocate(this%dq_dpsi(this%n))
    allocate(this%dR_dpsi(this%n))
    allocate(this%dZ_dpsi(this%n))
    allocate(this%d2R_dpsi_dtheta(this%n))
    allocate(this%d2Z_dpsi_dtheta(this%n))
    allocate(this%B0_phi(this%n))
    allocate(this%dB0R_dR(this%n))
    allocate(this%dB0R_dZ(this%n))
    allocate(this%dB0phi_dR(this%n))
    allocate(this%dB0phi_dZ(this%n))
    allocate(this%dB0Z_dR(this%n))
    allocate(this%dB0Z_dZ(this%n))
    allocate(this%B0_2(this%n))
  end subroutine coord_cache_ext_init

  subroutine coord_cache_ext_deinit(this)
    type(coord_cache_ext), intent(inout) :: this

    call coord_cache_deinit(this)
    this%nrad = 0
    this%npol = 0
    this%m = 0
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%q)) deallocate(this%q)
    if (allocated(this%dq_dpsi)) deallocate(this%dq_dpsi)
    if (allocated(this%dR_dpsi)) deallocate(this%dR_dpsi)
    if (allocated(this%dZ_dpsi)) deallocate(this%dZ_dpsi)
    if (allocated(this%d2R_dpsi_dtheta)) deallocate(this%d2R_dpsi_dtheta)
    if (allocated(this%d2Z_dpsi_dtheta)) deallocate(this%d2Z_dpsi_dtheta)
    if (allocated(this%B0_phi)) deallocate(this%B0_phi)
    if (allocated(this%dB0R_dR)) deallocate(this%dB0R_dR)
    if (allocated(this%dB0R_dZ)) deallocate(this%dB0R_dZ)
    if (allocated(this%dB0phi_dR)) deallocate(this%dB0phi_dR)
    if (allocated(this%dB0phi_dZ)) deallocate(this%dB0phi_dZ)
    if (allocated(this%dB0Z_dR)) deallocate(this%dB0Z_dR)
    if (allocated(this%dB0Z_dZ)) deallocate(this%dB0Z_dZ)
    if (allocated(this%B0_2)) deallocate(this%B0_2)
  end subroutine coord_cache_ext_deinit

  subroutine coord_cache_ext_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    type(coord_cache_ext), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call coord_cache_write(cache, file, dataset, comment)
    call h5_open_rw(file, h5id_root)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nrad', cache%nrad, &
         comment = 'number of radial divisions for ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/npol', cache%npol, &
         comment = 'number of poloidal divisions for ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m', cache%m, &
         comment = 'poloidal mode number of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad', cache%rad, &
         lbound(cache%rad), ubound(cache%rad), unit = 'cm', &
         comment = 'minor radius (on OX line) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/q', cache%q, &
         lbound(cache%q), ubound(cache%q), unit = '1', &
         comment = 'safety factor q at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dq_dpsi', cache%dq_dpsi, &
         lbound(cache%dq_dpsi), ubound(cache%dq_dpsi), unit = 'Mx^-1', &
         comment = 'q''(psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dR_dpsi', cache%dR_dpsi, &
         lbound(cache%dR_dpsi), ubound(cache%dR_dpsi), unit = 'cm Mx^-1', &
         comment = 'Jacobian element (R, psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dZ_dpsi', cache%dZ_dpsi, &
         lbound(cache%dZ_dpsi), ubound(cache%dZ_dpsi), unit = 'cm Mx^-1', &
         comment = 'Jacobian element (Z, psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2R_dpsi_dtheta', cache%d2R_dpsi_dtheta, &
         lbound(cache%d2R_dpsi_dtheta), ubound(cache%d2R_dpsi_dtheta), unit = 'cm Mx^-1 rad^-1', &
         comment = 'Hessian element (psi, theta) of R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2Z_dpsi_dtheta', cache%d2Z_dpsi_dtheta, &
         lbound(cache%d2Z_dpsi_dtheta), ubound(cache%d2Z_dpsi_dtheta), unit = 'cm Mx^-1 rad^-1', &
         comment = 'Hessian element (psi, theta) of Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache%B0_phi, &
         lbound(cache%B0_phi), ubound(cache%B0_phi), unit = 'G', &
         comment = 'physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache%dB0R_dR, &
         lbound(cache%dB0R_dR), ubound(cache%dB0R_dR), unit = 'G cm^-1', &
         comment = 'R derivative of B0_R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache%dB0R_dZ, &
         lbound(cache%dB0R_dZ), ubound(cache%dB0R_dZ), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache%dB0phi_dR, &
         lbound(cache%dB0phi_dR), ubound(cache%dB0phi_dR), unit = 'G cm^-1', &
         comment = 'R derivative of B0_phi at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache%dB0phi_dZ, &
         lbound(cache%dB0phi_dZ), ubound(cache%dB0phi_dZ), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_phi at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache%dB0Z_dR, &
         lbound(cache%dB0Z_dR), ubound(cache%dB0Z_dR), unit = 'G cm^-1', &
         comment = 'R derivative of B0_Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache%dB0Z_dZ, &
         lbound(cache%dB0Z_dZ), ubound(cache%dB0Z_dZ), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_2', cache%B0_2, &
         lbound(cache%B0_2), ubound(cache%B0_2), unit = 'G^2', &
         comment = 'square of equilibrium magnetic field  at ' // trim(adjustl(comment)))
    call h5_close(h5id_root)
  end subroutine coord_cache_ext_write

  subroutine coord_cache_ext_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(coord_cache_ext), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call coord_cache_read(cache, file, dataset)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nrad', cache%nrad)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/npol', cache%npol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m', cache%m)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/q', cache%q)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dq_dpsi', cache%dq_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dR_dpsi', cache%dR_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dZ_dpsi', cache%dZ_dpsi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/d2R_dpsi_dtheta', cache%d2R_dpsi_dtheta)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/d2Z_dpsi_dtheta', cache%d2Z_dpsi_dtheta)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache%B0_phi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache%dB0R_dR)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache%dB0R_dZ)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache%dB0phi_dR)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache%dB0phi_dZ)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache%dB0Z_dR)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache%dB0Z_dZ)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_2', cache%B0_2)
    call h5_close(h5id_root)
  end subroutine coord_cache_ext_read

  subroutine generate_mesh(unprocessed_geqdsk)
    use magdif_conf, only: conf, log
    use magdif_util, only: get_field_filenames, init_field

    character(len = *), intent(in) :: unprocessed_geqdsk
    character(len = 1024) :: gfile, pfile, convexfile

    call get_field_filenames(gfile, pfile, convexfile)
    log%msg = 'attempting to read unprocessed G EQDSK file ' // trim(unprocessed_geqdsk)
    if (log%info) call log%write
    call equil%read(trim(unprocessed_geqdsk))
    call equil%classify
    call equil%standardise
    if (conf%kilca_scale_factor /= 0) then
       call equil%scale(conf%kilca_scale_factor)
    end if
    log%msg = 'attempting to write processed G EQDSK file ' // trim(gfile)
    if (log%info) call log%write
    call equil%write(trim(gfile))
    call init_field(gfile, pfile, convexfile)

    if (conf%kilca_scale_factor /= 0) then
       mesh%n = conf%n * conf%kilca_scale_factor
    else
       mesh%n = conf%n
    end if
    call create_mesh_points(convexfile)
    call compare_gpec_coordinates
    call write_illustration_data(5, 8, 256, 256)
    call connect_mesh_points
    call write_meshfile_for_freefem
    call compute_sample_polmodes
    call compute_gpec_jacfac
    call cache_equilibrium_field
    call init_flux_variables
    call flux_func_cache_check
    call check_safety_factor
    call check_resonance_positions
    call compute_j0phi
    call check_curr0
  end subroutine generate_mesh

  subroutine write_mesh_cache
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    use magdif_conf, only: datafile
    integer(HID_T) :: h5id_root

    call mesh_write(mesh, datafile, 'mesh')
    call flux_func_cache_write(fs, datafile, 'cache/fs', 'on flux surfaces')
    call flux_func_cache_write(fs_half, datafile, 'cache/fs_half', 'between flux surfaces')
    call coord_cache_write(sample_polmodes, datafile, 'cache/sample_polmodes', &
         'poloidal mode sampling points')
    ! TODO: put in separate subroutine for edge_cache_type
    call h5_open_rw(datafile, h5id_root)
    ! TODO: revise naming and indexing when edge_cache type is working for GL quadrature in compute_currn
    call h5_add(h5id_root, 'cache/B0R_edge', B0R, lbound(B0R), ubound(B0R), &
         comment = 'R component of equilibrium magnetic field on triangle edge', unit = 'G')
    call h5_add(h5id_root, 'cache/B0phi_edge', B0phi, lbound(B0phi), ubound(B0phi), &
         comment = 'phi component of equilibrium magnetic field on triangle edge', unit = 'G')
    call h5_add(h5id_root, 'cache/B0Z_edge', B0Z, lbound(B0Z), ubound(B0Z), &
         comment = 'Z component of equilibrium magnetic field on triangle edge', unit = 'G')
    call h5_add(h5id_root, 'cache/B0_flux', B0flux, lbound(B0flux), ubound(B0flux), &
         comment = 'Equilibrium magnetic flux through triangle edge', unit = 'G cm^2')
    call h5_add(h5id_root, 'cache/B0R_centr', B0R_Omega, lbound(B0R_Omega), ubound(B0R_Omega), &
         comment = 'R component of equilibrium magnetic field on triangle ''centroid''', unit = 'G')
    call h5_add(h5id_root, 'cache/B0phi_centr', B0phi_Omega, lbound(B0phi_Omega), ubound(B0phi_Omega), &
         comment = 'phi component of equilibrium magnetic field on triangle ''centroid''', unit = 'G')
    call h5_add(h5id_root, 'cache/B0Z_centr', B0Z_Omega, lbound(B0Z_Omega), ubound(B0Z_Omega), &
         comment = 'Z component of equilibrium magnetic field on triangle ''centroid''', unit = 'G')
    call h5_add(h5id_root, 'cache/j0phi_edge', j0phi, lbound(j0phi), ubound(j0phi), &
         comment = 'phi component of equilibrium current density on triangle edge', unit = 'statA cm^-2')
    call h5_close(h5id_root)
  end subroutine write_mesh_cache

  subroutine read_mesh_cache
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use magdif_conf, only: datafile
    integer(HID_T) :: h5id_root

    call mesh_read(mesh, datafile, 'mesh')
    call flux_func_cache_init(fs, mesh%nflux, .false.)
    call flux_func_cache_init(fs_half, mesh%nflux, .true.)
    call coord_cache_init(sample_polmodes, mesh%ntri)
    call flux_func_cache_read(fs, datafile, 'cache/fs')
    call flux_func_cache_read(fs_half, datafile, 'cache/fs_half')
    call coord_cache_read(sample_polmodes, datafile, 'cache/sample_polmodes')
    ! TODO: revise naming and indexing when edge_cache type is working for GL quadrature in compute_currn
    allocate(B0R(3, mesh%ntri))
    allocate(B0phi(3, mesh%ntri))
    allocate(B0Z(3, mesh%ntri))
    allocate(B0R_Omega(mesh%ntri))
    allocate(B0phi_Omega(mesh%ntri))
    allocate(B0Z_Omega(mesh%ntri))
    allocate(B0flux(3, mesh%ntri))
    allocate(j0phi(3, mesh%ntri))
    call h5_open(datafile, h5id_root)
    call h5_get(h5id_root, 'cache/B0R_edge', B0R)
    call h5_get(h5id_root, 'cache/B0phi_edge', B0phi)
    call h5_get(h5id_root, 'cache/B0Z_edge', B0Z)
    call h5_get(h5id_root, 'cache/B0_flux', B0flux)
    call h5_get(h5id_root, 'cache/B0R_centr', B0R_Omega)
    call h5_get(h5id_root, 'cache/B0phi_centr', B0phi_Omega)
    call h5_get(h5id_root, 'cache/B0Z_centr', B0Z_Omega)
    call h5_get(h5id_root, 'cache/j0phi_edge', j0phi)
    call h5_close(h5id_root)
  end subroutine read_mesh_cache

  subroutine compute_resonance_positions(psi_sample, q_sample, psi2rho_norm)
    use magdif_conf, only: conf, log
    use magdif_util, only: flux_func
    use netlib_mod, only: zeroin
    real(dp), dimension(:), intent(in) :: psi_sample
    real(dp), dimension(:), intent(in) :: q_sample
    interface
       function psi2rho_norm(psi)
         import :: dp
         real(dp), intent(in) :: psi
         real(dp) :: psi2rho_norm
       end function psi2rho_norm
    end interface
    real(dp) :: psi_min, psi_max
    integer :: m
    type(flux_func) :: psi_eqd

    if (size(psi_sample) /= size(q_sample)) then
       call log%msg_arg_size('refine_resonant_surfaces', 'size(psi_sample)', &
            'size(q_sample)', size(psi_sample), size(q_sample))
       if (log%err) call log%write
       error stop
    end if
    psi_min = minval(psi_sample)
    psi_max = maxval(psi_sample)
    call psi_eqd%init(4, psi_sample)
    mesh%m_res_min = max(ceiling(minval(abs(q_sample)) * dble(mesh%n)), conf%n + 1)
    mesh%m_res_max = floor(maxval(abs(q_sample)) * dble(mesh%n))
    if (allocated(mesh%res_modes)) deallocate(mesh%res_modes)
    if (conf%kilca_scale_factor /= 0) then
       allocate(mesh%res_modes(1))
       mesh%res_modes(:) = [conf%kilca_pol_mode]
    else
       allocate(mesh%res_modes(mesh%m_res_max - mesh%m_res_min + 1))
       mesh%res_modes(:) = [(m, m = mesh%m_res_min, mesh%m_res_max)]
    end if
    if (allocated(mesh%psi_res)) deallocate(mesh%psi_res)
    if (allocated(mesh%rad_norm_res)) deallocate(mesh%rad_norm_res)
    allocate(mesh%psi_res(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%rad_norm_res(mesh%m_res_min:mesh%m_res_max))
    log%msg = 'resonance positions:'
    if (log%debug) call log%write
    do m = mesh%m_res_min, mesh%m_res_max
       mesh%psi_res(m) = zeroin(psi_min, psi_max, q_interp_resonant, 1d-9)
       mesh%rad_norm_res(m) = psi2rho_norm(mesh%psi_res(m))
       write (log%msg, '("m = ", i2, ", psi_m = ", es24.16e3, ", rho_m = ", f19.16)') &
            m, mesh%psi_res(m), mesh%rad_norm_res(m)
       if (log%debug) call log%write
    end do

  contains
    function q_interp_resonant(psi)
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant
      q_interp_resonant = psi_eqd%interp(abs(q_sample), psi) - dble(m) / dble(mesh%n)
    end function q_interp_resonant
  end subroutine compute_resonance_positions

  subroutine refine_eqd_partition(nref, deletions, additions, refinement, res, refined, &
       ref_ind)
    use magdif_conf, only: conf, log
    use magdif_util, only: linspace
    integer, intent(in) :: nref
    integer, dimension(:), intent(in) :: deletions, additions
    real(dp), dimension(:), intent(in) :: refinement, res
    real(dp), dimension(:), allocatable, intent(out) :: refined
    integer, dimension(:), intent(out) :: ref_ind
    integer :: kref, k
    integer, dimension(nref) :: coarse_lo, coarse_hi, fine_lo, fine_hi
    real(dp) :: coarse_sep
    real(dp), dimension(nref) :: fine_sep, factor
    real(dp), dimension(maxval(additions) + 1, nref) :: geom_ser

    if (allocated(refined)) deallocate(refined)
    if (nref < 1) then
       mesh%nflux = conf%nflux_unref
       allocate(refined(0:mesh%nflux))
       refined = linspace(0d0, 1d0, mesh%nflux + 1, 0, 0)
       return
    end if
    if (nref /= size(deletions)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(deletions)', nref, &
            size(deletions))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(additions)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(additions)', nref, &
            size(additions))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(refinement)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(refinement)', nref, &
            size(refinement))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(res)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(res)', nref, size(res))
       if (log%err) call log%write
       error stop
    end if
    if (nref /= size(ref_ind)) then
       call log%msg_arg_size('refine_eqd_partition', 'nref', 'size(ref_ind)', nref, &
            size(ref_ind))
       if (log%err) call log%write
       error stop
    end if
    mesh%nflux = conf%nflux_unref + 2 * sum(additions - deletions)
    allocate(refined(0:mesh%nflux))
    refined = 0d0
    coarse_lo = floor(res * (conf%nflux_unref + 1)) - deletions
    coarse_hi = ceiling(res * (conf%nflux_unref + 1)) + deletions
    ! compute upper and lower array indices of refined regions
    fine_lo = coarse_lo + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 0, nref - 1)]
    fine_hi = coarse_hi + [(2 * sum(additions(1:kref) - deletions(1:kref)), &
         kref = 1, nref)]
    ref_ind = fine_hi - additions
    ! compute separations between flux surfaces using geometric series
    coarse_sep = 1d0 / dble(conf%nflux_unref)
    fine_sep = coarse_sep * refinement
    factor = (coarse_sep / fine_sep) ** (1d0 / dble(additions + 1))
    geom_ser = reshape([(((factor(kref) ** k - factor(kref)) / (factor(kref) - 1d0), &
         k = 1, maxval(additions) + 1), kref = 1, nref)], [maxval(additions) + 1, nref])
    ! compute refined regions around resonant flux surfaces
    do kref = 1, nref
       refined(fine_lo(kref):fine_lo(kref)+additions(kref)) = res(kref) &
            - (geom_ser(additions(kref)+1:1:-1, kref) + 0.5d0) * fine_sep(kref)
       refined(fine_hi(kref)-additions(kref):fine_hi(kref)) = res(kref) &
            + (geom_ser(1:additions(kref)+1, kref) + 0.5d0) * fine_sep(kref)
    end do
    ! compute equidistant positions between refined regions
    refined(0:fine_lo(1)-1) = linspace(0d0, refined(fine_lo(1)), fine_lo(1), 0, 1)
    do kref = 2, nref
       refined(fine_hi(kref-1)+1:fine_lo(kref)-1) = linspace(refined(fine_hi(kref-1)), &
            refined(fine_lo(kref)), fine_lo(kref) - fine_hi(kref-1) - 1, 1, 1)
    end do
    refined(fine_hi(nref)+1:mesh%nflux) = linspace(refined(fine_hi(nref)), 1d0, &
         mesh%nflux - fine_hi(nref), 1, 0)
  end subroutine refine_eqd_partition

  subroutine refine_resonant_surfaces(rho_norm_ref)
    use magdif_conf, only: conf_arr, log
    real(dp), dimension(:), allocatable, intent(out) :: rho_norm_ref
    integer :: m, kref
    integer, dimension(:), allocatable :: ref_ind
    logical, dimension(:), allocatable :: mask

    if (allocated(mesh%refinement)) deallocate(mesh%refinement)
    if (allocated(mesh%deletions)) deallocate(mesh%deletions)
    if (allocated(mesh%additions)) deallocate(mesh%additions)
    allocate(mesh%refinement(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%deletions(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%additions(mesh%m_res_min:mesh%m_res_max))
    mesh%refinement(:) = conf_arr%refinement
    mesh%deletions(:) = conf_arr%deletions
    mesh%additions(:) = conf_arr%additions
    allocate(mask(mesh%m_res_min:mesh%m_res_max))
    mask = 0d0 < mesh%refinement .and. mesh%refinement < 1d0
    allocate(ref_ind(count(mask)))
    call refine_eqd_partition(count(mask), pack(mesh%deletions, mask), &
         pack(mesh%additions, mask), pack(mesh%refinement, mask), &
         pack(mesh%rad_norm_res, mask), rho_norm_ref, ref_ind)
    log%msg = 'refinement positions:'
    if (log%debug) call log%write
    kref = 0
    do m = mesh%m_res_min, mesh%m_res_max
       if (.not. mask(m)) cycle
       kref = kref + 1
       write (log%msg, '("m = ", i0, ", kf = ", i0, ' // &
            '", rho: ", f19.16, 2(" < ", f19.16))') m, ref_ind(kref), &
            rho_norm_ref(ref_ind(kref) - 1), mesh%rad_norm_res(m), rho_norm_ref(ref_ind(kref))
       if (log%debug) call log%write
    end do
    if (allocated(ref_ind)) deallocate(ref_ind)
    if (allocated(mask)) deallocate(mask)
  end subroutine refine_resonant_surfaces

  subroutine cache_resonance_positions
    use magdif_conf, only: log
    use magdif_util, only: binsearch
    integer :: m, kf_res

    log%msg = 'resonance positions:'
    if (log%debug) call log%write
    allocate(mesh%m_res(mesh%nflux))
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    mesh%m_res = 0
    mesh%res_ind = 0
    ! if more two or more resonance positions are within the same flux surface,
    ! assign the lowest mode number
    do m = mesh%m_res_max, mesh%m_res_min, -1
       call binsearch(fs%psi, 0, mesh%psi_res(m), kf_res)
       mesh%res_ind(m) = kf_res
       mesh%m_res(kf_res) = m
       write (log%msg, '("m = ", i2, ", kf = ", i3, ", rho: ", f19.16, 2(" < ", f19.16))') &
            m, kf_res, fs%rad(kf_res - 1) / fs%rad(mesh%nflux), mesh%rad_norm_res(m), &
            fs%rad(kf_res) / fs%rad(mesh%nflux)
       if (log%debug) call log%write
    end do
  end subroutine cache_resonance_positions

  subroutine write_kilca_convexfile(rho_max, convexfile)
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90

    integer, parameter :: nrz = 96  ! at most 100 values are read in by field_divB0.f90
    real(dp), intent(in) :: rho_max
    character(len = *), intent(in) :: convexfile
    real(dp) :: theta
    integer :: fid, kp

    open(newunit = fid, file = convexfile, status = 'replace')
    do kp = 1, nrz
       theta = dble(kp) / dble(nrz) * 2d0 * pi
       write (fid, '(2(1x, es24.16e3))') equil%rmaxis + rho_max * cos(theta), &
            equil%zmaxis + rho_max * sin(theta)
    end do
    close(fid)
  end subroutine write_kilca_convexfile

  subroutine create_mesh_points(convexfile)
    use constants, only: pi
    use magdif_conf, only: conf, conf_arr, log
    use magdif_util, only: interp_psi_pol, flux_func
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
         raxis, zaxis, load_magdata_in_symfluxcoord
    use field_line_integration_mod, only: circ_mesh_scale, o_point, x_point, &
         theta0_at_xpoint, theta_axis, theta0
    use points_2d, only: s_min, create_points_2d

    character(len = *), intent(in) :: convexfile
    integer :: kf, fid
    integer, dimension(:), allocatable :: n_theta
    real(dp), dimension(:), allocatable :: rho_norm_eqd, rho_norm_ref
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi
    type(flux_func) :: psi_interpolator
    real(dp) :: psi_axis, rho_max

    theta0_at_xpoint = .true.
    circ_mesh_scale = conf%kilca_scale_factor
    if (conf%kilca_scale_factor /= 0) then
       ! calculate maximal extent from magnetic axis
       rho_max = min(equil%rmaxis - equil%rleft, &
            equil%rleft + equil%rdim - equil%rmaxis, &
            equil%zdim * 0.5d0 + equil%zmid - equil%zmaxis, &
            equil%zdim * 0.5d0 - equil%zmid + equil%zmaxis)
       call write_kilca_convexfile(rho_max, convexfile)
       o_point = [equil%rmaxis, equil%zmaxis]
       x_point = o_point + [rho_max, 0d0]
    end if
    ! calculates points on a fine grid in the core region by integrating along field lines
    call preload_for_SYNCH
    ! loads points that are calculated in preload_for_SYNCH into module variables
    ! and spline interpolates them for use with magdata_in_symfluxcoord_ext
    call load_magdata_in_symfluxcoord
    mesh%R_O = raxis
    mesh%Z_O = zaxis
    mesh%R_X = X_point(1)
    mesh%Z_X = X_point(2)
    mesh%theta0 = theta0
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    call psi_interpolator%init(4, psisurf(1:) * psipol_max + psi_axis)
    ! interpolate between psi and rho
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd = rbeg / hypot(theta_axis(1), theta_axis(2))

    call compute_resonance_positions(psisurf(1:) * psipol_max + psi_axis, qsaf, psi2rho_norm)
    call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
    call refine_resonant_surfaces(rho_norm_ref)
    call fs%init(mesh%nflux, .false.)
    call fs_half%init(mesh%nflux, .true.)
    fs%rad = rho_norm_ref
    fs%psi = [(interp_psi_pol(raxis + fs%rad(kf) * theta_axis(1), &
         zaxis + fs%rad(kf) * theta_axis(2)), kf = 0, mesh%nflux)]
    fs_half%rad = 0.5d0 * (fs%rad(0:mesh%nflux-1) + fs%rad(1:mesh%nflux))
    fs_half%psi = [(interp_psi_pol(raxis + fs_half%rad(kf) * theta_axis(1), &
         zaxis + fs_half%rad(kf) * theta_axis(2)), kf = 1, mesh%nflux)]
    fs%rad = fs%rad * hypot(theta_axis(1), theta_axis(2))
    fs_half%rad = fs_half%rad * hypot(theta_axis(1), theta_axis(2))
    call flux_func_cache_check
    call cache_resonance_positions
    allocate(mesh%theta_offset(mesh%nflux))
    do kf = 1, mesh%nflux
       mesh%theta_offset(kf) = theta_offset(fs_half%psi(kf))
    end do
    ! dump presumedly optimal values for poloidal resolution
    open(newunit = fid, file = 'optpolres.dat', status = 'replace')
    do kf = 1, mesh%nflux - 1
       write (fid, '(2(1x, es24.16e3))') fs%rad(kf), 2d0 * pi * fs%rad(kf) / &
            (fs_half%rad(kf + 1) - fs_half%rad(kf))
    end do
    write (fid, '(2(1x, es24.16e3))') fs%rad(mesh%nflux), 2d0 * pi * fs%rad(mesh%nflux) / &
         (fs%rad(mesh%nflux) - fs%rad(mesh%nflux - 1))
    close(fid)

    call init_indices
    mesh%ntri = mesh%kt_low(mesh%nflux + 1)
    mesh%npoint = mesh%kp_low(mesh%nflux + 1)
    allocate(mesh%node_R(mesh%npoint))
    allocate(mesh%node_Z(mesh%npoint))
    allocate(mesh%node_theta_flux(mesh%npoint))
    allocate(mesh%node_theta_geom(mesh%npoint))

    allocate(n_theta(mesh%nflux))
    allocate(points(3, mesh%npoint))
    allocate(points_s_theta_phi(3, mesh%npoint))
    n_theta = conf%nkpol
    s_min = 1d-16
    ! inp_label 2 to use poloidal psi with magdata_in_symfluxcoord_ext
    ! psi is not normalized by psipol_max, but shifted by -psi_axis
    call create_points_2d(2, n_theta, points, points_s_theta_phi, r_scaling_func = psi_ref)
    mesh%node_R(1) = mesh%R_O
    mesh%node_R(2:) = points(1, 2:)
    mesh%node_Z(1) = mesh%Z_O
    mesh%node_Z(2:) = points(3, 2:)
    mesh%node_theta_flux(1) = 0d0
    mesh%node_theta_flux(2:) = points_s_theta_phi(2, 2:)
    mesh%node_theta_geom(1) = 0d0
    mesh%node_theta_geom(2:) = atan2(mesh%node_Z(2:) - mesh%Z_O, mesh%node_R(2:) - mesh%R_O)
    mesh%R_min = minval(mesh%node_R)
    mesh%R_max = maxval(mesh%node_R)
    mesh%Z_min = minval(mesh%node_Z)
    mesh%Z_max = maxval(mesh%node_Z)

    if (allocated(n_theta)) deallocate(n_theta)
    if (allocated(points)) deallocate(points)
    if (allocated(points_s_theta_phi)) deallocate(points_s_theta_phi)
    if (allocated(rho_norm_ref)) deallocate(rho_norm_ref)
    if (allocated(rho_norm_eqd)) deallocate(rho_norm_eqd)

  contains
    function psi2rho_norm(psi) result(rho_norm)
      real(dp), intent(in) :: psi
      real(dp) :: rho_norm
      rho_norm = psi_interpolator%interp(rho_norm_eqd, psi)
    end function psi2rho_norm
    function psi_ref(psi_eqd)
      real(dp), dimension(:), intent(in) :: psi_eqd
      real(dp), dimension(size(psi_eqd)) :: psi_ref
      integer :: kf
      if (mesh%nflux /= size(psi_eqd)) then
         call log%msg_arg_size('psi_ref', 'nflux', 'size(psi_eqd)', &
              mesh%nflux, size(psi_eqd))
         if (log%err) call log%write
         error stop
      end if
      psi_ref = [(interp_psi_pol(raxis + rho_norm_ref(kf) * theta_axis(1), &
           zaxis + rho_norm_ref(kf) * theta_axis(2)) - psi_axis, kf = 1, mesh%nflux)]
    end function psi_ref
  end subroutine create_mesh_points

  function theta_offset(psi)
    use netlib_mod, only: zeroin
    use constants, only: pi  ! orbit_mod.f90
    real(dp), intent(in) :: psi
    real(dp) :: theta_offset

    theta_offset = zeroin(0d0, pi, Z_offset, 1d-9)
  contains
    function Z_offset(theta)
      use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
      real(dp), intent(in) :: theta
      real(dp) :: Z_offset
      real(dp) :: dum
      call magdata_in_symfluxcoord_ext(2, dum, psi - fs%psi(0), theta, dum, dum, &
           dum, dum, dum, dum, dum, dum, Z_offset, dum, dum)
      Z_offset = Z_offset - mesh%Z_O
    end function Z_offset
  end function theta_offset

  !> Allocates and initializes #kp_low, #kp_max, #kt_low and #kt_max based on the values
  !> of #magdif_config::nflux and #magdif_config::nkpol. Deallocation is done in
  !> magdif_cleanup().
  subroutine init_indices
    use magdif_conf, only: conf, log
    integer :: kf

    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux + 1))
    allocate(mesh%kt_low(mesh%nflux + 1))

    mesh%kp_max = conf%nkpol
    mesh%kt_max = 2 * conf%nkpol
    mesh%kt_max(1) = conf%nkpol

    mesh%kp_low(1) = 1
    do kf = 2, mesh%nflux + 1
       mesh%kp_low(kf) = mesh%kp_low(kf-1) + mesh%kp_max(kf-1)
    end do
    mesh%kt_low(1) = 0
    do kf = 2, mesh%nflux + 1
       mesh%kt_low(kf) = mesh%kt_low(kf-1) + mesh%kt_max(kf-1)
    end do

    write (log%msg, '("Number of points up to LCFS: ", i0)') mesh%kp_low(mesh%nflux + 1)
    if (log%info) call log%write
    write (log%msg, '("Number of triangles up to LCFS: ", i0)') mesh%kt_low(mesh%nflux + 1)
    if (log%info) call log%write
  end subroutine init_indices

  function tri_area(ktri)
    integer, intent(in) :: ktri
    real(dp) :: tri_area, e1_R, e1_Z, e2_R, e2_Z
    e1_R = mesh%node_R(mesh%tri_node(1, ktri)) - mesh%node_R(mesh%tri_node(3, ktri))
    e1_Z = mesh%node_Z(mesh%tri_node(1, ktri)) - mesh%node_Z(mesh%tri_node(3, ktri))
    e2_R = mesh%node_R(mesh%tri_node(2, ktri)) - mesh%node_R(mesh%tri_node(3, ktri))
    e2_Z = mesh%node_Z(mesh%tri_node(2, ktri)) - mesh%node_Z(mesh%tri_node(3, ktri))
    tri_area = 0.5d0 * abs(e1_R * e2_Z - e1_Z * e2_R)
  end function tri_area

  !> Returns the indices of the two triangles sharing an edge.
  !>
  !> @param knot1 index of first knot of the edge
  !> @param knot2 index of second knot of the edge
  !> @param common_tri indices of the triangles sharing the given edge
  !>
  !> The program is halted if the input data is invalid, i.e. if more than two triangles
  !> appear to share the edge.
  subroutine common_triangles(knot1, knot2, common_tri)
    use magdif_conf, only: log
    integer, intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)
    logical :: tri_mask(mesh%ntri)

    common_tri = [0, 0]
    tri_mask = any(mesh%tri_node == knot1, 1) .and. any(mesh%tri_node == knot2, 1)
    select case (count(tri_mask))
    case (0)
       common_tri = [0, 0]
    case (1)
       ! TODO: use findloc intrinsic when implementation is readily available
       common_tri = [maxloc(merge(1, 0, tri_mask), 1), 0]
    case (2)
       ! TODO: use findloc intrinsic when implementation is readily available
       common_tri = [maxloc(merge(1, 0, tri_mask), 1), &
            mesh%ntri + 1 - maxloc(merge(1, 0, tri_mask(mesh%ntri:1:-1)), 1)]
    case default
       write (log%msg, '("More than two common triangles for knots ", ' // &
            'i0, " and ", i0)') knot1, knot2
       if (log%err) call log%write
       error stop
    end select
  end subroutine common_triangles

  subroutine connect_mesh_points
    integer :: kf, kp, kt, ktri, ktri_adj, common_tri(2), kedge, ke, ke_adj

    allocate(mesh%tri_node(3, mesh%ntri))
    allocate(mesh%tri_node_F(mesh%ntri))
    allocate(mesh%adj_tri(3, mesh%ntri))
    allocate(mesh%adj_edge(3, mesh%ntri))
    mesh%tri_node = 0
    mesh%tri_node_F = 0
    allocate(mesh%area(mesh%ntri))
    ktri = mesh%kt_low(1)
    ! define trianlges on innermost flux surface
    kf = 1
    do kp = 1, mesh%kp_max(kf)
       ktri = ktri + 1
       mesh%tri_node(1, ktri) = mesh%kp_low(kf) + kp
       mesh%tri_node(2, ktri) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
       mesh%tri_node(3, ktri) = mesh%kp_low(kf)
       mesh%tri_node_F(ktri) = 3
       mesh%area(ktri) = tri_area(ktri)
    end do
    ! define triangles on outer flux surfaces
    do kf = 2, mesh%nflux
       do kp = 1, mesh%kp_max(kf)
          ktri = ktri + 1
          mesh%tri_node(1, ktri) = mesh%kp_low(kf) + kp
          mesh%tri_node(2, ktri) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
          mesh%tri_node(3, ktri) = mesh%kp_low(kf - 1) + kp
          mesh%tri_node_F(ktri) = 3
          mesh%area(ktri) = tri_area(ktri)
          ktri = ktri + 1
          mesh%tri_node(1, ktri) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
          mesh%tri_node(2, ktri) = mesh%kp_low(kf - 1) + mod(kp, mesh%kp_max(kf - 1)) + 1
          mesh%tri_node(3, ktri) = mesh%kp_low(kf - 1) + kp
          mesh%tri_node_F(ktri) = 1
          mesh%area(ktri) = tri_area(ktri)
       end do
    end do
    ! set neighbours for edges i and o on innermost flux surface
    kf = 1
    do kt = 1, mesh%kt_max(kf)
       ktri = mesh%kt_low(kf) + kt
       ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
       mesh%adj_tri(2, ktri) = ktri_adj
       mesh%adj_edge(2, ktri) = 3
       mesh%adj_tri(3, ktri_adj) = ktri
       mesh%adj_edge(3, ktri_adj) = 2
    end do
    ! set neighbours for edges i and o on outer flux surfaces
    do kf = 2, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
          if (mod(kt, 2) == 1) then
             mesh%adj_tri(2, ktri) = ktri_adj
             mesh%adj_edge(2, ktri) = 3
             mesh%adj_tri(3, ktri_adj) = ktri
             mesh%adj_edge(3, ktri_adj) = 2
          else
             mesh%adj_tri(1, ktri) = ktri_adj
             mesh%adj_edge(1, ktri) = 3
             mesh%adj_tri(3, ktri_adj) = ktri
             mesh%adj_edge(3, ktri_adj) = 1
          end if
       end do
    end do
    ! set neighbours for edges f
    do kf = 1, mesh%nflux - 1
       do kp = 1, mesh%kp_max(kf)
          call common_triangles(mesh%kp_low(kf) + kp, &
               mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1, common_tri)
          ktri = minval(common_tri)
          ktri_adj = maxval(common_tri)
          mesh%adj_tri(1, ktri) = ktri_adj
          mesh%adj_edge(1, ktri) = 2
          mesh%adj_tri(2, ktri_adj) = ktri
          mesh%adj_edge(2, ktri_adj) = 1
       end do
    end do
    ! set dummy values for edges on boundary
    kf = mesh%nflux
    do kt = 1, mesh%kt_max(kf)
       if (mod(kt, 2) == 1) then
          ktri = mesh%kt_low(kf) + kt
          ktri_adj = ktri + mesh%kt_max(kf) + 1
          mesh%adj_tri(1, ktri) = ktri_adj
          mesh%adj_edge(1, ktri) = 2
       end if
    end do
    ! define edges and 'centroids'
    allocate(mesh%li(2, mesh%ntri))
    allocate(mesh%lo(2, mesh%ntri))
    allocate(mesh%lf(2, mesh%ntri))
    allocate(mesh%ei(mesh%ntri))
    allocate(mesh%eo(mesh%ntri))
    allocate(mesh%ef(mesh%ntri))
    allocate(mesh%orient(mesh%ntri))
    allocate(mesh%R_Omega(mesh%ntri))
    allocate(mesh%Z_Omega(mesh%ntri))
    mesh%nedge = (3 * mesh%kt_low(mesh%nflux + 1) + mesh%kp_max(mesh%nflux)) / 2
    allocate(mesh%edge_map2global(3, mesh%ntri))
    allocate(mesh%edge_map2ktri(2, mesh%nedge))
    allocate(mesh%edge_map2ke(2, mesh%nedge))
    mesh%edge_map2global = 0
    mesh%edge_map2ktri = 0
    mesh%edge_map2ke = 0
    kedge = 1
    do ktri = 1, mesh%ntri
       call get_labeled_edges(ktri, mesh%li(:, ktri), mesh%lo(:, ktri), mesh%lf(:, ktri), &
            mesh%ei(ktri), mesh%eo(ktri), mesh%ef(ktri), mesh%orient(ktri))
       call ring_centered_avg_coord(ktri, mesh%R_Omega(ktri), mesh%Z_Omega(ktri))
       do ke = 1, 3
          if (mesh%edge_map2global(ke, ktri) == 0) then
             ktri_adj = mesh%adj_tri(ke, ktri)
             ke_adj = mesh%adj_edge(ke, ktri)
             if (ktri_adj > mesh%ntri) then
                mesh%edge_map2global(ke, ktri) = kedge
                mesh%edge_map2ktri(:, kedge) = [ktri, -1]
                mesh%edge_map2ke(:, kedge) = [ke, -1]
                kedge = kedge + 1
             else
                mesh%edge_map2global(ke, ktri) = kedge
                mesh%edge_map2global(ke_adj, ktri_adj) = kedge
                mesh%edge_map2ktri(:, kedge) = [ktri, ktri_adj]
                mesh%edge_map2ke(:, kedge) = [ke, ke_adj]
                kedge = kedge + 1
             end if
          end if
       end do
    end do
  end subroutine connect_mesh_points

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

  subroutine get_labeled_edges(ktri, li, lo, lf, ei, eo, ef, orient)
    use magdif_conf, only: log
    integer, intent(in) :: ktri
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
    knot_f = mesh%tri_node_F(ktri)
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
    if (mesh%tri_node(i1, ktri) == mesh%tri_node(i2, ktri)) then
       if (log%err) call log%write
       error stop
    end if
    ! last triangle in strip if indices not next to each other
    closing_loop = abs(mesh%tri_node(i1, ktri) - mesh%tri_node(i2, ktri)) /= 1
    i_knot_diff = mesh%tri_node(:, ktri) - mesh%tri_node(knot_f, ktri)
    if (all(i_knot_diff >= 0)) then
       ! knot_f lies on inner surface
       orient = .true.
       if ((mesh%tri_node(i1, ktri) < mesh%tri_node(i2, ktri)) .neqv. closing_loop) then
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
       li = [mesh%tri_node(knot_f, ktri), mesh%tri_node(knot_o, ktri)]
       lo = [mesh%tri_node(knot_i, ktri), mesh%tri_node(knot_f, ktri)]
       lf = [mesh%tri_node(knot_o, ktri), mesh%tri_node(knot_i, ktri)]
    else if (all(i_knot_diff <= 0)) then
       ! knot_f lies on outer surface
       orient = .false.
       if ((mesh%tri_node(i1, ktri) > mesh%tri_node(i2, ktri)) .neqv. closing_loop) then
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
       li = [mesh%tri_node(knot_o, ktri), mesh%tri_node(knot_f, ktri)]
       lo = [mesh%tri_node(knot_f, ktri), mesh%tri_node(knot_i, ktri)]
       lf = [mesh%tri_node(knot_i, ktri), mesh%tri_node(knot_o, ktri)]
    else
       if (log%err) call log%write
       error stop
    end if
  end subroutine get_labeled_edges

  !> Computes the "weighted" centroid for a triangle so that it is approximately
  !> equidistant between the enclosing flux surfaces, independent of triangle orientation.
  !>
  !> @param elem the triangle for which the centroid is to be computed
  !> @param r radial cylindrical coordinate of the centroid
  !> @param z axial cylindrical coordinate of the centroid
  !>
  !> Depending on the orientation of the triangle (see also \p orient of
  !> get_labeled_edges()), two knots lie on the inner flux surface and one on the outer
  !> one, or vice versa. A simple arithmetic mean of the three knots' coordinates would
  !> place the centroid closer to the inner flux surface for one orientation and closer
  !> to the outer one for the other. To counteract this, the "lonely" knot is counted
  !> twice in the averaging procedure, i.e. with double weighting.
  pure subroutine ring_centered_avg_coord(ktri, R, Z)
    integer, intent(in) :: ktri
    real(dp), intent(out) :: R, Z
    integer :: nodes(4)

    nodes = [mesh%tri_node(:, ktri), mesh%tri_node(mesh%tri_node_F(ktri), ktri)]
    R = sum(mesh%node_R(nodes)) * 0.25d0
    Z = sum(mesh%node_Z(nodes)) * 0.25d0
  end subroutine ring_centered_avg_coord

  subroutine compute_gpec_jacfac
    use magdif_util, only: imun
    integer, parameter :: m_max = 16
    integer :: kf, kt, ktri, m
    complex(dp) :: fourier_basis(-m_max:m_max)

    allocate(mesh%gpec_jacfac(-m_max:m_max, mesh%nflux))
    mesh%gpec_jacfac(:, :) = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          associate (s => sample_polmodes)
            fourier_basis = [(exp(-imun * m * s%theta(ktri)), m = -m_max, m_max)]
            mesh%gpec_jacfac(:, kf) = mesh%gpec_jacfac(:, kf) + s%sqrt_g(ktri) * &
                 s%R(ktri) * hypot(s%B0_Z(ktri), -s%B0_Z(ktri)) * fourier_basis
          end associate
       end do
       mesh%gpec_jacfac(:, kf) = mesh%gpec_jacfac(:, kf) / mesh%kt_max(kf)
    end do
  end subroutine compute_gpec_jacfac

  !> Compute coarse grid for poloidal mode sampling points - one point per triangle.
  subroutine compute_sample_polmodes
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use magdif_conf, only: conf
    integer :: ktri, kf, kt
    real(dp) :: dum, q

    call coord_cache_init(sample_polmodes, mesh%ntri)
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          associate (s => sample_polmodes)
            s%psi(ktri) = fs_half%psi(kf)
            s%theta(ktri) = 2d0 * pi * dble(kt - 1) / dble(mesh%kt_max(kf))
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%R(ktri) = mesh%R_O + fs_half%rad(kf) * cos(s%theta(ktri))
               s%Z(ktri) = mesh%Z_O + fs_half%rad(kf) * sin(s%theta(ktri))
               s%dR_dtheta(ktri) = -fs_half%rad(kf) * sin(s%theta(ktri))
               s%dZ_dtheta(ktri) =  fs_half%rad(kf) * cos(s%theta(ktri))
            else
               ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
               call magdata_in_symfluxcoord_ext(2, dum, s%psi(ktri) - fs%psi(0), &
                    s%theta(ktri) + mesh%theta_offset(kf), q, dum, s%sqrt_g(ktri), dum, dum, &
                    s%R(ktri), dum, s%dR_dtheta(ktri), s%Z(ktri), dum, s%dZ_dtheta(ktri))
            end if
            s%ktri(ktri) = point_location(s%R(ktri), s%Z(ktri))
            call field(s%R(ktri), 0d0, s%Z(ktri), s%B0_R(ktri), dum, s%B0_Z(ktri), &
                 dum, dum, dum, dum, dum, dum, dum, dum, dum)
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%sqrt_g(ktri) = equil%cocos%sgn_dpsi * fs_half%rad(kf) / &
                    (-s%B0_R(ktri) * sin(s%theta(ktri)) + s%B0_Z(ktri) * cos(s%theta(ktri)))
            else
               ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
               ! taken together, these three always yield a positive sign in COCOS 3
               s%sqrt_g(ktri) = s%sqrt_g(ktri) * abs(q)
            end if
          end associate
       end do
    end do
  end subroutine compute_sample_polmodes

  !> Compute fine grid for parallel current sampling points.
  subroutine compute_sample_Ipar(sample_Ipar, m)
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: psipol_max, magdata_in_symfluxcoord_ext
    use magdif_conf, only: conf
    use magdif_util, only: linspace, interp_psi_pol
    type(coord_cache_ext), intent(inout) :: sample_Ipar
    integer, intent(in) :: m
    integer :: kf_min, kf_max, krad, kpol, k
    real(dp) :: theta(sample_Ipar%npol), rad_eqd(equil%nw), drad_dpsi, dum
    real(dp), dimension(sample_Ipar%nrad) :: rad, psi

    sample_Ipar%m = m
    rad_eqd(:) = linspace(fs%rad(0), fs%rad(mesh%nflux), equil%nw, 0, 0)
    ! res_ind is half-grid, but evaluation limits refer to full-grid indices
    ! deletions is used to estimate the doubled width of the refined interval
    kf_min = mesh%res_ind(m) - mesh%additions(m) - mesh%deletions(m) - 1
    if (kf_min < 0) kf_min = 0
    kf_max = mesh%res_ind(m) + mesh%additions(m) + mesh%deletions(m)
    if (kf_max > mesh%nflux) kf_max = mesh%nflux
    rad(:) = linspace(fs%rad(kf_min), fs%rad(kf_max), sample_Ipar%nrad, 0, 0)
    ! TODO: replace by dedicated interpolation function
    if (conf%kilca_scale_factor /= 0) then
       psi(:) = [(interp_psi_pol(mesh%R_O, mesh%Z_O + rad(krad)), krad = 1, sample_Ipar%nrad)]
    else
       psi(:) = [(interp_psi_pol(mesh%R_O + rad(krad) * cos(mesh%theta0), &
            mesh%Z_O + rad(krad) * sin(mesh%theta0)), krad = 1, sample_Ipar%nrad)]
    end if
    theta(:) = 2d0 * pi * [(dble(kpol) - 0.5d0, kpol = 1, sample_Ipar%npol)] / dble(sample_Ipar%npol)
    do krad = 1, sample_Ipar%nrad
       do kpol = 1, sample_Ipar%npol
          associate (s => sample_Ipar)
            k = (krad - 1) * s%npol + kpol
            s%rad(k) = rad(krad)
            s%psi(k) = psi(krad)
            s%theta(k) = theta(kpol)
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%R(k) = mesh%R_O + rad(krad) * cos(theta(kpol))
               s%Z(k) = mesh%Z_O + rad(krad) * sin(theta(kpol))
               s%dR_dtheta(k) = -rad(krad) * sin(theta(kpol))
               s%dZ_dtheta(k) =  rad(krad) * cos(theta(kpol))
               s%dq_dpsi(k) = fluxvar%interp(equil%qpsi, psi(krad))
               drad_dpsi = 1d0 / fluxvar%interp(rad_eqd, psi(krad))
               s%dR_dpsi(k) = drad_dpsi * cos(theta(kpol))
               s%dZ_dpsi(k) = drad_dpsi * sin(theta(kpol))
               s%d2R_dpsi_dtheta(k) = -drad_dpsi * sin(theta(kpol))
               s%d2Z_dpsi_dtheta(k) =  drad_dpsi * cos(theta(kpol))
            else
               ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
               call magdata_in_symfluxcoord_ext(2, dum, psi(krad) - fs%psi(0), theta(kpol), &
                    s%q(k), s%dq_dpsi(k), s%sqrt_g(k), dum, dum, &
                    s%R(k), s%dR_dpsi(k), s%dR_dtheta(k), &
                    s%Z(k), s%dZ_dpsi(k), s%dZ_dtheta(k), &
                    s%d2R_dpsi_dtheta(k), s%d2Z_dpsi_dtheta(k))
               ! psi is normalized in derivatives - rescale
               s%dq_dpsi(k) = s%dq_dpsi(k) / psipol_max
               s%dR_dpsi(k) = s%dR_dpsi(k) / psipol_max
               s%dZ_dpsi(k) = s%dZ_dpsi(k) / psipol_max
               s%d2R_dpsi_dtheta(k) = s%d2R_dpsi_dtheta(k) / psipol_max
               s%d2Z_dpsi_dtheta(k) = s%d2Z_dpsi_dtheta(k) / psipol_max
            end if
            s%ktri(k) = point_location(s%R(k), s%Z(k))
            call field(s%R(k), 0d0, s%Z(k), s%B0_R(k), s%B0_phi(k), s%B0_Z(k), &
                 s%dB0R_dR(k), dum, s%dB0R_dZ(k), s%dB0phi_dR(k), dum, s%dB0phi_dZ(k), &
                 s%dB0Z_dR(k), dum, s%dB0Z_dZ(k))
            s%B0_2(k) = s%B0_R(k) ** 2 + s%B0_Z(k) ** 2 + s%B0_phi(k) ** 2
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%sqrt_g(k) = equil%cocos%sgn_dpsi * s%rad(k) / &
                    (-s%B0_R(k) * sin(s%theta(k)) + s%B0_Z(k) * cos(s%theta(k)))
            else
               ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
               ! taken together, these three always yield a positive sign in COCOS 3
               s%sqrt_g(k) = s%sqrt_g(k) * abs(s%q(k))
            end if
          end associate
       end do
    end do
  end subroutine compute_sample_Ipar

  function point_location(r, z) result(location)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf
    use magdif_util, only: interp_psi_pol, binsearch
    real(dp), intent(in) :: r, z
    integer :: location

    integer :: kf, kq, k, k_max, ktri, candidates(6)
    real(dp) :: psi, pol_frac, pol_offset, pol_interval(0:conf%nkpol), thickness

    location = -1
    if (R < mesh%R_min .or. R > mesh%R_max .or. Z < mesh%Z_min .or. Z > mesh%Z_max) return
    location = -2
    psi = interp_psi_pol(r, z)
    if (equil%cocos%sgn_dpsi == +1) then
       if (psi > fs%psi(mesh%nflux)) return
    else
       if (psi < fs%psi(mesh%nflux)) return
    end if
    call binsearch(fs%psi, lbound(fs%psi, 1), psi, kf)
    location = -3

    pol_interval = 0d0
    pol_interval(0:mesh%kp_max(kf)-1) = 0.5d0 / pi * atan2( &
         mesh%node_Z((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf))) - mesh%Z_O, &
         mesh%node_R((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf))) - mesh%R_O)
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
    elseif (kf == mesh%nflux) then
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

    thickness = fs%rad(mesh%nflux) * sqrt(epsilon(1d0)) * 8d0
    do k = 1, k_max
       ktri = candidates(k)
       if (point_in_triangle(ktri, R, Z, thickness)) then
          location = ktri
          exit
       end if
    end do
  end function point_location

  ! based on http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
  function point_in_triangle(ktri, R, Z, thickness) result(probably)
    integer, intent(in) :: ktri
    real(dp), intent(in) :: R, Z
    real(dp), intent(in), optional :: thickness
    logical :: probably
    real(dp), dimension(1:4) :: node_R, node_Z, dist_R, dist_Z
    real(dp), dimension(1:3) :: edge_R, edge_Z, edge_2, dist_2, dotprod

    probably = .false.
    if (ktri <= 0) return
    node_R = mesh%node_R([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
    node_Z = mesh%node_Z([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
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

  subroutine mesh_write(mesh, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(mesh_t), intent(in) :: mesh
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer :: orient(mesh%ntri)

    where (mesh%orient)
       orient = 1
    elsewhere
       orient = 0
    end where
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_O', mesh%R_O, &
         comment = 'R coordinate of O point', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_O', mesh%Z_O, &
         comment = 'Z coordinate of O point', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_X', mesh%R_X, &
         comment = 'R coordinate of X point', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_X', mesh%Z_X, &
         comment = 'Z coordinate of X point', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/theta0', mesh%theta0, unit = 'rad', &
         comment = 'geometric poloidal angle of symmetry flux poloidal angle origin')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_min', mesh%R_min, &
         comment = 'minimal R coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_min', mesh%Z_min, &
         comment = 'minimal Z coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_max', mesh%R_max, &
         comment = 'maximal R coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_max', mesh%Z_max, &
         comment = 'maximal Z coordinate of computational grid', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nflux', mesh%nflux, &
         comment = 'number of flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/npoint', mesh%npoint, &
         comment = 'number of points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ntri', mesh%ntri, &
         comment = 'number of triangles')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nedge', mesh%nedge, &
         comment = 'number of edges')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n', mesh%n, &
         comment = 'toroidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_res_min', mesh%m_res_min, &
         comment = 'minimal absolute poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_res_max', mesh%m_res_max, &
         comment = 'maximal absolute poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/deletions', mesh%deletions, &
         lbound(mesh%deletions), ubound(mesh%deletions), &
         comment = 'number of unrefined flux surfaces to be replaced by refined ones')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/additions', mesh%additions, &
         lbound(mesh%additions), ubound(mesh%additions), &
         comment = 'number of refined flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/refinement', mesh%refinement, &
         lbound(mesh%refinement), ubound(mesh%refinement), &
         comment = 'relative size of most refined flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_res', mesh%m_res, &
         lbound(mesh%m_res), ubound(mesh%m_res), &
         comment = 'poloidal mode number m in resonance at given flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/res_ind', mesh%res_ind, &
         lbound(mesh%res_ind), ubound(mesh%res_ind), &
         comment = 'flux surface index in resonance with given poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi_res', mesh%psi_res, &
         lbound(mesh%psi_res), ubound(mesh%psi_res), &
         comment = 'poloidal flux in resonance with given poloidal mode number', unit = 'Mx')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/res_modes', mesh%res_modes, &
         lbound(mesh%res_modes), ubound(mesh%res_modes), &
         comment = 'poloidal modes that are expected to actually be in resonance')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad_norm_res', mesh%rad_norm_res, &
         lbound(mesh%rad_norm_res), ubound(mesh%rad_norm_res), unit = '1', &
         comment = 'normalized minor radius (along X-O line) in resonance with given poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kp_max', mesh%kp_max, &
         lbound(mesh%kp_max), ubound(mesh%kp_max), &
         comment = 'number of nodes on flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kp_low', mesh%kp_low, &
         lbound(mesh%kp_low), ubound(mesh%kp_low), &
         comment = 'index of last node on previous flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kt_max', mesh%kt_max, &
         lbound(mesh%kt_max), ubound(mesh%kt_max), &
         comment = 'number of triangles on flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/kt_low', mesh%kt_low, &
         lbound(mesh%kt_low), ubound(mesh%kt_low), &
         comment = 'index of last triangle on previous flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/theta_offset', mesh%theta_offset, &
         lbound(mesh%theta_offset), ubound(mesh%theta_offset), unit = 'rad', &
         comment = 'flux poloidal angle at zero geometric poloidal angle between flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/node_R', mesh%node_R, &
         lbound(mesh%node_R), ubound(mesh%node_R), &
         comment = 'R coordinates of mesh points', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/node_Z', mesh%node_Z, &
         lbound(mesh%node_Z), ubound(mesh%node_Z), &
         comment = 'Z coordinates of mesh points', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/node_theta_flux', mesh%node_theta_flux, &
         lbound(mesh%node_theta_flux), ubound(mesh%node_theta_flux), &
         comment = 'flux poloidal angles of mesh points', unit = 'rad')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/node_theta_geom', mesh%node_theta_geom, &
         lbound(mesh%node_theta_geom), ubound(mesh%node_theta_geom), &
         comment = 'geometric poloidal angles of mesh points', unit = 'rad')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_node', mesh%tri_node, &
         lbound(mesh%tri_node), ubound(mesh%tri_node), &
         comment = 'triangle node indices')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_node_F', mesh%tri_node_F, &
         lbound(mesh%tri_node_F), ubound(mesh%tri_node_F), &
         comment = 'local node index of node F')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/li', mesh%li, &
         lbound(mesh%li), ubound(mesh%li), &
         comment = 'local node indices of edge i, counter-clockwise')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/lo', mesh%lo, &
         lbound(mesh%lo), ubound(mesh%lo), &
         comment = 'local node indices of edge o, counter-clockwise')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/lf', mesh%lf, &
         lbound(mesh%lf), ubound(mesh%lf), &
         comment = 'local node indices of edge f, counter-clockwise')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ei', mesh%ei, &
         lbound(mesh%ei), ubound(mesh%ei), &
         comment = 'local edge index of edge i')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/eo', mesh%eo, &
         lbound(mesh%eo), ubound(mesh%eo), &
         comment = 'local edge index of edge o')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ef', mesh%ef, &
         lbound(mesh%ef), ubound(mesh%ef), &
         comment = 'local edge index of edge f')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/orient', orient, &
         lbound(orient), ubound(orient), &
         comment = 'triangle orientation: true (1) if edge f lies on outer flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/adj_tri', mesh%adj_tri, &
         lbound(mesh%adj_tri), ubound(mesh%adj_tri), &
         comment = 'adjacent triangle indices')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/adj_edge', mesh%adj_edge, &
         lbound(mesh%adj_edge), ubound(mesh%adj_edge), &
         comment = 'adjacent triangle edge indices')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_map2kedge', mesh%edge_map2global, &
         lbound(mesh%edge_map2global), ubound(mesh%edge_map2global), &
         comment = 'mapping of triangle & local edge index to global edge index')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_map2ktri', mesh%edge_map2ktri, &
         lbound(mesh%edge_map2ktri), ubound(mesh%edge_map2ktri), &
         comment = 'mapping of global edge index to triangle index')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_map2ke', mesh%edge_map2ke, &
         lbound(mesh%edge_map2ke), ubound(mesh%edge_map2ke), &
         comment = 'mapping of global edge index to local edge index')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R_Omega', mesh%R_Omega, &
         lbound(mesh%R_Omega), ubound(mesh%R_Omega), &
         comment = 'R coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z_Omega', mesh%Z_Omega, &
         lbound(mesh%Z_Omega), ubound(mesh%Z_Omega), &
         comment = 'Z coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/area', mesh%area, &
         lbound(mesh%area), ubound(mesh%area), &
         comment = 'triangle area', unit = 'cm^2')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac, &
         lbound(mesh%gpec_jacfac), ubound(mesh%gpec_jacfac), unit = 'cm^2', &
         comment = 'Jacobian surface factor between flux surfaces')
    call h5_close(h5id_root)
  end subroutine mesh_write

  subroutine write_meshfile_for_freefem
    integer :: fid, kpoi, ktri, kp

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace')
    write (fid, '(3(1x, i0))') mesh%npoint, mesh%ntri, mesh%kp_max(mesh%nflux) - 1
    do kpoi = 1, mesh%kp_low(mesh%nflux + 1)
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh%node_R(kpoi), mesh%node_Z(kpoi), 0
    end do
    do kpoi = mesh%kp_low(mesh%nflux + 1) + 1, mesh%npoint
       write (fid, '(2(1x, es23.15e3), 1x, i0)') &
            mesh%node_R(kpoi), mesh%node_Z(kpoi), 1
    end do
    do ktri = 1, mesh%ntri
       write (fid, '(4(1x, i0))') mesh%tri_node(:, ktri), 0
    end do
    do kp = 1, mesh%kp_max(mesh%nflux) - 1
       write (fid, '(4(1x, i0))') mesh%kp_low(mesh%nflux) + kp, mesh%kp_low(mesh%nflux) + kp + 1, 1
    end do
    close(fid)
  end subroutine write_meshfile_for_freefem

  subroutine mesh_read(mesh, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use magdif_conf, only: conf
    type(mesh_t), intent(inout) :: mesh
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer, allocatable :: orient(:)

    call magdif_mesh_destructor(mesh)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_O', mesh%R_O)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_O', mesh%Z_O)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_X', mesh%R_X)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_X', mesh%Z_X)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/theta0', mesh%theta0)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_min', mesh%R_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_min', mesh%Z_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_max', mesh%R_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_max', mesh%Z_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nflux', mesh%nflux)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/npoint', mesh%npoint)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/ntri', mesh%ntri)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nedge', mesh%nedge)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/n', mesh%n)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_res_min', mesh%m_res_min)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_res_max', mesh%m_res_max)
    ! TODO: allocate deferred-shape arrays in hdf5_tools and skip allocation here
    allocate(mesh%deletions(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%additions(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%refinement(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%m_res(mesh%nflux))
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%psi_res(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%rad_norm_res(mesh%m_res_min:mesh%m_res_max))
    if (conf%kilca_scale_factor /= 0) then
       allocate(mesh%res_modes(1))
    else
       allocate(mesh%res_modes(1:mesh%m_res_max-mesh%m_res_min))
    end if
    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux + 1))
    allocate(mesh%kt_low(mesh%nflux + 1))
    allocate(mesh%theta_offset(mesh%nflux))
    allocate(mesh%node_R(mesh%npoint))
    allocate(mesh%node_Z(mesh%npoint))
    allocate(mesh%node_theta_flux(mesh%npoint))
    allocate(mesh%node_theta_geom(mesh%npoint))
    allocate(mesh%tri_node(3, mesh%ntri))
    allocate(mesh%tri_node_F(mesh%ntri))
    allocate(mesh%li(2, mesh%ntri))
    allocate(mesh%lo(2, mesh%ntri))
    allocate(mesh%lf(2, mesh%ntri))
    allocate(mesh%ei(mesh%ntri))
    allocate(mesh%eo(mesh%ntri))
    allocate(mesh%ef(mesh%ntri))
    allocate(mesh%orient(mesh%ntri))
    allocate(orient(mesh%ntri))
    allocate(mesh%adj_tri(3, mesh%ntri))
    allocate(mesh%adj_edge(3, mesh%ntri))
    allocate(mesh%edge_map2global(3, mesh%ntri))
    allocate(mesh%edge_map2ktri(2, mesh%nedge))
    allocate(mesh%edge_map2ke(2, mesh%nedge))
    allocate(mesh%gpec_jacfac(-16:16, mesh%nflux))
    allocate(mesh%area(mesh%ntri))
    allocate(mesh%R_Omega(mesh%ntri))
    allocate(mesh%Z_Omega(mesh%ntri))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/deletions', mesh%deletions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/additions', mesh%additions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/refinement', mesh%refinement)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/m_res', mesh%m_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_ind', mesh%res_ind)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi_res', mesh%psi_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rad_norm_res', mesh%rad_norm_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_modes', mesh%res_modes)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kp_max', mesh%kp_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kp_low', mesh%kp_low)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kt_max', mesh%kt_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kt_low', mesh%kt_low)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/theta_offset', mesh%theta_offset)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_R', mesh%node_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_Z', mesh%node_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_theta_flux', mesh%node_theta_flux)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_theta_geom', mesh%node_theta_geom)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_node', mesh%tri_node)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_node_F', mesh%tri_node_F)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/li', mesh%li)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/lo', mesh%lo)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/lf', mesh%lf)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/ei', mesh%ei)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/eo', mesh%eo)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/ef', mesh%ef)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/orient', orient)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/adj_tri', mesh%adj_tri)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/adj_edge', mesh%adj_edge)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_map2kedge', mesh%edge_map2global)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_map2ktri', mesh%edge_map2ktri)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_map2ke', mesh%edge_map2ke)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_Omega', mesh%R_Omega)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_Omega', mesh%Z_Omega)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/area', mesh%area)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac)
    call h5_close(h5id_root)
    where (orient == 1)
       mesh%orient = .true.
    elsewhere
       mesh%orient = .false.
    end where
    deallocate(orient)
  end subroutine mesh_read

  subroutine magdif_mesh_destructor(this)
    type(mesh_t), intent(inout) :: this
    if (allocated(this%deletions)) deallocate(this%deletions)
    if (allocated(this%additions)) deallocate(this%additions)
    if (allocated(this%refinement)) deallocate(this%refinement)
    if (allocated(this%m_res)) deallocate(this%m_res)
    if (allocated(this%res_ind)) deallocate(this%res_ind)
    if (allocated(this%psi_res)) deallocate(this%psi_res)
    if (allocated(this%rad_norm_res)) deallocate(this%rad_norm_res)
    if (allocated(this%res_modes)) deallocate(this%res_modes)
    if (allocated(this%kp_max)) deallocate(this%kp_max)
    if (allocated(this%kt_max)) deallocate(this%kt_max)
    if (allocated(this%kp_low)) deallocate(this%kp_low)
    if (allocated(this%kt_low)) deallocate(this%kt_low)
    if (allocated(this%theta_offset)) deallocate(this%theta_offset)
    if (allocated(this%node_R)) deallocate(this%node_R)
    if (allocated(this%node_Z)) deallocate(this%node_Z)
    if (allocated(this%node_theta_flux)) deallocate(this%node_theta_flux)
    if (allocated(this%node_theta_geom)) deallocate(this%node_theta_geom)
    if (allocated(this%tri_node)) deallocate(this%tri_node)
    if (allocated(this%tri_node_F)) deallocate(this%tri_node_F)
    if (allocated(this%li)) deallocate(this%li)
    if (allocated(this%lo)) deallocate(this%lo)
    if (allocated(this%lf)) deallocate(this%lf)
    if (allocated(this%ei)) deallocate(this%ei)
    if (allocated(this%eo)) deallocate(this%eo)
    if (allocated(this%ef)) deallocate(this%ef)
    if (allocated(this%orient)) deallocate(this%orient)
    if (allocated(this%adj_tri)) deallocate(this%adj_tri)
    if (allocated(this%adj_edge)) deallocate(this%adj_edge)
    if (allocated(this%edge_map2global)) deallocate(this%edge_map2global)
    if (allocated(this%edge_map2ktri)) deallocate(this%edge_map2ktri)
    if (allocated(this%edge_map2ke)) deallocate(this%edge_map2ke)
    if (allocated(this%gpec_jacfac)) deallocate(this%gpec_jacfac)
    if (allocated(this%area)) deallocate(this%area)
    if (allocated(this%R_Omega)) deallocate(this%R_Omega)
    if (allocated(this%Z_Omega)) deallocate(this%Z_Omega)
  end subroutine magdif_mesh_destructor

  subroutine compare_gpec_coordinates
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, nf90_inq_varid, &
         nf90_inquire_dimension, nf90_get_var, nf90_close, nf90_global, nf90_get_att
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext, psipol_max
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: conf, log, datafile
    character(len = *), parameter :: dataset = 'debug_GPEC'
    character(len = 1024) :: filename
    logical :: file_exists
    integer(HID_T) :: h5id_root
    integer :: ncid_file, ncid, nrad, npol, krad, kpol, idum, fid
    real(dp) :: dum, B0_R, B0_Z, unit_normal(0:1), q
    real(dp), allocatable :: psi(:), theta(:), R(:, :), Z(:, :), theta_shift(:), &
         xi_n(:, :, :), jac(:, :), sqrt_g(:, :)
    complex(dp), allocatable :: xi_n_R(:), xi_n_Z(:)

    write (filename, '("gpec_profile_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (log%msg, '("File ", a, " found, performing GPEC coordinate comparison.")') &
         trim(filename)
    if (log%info) call log%write
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "psi_n", ncid))
    call check_error("nf90_inquire_dimension", &
         nf90_inquire_dimension(ncid_file, ncid, len = nrad))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "theta_dcon", ncid))
    call check_error("nf90_inquire_dimension", &
         nf90_inquire_dimension(ncid_file, ncid, len = npol))
    allocate(psi(nrad), theta(npol), R(nrad, npol), Z(nrad, npol), xi_n(nrad, npol, 0:1))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "psi_n", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, psi))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "theta_dcon", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, theta))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "R", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, R))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "z", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Z))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "xi_n_fun", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, xi_n))
    call check_error("nf90_close", nf90_close(ncid_file))
    psi(:) = psipol_max * psi
    theta(:) = 2d0 * pi * theta
    R(:, :) = R * 1d2
    Z(:, :) = Z * 1d2
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/psi', psi, lbound(psi), ubound(psi), &
         unit = 'Mx', comment = 'poloidal flux (shifted to zero at axis)')
    call h5_add(h5id_root, dataset // '/theta', theta, lbound(theta), ubound(theta), &
         unit = 'rad', comment = 'flux poloidal angle')
    call h5_add(h5id_root, dataset // '/R_GPEC', R, lbound(R), ubound(R), &
         unit = 'cm', comment = 'R(psi, theta) from GPEC')
    call h5_add(h5id_root, dataset // '/Z_GPEC', Z, lbound(Z), ubound(Z), &
         unit = 'cm', comment = 'Z(psi, theta) from GPEC')
    allocate(theta_shift(nrad))
    do krad = 1, nrad
       theta_shift(krad) = theta_offset(psi(krad) + fs%psi(0))
    end do
    do kpol = 1, npol
       do krad = 1, nrad
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol) + theta_shift(krad), &
               dum, dum, dum, dum, dum, R(krad, kpol), dum, dum, Z(krad, kpol), dum, dum)
       end do
    end do
    allocate(xi_n_R(npol), xi_n_Z(npol))
    krad = nrad
    do kpol = 1, npol
       call field(R(krad, kpol), 0d0, Z(krad, kpol), B0_R, dum, B0_Z, &
            dum, dum, dum, dum, dum, dum, dum, dum, dum)
       unit_normal = [R(krad, kpol) * B0_Z, -R(krad, kpol) * B0_R] * equil%cocos%sgn_dpsi
       unit_normal = unit_normal / norm2(unit_normal)
       xi_n_R(kpol) = unit_normal(0) * 1d2 * cmplx(xi_n(krad, kpol, 0), xi_n(krad, kpol, 1), dp)
       xi_n_Z(kpol) = unit_normal(1) * 1d2 * cmplx(xi_n(krad, kpol, 0), xi_n(krad, kpol, 1), dp)
    end do
    if (allocated(theta_shift)) deallocate(theta_shift)
    call h5_add(h5id_root, dataset // '/R', R, lbound(R), ubound(R), &
         unit = 'cm', comment = 'R(psi, theta)')
    call h5_add(h5id_root, dataset // '/Z', Z, lbound(Z), ubound(Z), &
         unit = 'cm', comment = 'Z(psi, theta)')
    call h5_add(h5id_root, dataset // '/xi_n_R', xi_n_R, lbound(xi_n_R), ubound(xi_n_R), &
         unit = 'cm', comment = 'Radial component of normal displacement xi_n(theta) at last flux surface')
    call h5_add(h5id_root, dataset // '/xi_n_Z', xi_n_Z, lbound(xi_n_Z), ubound(xi_n_Z), &
         unit = 'cm', comment = 'Axial component of normal displacement xi_n(theta) at last flux surface')
    call h5_close(h5id_root)
    if (allocated(psi)) deallocate(psi)
    if (allocated(theta)) deallocate(theta)
    if (allocated(R)) deallocate(R)
    if (allocated(Z)) deallocate(Z)
    if (allocated(xi_n)) deallocate(xi_n)
    if (allocated(xi_n_R)) deallocate(xi_n_R)
    if (allocated(xi_n_Z)) deallocate(xi_n_Z)
    filename = '2d.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (filename, '("dcon_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (log%msg, '("Files ", a, " and 2d.out found, performing GPEC Jacobian comparison.")') &
         trim(filename)
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_get_att", nf90_get_att(ncid_file, nf90_global, 'mpsi', nrad))
    call check_error("nf90_get_att", nf90_get_att(ncid_file, nf90_global, 'mtheta', npol))
    call check_error("nf90_close", nf90_close(ncid_file))
    nrad = nrad + 1
    allocate(psi(nrad), theta(npol), theta_shift(nrad), jac(nrad, npol), sqrt_g(nrad, npol), &
         R(nrad, npol), Z(nrad, npol))
    open(newunit = fid, file = '2d.out', status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
       read (fid, '(6x, i3, 6x, e24.16)') idum, psi(krad)
       psi(krad) = psi(krad) * psipol_max
       theta_shift(krad) = theta_offset(psi(krad) + fs%psi(0))
       read (fid, *)
       read (fid, *)
       read (fid, *)
       do kpol = 1, npol
          read (fid, '(i6, 1p, 8e24.16)') idum, theta(kpol), dum, dum, dum, dum, &
               R(krad, kpol), Z(krad, kpol), jac(krad, kpol)
          theta(kpol) = theta(kpol) * 2d0 * pi
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol) + theta_shift(krad), &
               q, dum, sqrt_g(krad, kpol), dum, dum, dum, dum, dum, dum, dum, dum)
          sqrt_g(krad, kpol) = sqrt_g(krad, kpol) * abs(q)
       end do
       read (fid, *)  ! skip theta = 2 pi
       read (fid, *)
       read (fid, *)
       read (fid, *)
    end do
    R(:, :) = R * 1d2  ! m to cm
    Z(:, :) = Z * 1d2  ! m to cm
    jac(:, :) = jac * 1d-2  ! m per T to cm per G
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/jac/')
    call h5_add(h5id_root, dataset // '/jac/psi', psi, lbound(psi), ubound(psi), &
         unit = 'Mx', comment = 'poloidal flux (shifted to zero at axis)')
    call h5_add(h5id_root, dataset // '/jac/theta', theta, lbound(theta), ubound(theta), &
         unit = 'rad', comment = 'flux poloidal angle')
    call h5_add(h5id_root, dataset // '/jac/R', R, lbound(R), ubound(R), &
         unit = 'cm', comment = 'R(psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/Z', Z, lbound(Z), ubound(Z), &
         unit = 'cm', comment = 'Z(psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/sqrt_g', sqrt_g, lbound(sqrt_g), ubound(sqrt_g), &
         unit = 'cm', comment = 'MEPHIT Jacobian at (psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/jac', jac, lbound(jac), ubound(jac), &
         unit = 'cm', comment = 'GPEC Jacobian at (psi, theta)')
    call h5_close(h5id_root)
    deallocate(psi, theta, theta_shift, jac, sqrt_g, R, Z)

  contains
    subroutine check_error(funcname, status)
      character(len = *), intent(in) :: funcname
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
         write (log%msg, '(a, " returned error ", i0)') funcname, status
         if (log%err) call log%write
         error stop
      end if
    end subroutine check_error
  end subroutine compare_gpec_coordinates

  subroutine write_illustration_data(npsi, ntheta, nrad, npol)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext, psipol_max
    use constants, only: pi  ! orbit_mod.f90
    use magdif_util, only: linspace
    integer, intent(in) :: npsi, ntheta, nrad, npol
    integer :: fid, krad, kpol
    real(dp) :: dum
    real(dp), allocatable :: psi(:), theta(:), theta_shift(:), R(:, :), Z(:, :)

    open(newunit = fid, file = 'illustration.asc', status = 'replace', form = 'formatted')
    write (fid, '(6(1x, i0))') npsi, ntheta, nrad, npol, equil%nbbbs, equil%limitr
    allocate(psi(npsi), theta(npol), theta_shift(npsi), R(npol, npsi), Z(npol, npsi))
    psi(:) = linspace(0d0, psipol_max, npsi, 1, 1)
    theta(:) = linspace(0d0, 2d0 * pi, npol, 0, 1)
    do krad = 1, npsi
       theta_shift(krad) = theta_offset(psi(krad) + fs%psi(0))
    end do
    do krad = 1, npsi
       do kpol = 1, npol
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol) + theta_shift(krad), &
               dum, dum, dum, dum, dum, R(kpol, krad), dum, dum, Z(kpol, krad), dum, dum)
       end do
    end do
    do krad = 1, npsi
       do kpol = 1, npol
          write (fid, '(es24.16e3, 1x, es24.16e3)') R(kpol, krad), Z(kpol, krad)
       end do
    end do
    deallocate(psi, theta, theta_shift, R, Z)
    allocate(psi(nrad), theta(ntheta), theta_shift(nrad), R(nrad, ntheta), Z(nrad, ntheta))
    psi(:) = linspace(0d0, psipol_max, nrad, 0, 0)
    theta(:) = linspace(0d0, 2d0 * pi, ntheta, 0, 1)
    do krad = 1, nrad
       theta_shift(krad) = theta_offset(psi(krad) + fs%psi(0))
    end do
    do kpol = 1, ntheta
       do krad = 2, nrad
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol) + theta_shift(krad), &
               dum, dum, dum, dum, dum, R(krad, kpol), dum, dum, Z(krad, kpol), dum, dum)
       end do
    end do
    do kpol = 1, ntheta
       write (fid, '(es24.16e3, 1x, es24.16e3)') mesh%R_O, mesh%Z_O
       do krad = 2, nrad
          write (fid, '(es24.16e3, 1x, es24.16e3)') R(krad, kpol), Z(krad, kpol)
       end do
    end do
    deallocate(psi, theta, theta_shift, R, Z)
    do kpol = 1, equil%nbbbs
       write (fid, '(es24.16e3, 1x, es24.16e3)') equil%rbbbs(kpol), equil%zbbbs(kpol)
    end do
    do kpol = 1, equil%limitr
       write (fid, '(es24.16e3, 1x, es24.16e3)') equil%rlim(kpol), equil%zlim(kpol)
    end do
    close(fid)
  end subroutine write_illustration_data

  subroutine init_flux_variables
    use magdif_conf, only: conf, log, pres_prof_eps, pres_prof_par, pres_prof_geqdsk, &
         q_prof_flux, q_prof_rot, q_prof_geqdsk
    integer :: kf

    ! initialize fluxvar with equidistant psi values
    call fluxvar%init(4, equil%psi_eqd)

    select case (conf%pres_prof)
    case (pres_prof_eps)
       call compute_pres_prof_eps
    case (pres_prof_par)
       call compute_pres_prof_par
    case (pres_prof_geqdsk)
       call compute_pres_prof_geqdsk
    case default
       write (log%msg, '("unknown pressure profile selection", i0)') conf%pres_prof
       if (log%err) call log%write
       error stop
    end select

    select case (conf%q_prof)
    case (q_prof_flux)
       call compute_safety_factor_flux
    case (q_prof_rot)
       call compute_safety_factor_rot
    case (q_prof_geqdsk)
       call compute_safety_factor_geqdsk
    case default
       write (log%msg, '("unknown q profile selection: ", i0)') conf%q_prof
       if (log%err) call log%write
       error stop
    end select

    fs%F = [(fluxvar%interp(equil%fpol, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%F = [(fluxvar%interp(equil%fpol, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    fs_half%FdF_dpsi = [(fluxvar%interp(equil%ffprim, fs_half%psi(kf)), &
         kf = 1, mesh%nflux)]
  end subroutine init_flux_variables

  subroutine compute_pres_prof_eps
    use constants, only: ev2erg  ! orbit_mod.f90
    use magdif_conf, only: conf, log
    real(dp) :: ddens_dpsi, dtemp_dpsi, psi_int, psi_ext

    ! Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
    real(dp) :: dens(0:mesh%nflux)

    ! Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
    real(dp) :: temp(0:mesh%nflux)

    dens = 0d0
    temp = 0d0
    if (equil%cocos%sgn_dpsi == -1) then
       psi_ext = minval(equil%psirz)
       psi_int = maxval(equil%psirz)
    else
       psi_ext = maxval(equil%psirz)
       psi_int = minval(equil%psirz)
    end if
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
  end subroutine compute_pres_prof_eps

  subroutine compute_pres_prof_par
    use constants, only: ev2erg  ! orbit_mod.f90
    use magdif_conf, only: conf
    real(dp) :: ddens_dpsi, dtemp_dpsi, psi_int, psi_ext

    ! Density \f$ \frac{N}{V} \f$ on flux surface in cm^-3.
    real(dp) :: dens(0:mesh%nflux)

    ! Temperature \f$ T \f$ on flux surface with \f$ k_{\mathrm{B}} T \f$ in eV.
    real(dp) :: temp(0:mesh%nflux)

    dens = 0d0
    temp = 0d0
    if (equil%cocos%sgn_dpsi == -1) then
       psi_ext = minval(equil%psirz)
       psi_int = maxval(equil%psirz)
    else
       psi_ext = maxval(equil%psirz)
       psi_int = minval(equil%psirz)
    end if
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
  end subroutine compute_pres_prof_par

  subroutine compute_pres_prof_geqdsk
    integer :: kf
    fs%p = [(fluxvar%interp(equil%pres, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs%dp_dpsi = [(fluxvar%interp(equil%pprime, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%p = [(fluxvar%interp(equil%pres, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    fs_half%dp_dpsi = [(fluxvar%interp(equil%pprime, fs_half%psi(kf)), kf = 1, mesh%nflux)]
  end subroutine compute_pres_prof_geqdsk

  subroutine compute_safety_factor_flux
    use constants, only: pi  ! orbit_mod.f90
    use magdif_util, only: flux_func
    integer :: kf, kt, ktri
    type(flux_func) :: psi_interpolator

    fs_half%q = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          fs_half%q(kf) = fs_half%q(kf) + B0phi_Omega(ktri) * mesh%area(ktri)
       end do
       fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    call psi_interpolator%init(4, fs_half%psi)
    ! Lagrange polynomial extrapolation for values at separatrix and magnetic axis
    fs%q = [(psi_interpolator%interp(fs_half%q, fs%psi(kf)), kf = 0, mesh%nflux)]
  end subroutine compute_safety_factor_flux

  subroutine compute_safety_factor_rot
    use magdif_util, only: flux_func
    use magdata_in_symfluxcoor_mod, only: psipol_max, psisurf, qsaf
    integer :: kf
    type(flux_func) :: psi_interpolator

    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    call psi_interpolator%init(4, psisurf(1:) * psipol_max + fs%psi(0))
    ! Lagrange polynomial extrapolation for value at magnetic axis
    fs%q = [(psi_interpolator%interp(qsaf, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%q = [(psi_interpolator%interp(qsaf, fs_half%psi(kf)), kf = 1, mesh%nflux)]
  end subroutine compute_safety_factor_rot

  subroutine compute_safety_factor_geqdsk
    integer :: kf

    fs%q = [(fluxvar%interp(equil%qpsi, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%q = [(fluxvar%interp(equil%qpsi, fs_half%psi(kf)), kf = 1, mesh%nflux)]
  end subroutine compute_safety_factor_geqdsk

  subroutine check_resonance_positions
    use magdif_conf, only: log
    integer :: m, kf
    real(dp), dimension(mesh%nflux) :: abs_err

    do m = mesh%m_res_max, mesh%m_res_min, -1
       abs_err = [(abs(abs(fs_half%q(kf)) - dble(m) / dble(mesh%n)), kf = 1, mesh%nflux)]
       kf = minloc(abs_err, 1)
       if (kf /= mesh%res_ind(m)) then
          write (log%msg, '("m = ", i0, ": q is resonant at index ", i0, ' // &
               '", but psi is resonant at index ", i0)') m, kf, mesh%res_ind(m)
          if (log%warn) call log%write
       end if
    end do
  end subroutine check_resonance_positions

  subroutine check_safety_factor
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: psisurf, qsaf
    use magdif_conf, only: datafile
    character(len = *), parameter :: grp = 'debug_q'
    integer(HID_T) :: h5id_root

    integer :: kf, kt, ktri
    real(dp) :: step_q(mesh%nflux), step_psi_norm(mesh%nflux), eqd_psi_norm(equil%nw)

    step_q(:) = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          step_q(kf) = step_q(kf) + B0phi_Omega(ktri) * mesh%area(ktri)
       end do
       step_q(kf) = step_q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    step_psi_norm(:) = (fs_half%psi - fs%psi(0)) / (fs%psi(mesh%nflux) - fs%psi(0))
    eqd_psi_norm(:) = (equil%psi_eqd - equil%simag) / (equil%sibry - equil%simag)
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/step_psi_norm', step_psi_norm, &
         lbound(step_psi_norm), ubound(step_psi_norm), &
         comment = 'normalized poloidal flux (stepped)')
    call h5_add(h5id_root, grp // '/step_q', step_q, lbound(step_q), ubound(step_q), &
         comment = 'safety factor (stepped)')
    call h5_add(h5id_root, grp // '/GEQDSK_psi_norm', eqd_psi_norm, &
         lbound(eqd_psi_norm), ubound(eqd_psi_norm), &
         comment = 'normalized poloidal flux (GEQDSK)')
    call h5_add(h5id_root, grp // '/GEQDSK_q', equil%qpsi, lbound(equil%qpsi), ubound(equil%qpsi), &
         comment = 'safety factor (GEQDSK)')
    call h5_add(h5id_root, grp // '/RK_psi_norm', psisurf(1:), &
         lbound(psisurf(1:)), ubound(psisurf(1:)), &
         comment = 'normalized poloidal flux (Runge-Kutta field line integration)')
    call h5_add(h5id_root, grp // '/RK_q', qsaf, lbound(qsaf), ubound(qsaf), &
         comment = 'safety factor (Runge-Kutta field line integration)')
    call h5_close(h5id_root)
  end subroutine check_safety_factor

  subroutine cache_equilibrium_field
    real(dp) :: R, Z, B0_R, B0_phi, B0_Z, dum
    integer :: kf, kt, ktri, ke
    integer :: base, tip
    real(dp) :: n_r, n_z

    allocate(B0R(3, mesh%ntri))
    allocate(B0phi(3, mesh%ntri))
    allocate(B0Z(3, mesh%ntri))
    allocate(B0R_Omega(mesh%ntri))
    allocate(B0phi_Omega(mesh%ntri))
    allocate(B0Z_Omega(mesh%ntri))
    allocate(B0flux(3, mesh%ntri))
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          do ke = 1, 3
             base = mesh%tri_node(ke, ktri)
             tip = mesh%tri_node(mod(ke, 3) + 1, ktri)
             R = (mesh%node_R(base) + mesh%node_R(tip)) * 0.5d0
             Z = (mesh%node_Z(base) + mesh%node_Z(tip)) * 0.5d0
             call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dum, dum, dum, &
                  dum, dum, dum, dum, dum, dum)
             B0R(ke, ktri) = B0_R
             B0phi(ke, ktri) = B0_phi
             B0Z(ke, ktri) = B0_Z
             n_R = mesh%node_Z(tip) - mesh%node_Z(base)
             n_Z = mesh%node_R(base) - mesh%node_R(tip)
             B0flux(ke, ktri) = R * (B0_R * n_R + B0_Z * n_Z)
          end do
          call ring_centered_avg_coord(ktri, R, Z)
          call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dum, dum, dum, &
               dum, dum, dum, dum, dum, dum)
          B0R_Omega(ktri) = B0_R
          B0phi_Omega(ktri) = B0_phi
          B0Z_Omega(ktri) = B0_Z
       end do
    end do
  end subroutine cache_equilibrium_field

  subroutine compute_j0phi
    use magdif_conf, only: conf, log, curr_prof_ps, curr_prof_rot, curr_prof_geqdsk
    real(dp) :: plot_j0phi(mesh%ntri)

    allocate(j0phi(3, mesh%ntri))
    j0phi = 0d0
    select case (conf%curr_prof)
    case (curr_prof_ps)
       call compute_j0phi_ps(plot_j0phi)
    case (curr_prof_rot)
       call compute_j0phi_rot(plot_j0phi)
    case (curr_prof_geqdsk)
       call compute_j0phi_geqdsk(plot_j0phi)
    case default
       write (log%msg, '("unknown current profile selection: ", i0)') conf%curr_prof
       if (log%err) call log%write
       error stop
    end select
    ! TODO: write out plot_j0phi when edge_cache type is working
    call check_redundant_edges(j0phi, 'j0phi')
  end subroutine compute_j0phi

  subroutine compute_j0phi_ps(plot_j0phi)
    use magdif_util, only: clight
    real(dp), intent(out) :: plot_j0phi(:)
    integer :: kf, kt, ktri
    real(dp) :: R
    real(dp) :: Btor2
    real(dp), dimension(mesh%nflux) :: B2avg, B2avg_half

    B2avg = 0d0
    B2avg_half = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          R = sum(mesh%node_R(mesh%lf(:, ktri))) * 0.5d0
          B2avg(kf) = B2avg(kf) + B0R(mesh%ef(ktri), ktri) ** 2 + &
               B0phi(mesh%ef(ktri), ktri) ** 2 + B0Z(mesh%ef(ktri), ktri) ** 2
          R = sum(mesh%node_R(mesh%li(:, ktri))) * 0.5d0
          B2avg_half(kf) = B2avg_half(kf) + B0R(mesh%ei(ktri), ktri) ** 2 + &
               B0phi(mesh%ei(ktri), ktri) ** 2 + B0Z(mesh%ei(ktri), ktri) ** 2
       end do
       B2avg(kf) = B2avg(kf) / mesh%kt_max(kf)
       B2avg_half(kf) = B2avg_half(kf) / mesh%kt_max(kf)
    end do
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ! edge f
          R = sum(mesh%node_R(mesh%lf(:, ktri))) * 0.5d0
          Btor2 = B0phi(ktri, mesh%ef(ktri)) ** 2
          if (kf == 1 .or. mesh%orient(ktri)) then
             j0phi(mesh%ef(ktri), ktri) = clight * R * fs%dp_dpsi(kf) * (1d0 - &
                  Btor2 / B2avg(kf))
          else
             j0phi(mesh%ef(ktri), ktri) = clight * R * fs%dp_dpsi(kf-1) * (1d0 - &
                  Btor2 / B2avg(kf-1))
          end if
          ! edge i
          R = sum(mesh%node_R(mesh%li(:, ktri))) * 0.5d0
          Btor2 = B0phi(mesh%ei(ktri), ktri) ** 2
          j0phi(mesh%ei(ktri), ktri) = clight * R * fs_half%dp_dpsi(kf) * (1d0 - &
               Btor2 / B2avg_half(kf))
          ! edge o
          R = sum(mesh%node_R(mesh%lo(:, ktri))) * 0.5d0
          Btor2 = B0phi(mesh%eo(ktri), ktri) ** 2
          j0phi(mesh%eo(ktri), ktri) = clight * R * fs_half%dp_dpsi(kf) * (1d0 - &
               Btor2 / B2avg_half(kf))
          ! centroid
          Btor2 = B0phi_Omega(ktri) ** 2
          plot_j0phi(ktri) = clight * mesh%R_Omega(ktri) * fs_half%dp_dpsi(kf) * (1d0 - &
               Btor2 / B2avg_half(kf))
       end do
    end do
  end subroutine compute_j0phi_ps

  subroutine compute_j0phi_rot(plot_j0phi)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_util, only: clight
    real(dp), intent(out) :: plot_j0phi(:)
    integer :: kf, kt, ktri
    real(dp) :: R, Z, dum, dB0R_dZ, dB0Z_dR

    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ! edge f
          R = sum(mesh%node_R(mesh%lf(:, ktri))) * 0.5d0
          Z = sum(mesh%node_Z(mesh%lf(:, ktri))) * 0.5d0
          call field(R, 0d0, Z, dum, dum, dum, dum, dum, dB0R_dZ, &
               dum, dum, dum, dB0Z_dR, dum, dum)
          j0phi(mesh%ef(ktri), ktri) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          ! edge i
          R = sum(mesh%node_R(mesh%li(:, ktri))) * 0.5d0
          Z = sum(mesh%node_Z(mesh%li(:, ktri))) * 0.5d0
          call field(R, 0d0, Z, dum, dum, dum, dum, dum, dB0R_dZ, &
               dum, dum, dum, dB0Z_dR, dum, dum)
          j0phi(mesh%ei(ktri), ktri) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          ! edge o
          R = sum(mesh%node_R(mesh%lo(:, ktri))) * 0.5d0
          Z = sum(mesh%node_Z(mesh%lo(:, ktri))) * 0.5d0
          call field(R, 0d0, Z, dum, dum, dum, dum, dum, dB0R_dZ, &
               dum, dum, dum, dB0Z_dR, dum, dum)
          j0phi(mesh%eo(ktri), ktri) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          ! centroid
          call field(mesh%R_Omega(ktri), 0d0, mesh%Z_Omega(ktri), dum, dum, dum, &
               dum, dum, dB0R_dZ, dum, dum, dum, dB0Z_dR, dum, dum)
          plot_j0phi(ktri) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
       end do
    end do
  end subroutine compute_j0phi_rot

  subroutine compute_j0phi_geqdsk(plot_j0phi)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_util, only: clight
    real(dp), intent(out) :: plot_j0phi(:)
    integer :: kf, kt, ktri
    real(dp) :: R

    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ! edge f
          R = sum(mesh%node_R(mesh%lf(:, ktri))) * 0.5d0
          if (kf == 1 .or. mesh%orient(ktri)) then
             j0phi(mesh%ef(ktri), ktri) = clight * (fs%dp_dpsi(kf) * R + &
                  0.25d0 / pi * fs%FdF_dpsi(kf) / R)
          else
             j0phi(mesh%ef(ktri), ktri) = clight * (fs%dp_dpsi(kf-1) * R + &
                  0.25d0 / pi * fs%FdF_dpsi(kf-1) / R)
          end if
          ! edge i
          R = sum(mesh%node_R(mesh%li(:, ktri))) * 0.5d0
          j0phi(mesh%ei(ktri), ktri) = clight * (fs_half%dp_dpsi(kf) * R + &
               0.25d0 / pi * fs_half%FdF_dpsi(kf) / R)
          ! edge o
          R = sum(mesh%node_R(mesh%lo(:, ktri))) * 0.5d0
          j0phi(mesh%eo(ktri), ktri) = clight * (fs_half%dp_dpsi(kf) * R + &
               0.25d0 / pi * fs_half%FdF_dpsi(kf) / R)
          ! centroid
          plot_j0phi(ktri) = clight * (fs_half%dp_dpsi(kf) * mesh%R_Omega(ktri) + &
               0.25d0 / pi * fs_half%FdF_dpsi(kf) / mesh%R_Omega(ktri))
       end do
    end do
  end subroutine compute_j0phi_geqdsk

  subroutine check_curr0
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use constants, only: pi  ! orbit_mod.f90
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use magdif_conf, only: datafile
    use magdif_util, only: clight, linspace
    character(len = *), parameter :: grp = 'debug_equil'
    integer, parameter :: ntheta = 512
    integer(HID_T) :: h5id_root
    integer :: kf, ktheta
    real(dp) :: R, Z, dum, B0_R, B0_phi, B0_Z, dB0R_dZ, dB0phi_dR, dB0phi_dZ, dB0Z_dR, &
         grad_psi(2), psi(2:equil%nw - 1), theta(ntheta)
    real(dp), dimension(ntheta, 2:equil%nw - 1) :: grad_p0, amp_lorentz, gs_lorentz, &
         amp_J0_R, amp_J0_phi, amp_J0_Z, gs_J0_R, gs_J0_phi, gs_J0_Z

    psi(:) = equil%psi_eqd(2:equil%nw-1) - equil%simag
    theta(:) = linspace(0d0, 2d0 * pi, ntheta, 0, 1)
    do kf = 2, equil%nw - 1
       do ktheta = 1, ntheta
          ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
          call magdata_in_symfluxcoord_ext(2, dum, psi(kf), theta(ktheta), &
               dum, dum, dum, dum, dum, R, dum, dum, Z, dum, dum)
          call field(R, 0d0, Z, B0_R, B0_phi, B0_Z, dum, dum, dB0R_dZ, dB0phi_dR, &
               dum, dB0phi_dZ, dB0Z_dR, dum, dum)
          ! left-hand side of iMHD force balance
          grad_psi = [R * B0_Z, -R * B0_R]
          grad_p0(ktheta, kf) = equil%pprime(kf) * dot_product(grad_psi, grad_psi) / &
               norm2(grad_psi)
          ! current density via Grad-Shafranov equation
          gs_J0_R(ktheta, kf) = 0.25d0 / pi * clight * equil%ffprim(kf) / equil%fpol(kf) * B0_R
          gs_J0_Z(ktheta, kf) = 0.25d0 / pi * clight * equil%ffprim(kf) / equil%fpol(kf) * B0_Z
          gs_J0_phi(ktheta, kf) = clight * (equil%pprime(kf) * R + &
               0.25d0 / pi * equil%ffprim(kf) / R)
          gs_lorentz(ktheta, kf) = 1d0 / norm2(grad_psi) / clight * dot_product(grad_psi, &
               [gs_J0_phi(ktheta, kf) * B0_Z - gs_J0_Z(ktheta, kf) * B0_phi, &
               gs_J0_R(ktheta, kf) * B0_phi - gs_J0_phi(ktheta, kf) * B0_R])
          ! current density via Ampere's equation
          amp_J0_R(ktheta, kf) = 0.25d0 / pi * clight * (-dB0phi_dZ)
          amp_J0_phi(ktheta, kf) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          amp_J0_Z(ktheta, kf) = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
          amp_lorentz(ktheta, kf) = 1d0 / norm2(grad_psi) / clight * dot_product(grad_psi, &
               [amp_J0_phi(ktheta, kf) * B0_Z - amp_J0_Z(ktheta, kf) * B0_phi, &
               amp_J0_R(ktheta, kf) * B0_phi - amp_J0_phi(ktheta, kf) * B0_R])
       end do
    end do
    psi(:) = psi / (equil%sibry - equil%simag)
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    call h5_add(h5id_root, grp // '/psi', psi, lbound(psi), ubound(psi), &
         unit = 'Mx', comment = 'poloidal flux')
    call h5_add(h5id_root, grp // '/theta', theta, lbound(theta), ubound(theta), &
         unit = 'rad', comment = 'flux poloidal angle')
    call h5_add(h5id_root, grp // '/grad_p0', grad_p0, lbound(grad_p0), ubound(grad_p0), &
         unit = 'dyn cm^-3', comment = 'pressure gradient')
    call h5_add(h5id_root, grp // '/GS_Lorentz', gs_lorentz, lbound(gs_lorentz), ubound(gs_lorentz), &
         unit = 'dyn cm^-3', comment = 'Lorentz force density via Grad-Shafranov equation')
    call h5_add(h5id_root, grp // '/Ampere_Lorentz', amp_lorentz, lbound(amp_lorentz), ubound(amp_lorentz), &
         unit = 'dyn cm^-3', comment = 'Lorentz force density via Ampere equation')
    call h5_add(h5id_root, grp // '/GS_j0_R', gs_J0_R, &
         lbound(gs_j0_R), ubound(gs_J0_R), unit = 'statA cm^-2', &
         comment = 'R component of equilibrium current density via Grad-Shafranov equation')
    call h5_add(h5id_root, grp // '/GS_j0_phi', gs_J0_phi, &
         lbound(gs_j0_phi), ubound(gs_J0_phi), unit = 'statA cm^-2', &
         comment = 'physical phi component of equilibrium current density via Grad-Shafranov equation')
    call h5_add(h5id_root, grp // '/GS_j0_Z', gs_J0_Z, &
         lbound(gs_j0_Z), ubound(gs_J0_Z), unit = 'statA cm^-2', &
         comment = 'Z component of equilibrium current density via Grad-Shafranov equation')
    call h5_add(h5id_root, grp // '/Ampere_j0_R', amp_J0_R, &
         lbound(amp_j0_R), ubound(amp_J0_R), unit = 'statA cm^-2', &
         comment = 'R component of equilibrium current density via Ampere equation')
    call h5_add(h5id_root, grp // '/Ampere_j0_phi', amp_J0_phi, &
         lbound(amp_j0_phi), ubound(amp_J0_phi), unit = 'statA cm^-2', &
         comment = 'physical phi component of equilibrium current density via Ampere equation')
    call h5_add(h5id_root, grp // '/Ampere_j0_Z', amp_J0_Z, &
         lbound(amp_j0_Z), ubound(amp_J0_Z), unit = 'statA cm^-2', &
         comment = 'Z component of equilibrium current density via Ampere equation')
    call h5_close(h5id_root)
  end subroutine check_curr0

  subroutine check_redundant_edges(cache, name)
    use magdif_conf, only: log
    real(dp), intent(in) :: cache(:, :)
    character(len = *), intent(in) :: name
    integer :: kedge, ktri, ktri_adj, ke, ke_adj
    logical :: inconsistent
    real(dp), parameter :: eps = epsilon(1d0), small = tiny(0d0)

    do kedge = 1, mesh%nedge
       ktri = mesh%edge_map2ktri(1, kedge)
       ktri_adj = mesh%edge_map2ktri(2, kedge)
       ke = mesh%edge_map2ke(1, kedge)
       ke_adj = mesh%edge_map2ke(2, kedge)
       if (ktri_adj <= 0) cycle
       inconsistent = .false.
       if (abs(cache(ke, ktri)) < small) then
          inconsistent = inconsistent .or. abs(cache(ke_adj, ktri_adj)) >= small
       else
          inconsistent = inconsistent .or. eps < abs(1d0 - &
               cache(ke_adj, ktri_adj) / cache(ke, ktri))
       end if
       if (inconsistent) then
          write (log%msg, '("inconsistent redundant edges: ", ' // &
               'a, "(", i0, ", ", i0, ") = ", es24.16e3, ", ", ' // &
               'a, "(", i0, ", ", i0, ") = ", es24.16e3)') &
               trim(name), ke, ktri, cache(ke, ktri), &
               trim(name), ke_adj, ktri_adj, cache(ke_adj, ktri_adj)
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

end module magdif_mesh
