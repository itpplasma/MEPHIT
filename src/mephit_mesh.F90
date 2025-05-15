module mephit_mesh

  use iso_fortran_env, only: dp => real64
  use geqdsk_tools, only: geqdsk_t
  use mephit_util, only: interp1d, func1d_t

  implicit none

  private

  ! types and associated procedures
  public :: flux_func_cache, flux_func_cache_init, flux_func_cache_deinit
  public :: mesh_t, mesh_write, mesh_read, mesh_deinit, &
    generate_mesh, write_cache, read_cache, point_location, mesh_interp_theta_flux
  public :: coord_cache_t, field_cache_t, shielding_t, cache_t, &
    cache_write, cache_read, cache_init, cache_deinit
  public :: read_profiles, compute_auxiliary_profiles, resample_profiles, &
    write_profiles_hdf5, read_profiles_hdf5, deinit_profiles

  ! testing and debugging procedures
  public :: check_mesh, write_illustration_data, flux_func_cache_check, &
    check_safety_factor, check_curr0, equilibrium_field, curr0_geqdsk

  ! module variables
  public :: equil, psi_fine, fs, fs_half, mesh, cache
  public :: dens_e, temp_e, temp_i, E_r, Phi0, dPhi0_dpsi, nu_e, nu_i

  !> Electron density profile \f$ n \f$ in cm^-3.
  type(func1d_t) :: dens_e

  !> Electron temperature profile \f$ T_{e} \f$ in eV.
  type(func1d_t) :: temp_e

  !> Ion temperature profile \f$ T_{i} \f$ in eV.
  type(func1d_t) :: temp_i

  !> Radial electric field profile \f$ E_{r} \f$ in statV cm^-1.
  type(func1d_t) :: E_r

  !> Electric potential profile \f$ \Phi_{0} \f$ in statV.
  type(func1d_t) :: Phi0

  !> psi derivative of electric potential profile \f$ \Phi_{0} \f$ in statV Mx^-1.
  type(func1d_t) :: dPhi0_dpsi

  !> Electron collision frequency profile \f$ \nu_{e} \f$ in s^-1.
  type(func1d_t) :: nu_e

  !> Ion collision frequency profile \f$ \nu_{e} \f$ in s^-1.
  type(func1d_t) :: nu_i

  type(geqdsk_t) :: equil

  real(dp), dimension(:), allocatable :: psi_fine  ! intermediate step, to be removed

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

    !> Minor radius \f$ r \f$ in centimeter, measured outward from magnetic axis.
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

    !> Equivalent radius of poloidal cross-section area, in cm
    real(dp), dimension(:), allocatable, public :: rsmall

    !> Perimeter of poloidal cross-section, in cm
    real(dp), dimension(:), allocatable, public :: perimeter
  contains
    procedure :: init => flux_func_cache_init
    procedure :: deinit => flux_func_cache_deinit
  end type flux_func_cache

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

    !> Minimal R value on computational grid in cm.
    real(dp) :: R_min

    !> Maximal Z value on computational grid in cm.
    real(dp) :: R_max

    !> Minimal R value on computational grid in cm.
    real(dp) :: Z_min

    !> Maximal Z value on computational grid in cm.
    real(dp) :: Z_max

    !> Number of flux surfaces. May differ from #mephit_conf::mephit_conf::nflux due
    !> to refinement of flux surfaces.
    integer :: nflux

    !> Number of triangle nodes.
    integer :: npoint

    !> Number of triangles.
    integer :: ntri

    !> Number of triangle edges.
    integer :: nedge

    !> Toroidal mode number. May differ from #mephit_conf::config_t#n due to large
    !> aspect ratio scaling.
    integer :: n

    !> Minimal poloidal mode number in resonance.
    integer :: m_res_min

    !> Maximal poloidal mode number in resonance
    integer :: m_res_max

    !> Indices of flux surfaces where resonance corresponding to a poloidal mode (given as
    !> array index) occurs.
    integer, allocatable :: res_ind(:)

    !> Poloidal flux at resonance position corresponding to a poloidal mode (given as
    !> array index).
    real(dp), allocatable :: psi_res(:)

    !> Normalized small radius (outboard from O point, or towards X point) at resonance
    !> position corresponding to a poloidal mode (given as array index).
    real(dp), allocatable :: rad_norm_res(:)

    !> Small radius (of equivalent-area circle) at resonance
    !> position corresponding to a poloidal mode (given as array index).
    real(dp), allocatable :: rsmall_res(:)

    !> Poloidal modes that are expected to be in resonance. This might be different from
    !> m_res_min:m_res_max for specially constructed vacuum perturbation fields.
    integer, allocatable :: res_modes(:)

    !> Estimated resonant layer width in cm, computed from kinetic profiles.
    real(dp), allocatable :: delta_mn(:)

    !> Estimated resonant layer width in Mx (of psi), computed from kinetic profiles.
    real(dp), allocatable :: delta_psi_mn(:)

    !> Estimated resonant layer width in cm (of rad), computed from kinetic profiles.
    real(dp), allocatable :: delta_rad_mn(:)

    !> Damping factor for MDE solutions when "KiLCA" currents are used.
    real(dp), allocatable :: damping(:)

    !> Number of knots on the flux surface given by the array index.
    !>
    !> The array index ranges from 1 for the innermost flux surface to
    !> #mephit_conf::nflux for the last closed flux surface.
    integer, allocatable :: kp_max(:)

    !> Number of triangles inside the flux surface given by the array index.
    !>
    !> The array index ranges from 1 for the innermost flux surface to
    !> #mephit_conf::nflux for the last closed flux surface.
    integer, allocatable :: kt_max(:)

    !> Global index of the last knot of the previous flux surface given by the array index.
    !>
    !> The global index of knots in #mesh_mod::mesh_point on the flux surface kf runs from
    !> #kp_low (kf)+1 to #kp_low (kf)+#kp_max (kf), so #kp_low is determined by cumulatively
    !> adding consecutive values of #kp_max. The array index ranges from 1, giving the
    !> global index of the knot on the magnetic axis (which has to be 1), to
    !> #mephit_conf::nflux+1, effectively giving the last knot on the last closed
    !> flux surface.
    integer, allocatable :: kp_low(:)

    !> Global index of the last triangle of the previous flux surface given by the array
    !> index.
    !>
    !> The global index of triangles in #mesh_mod::mesh_element inside the flux surface kf
    !> runs from #kt_low (kf)+1 to #kt_low (kf)+#kt_max (kf), so #kt_low is determined by
    !> cumulatively adding consecutive values of #kt_max. The array index ranges from 1,
    !> giving the global index of the non-existent triangle on the magnetic axis (which is
    !> therefore 0), to #mephit_conf::nflux+1, giving the last triangle inside the last
    !> closed flux surface.
    integer, allocatable :: kt_low(:)

    real(dp), allocatable :: node_R(:)
    real(dp), allocatable :: node_Z(:)
    real(dp), allocatable :: node_theta_flux(:)
    real(dp), allocatable :: node_theta_geom(:)

    integer, allocatable :: tri_node(:, :)
    integer, allocatable :: tri_node_F(:)
    logical, allocatable :: tri_theta2pi(:)
    real(dp), allocatable :: tri_theta_extent(:, :)
    real(dp), allocatable :: tri_RZ_extent(:, :, :)

    !> true if edge f of given triangle lies on the outer flux surface, false otherwise
    logical, allocatable :: orient(:)
    integer, allocatable :: edge_node(:, :)
    integer, allocatable :: edge_tri(:, :)
    !> Edges [f, o, i] for a given triangle
    integer, allocatable :: tri_edge(:, :)

    real(dp), allocatable :: mid_R(:)
    real(dp), allocatable :: mid_Z(:)
    real(dp), allocatable :: edge_R(:)
    real(dp), allocatable :: edge_Z(:)
    real(dp), allocatable :: area(:)
    real(dp), allocatable :: cntr_R(:)
    real(dp), allocatable :: cntr_Z(:)

    !> Order of Gauss-Legendre quadrature
    integer :: GL_order = 2
    !> Weights of Gauss-Legendre quadrature
    real(dp), allocatable :: GL_weights(:)
    !> R coordinate of Gauss-Legendre quadrature point on edge
    real(dp), allocatable :: GL_R(:, :)
    !> Z coordinate of Gauss-Legendre quadrature point on edge
    real(dp), allocatable :: GL_Z(:, :)

    !> Number of 2D Gauss-Legendre quadrature points on triangle
    integer :: GL2_order = 6
    !> Weights of Gauss-Legendre quadrature
    real(dp), allocatable :: GL2_weights(:)
    !> R coordinate of Gauss-Legendre quadrature point on triangle
    real(dp), allocatable :: GL2_R(:, :)
    !> Z coordinate of Gauss-Legendre quadrature point on triangle
    real(dp), allocatable :: GL2_Z(:, :)

    !> Surface integral Jacobian used for normalization of poloidal modes in GPEC
    real(dp), allocatable :: gpec_jacfac(:)

    !> Flux surface average \f$ \langle R^{2} \lVert \nabla \psi \rVert^{2} \rangle \f$.
    real(dp), allocatable :: avg_R2gradpsi2(:)

  end type mesh_t

  interface
    subroutine gauss_legendre_unit_interval(order, points, weights) &
      bind(C, name = 'gauss_legendre_unit_interval')
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in), value :: order
      real(c_double), intent(out), dimension(1:order) :: points, weights
    end subroutine gauss_legendre_unit_interval

    subroutine FEM_triangulate_external(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname) &
      bind(C, name = 'FEM_triangulate_external')
      use iso_c_binding, only: c_char, c_int, c_double
      integer(c_int), intent(in), value :: npt_inner, npt_outer
      real(c_double), intent(in), dimension(1:npt_inner + npt_outer) :: node_R, node_Z
      real(c_double), intent(in), value :: R_O, Z_O
      character(c_char), intent(in) :: fname(*)
    end subroutine FEM_triangulate_external

    subroutine Rtree_init(ntri, tri_bb) bind(C, name = 'Rtree_init')
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in), value :: ntri
      real(c_double), intent(in), dimension(2, 2, ntri) :: tri_bb
    end subroutine rtree_init

    subroutine Rtree_query(R, Z, result_size, results) bind(C, name = "Rtree_query")
      use iso_c_binding, only: c_double, c_int, c_ptr
      real(c_double), intent(in), value :: R, Z
      integer(c_int), intent(out) :: result_size
      type(c_ptr), intent(out) :: results
    end subroutine rtree_query
  end interface

  type(mesh_t) :: mesh

  type :: coord_cache_t
    integer :: ktri
    real(dp) :: R, Z, psi, theta, sqrt_g, B0_R, B0_phi, B0_Z, dR_dtheta, dZ_dtheta
  end type coord_cache_t

  interface coord_cache_write
    module procedure coord_cache_write_1
    module procedure coord_cache_write_2
  end interface coord_cache_write

  interface coord_cache_read
    module procedure coord_cache_read_1
    module procedure coord_cache_read_2
  end interface coord_cache_read

  type :: field_cache_t
    real(dp) :: psi, theta, B0(3), j0(3), Bmod, dBmod_dR, dBmod_dZ, &
      dB0_dR(3), dB0_dZ(3), dj0_dR(3), dj0_dZ(3)
  end type field_cache_t

  interface field_cache_write
    module procedure field_cache_write_1
    module procedure field_cache_write_2
  end interface field_cache_write

  interface field_cache_read
    module procedure field_cache_read_1
    module procedure field_cache_read_2
  end interface field_cache_read

  type :: shielding_t
    real(dp), allocatable :: GL_weights(:)
    type(coord_cache_t), allocatable :: sample_Ires(:, :)

    !> Free parameter in the compensated scheme for shielding
    real(dp) :: coeff

    real(dp), allocatable :: cross_fade(:)
  end type shielding_t

  type :: cache_t
    integer :: GL_order
    type(shielding_t), allocatable :: shielding(:)
    type(coord_cache_t), allocatable :: sample_polmodes_half(:), sample_polmodes(:)
    type(coord_cache_t), allocatable :: sample_jnperp(:)
    type(field_cache_t), allocatable :: edge_fields(:, :), area_fields(:, :)
    type(field_cache_t), allocatable :: mid_fields(:), cntr_fields(:)
  end type cache_t

  type(cache_t) :: cache

contains

  !> Set up arrays of cached values of flux functions.
  !>
  !> @nflux number of flux surfaces
  !> @half_step values are taken at flux surfaces (false) or between flux surfaces (true)
  !>
  !> For full-grid quantities, values are taken on flux surfaces with indices running
  !> from 0 to \p nflux, i.e. from the magnetic axis to the separatrix. An exception is
  !> made for \psi, where the index runs up to \p nflux +1. This value is extrapolated for
  !> finite differences in mephit_iter::compute_presn() and mephit_iter::compute_bn_nonres().
  !> For half-grid quantities, values are taken between two flux surfaces with indices
  !> running from 1 to \p nflux, i.e. from the triangle strip surrounding the magnetic
  !> axis to the triangle strip just inside the separatrix.
  subroutine flux_func_cache_init(this, nflux, half_step)
    class(flux_func_cache), intent(inout) :: this
    integer, intent(in) :: nflux
    logical, intent(in) :: half_step

    call flux_func_cache_deinit(this)
    if (half_step) then
      allocate(this%psi(nflux))
      allocate(this%rad(nflux))
      allocate(this%F(nflux))
      allocate(this%p(nflux))
      allocate(this%FdF_dpsi(nflux))
      allocate(this%dp_dpsi(nflux))
      allocate(this%q(nflux))
      allocate(this%rsmall(nflux))
      allocate(this%perimeter(nflux))
    else
      allocate(this%psi(0:nflux))
      allocate(this%rad(0:nflux))
      allocate(this%F(0:nflux))
      allocate(this%p(0:nflux))
      allocate(this%FdF_dpsi(0:nflux))
      allocate(this%dp_dpsi(0:nflux))
      allocate(this%q(0:nflux))
      allocate(this%rsmall(0:nflux))
      allocate(this%perimeter(0:nflux))
    end if
    this%psi = 0d0
    this%rad = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
    this%rsmall = 0d0
    this%perimeter = 0d0
  end subroutine flux_func_cache_init

  subroutine flux_func_cache_deinit(this)
    class(flux_func_cache), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%F)) deallocate(this%F)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%FdF_dpsi)) deallocate(this%FdF_dpsi)
    if (allocated(this%dp_dpsi)) deallocate(this%dp_dpsi)
    if (allocated(this%q)) deallocate(this%q)
    if (allocated(this%rsmall)) deallocate(this%rsmall)
    if (allocated(this%perimeter)) deallocate(this%perimeter)
  end subroutine flux_func_cache_deinit

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
      comment = 'radial position outward from magnetic axis ' // trim(adjustl(comment)))
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rsmall', &
      cache%rsmall, lbound(cache%rsmall), ubound(cache%rsmall), unit = 'cm', &
      comment = 'equivalent radius of poloidal cross-section area ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/perimeter', &
      cache%perimeter, lbound(cache%perimeter), ubound(cache%perimeter), unit = 'cm', &
      comment = 'poloidal cross-section perimeter ' // trim(adjustl(comment)))
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rsmall', cache%rsmall)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/perimeter', cache%perimeter)
    call h5_close(h5id_root)
  end subroutine flux_func_cache_read

  subroutine coord_cache_write_1(cache, file, group, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(coord_cache_t), dimension(:), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: comment
    character(len = len_trim(group)) :: grp
    character(len = len_trim(comment)) :: cmnt
    integer(HID_T) :: h5id_root
    integer, dimension(:), allocatable :: itemp
    real(dp), dimension(:), allocatable :: rtemp

    grp = trim(group)
    cmnt = trim(comment)
    allocate(itemp(size(cache)), rtemp(size(cache)))
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    itemp(:) = cache%ktri
    call h5_add(h5id_root, grp // '/ktri', itemp, &
      lbound(cache), ubound(cache), &
      comment = 'triangle index of ' // cmnt)
    rtemp(:) = cache%R
    call h5_add(h5id_root, grp // '/R', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm', &
      comment = 'R coordinate of ' // cmnt)
    rtemp(:) = cache%Z
    call h5_add(h5id_root, grp // '/Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm', &
      comment = 'Z coordinate of ' // cmnt)
    rtemp(:) = cache%psi
    call h5_add(h5id_root, grp // '/psi', rtemp, &
      lbound(cache), ubound(cache), unit = 'Mx', &
      comment = 'poloidal flux at ' // cmnt)
    rtemp(:) = cache%theta
    call h5_add(h5id_root, grp // '/theta', rtemp, &
      lbound(cache), ubound(cache), unit = 'rad', &
      comment = 'flux poloidal angle at ' // cmnt)
    rtemp(:) = cache%sqrt_g
    call h5_add(h5id_root, grp // '/sqrt_g', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm G^-1', &
      comment = 'Jacobian at ' // cmnt)
    rtemp(:) = cache%B0_R
    call h5_add(h5id_root, grp // '/B0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%B0_phi
    call h5_add(h5id_root, grp // '/B0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%B0_Z
    call h5_add(h5id_root, grp // '/B0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dR_dtheta
    call h5_add(h5id_root, grp // '/dR_dtheta', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm rad^-1', &
      comment = 'Jacobian element (R, theta) at ' // cmnt)
    rtemp(:) = cache%dZ_dtheta
    call h5_add(h5id_root, grp // '/dZ_dtheta', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm rad^-1', &
      comment = 'Jacobian element (Z, theta) at ' // cmnt)
    call h5_close(h5id_root)
    deallocate(itemp, rtemp)
  end subroutine coord_cache_write_1

  subroutine coord_cache_write_2(cache, file, group, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(coord_cache_t), dimension(:, :), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: comment
    character(len = len_trim(group)) :: grp
    character(len = len_trim(comment)) :: cmnt
    integer(HID_T) :: h5id_root
    integer, dimension(:, :), allocatable :: itemp
    real(dp), dimension(:, :), allocatable :: rtemp

    grp = trim(group)
    cmnt = trim(comment)
    allocate(itemp(size(cache, 1), size(cache, 2)), rtemp(size(cache, 1), size(cache, 2)))
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    itemp(:, :) = cache%ktri
    call h5_add(h5id_root, grp // '/ktri', itemp, &
      lbound(cache), ubound(cache), &
      comment = 'triangle index of ' // cmnt)
    rtemp(:, :) = cache%R
    call h5_add(h5id_root, grp // '/R', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm', &
      comment = 'R coordinate of ' // cmnt)
    rtemp(:, :) = cache%Z
    call h5_add(h5id_root, grp // '/Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm', &
      comment = 'Z coordinate of ' // cmnt)
    rtemp(:, :) = cache%psi
    call h5_add(h5id_root, grp // '/psi', rtemp, &
      lbound(cache), ubound(cache), unit = 'Mx', &
      comment = 'poloidal flux at ' // cmnt)
    rtemp(:, :) = cache%theta
    call h5_add(h5id_root, grp // '/theta', rtemp, &
      lbound(cache), ubound(cache), unit = 'rad', &
      comment = 'flux poloidal angle at ' // cmnt)
    rtemp(:, :) = cache%sqrt_g
    call h5_add(h5id_root, grp // '/sqrt_g', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm G^-1', &
      comment = 'Jacobian at ' // cmnt)
    rtemp(:, :) = cache%B0_R
    call h5_add(h5id_root, grp // '/B0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%B0_phi
    call h5_add(h5id_root, grp // '/B0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%B0_Z
    call h5_add(h5id_root, grp // '/B0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dR_dtheta
    call h5_add(h5id_root, grp // '/dR_dtheta', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm rad^-1', &
      comment = 'Jacobian element (R, theta) at ' // cmnt)
    rtemp(:, :) = cache%dZ_dtheta
    call h5_add(h5id_root, grp // '/dZ_dtheta', rtemp, &
      lbound(cache), ubound(cache), unit = 'cm rad^-1', &
      comment = 'Jacobian element (Z, theta) at ' // cmnt)
    call h5_close(h5id_root)
    deallocate(itemp, rtemp)
  end subroutine coord_cache_write_2

  subroutine coord_cache_read_1(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(coord_cache_t), dimension(:), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer, dimension(:), allocatable :: itemp
    real(dp), dimension(:), allocatable :: rtemp

    grp = trim(group)
    allocate(itemp(size(cache)), rtemp(size(cache)))
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/ktri', itemp)
    cache%ktri = itemp
    call h5_get(h5id_root, grp // '/R', rtemp)
    cache%R = rtemp
    call h5_get(h5id_root, grp // '/Z', rtemp)
    cache%Z = rtemp
    call h5_get(h5id_root, grp // '/psi', rtemp)
    cache%psi = rtemp
    call h5_get(h5id_root, grp // '/theta', rtemp)
    cache%theta = rtemp
    call h5_get(h5id_root, grp // '/sqrt_g', rtemp)
    cache%sqrt_g = rtemp
    call h5_get(h5id_root, grp // '/B0_R', rtemp)
    cache%B0_R = rtemp
    call h5_get(h5id_root, grp // '/B0_phi', rtemp)
    cache%B0_phi = rtemp
    call h5_get(h5id_root, grp // '/B0_Z', rtemp)
    cache%B0_Z = rtemp
    call h5_get(h5id_root, grp // '/dR_dtheta', rtemp)
    cache%dR_dtheta = rtemp
    call h5_get(h5id_root, grp // '/dZ_dtheta', rtemp)
    cache%dZ_dtheta = rtemp
    call h5_close(h5id_root)
    deallocate(itemp, rtemp)
  end subroutine coord_cache_read_1

  subroutine coord_cache_read_2(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(coord_cache_t), dimension(:, :), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer, dimension(:, :), allocatable :: itemp
    real(dp), dimension(:, :), allocatable :: rtemp

    grp = trim(group)
    allocate(itemp(size(cache, 1), size(cache, 2)), rtemp(size(cache, 1), size(cache, 2)))
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/ktri', itemp)
    cache%ktri = itemp
    call h5_get(h5id_root, grp // '/R', rtemp)
    cache%R = rtemp
    call h5_get(h5id_root, grp // '/Z', rtemp)
    cache%Z = rtemp
    call h5_get(h5id_root, grp // '/psi', rtemp)
    cache%psi = rtemp
    call h5_get(h5id_root, grp // '/theta', rtemp)
    cache%theta = rtemp
    call h5_get(h5id_root, grp // '/sqrt_g', rtemp)
    cache%sqrt_g = rtemp
    call h5_get(h5id_root, grp // '/B0_R', rtemp)
    cache%B0_R = rtemp
    call h5_get(h5id_root, grp // '/B0_phi', rtemp)
    cache%B0_phi = rtemp
    call h5_get(h5id_root, grp // '/B0_Z', rtemp)
    cache%B0_Z = rtemp
    call h5_get(h5id_root, grp // '/dR_dtheta', rtemp)
    cache%dR_dtheta = rtemp
    call h5_get(h5id_root, grp // '/dZ_dtheta', rtemp)
    cache%dZ_dtheta = rtemp
    call h5_close(h5id_root)
    deallocate(itemp, rtemp)
  end subroutine coord_cache_read_2

  subroutine field_cache_write_1(cache, file, group, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(field_cache_t), dimension(:), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: comment
    character(len = len_trim(group)) :: grp
    character(len = len_trim(comment)) :: cmnt
    integer(HID_T) :: h5id_root
    real(dp), dimension(:), allocatable :: rtemp

    grp = trim(group)
    cmnt = trim(comment)
    allocate(rtemp(size(cache)))
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    rtemp(:) = cache%psi
    call h5_add(h5id_root, grp // '/psi', rtemp, &
      lbound(cache), ubound(cache), unit = 'Mx', &
      comment = 'poloidal flux at ' // cmnt)
    rtemp(:) = cache%theta
    call h5_add(h5id_root, grp // '/theta', rtemp, &
      lbound(cache), ubound(cache), unit = 'rad', &
      comment = 'flux poloidal angle at ' // cmnt)
    rtemp(:) = cache%B0(1)
    call h5_add(h5id_root, grp // '/B0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%B0(3)
    call h5_add(h5id_root, grp // '/B0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%B0(2)
    call h5_add(h5id_root, grp // '/B0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%j0(1)
    call h5_add(h5id_root, grp // '/j0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'R component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%j0(3)
    call h5_add(h5id_root, grp // '/j0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'Z component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%j0(2)
    call h5_add(h5id_root, grp // '/j0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%Bmod
    call h5_add(h5id_root, grp // '/B0', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dBmod_dR
    call h5_add(h5id_root, grp // '/dB0_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dBmod_dZ
    call h5_add(h5id_root, grp // '/dB0_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dR(1)
    call h5_add(h5id_root, grp // '/dB0R_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dZ(1)
    call h5_add(h5id_root, grp // '/dB0R_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dR(2)
    call h5_add(h5id_root, grp // '/dB0phi_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dZ(2)
    call h5_add(h5id_root, grp // '/dB0phi_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dR(3)
    call h5_add(h5id_root, grp // '/dB0Z_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dB0_dZ(3)
    call h5_add(h5id_root, grp // '/dB0Z_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:) = cache%dj0_dR(1)
    call h5_add(h5id_root, grp // '/dj0R_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of R component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%dj0_dZ(1)
    call h5_add(h5id_root, grp // '/dj0R_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of R component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%dj0_dR(2)
    call h5_add(h5id_root, grp // '/dj0phi_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%dj0_dZ(2)
    call h5_add(h5id_root, grp // '/dj0phi_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%dj0_dR(3)
    call h5_add(h5id_root, grp // '/dj0Z_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of Z component of equilibrium current density at ' // cmnt)
    rtemp(:) = cache%dj0_dZ(3)
    call h5_add(h5id_root, grp // '/dj0Z_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of Z component of equilibrium current density at ' // cmnt)
    call h5_close(h5id_root)
    deallocate(rtemp)
  end subroutine field_cache_write_1

  subroutine field_cache_write_2(cache, file, group, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(field_cache_t), dimension(:, :), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = *), intent(in) :: comment
    character(len = len_trim(group)) :: grp
    character(len = len_trim(comment)) :: cmnt
    integer(HID_T) :: h5id_root
    real(dp), dimension(:, :), allocatable :: rtemp

    grp = trim(group)
    cmnt = trim(comment)
    allocate(rtemp(size(cache, 1), size(cache, 2)))
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, grp // '/')
    rtemp(:, :) = cache%psi
    call h5_add(h5id_root, grp // '/psi', rtemp, &
      lbound(cache), ubound(cache), unit = 'Mx', &
      comment = 'poloidal flux at ' // cmnt)
    rtemp(:, :) = cache%theta
    call h5_add(h5id_root, grp // '/theta', rtemp, &
      lbound(cache), ubound(cache), unit = 'rad', &
      comment = 'flux poloidal angle at ' // cmnt)
    rtemp(:, :) = cache%B0(1)
    call h5_add(h5id_root, grp // '/B0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%B0(3)
    call h5_add(h5id_root, grp // '/B0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%B0(2)
    call h5_add(h5id_root, grp // '/B0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%j0(1)
    call h5_add(h5id_root, grp // '/j0_R', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'R component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%j0(3)
    call h5_add(h5id_root, grp // '/j0_Z', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'Z component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%j0(2)
    call h5_add(h5id_root, grp // '/j0_phi', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA', &
      comment = 'physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%Bmod
    call h5_add(h5id_root, grp // '/B0', rtemp, &
      lbound(cache), ubound(cache), unit = 'G', &
      comment = 'magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dBmod_dR
    call h5_add(h5id_root, grp // '/dB0_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dBmod_dZ
    call h5_add(h5id_root, grp // '/dB0_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of magnitude of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dR(1)
    call h5_add(h5id_root, grp // '/dB0R_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dZ(1)
    call h5_add(h5id_root, grp // '/dB0R_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of R component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dR(2)
    call h5_add(h5id_root, grp // '/dB0phi_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dZ(2)
    call h5_add(h5id_root, grp // '/dB0phi_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of physical phi component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dR(3)
    call h5_add(h5id_root, grp // '/dB0Z_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'R derivative of Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dB0_dZ(3)
    call h5_add(h5id_root, grp // '/dB0Z_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'G cm^-1', &
      comment = 'Z derivative of Z component of equilibrium magnetic field at ' // cmnt)
    rtemp(:, :) = cache%dj0_dR(1)
    call h5_add(h5id_root, grp // '/dj0R_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of R component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%dj0_dZ(1)
    call h5_add(h5id_root, grp // '/dj0R_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of R component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%dj0_dR(2)
    call h5_add(h5id_root, grp // '/dj0phi_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%dj0_dZ(2)
    call h5_add(h5id_root, grp // '/dj0phi_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of physical phi component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%dj0_dR(3)
    call h5_add(h5id_root, grp // '/dj0Z_dR', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'R derivative of Z component of equilibrium current density at ' // cmnt)
    rtemp(:, :) = cache%dj0_dZ(3)
    call h5_add(h5id_root, grp // '/dj0Z_dZ', rtemp, &
      lbound(cache), ubound(cache), unit = 'statA cm^-1', &
      comment = 'Z derivative of Z component of equilibrium current density at ' // cmnt)
    call h5_close(h5id_root)
    deallocate(rtemp)
  end subroutine field_cache_write_2

  subroutine field_cache_read_1(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(field_cache_t), dimension(:), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    real(dp), dimension(:), allocatable :: rtemp

    grp = trim(group)
    allocate(rtemp(size(cache)))
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/psi', rtemp)
    cache%psi = rtemp
    call h5_get(h5id_root, grp // '/theta', rtemp)
    cache%theta = rtemp
    call h5_get(h5id_root, grp // '/B0_R', rtemp)
    cache%B0(1) = rtemp
    call h5_get(h5id_root, grp // '/B0_Z', rtemp)
    cache%B0(3) = rtemp
    call h5_get(h5id_root, grp // '/B0_phi', rtemp)
    cache%B0(2) = rtemp
    call h5_get(h5id_root, grp // '/j0_R', rtemp)
    cache%j0(1) = rtemp
    call h5_get(h5id_root, grp // '/j0_Z', rtemp)
    cache%j0(3) = rtemp
    call h5_get(h5id_root, grp // '/j0_phi', rtemp)
    cache%j0(2) = rtemp
    call h5_get(h5id_root, grp // '/B0', rtemp)
    cache%Bmod = rtemp
    call h5_get(h5id_root, grp // '/dB0_dR', rtemp)
    cache%dBmod_dR = rtemp
    call h5_get(h5id_root, grp // '/dB0_dZ', rtemp)
    cache%dBmod_dZ = rtemp
    call h5_get(h5id_root, grp // '/dB0R_dR', rtemp)
    cache%dB0_dR(1) = rtemp
    call h5_get(h5id_root, grp // '/dB0R_dZ', rtemp)
    cache%dB0_dZ(1) = rtemp
    call h5_get(h5id_root, grp // '/dB0phi_dR', rtemp)
    cache%dB0_dR(2) = rtemp
    call h5_get(h5id_root, grp // '/dB0phi_dZ', rtemp)
    cache%dB0_dZ(2) = rtemp
    call h5_get(h5id_root, grp // '/dB0Z_dR', rtemp)
    cache%dB0_dR(3) = rtemp
    call h5_get(h5id_root, grp // '/dB0Z_dZ', rtemp)
    cache%dB0_dZ(3) = rtemp
    call h5_get(h5id_root, grp // '/dj0R_dR', rtemp)
    cache%dj0_dR(1) = rtemp
    call h5_get(h5id_root, grp // '/dj0R_dZ', rtemp)
    cache%dj0_dZ(1) = rtemp
    call h5_get(h5id_root, grp // '/dj0phi_dR', rtemp)
    cache%dj0_dR(2) = rtemp
    call h5_get(h5id_root, grp // '/dj0phi_dZ', rtemp)
    cache%dj0_dZ(2) = rtemp
    call h5_get(h5id_root, grp // '/dj0Z_dR', rtemp)
    cache%dj0_dR(3) = rtemp
    call h5_get(h5id_root, grp // '/dj0Z_dZ', rtemp)
    cache%dj0_dZ(3) = rtemp
    call h5_close(h5id_root)
    deallocate(rtemp)
  end subroutine field_cache_read_1

  subroutine field_cache_read_2(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(field_cache_t), dimension(:, :), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    real(dp), dimension(:, :), allocatable :: rtemp

    grp = trim(group)
    allocate(rtemp(size(cache, 1), size(cache, 2)))
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/psi', rtemp)
    cache%psi = rtemp
    call h5_get(h5id_root, grp // '/theta', rtemp)
    cache%theta = rtemp
    call h5_get(h5id_root, grp // '/B0_R', rtemp)
    cache%B0(1) = rtemp
    call h5_get(h5id_root, grp // '/B0_Z', rtemp)
    cache%B0(3) = rtemp
    call h5_get(h5id_root, grp // '/B0_phi', rtemp)
    cache%B0(2) = rtemp
    call h5_get(h5id_root, grp // '/j0_R', rtemp)
    cache%j0(1) = rtemp
    call h5_get(h5id_root, grp // '/j0_Z', rtemp)
    cache%j0(3) = rtemp
    call h5_get(h5id_root, grp // '/j0_phi', rtemp)
    cache%j0(2) = rtemp
    call h5_get(h5id_root, grp // '/B0', rtemp)
    cache%Bmod = rtemp
    call h5_get(h5id_root, grp // '/dB0_dR', rtemp)
    cache%dBmod_dR = rtemp
    call h5_get(h5id_root, grp // '/dB0_dZ', rtemp)
    cache%dBmod_dZ = rtemp
    call h5_get(h5id_root, grp // '/dB0R_dR', rtemp)
    cache%dB0_dR(1) = rtemp
    call h5_get(h5id_root, grp // '/dB0R_dZ', rtemp)
    cache%dB0_dZ(1) = rtemp
    call h5_get(h5id_root, grp // '/dB0phi_dR', rtemp)
    cache%dB0_dR(2) = rtemp
    call h5_get(h5id_root, grp // '/dB0phi_dZ', rtemp)
    cache%dB0_dZ(2) = rtemp
    call h5_get(h5id_root, grp // '/dB0Z_dR', rtemp)
    cache%dB0_dR(3) = rtemp
    call h5_get(h5id_root, grp // '/dB0Z_dZ', rtemp)
    cache%dB0_dZ(3) = rtemp
    call h5_get(h5id_root, grp // '/dj0R_dR', rtemp)
    cache%dj0_dR(1) = rtemp
    call h5_get(h5id_root, grp // '/dj0R_dZ', rtemp)
    cache%dj0_dZ(1) = rtemp
    call h5_get(h5id_root, grp // '/dj0phi_dR', rtemp)
    cache%dj0_dR(2) = rtemp
    call h5_get(h5id_root, grp // '/dj0phi_dZ', rtemp)
    cache%dj0_dZ(2) = rtemp
    call h5_get(h5id_root, grp // '/dj0Z_dR', rtemp)
    cache%dj0_dR(3) = rtemp
    call h5_get(h5id_root, grp // '/dj0Z_dZ', rtemp)
    cache%dj0_dZ(3) = rtemp
    call h5_close(h5id_root)
    deallocate(rtemp)
  end subroutine field_cache_read_2

  elemental subroutine shielding_init(s)
    type(shielding_t), intent(inout) :: s

    call shielding_deinit(s)
    allocate(s%cross_fade(0:mesh%nflux))
  end subroutine shielding_init

  elemental subroutine shielding_deinit(s)
    type(shielding_t), intent(inout) :: s

    s%coeff = 0d0
    if (allocated(s%GL_weights)) deallocate(s%GL_weights)
    if (allocated(s%sample_Ires)) deallocate(s%sample_Ires)
    if (allocated(s%cross_fade)) deallocate(s%cross_fade)
  end subroutine shielding_deinit

  subroutine shielding_write(shielding, file, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    type(shielding_t), intent(in) :: shielding
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root

    grp = trim(group)
    call coord_cache_write(shielding%sample_Ires, file, grp // &
      '/sample_Ires', 'resonant current sampling points')
    call h5_open_rw(file, h5id_root)
    call h5_add(h5id_root, grp // '/GL_weights', shielding%GL_weights, &
      lbound(shielding%GL_weights), ubound(shielding%GL_weights), &
      unit = 'Mx', comment = 'G-L quadrature weights including psi interval')
    call h5_add(h5id_root, grp // '/coeff', shielding%coeff, unit = 'g^-1 cm s', &
      comment = 'coefficient in the compensated scheme for shielding')
    call h5_add(h5id_root, grp // '/cross_fade', shielding%cross_fade, &
      lbound(shielding%cross_fade), ubound(shielding%cross_fade), &
      unit = '1', comment = 'Transition function between current components')
    call h5_close(h5id_root)
  end subroutine shielding_write

  subroutine shielding_read(shielding, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_get_bounds, h5_close
    type(shielding_t), intent(inout) :: shielding
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    integer(HID_T) :: h5id_root
    integer :: lb1, lb2, ub1, ub2

    grp = trim(group)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/coeff', shielding%coeff)
    call h5_get(h5id_root, grp // '/cross_fade', shielding%cross_fade)
    call h5_get_bounds(h5id_root, grp // '/sample_Ires/R', lb1, lb2, ub1, ub2)
    if (allocated(shielding%GL_weights)) deallocate(shielding%GL_weights)
    allocate(shielding%GL_weights(lb2:ub2))
    call h5_get(h5id_root, grp // '/GL_weights', shielding%GL_weights)
    call h5_close(h5id_root)
    if (allocated(shielding%sample_Ires)) deallocate(shielding%sample_Ires)
    allocate(shielding%sample_Ires(lb1:ub1, lb2:ub2))
    call coord_cache_read(shielding%sample_Ires, file, grp // '/sample_Ires')
  end subroutine shielding_read

  subroutine cache_init(cache, GL_order)
    type(cache_t), intent(inout) :: cache
    integer, intent(in) :: GL_order

    call cache_deinit(cache)
    cache%GL_order = GL_order
    allocate(cache%sample_polmodes(mesh%npoint))
    allocate(cache%sample_polmodes_half(mesh%ntri))
    allocate(cache%sample_jnperp(mesh%npoint - 1))
    allocate(cache%shielding(mesh%m_res_min:mesh%m_res_max))
    call shielding_init(cache%shielding)
    allocate(cache%edge_fields(mesh%GL_order, mesh%nedge), cache%area_fields(mesh%GL2_order, mesh%ntri))
    allocate(cache%mid_fields(mesh%nedge), cache%cntr_fields(mesh%ntri))
  end subroutine cache_init

  subroutine cache_deinit(cache)
    type(cache_t), intent(inout) :: cache

    cache%GL_order = 0
    if (allocated(cache%sample_polmodes)) deallocate(cache%sample_polmodes)
    if (allocated(cache%sample_polmodes_half)) deallocate(cache%sample_polmodes_half)
    if (allocated(cache%sample_jnperp)) deallocate(cache%sample_jnperp)
    if (allocated(cache%shielding)) then
      call shielding_deinit(cache%shielding)
      deallocate(cache%shielding)
    end if
    if (allocated(cache%edge_fields)) deallocate(cache%edge_fields)
    if (allocated(cache%area_fields)) deallocate(cache%area_fields)
    if (allocated(cache%mid_fields)) deallocate(cache%mid_fields)
    if (allocated(cache%cntr_fields)) deallocate(cache%cntr_fields)
  end subroutine cache_deinit

  subroutine cache_write(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    type(cache_t), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    character(len = *), parameter :: fmt = '("_", i2.2)'
    character(len = 3) :: suffix
    integer(HID_T) :: h5id_root
    integer :: m

    grp = trim(group)
    call coord_cache_write(cache%sample_polmodes_half, file, grp // &
      '/sample_polmodes_half', 'poloidal mode sampling points between flux surfaces')
    call coord_cache_write(cache%sample_polmodes, file, grp // &
      '/sample_polmodes', 'poloidal mode sampling points on flux surfaces')
    call coord_cache_write(cache%sample_jnperp, file, grp // &
      '/sample_jnperp', 'sampling points for FDM computation of perpendicular current')
    call field_cache_write(cache%edge_fields, file, grp // &
      '/edge_fields', 'GL quadrature points on triangle edges')
    call field_cache_write(cache%area_fields, file, grp // &
      '/area_fields', 'GL quadrature points on triangle areas')
    call field_cache_write(cache%mid_fields, file, grp // &
      '/mid_fields', 'triangle edge midpoints')
    call field_cache_write(cache%cntr_fields, file, grp // &
      '/cntr_fields', 'weighted triangle centroids')
    call h5_open_rw(file, h5id_root)
    call h5_add(h5id_root, grp // '/GL_order', cache%GL_order, &
      comment = 'order of G-L quadrature for resonant current sampling points')
    call h5_close(h5id_root)
    do m = mesh%m_res_min, mesh%m_res_max
      write (suffix, fmt) m
      call shielding_write(cache%shielding(m), file, grp // &
        '/shielding' // suffix)
    end do
  end subroutine cache_write

  subroutine cache_read(cache, file, group)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(cache_t), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp
    character(len = *), parameter :: fmt = '("_", i2.2)'
    character(len = 3) :: suffix
    integer(HID_T) :: h5id_root
    integer :: m, GL_order

    grp = trim(group)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, grp // '/GL_order', GL_order)
    call h5_close(h5id_root)
    call cache_init(cache, GL_order)
    call coord_cache_read(cache%sample_polmodes_half, file, &
      grp // '/sample_polmodes_half')
    call coord_cache_read(cache%sample_polmodes, file, &
      grp // '/sample_polmodes')
    call coord_cache_read(cache%sample_jnperp, file, &
      grp // '/sample_jnperp')
    do m = mesh%m_res_min, mesh%m_res_max
      write (suffix, fmt) m
      call shielding_read(cache%shielding(m), file, grp // '/shielding' // suffix)
    end do
    call field_cache_read(cache%edge_fields, file, grp // '/edge_fields')
    call field_cache_read(cache%area_fields, file, grp // '/area_fields')
    call field_cache_read(cache%mid_fields, file, grp // '/mid_fields')
    call field_cache_read(cache%cntr_fields, file, grp // '/cntr_fields')
  end subroutine cache_read

  subroutine generate_mesh
    use mephit_conf, only: conf
    integer :: m

    if (conf%kilca_scale_factor /= 0) then
      mesh%n = conf%n * conf%kilca_scale_factor
    else
      mesh%n = conf%n
    end if
    call create_mesh_points
    call init_flux_variables
    call compare_gpec_coordinates
    call connect_mesh_points
#ifdef USE_MFEM
    call mesh_write_MFEM
#endif
    call write_FreeFem_mesh
    call cache_init(cache, 4)
    call cache_equilibrium_field
    call compute_sample_polmodes(cache%sample_polmodes_half, .true.)
    call compute_sample_polmodes(cache%sample_polmodes, .false.)
    call compute_sample_jnperp(cache%sample_jnperp)
    do m = mesh%m_res_min, mesh%m_res_max
      call compute_shielding_auxiliaries(cache%shielding(m), m)
      call compute_sample_Ires(cache%shielding(m)%sample_Ires, &
        cache%shielding(m)%GL_weights, cache%GL_order, m)
    end do
    call compute_kilca_auxiliaries
    call compute_gpec_jacfac
    call check_resonance_positions
    call compute_curr0
  end subroutine generate_mesh

  subroutine write_cache
    use mephit_conf, only: conf_arr, datafile

    call conf_arr%export_hdf5(datafile, 'config')
    call flux_func_cache_write(fs, datafile, 'cache/fs', 'on flux surfaces')
    call flux_func_cache_write(fs_half, datafile, 'cache/fs_half', 'between flux surfaces')
    call cache_write(cache, datafile, 'cache')
  end subroutine write_cache

  subroutine read_cache
    use mephit_conf, only: datafile

    call flux_func_cache_init(fs, mesh%nflux, .false.)
    call flux_func_cache_init(fs_half, mesh%nflux, .true.)
    call flux_func_cache_read(fs, datafile, 'cache/fs')
    call flux_func_cache_read(fs_half, datafile, 'cache/fs_half')
    call cache_read(cache, datafile, 'cache')
  end subroutine read_cache

  subroutine read_profiles
    use mephit_conf, only: conf
    use mephit_util, only: func1d_read_formatted

    call func1d_read_formatted(dens_e, conf%dens_file)
    if (abs(dens_e%x(ubound(dens_e%x, 1)) - 1d0) <= 0.05d0) then
      dens_e%y(:) = dens_e%y * 1d-6  ! SI units
    end if
    call func1d_read_formatted(temp_e, conf%temp_e_file)
    call func1d_read_formatted(temp_i, conf%temp_i_file)
    call func1d_read_formatted(E_r, conf%E_r_file)
  end subroutine read_profiles

  subroutine compute_auxiliary_profiles
    use magdata_in_symfluxcoor_mod, only: rsmall, psisurf, psipol_max
    use mephit_conf, only: conf, logger
    use mephit_util, only: func1d_init, resample1d
    integer :: nrad, krad
    real(dp), allocatable :: rsmall_interp(:), psi_interp(:)

    if (size(dens_e%x) /= size(temp_e%x)) then
      call logger%msg_arg_size('calculate_auxiliary_profiles', &
        'size(dens_e)', 'size(temp_e)', size(dens_e%x), size(temp_e%x))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (size(dens_e%x) /= size(temp_i%x)) then
      call logger%msg_arg_size('calculate_auxiliary_profiles', &
        'size(dens_e)', 'size(temp_i)', size(dens_e%x), size(temp_i%x))
      if (logger%err) call logger%write_msg
      error stop
    end if
    nrad = size(dens_e%x)
    ! transverse diffusion rate of fast electrons in ion background
    ! (NRL Plasma Formulary 2016, p. 32)
    call func1d_init(nu_e, 1, nrad)
    nu_e%x(:) = dens_e%x
    nu_e%y(:) = 7.7d-6 * (1d0 + conf%Z_i) * dens_e%y / temp_e%y ** 1.5d0 &
      * (24.0d0 - log(sqrt(dens_e%y) / temp_e%y))
    ! transverse diffusion rate of fast ions in ion background
    ! (NRL Plasma Formulary 2016, p. 32)
    call func1d_init(nu_i, 1, nrad)
    nu_i%x(:) = dens_e%x
    nu_i%y(:) = 1.8d-7 * conf%Z_i ** 3 / sqrt(conf%m_i) * dens_e%y / temp_i%y ** 1.5d0 &
      * (23.0d0 - log(conf%Z_i ** 2.5d0 * sqrt(2.0d0 * dens_e%y) / temp_i%y ** 1.5d0))
    ! add offset to electrical field
    E_r%y = E_r%y + conf%offset_E_r
    ! integrate electric field over rsmall to yield the electric potential
    nrad = size(E_r%x)
    allocate(rsmall_interp(nrad), psi_interp(nrad))
    if (abs(E_r%x(nrad) - 1d0) <= 0.05d0) then
      psi_interp(:) = E_r%x ** 2
      call resample1d(psisurf(1:), rsmall, psi_interp, rsmall_interp, 3)
    else
      rsmall_interp(:) = E_r%x
      call resample1d(rsmall, psisurf(1:), rsmall_interp, psi_interp, 3)
    end if
    psi_interp(:) = psi_interp * psipol_max  ! for correct scaling of derivative
    call func1d_init(Phi0, 1, nrad)
    Phi0%y(:) = 0d0
    do krad = 2, nrad
      Phi0%y(krad) = Phi0%y(krad - 1) - (rsmall_interp(krad) - rsmall_interp(krad - 1)) &
        * 0.5d0 * (E_r%y(krad) + E_r%y(krad - 1))
    end do
    Phi0%x(:) = E_r%x
    ! psi derivative of Phi0
    call func1d_init(dPhi0_dpsi, 1, nrad)
    call resample1d(psi_interp, Phi0%y, psi_interp, dPhi0_dpsi%y, 3, .true.)
    dPhi0_dpsi%x(:) = E_r%x
    deallocate(rsmall_interp, psi_interp)
  end subroutine compute_auxiliary_profiles

  subroutine resample_profiles
    use mephit_util, only: func1d_init, resample1d
    real(dp), dimension(0:mesh%nflux) :: rho_pol, resampled

    rho_pol = sqrt((fs%psi - fs%psi(0)) / (fs%psi(mesh%nflux) - fs%psi(0)))
    call resample_profile(dens_e)
    call resample_profile(temp_e)
    call resample_profile(temp_i)
    call resample_profile(E_r)
    call resample_profile(Phi0)
    call resample_profile(dPhi0_dpsi)
    call resample_profile(nu_e)
    call resample_profile(nu_i)

  contains
    subroutine resample_profile(profile)
      type(func1d_t), intent(inout) :: profile

      if (abs(profile%x(ubound(profile%x, 1)) - 1d0) <= 0.05d0) then
        call resample1d(profile%x, profile%y, rho_pol, resampled, 3)
      else
        call resample1d(profile%x, profile%y, fs%rsmall, resampled, 3)
      end if
      call func1d_init(profile, 0, mesh%nflux)
      profile%y(:) = resampled
      profile%x(:) = fs%psi
    end subroutine resample_profile
  end subroutine resample_profiles

  subroutine write_profiles_hdf5(file, group)
    use mephit_util, only: func1d_write
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_write(dens_e, file, grp // '/dens_e', &
      'poloidal flux', 'Mx', 'electron density', 'cm^-3')
    call func1d_write(temp_e, file, grp // '/temp_e', &
      'poloidal flux', 'Mx', 'electron temperature', 'eV')
    call func1d_write(temp_i, file, grp // '/temp_i', &
      'poloidal flux', 'Mx', 'ion temperature', 'eV')
    call func1d_write(E_r, file, grp // '/E_r', &
      'poloidal flux', 'Mx', 'radial electric field', 'statV cm^-1')
    call func1d_write(Phi0, file, grp // '/Phi0', &
      'poloidal flux', 'Mx', 'electric potential', 'statV')
    call func1d_write(dPhi0_dpsi, file, grp // '/dPhi0_dpsi', &
      'poloidal flux', 'Mx', 'psi derivative of electric potential', 'statV Mx^-1')
    call func1d_write(nu_e, file, grp // '/nu_e', &
      'poloidal flux', 'Mx', 'electron collision frequency', 's^-1')
    call func1d_write(nu_i, file, grp // '/nu_i', &
      'poloidal flux', 'Mx', 'ion collision frequency', 's^-1')
  end subroutine write_profiles_hdf5

  subroutine read_profiles_hdf5(file, group)
    use mephit_util, only: func1d_read
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: group
    character(len = len_trim(group)) :: grp

    grp = trim(group)
    call func1d_read(dens_e, file, grp // '/dens_e')
    call func1d_read(temp_e, file, grp // '/temp_e')
    call func1d_read(temp_i, file, grp // '/temp_i')
    call func1d_read(E_r, file, grp // '/E_r')
    call func1d_read(Phi0, file, grp // '/Phi0')
    call func1d_read(dPhi0_dpsi, file, grp // '/dPhi0_dpsi')
    call func1d_read(nu_e, file, grp // '/nu_e')
    call func1d_read(nu_i, file, grp // '/nu_i')
  end subroutine read_profiles_hdf5

  subroutine deinit_profiles
    use mephit_util, only: func1d_deinit

    call func1d_deinit(dens_e)
    call func1d_deinit(temp_e)
    call func1d_deinit(temp_i)
    call func1d_deinit(E_r)
    call func1d_deinit(Phi0)
    call func1d_deinit(dPhi0_dpsi)
    call func1d_deinit(nu_e)
    call func1d_deinit(nu_i)
  end subroutine deinit_profiles

  subroutine compute_resonance_positions(psi_sample, q_sample, rad_norm_sample)
    use mephit_conf, only: conf, logger
    use mephit_util, only: interp1d
    use netlib_mod, only: zeroin
    real(dp), dimension(:), intent(in) :: psi_sample
    real(dp), dimension(:), intent(in) :: q_sample
    real(dp), dimension(:), intent(in) :: rad_norm_sample
    real(dp) :: psi_min, psi_max
    integer :: m

    if (size(psi_sample) /= size(q_sample)) then
      call logger%msg_arg_size('compute_resonance_positions', 'size(psi_sample)', &
        'size(q_sample)', size(psi_sample), size(q_sample))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (size(psi_sample) /= size(rad_norm_sample)) then
      call logger%msg_arg_size('compute_resonance_positions', 'size(psi_sample)', &
        'size(rad_norm_sample)', size(psi_sample), size(rad_norm_sample))
      if (logger%err) call logger%write_msg
      error stop
    end if
    psi_min = minval(psi_sample)
    psi_max = maxval(psi_sample)
    mesh%m_res_min = ceiling(minval(abs(q_sample)) * dble(mesh%n))
    if (conf%ignore_q1_res) then
      mesh%m_res_min = max(mesh%m_res_min, conf%n + 1)
    else
      mesh%m_res_min = max(mesh%m_res_min, 1)
    end if
    mesh%m_res_max = floor(maxval(abs(q_sample)) * dble(mesh%n))
    if (allocated(mesh%res_modes)) deallocate(mesh%res_modes)
    if (conf%kilca_scale_factor /= 0 .and. conf%kilca_pol_mode /= 0) then
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
    logger%msg = 'resonance positions:'
    if (logger%debug) call logger%write_msg
    do m = mesh%m_res_min, mesh%m_res_max
      mesh%psi_res(m) = zeroin(psi_min, psi_max, q_interp_resonant, 1d-9)
      mesh%rad_norm_res(m) = interp1d(psi_sample, rad_norm_sample, mesh%psi_res(m), 3)
      write (logger%msg, '("m = ", i2, ", psi = ", es24.16e3, ", rad = ", f19.16)') &
        m, mesh%psi_res(m), mesh%rad_norm_res(m)
      if (logger%debug) call logger%write_msg
    end do

  contains
    function q_interp_resonant(psi)
      use mephit_util, only: interp1d
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant

      q_interp_resonant = interp1d(psi_sample, abs(q_sample), psi, 3) - dble(m) / dble(mesh%n)
    end function q_interp_resonant
  end subroutine compute_resonance_positions

  subroutine compute_resonant_layer_widths
    use magdata_in_symfluxcoor_mod, only: qsaf, rsmall, psisurf, psipol_max, rbeg
    use mephit_conf, only: logger
    use mephit_util, only: clight, resample1d
    integer :: m
    real(dp) :: delta(2)
    real(dp), dimension(mesh%m_res_min:mesh%m_res_max) :: rho_pol_res, &
      q, V_ExB, k_perp, v_Te, nu_e_interp

    if (allocated(mesh%rsmall_res)) deallocate(mesh%rsmall_res)
    allocate(mesh%rsmall_res(mesh%m_res_min:mesh%m_res_max))
    call resample1d(psi_fine, rsmall, mesh%psi_res, mesh%rsmall_res, 3)
    call resample1d(psi_fine, sqrt(psisurf(1:)), mesh%psi_res, rho_pol_res, 3)
    call resample1d(psi_fine, abs(qsaf), mesh%psi_res, q, 3)
    call resample1d(psi_fine, rsmall, mesh%psi_res, k_perp, 3)
    k_perp = [(dble(m), m = mesh%m_res_min, mesh%m_res_max)] / k_perp  ! k_perp => rsmall
    if (abs(E_r%x(ubound(E_r%x, 1)) - 1d0) <= 0.05d0) then
      call resample1d(E_r%x, E_r%y, rho_pol_res, V_ExB, 3)
    else
      call resample1d(E_r%x, E_r%y, mesh%rsmall_res, V_ExB, 3)
    end if
    V_ExB = clight * abs(V_ExB / equil%bcentr)  ! V_ExB => E_r
    if (abs(temp_e%x(ubound(temp_e%x, 1)) - 1d0) <= 0.05d0) then
      call resample1d(temp_e%x, temp_e%y, rho_pol_res, v_Te, 3)
    else
      call resample1d(temp_e%x, temp_e%y, mesh%rsmall_res, v_Te, 3)
    end if
    v_Te = 4.19e7 * sqrt(abs(v_Te))  ! v_Te => temp_e
    if (abs(nu_e%x(ubound(nu_e%x, 1)) - 1d0) <= 0.05d0) then
      call resample1d(nu_e%x, nu_e%y, rho_pol_res, nu_e_interp, 3)
    else
      call resample1d(nu_e%x, nu_e%y, mesh%rsmall_res, nu_e_interp, 3)
    end if
    if (allocated(mesh%delta_mn)) deallocate(mesh%delta_mn)
    allocate(mesh%delta_mn(mesh%m_res_min:mesh%m_res_max))
    if (allocated(mesh%delta_psi_mn)) deallocate(mesh%delta_psi_mn)
    allocate(mesh%delta_psi_mn(mesh%m_res_min:mesh%m_res_max))
    if (allocated(mesh%delta_rad_mn)) deallocate(mesh%delta_rad_mn)
    allocate(mesh%delta_rad_mn(mesh%m_res_min:mesh%m_res_max))
    logger%msg = 'resonant layer widths:'
    if (logger%debug) call logger%write_msg
    do m = mesh%m_res_min, mesh%m_res_max
      mesh%delta_mn(m) = q(m) * mesh%R_O * V_ExB(m) / v_Te(m) * &
        max(1d0, sqrt(nu_e_interp(m) / (k_perp(m) * V_ExB(m))))
      write (logger%msg, '("m = ", i2, ", rsmall_mn = ", es24.16e3, ", delta_mn = ", es24.16e3)') &
        m, mesh%rsmall_res(m), mesh%delta_mn(m)
      if (logger%debug) call logger%write_msg
      call resample1d(rsmall, psisurf(1:) * psipol_max, &
        mesh%rsmall_res(m) + [-0.5d0, 0.5d0] * mesh%delta_mn(m), delta, 3)
      mesh%delta_psi_mn(m) = abs(delta(2) - delta(1))
      call resample1d(rsmall, rbeg, &
        mesh%rsmall_res(m) + [-0.5d0, 0.5d0] * mesh%delta_mn(m), delta, 3)
      mesh%delta_rad_mn(m) = delta(2) - delta(1)
    end do
  end subroutine compute_resonant_layer_widths

  subroutine refine_unit_partition_gaussian(nref, coarse_sep, refinement, resonances, widths, nflux, partition)
    use mephit_conf, only: logger
    use mephit_util, only: linspace
    integer, intent(in) :: nref
    real(dp), intent(in) :: coarse_sep
    real(dp), intent(in) :: refinement(:)
    real(dp), intent(in) :: resonances(:)
    real(dp), intent(in) :: widths(:)
    integer, intent(out) :: nflux
    real(dp), allocatable, intent(out) :: partition(:)
    real(dp) :: pos_curr, pos_next, pos_max
    integer :: k

    if (allocated(partition)) deallocate(partition)
    if (nref < 1) then
      nflux = ceiling(1d0 / coarse_sep)
      allocate(partition(0:nflux))
      partition(:) = linspace(0d0, 1d0, nflux + 1, 0, 0)
      return
    end if
    if (nref /= size(refinement)) then
      call logger%msg_arg_size('refine_unit_partition_KilCA', 'nref', 'size(refinement)', nref, &
        size(refinement))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(resonances)) then
      call logger%msg_arg_size('refine_unit_partition_KiLCA', 'nref', 'size(resonances)', nref, &
        size(resonances))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(resonances)) then
      call logger%msg_arg_size('refine_unit_partition_KiLCA', 'nref', 'size(resonances)', nref, &
        size(resonances))
      if (logger%err) call logger%write_msg
      error stop
    end if
    nflux = 0
    pos_curr = 0d0
    do while (pos_curr < 1d0)
      pos_next = pos_curr + coarse_sep / recnsplit(pos_curr)
      pos_curr = 0.5d0 * (pos_curr + pos_next + coarse_sep / recnsplit(pos_next))
      nflux = nflux + 1
    end do
    pos_max = pos_curr
    allocate(partition(0:nflux))
    partition(:) = 0d0
    pos_curr = 0d0
    do k = 1, nflux
      pos_next = pos_curr + coarse_sep / recnsplit(pos_curr)
      pos_curr = 0.5d0 * (pos_curr + pos_next + coarse_sep / recnsplit(pos_next))
      partition(k) = pos_curr / pos_max
    end do

  contains
    pure function recnsplit(pos)
      real(dp), intent(in) :: pos
      real(dp) :: recnsplit

      recnsplit = 1d0 + sum(refinement * exp(-(pos - resonances) ** 2 / widths ** 2))
    end function recnsplit
  end subroutine refine_unit_partition_gaussian

  subroutine refine_unit_partition(nref, coarse_sep, fine_sep, add_fine, refinement, resonances, diverging_q, &
    nflux, partition, kf_ref)
    use mephit_conf, only: logger
    use mephit_util, only: linspace
    integer, intent(in) :: nref
    real(dp), intent(in) :: coarse_sep
    real(dp), intent(in) :: fine_sep(:)
    integer, intent(in) :: add_fine(:)
    real(dp), intent(in) :: refinement(:)
    real(dp), intent(in) :: resonances(:)
    logical, intent(in) :: diverging_q
    integer, intent(out) :: nflux
    real(dp), allocatable, intent(out) :: partition(:)
    integer, intent(out) :: kf_ref(:)
    integer :: kref, k, max_extent, fine_lo(nref), fine_hi(nref), grow_lo(nref), grow_hi(nref), &
      inter(nref + 1), kf_lo, kf_hi
    real(dp) :: sep_lo, sep_hi, sep_inter, dist_inter
    real(dp), allocatable :: pos_lo(:, :), pos_hi(:, :)

    if (allocated(partition)) deallocate(partition)
    if (nref < 1) then
      nflux = ceiling(1d0 / coarse_sep)
      allocate(partition(0:nflux))
      partition(:) = linspace(0d0, 1d0, nflux + 1, 0, 0)
      return
    end if
    if (nref /= size(fine_sep)) then
      call logger%msg_arg_size('refine_unit_partition', 'nref', 'size(fine_sep)', nref, &
        size(fine_sep))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(add_fine)) then
      call logger%msg_arg_size('refine_unit_partition', 'nref', 'size(add_fine)', nref, &
        size(add_fine))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(refinement)) then
      call logger%msg_arg_size('refine_unit_partition', 'nref', 'size(refinement)', nref, &
        size(refinement))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(resonances)) then
      call logger%msg_arg_size('refine_unit_partition', 'nref', 'size(resonances)', nref, &
        size(resonances))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (nref /= size(kf_ref)) then
      call logger%msg_arg_size('refine_unit_partition', 'nref', 'size(kf_ref)', nref, &
        size(kf_ref))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (any(refinement <= 1d0)) then
      logger%msg = 'refine_unit_partition: refinement factors need to be larger than one'
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (any(fine_sep >= coarse_sep)) then
      logger%msg = 'refine_unit_partition: all elements of fine_sep need to be less than coarse_sep'
      if (logger%err) call logger%write_msg
      error stop
    end if
    ! construct refined regions via geometric series
    grow_lo = floor(log(coarse_sep / fine_sep) / log(refinement))
    grow_hi = grow_lo
    if (diverging_q) then
      ! continue with fine separation outside last refined region
      grow_hi(nref) = 0
    end if
    fine_lo = max(0, add_fine)
    fine_hi = fine_lo
    max_extent = maxval(grow_lo + fine_lo)
    allocate(pos_lo(0:max_extent, 1:nref), pos_hi(0:max_extent, 1:nref))
    do kref = 1, nref
      pos_lo(:fine_lo(kref), kref) = resonances(kref) - fine_sep(kref) * [(0.5d0 + dble(k), k = 0, fine_lo(kref))]
      pos_lo(fine_lo(kref):, kref) = pos_lo(fine_lo(kref), kref) - fine_sep(kref) * (refinement(kref) ** &
        [(k, k = 1, max_extent - fine_lo(kref) + 1)] - refinement(kref) ** 1) / (refinement(kref) - 1d0)
      pos_hi(:fine_hi(kref), kref) = resonances(kref) + fine_sep(kref) * [(0.5d0 + dble(k), k = 0, fine_hi(kref))]
      pos_hi(fine_hi(kref):, kref) = pos_hi(fine_hi(kref), kref) + fine_sep(kref) * (refinement(kref) ** &
        [(k, k = 1, max_extent - fine_hi(kref) + 1)] - refinement(kref) ** 1) / (refinement(kref) - 1d0)
    end do
    ! reduce number of points in geometric series if refined regions overlap
    ! and if transition between refined and unrefined regions is not concave
    inter = 0
    do kref = 2, nref
      inner: do
        sep_lo = fine_sep(kref) * refinement(kref) ** grow_lo(kref)
        sep_hi = fine_sep(kref - 1) * refinement(kref - 1) ** grow_hi(kref - 1)
        dist_inter = pos_lo(fine_lo(kref) + grow_lo(kref), kref) - &
          pos_hi(fine_hi(kref - 1) + grow_hi(kref - 1), kref - 1)
        if (dist_inter > 0d0) then  ! refined regions do not overlap
          inter(kref) = ceiling(dist_inter / min(coarse_sep, &
            sep_lo * refinement(kref), sep_hi * refinement(kref - 1)))
          sep_inter = dist_inter / dble(inter(kref))
          if (sep_lo <= sep_inter .and. sep_hi <= sep_inter) then
            exit inner
          end if
        end if
        if (grow_lo(kref) > 0 .and. grow_hi(kref - 1) > 0) then
          if (sep_lo <= sep_hi) then
            grow_hi(kref - 1) = grow_hi(kref - 1) - 1
          else
            grow_lo(kref) = grow_lo(kref) - 1
          end if
        elseif (grow_hi(kref - 1) > 0) then
          grow_hi(kref - 1) = grow_hi(kref - 1) - 1
        elseif (grow_lo(kref) > 0) then
          grow_lo(kref) = grow_lo(kref) - 1
        else
          write (logger%msg, '("refine_unit_partition: refinement failed at interval ", i0)') kref
          if (logger%err) call logger%write_msg
          error stop
        end if
      end do inner
    end do
    do
      sep_lo = fine_sep(1) * refinement(1) ** grow_lo(1)
      dist_inter = pos_lo(fine_lo(1) + grow_lo(1), 1)
      if (dist_inter > 0d0) then  ! refined region does not extend beyond axis
        inter(1) = ceiling(dist_inter / min(coarse_sep, sep_lo * refinement(1)))
        sep_inter = dist_inter / dble(inter(1))
        if (sep_lo <= sep_inter) then
          exit
        end if
      end if
      if (grow_lo(1) > 0) then
        grow_lo(1) = grow_lo(1) - 1
      else
        logger%msg = 'refine_unit_partition: refinement failed at interval 1'
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do
    do
      sep_hi = fine_sep(nref) * refinement(nref) ** grow_hi(nref)
      dist_inter = 1d0 - pos_hi(fine_hi(nref) + grow_hi(nref), nref)
      if (dist_inter > 0d0) then  ! refined region does not extend beyond separatrix
        inter(nref + 1) = ceiling(dist_inter / min(coarse_sep, sep_hi * refinement(nref)))
        sep_inter = dist_inter / dble(inter(nref + 1))
        if (sep_hi <= sep_inter) then
          exit
        end if
      end if
      if (grow_hi(nref) > 0) then
        grow_hi(nref) = grow_hi(nref) - 1
      else
        write (logger%msg, '("refine_unit_partition: refinement failed at interval ", i0)') nref + 1
        if (logger%err) call logger%write_msg
        error stop
      end if
    end do
    ! fill unrefined regions with equidistant intervals
    nflux = sum(fine_lo) + sum(fine_hi) + sum(grow_lo) + sum(grow_hi) + sum(inter) + nref
    allocate(partition(0:nflux))
    partition(:) = 0d0
    kf_hi = 0
    do kref = 1, nref
      kf_lo = kf_hi + inter(kref)
      partition(kf_hi + 1:kf_lo) = linspace(partition(kf_hi), pos_lo(grow_lo(kref) + fine_lo(kref), kref), inter(kref), 1, 0)
      kf_ref(kref) = kf_lo + grow_lo(kref) + fine_lo(kref) + 1
      ! kf_lo + 1 -> grow_lo + fine_lo - 1
      partition(kf_lo + 1:kf_ref(kref) - 1) = pos_lo(grow_lo(kref) + fine_lo(kref) - 1:0:-1, kref)
      kf_hi = kf_ref(kref) + grow_hi(kref) + fine_hi(kref)
      partition(kf_ref(kref):kf_hi) = pos_hi(0:fine_hi(kref) + grow_hi(kref), kref)
    end do
    partition(kf_hi + 1:nflux) = linspace(partition(kf_hi), 1d0, inter(nref + 1), 1, 0)
    deallocate(pos_lo, pos_hi)
  end subroutine refine_unit_partition

  subroutine refine_resonant_surfaces(coarse_sep, fine_sep, add_fine, refinement, widths, rho_norm_ref)
    use mephit_conf, only: conf, logger, refinement_scheme_geometric, refinement_scheme_gaussian
    real(dp), intent(in) :: coarse_sep
    real(dp), dimension(mesh%m_res_min:), intent(in) :: fine_sep
    integer, dimension(mesh%m_res_min:), intent(in) :: add_fine
    real(dp), dimension(mesh%m_res_min:), intent(in) :: refinement
    real(dp), dimension(mesh%m_res_min:), intent(in) :: widths
    real(dp), dimension(:), allocatable, intent(out) :: rho_norm_ref
    logical :: diverging_q
    integer :: m, m_dense, kref
    integer, dimension(:), allocatable :: ref_ind
    logical, dimension(:), allocatable :: mask

    select case (conf%refinement_scheme)
    case (refinement_scheme_geometric)
      allocate(mask(mesh%m_res_min:mesh%m_res_max))
      mask(:) = 1d0 < refinement .and. 0d0 < fine_sep .and. fine_sep < coarse_sep
      if (conf%kilca_scale_factor /= 0) then
        diverging_q = .false.
      else
        ! heuristic: if distance between resonances is too low,
        ! take inner resonance as last to be refined; outside, only the fine separation is used
        m_dense = mesh%m_res_min + 1
        do while (m_dense <= mesh%m_res_max)
          if (mesh%rad_norm_res(m_dense) - mesh%rad_norm_res(m_dense - 1) < &
            sum((0.5d0 + add_fine(m_dense - 1:m_dense) + refinement(m_dense - 1:m_dense)) * &
            fine_sep(m_dense - 1:m_dense))) then
            exit
          end if
          m_dense = m_dense + 1
        end do
        diverging_q = m_dense > mesh%m_res_max
        mask(m_dense:mesh%m_res_max) = .false.
      end if
      allocate(ref_ind(count(mask)))
      call refine_unit_partition(count(mask), coarse_sep, pack(fine_sep, mask), pack(add_fine, mask), &
        pack(refinement, mask), pack(mesh%rad_norm_res, mask), diverging_q, &
        mesh%nflux, rho_norm_ref, ref_ind)
      logger%msg = 'refinement positions:'
      if (logger%debug) call logger%write_msg
      kref = 0
      do m = mesh%m_res_min, mesh%m_res_max
        if (.not. mask(m)) cycle
        kref = kref + 1
        write (logger%msg, '("m = ", i2, ", kf = ", i3, ' // &
          '", rho: ", f19.16, 2(" < ", f19.16))') m, ref_ind(kref), &
          rho_norm_ref(ref_ind(kref) - 1), mesh%rad_norm_res(m), rho_norm_ref(ref_ind(kref))
        if (logger%debug) call logger%write_msg
      end do
      deallocate(ref_ind, mask)
    case (refinement_scheme_gaussian)
      call refine_unit_partition_gaussian(mesh%m_res_max - mesh%m_res_min + 1, coarse_sep, refinement, &
        mesh%rad_norm_res, widths, mesh%nflux, rho_norm_ref)
    case default
      write (logger%msg, '("unknown refinement scheme selection", i0)') conf%refinement_scheme
      if (logger%err) call logger%write_msg
      error stop
    end select
  end subroutine refine_resonant_surfaces

  subroutine cache_resonance_positions
    use mephit_conf, only: logger
    use mephit_util, only: binsearch
    integer :: m, kf_res

    logger%msg = 'resonance positions:'
    if (logger%debug) call logger%write_msg
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    mesh%res_ind = 0
    do m = mesh%m_res_min, mesh%m_res_max
      call binsearch(fs%psi, 0, mesh%psi_res(m), kf_res)
      if (kf_res <= 1) then
        write (logger%msg, '("Warning: resonance for m = ", i0, " occurs at flux surface index ", i0)') m, kf_res
        if (logger%warn) call logger%write_msg
      end if
      mesh%res_ind(m) = kf_res
      write (logger%msg, '("m = ", i2, ", kf = ", i3, ", rho: ", f19.16, 2(" < ", f19.16))') &
        m, kf_res, fs%rad(kf_res - 1) / fs%rad(mesh%nflux), mesh%rad_norm_res(m), &
        fs%rad(kf_res) / fs%rad(mesh%nflux)
      if (logger%debug) call logger%write_msg
    end do
  end subroutine cache_resonance_positions

  subroutine create_mesh_points
    use mephit_conf, only: conf, conf_arr, logger
    use mephit_util, only: interp_psi_pol, resample1d, pos_angle, generate_symfluxcoord
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
      rsmall, circumf, raxis, zaxis
    use field_line_integration_mod, only: circ_mesh_scale, o_point, x_point, theta0_at_xpoint
    use points_2d, only: s_min, create_points_2d

    integer :: kf, kp
    real(dp), dimension(:), allocatable :: rho_norm_eqd, rho_norm_ref, opt_pol_edge_len
    real(dp), dimension(:, :), allocatable :: points, points_s_theta_phi
    real(dp) :: psi_axis, rad_max

    theta0_at_xpoint = .false.
    circ_mesh_scale = conf%kilca_scale_factor
    if (conf%kilca_scale_factor /= 0) then
      ! calculate maximal extent from magnetic axis
      rad_max = min(equil%rmaxis - equil%rleft, &
        equil%rleft + equil%rdim - equil%rmaxis, &
        equil%zdim * 0.5d0 + equil%zmid - equil%zmaxis, &
        equil%zdim * 0.5d0 - equil%zmid + equil%zmaxis)
      o_point = [equil%rmaxis, equil%zmaxis]
      x_point = o_point + [rad_max, 0d0]
    end if
    call generate_symfluxcoord
    mesh%R_O = raxis
    mesh%Z_O = zaxis
    mesh%R_X = X_point(1)
    mesh%Z_X = X_point(2)
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    allocate(psi_fine(nlabel))  ! intermediate step, to be removed
    psi_fine(:) = psisurf(1:) * psipol_max + psi_axis
    ! interpolate between psi and rho
    rad_max = rbeg(nlabel)
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd(:) = rbeg / rad_max

    call compute_resonance_positions(psi_fine, qsaf, rho_norm_eqd)
    call read_profiles
    call compute_auxiliary_profiles
    call compute_resonant_layer_widths
    call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
    call refine_resonant_surfaces(conf%max_Delta_rad / rad_max, conf_arr%Delta_rad_res / rad_max, &
      conf_arr%add_fine, conf_arr%refinement, mesh%delta_rad_mn / rad_max, rho_norm_ref)
    call fs%init(mesh%nflux, .false.)
    call fs_half%init(mesh%nflux, .true.)
    fs%rad(:) = rho_norm_ref
    fs%psi(:) = [(interp_psi_pol(raxis + fs%rad(kf) * rad_max, zaxis), kf = 0, mesh%nflux)]
    fs_half%rad(:) = 0.5d0 * (fs%rad(0:mesh%nflux-1) + fs%rad(1:mesh%nflux))
    fs_half%psi(:) = [(interp_psi_pol(raxis + fs_half%rad(kf) * rad_max, zaxis), kf = 1, mesh%nflux)]
    fs%rad(:) = fs%rad * rad_max
    fs_half%rad(:) = fs_half%rad * rad_max
    call flux_func_cache_check
    call cache_resonance_positions

    call resample1d(psi_fine, rsmall, fs%psi, fs%rsmall, 3)
    call resample1d(psi_fine, circumf, fs%psi, fs%perimeter, 3)
    call resample1d(psi_fine, rsmall, fs_half%psi, fs_half%rsmall, 3)
    call resample1d(psi_fine, circumf, fs_half%psi, fs_half%perimeter, 3)
    allocate(opt_pol_edge_len(mesh%nflux))
    opt_pol_edge_len(:mesh%nflux - 1) = fs_half%rsmall(2:) - fs_half%rsmall(:mesh%nflux - 1)
    opt_pol_edge_len(mesh%nflux) = fs%rsmall(mesh%nflux) - fs%rsmall(mesh%nflux - 1)
    ! averaged radius corresponds to altitude - factor for edge length of equilateral triangle
    ! use geometric mean to reduce effect of radial refinement
    opt_pol_edge_len(:) = 2d0 / sqrt(3d0) * sqrt(opt_pol_edge_len * conf%max_Delta_rad)
    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux))
    allocate(mesh%kt_low(mesh%nflux))
    ! round to even numbers
    mesh%kp_max(:) = 2 * nint(0.5d0 * fs%perimeter(1:) / opt_pol_edge_len)
    if (conf%pol_max < conf%pol_min .and. conf%pol_max > 0) then
      write (logger%msg, '("config%pol_max = ", i0, " < config%pol_min = ", i0, ' // &
        '", number of points per flux surface is fixed at the former.")') conf%pol_max, conf%pol_min
      if (logger%warn) call logger%write_msg
    end if
    if (conf%pol_min > 0) then
      mesh%kp_max(:) = max(mesh%kp_max, 2 * nint(0.5d0 * dble(conf%pol_min)))
    end if
    if (conf%pol_max > 0) then
      mesh%kp_max(:) = min(mesh%kp_max, 2 * nint(0.5d0 * dble(conf%pol_max)))
    end if
    mesh%kp_low(1) = 1
    do kf = 2, mesh%nflux
      mesh%kp_low(kf) = mesh%kp_low(kf-1) + mesh%kp_max(kf-1)
    end do
    mesh%kt_max(1) = mesh%kp_max(1)
    mesh%kt_max(2:) = mesh%kp_max(:(mesh%nflux - 1)) + mesh%kp_max(2:)
    mesh%kt_low(1) = 0
    do kf = 2, mesh%nflux
      mesh%kt_low(kf) = mesh%kt_low(kf-1) + mesh%kt_max(kf-1)
    end do
    mesh%ntri = mesh%kt_low(mesh%nflux) + mesh%kt_max(mesh%nflux)
    mesh%npoint = mesh%kp_low(mesh%nflux) + mesh%kp_max(mesh%nflux)
    mesh%nedge = mesh%ntri + mesh%npoint - 1
    deallocate(opt_pol_edge_len)

    allocate(mesh%node_R(mesh%npoint))
    allocate(mesh%node_Z(mesh%npoint))
    allocate(mesh%node_theta_flux(mesh%npoint))
    allocate(mesh%node_theta_geom(mesh%npoint))
    allocate(points(3, mesh%npoint), points_s_theta_phi(3, mesh%npoint))
    s_min = 1d-16
    ! inp_label 2 to use poloidal psi with magdata_in_symfluxcoord_ext
    ! psi is not normalized by psipol_max, but shifted by -psi_axis
    call create_points_2d(2, mesh%kp_max, points, points_s_theta_phi, r_scaling_func = psi_ref)
    mesh%node_R(1) = mesh%R_O
    mesh%node_R(2:) = points(1, 2:)
    mesh%node_Z(1) = mesh%Z_O
    mesh%node_Z(2:) = points(3, 2:)
    mesh%node_theta_flux(1) = 0d0
    mesh%node_theta_flux(2:) = points_s_theta_phi(2, 2:)
    mesh%node_theta_geom(1) = 0d0
    mesh%node_theta_geom(2:) = pos_angle(atan2(mesh%node_Z(2:) - mesh%Z_O, mesh%node_R(2:) - mesh%R_O))
    ! set theta = 0 exactly
    do kf = 1, mesh%nflux
      mesh%node_theta_flux(mesh%kp_low(kf) + 1) = 0d0
      mesh%node_theta_geom(mesh%kp_low(kf) + 1) = 0d0
    end do
    ! reposition closest point exactly to X point
    kp = minloc(hypot(mesh%node_R((mesh%kp_low(mesh%nflux) + 1):) - mesh%R_X, &
      mesh%node_Z((mesh%kp_low(mesh%nflux) + 1):) - mesh%Z_X), 1)
    mesh%node_R(mesh%kp_low(mesh%nflux) + kp) = mesh%R_X
    mesh%node_Z(mesh%kp_low(mesh%nflux) + kp) = mesh%Z_X
    mesh%R_min = minval(mesh%node_R)
    mesh%R_max = maxval(mesh%node_R)
    mesh%Z_min = minval(mesh%node_Z)
    mesh%Z_max = maxval(mesh%node_Z)
    deallocate(points, points_s_theta_phi)
    deallocate(rho_norm_ref, rho_norm_eqd)

  contains
    function psi_ref(psi_eqd)
      real(dp), dimension(:), intent(in) :: psi_eqd
      real(dp), dimension(size(psi_eqd)) :: psi_ref
      integer :: kf

      if (mesh%nflux /= size(psi_eqd)) then
        call logger%msg_arg_size('psi_ref', 'nflux', 'size(psi_eqd)', &
          mesh%nflux, size(psi_eqd))
        if (logger%err) call logger%write_msg
        error stop
      end if
      psi_ref = [(interp_psi_pol(raxis + rho_norm_ref(kf) * rad_max, zaxis) - &
        psi_axis, kf = 1, mesh%nflux)]
    end function psi_ref
  end subroutine create_mesh_points

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
    use mephit_conf, only: logger
    integer, intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)
    logical :: tri_mask(mesh%ntri)

    common_tri = [0, 0]
    tri_mask = any(mesh%tri_node == knot1, 1) .and. any(mesh%tri_node == knot2, 1)
    select case (count(tri_mask))
    case (0)
      common_tri = [0, 0]
    case (1)
      common_tri = [findloc(tri_mask, .true., 1), 0]
    case (2)
      common_tri = [findloc(tri_mask, .true., 1), findloc(tri_mask, .true., 1, back = .true.)]
      if (any((common_tri(1) == mesh%kt_low + 1) .and. &
        (common_tri(2) == mesh%kt_low + mesh%kt_max))) then
        common_tri = common_tri([2, 1])
      end if
    case default
      write (logger%msg, '("More than two common triangles for knots ", ' // &
        'i0, " and ", i0)') knot1, knot2
      if (logger%debug) call logger%write_msg
    end select
  end subroutine common_triangles

  subroutine connect_mesh_points
    use mephit_conf, only: logger
    use mephit_util, only: pi
    integer :: kf, kp, kp_lo, kp_hi, kt, ktri, ktri_adj, kedge, nodes(4), k
    real(dp) :: mat(3, 3), points(mesh%GL_order), points2(3, mesh%GL2_order)
    logical :: theta2pi

    allocate(mesh%tri_node(3, mesh%ntri))
    allocate(mesh%tri_node_F(mesh%ntri))
    allocate(mesh%tri_theta2pi(mesh%ntri))
    allocate(mesh%tri_theta_extent(2, mesh%ntri))
    allocate(mesh%orient(mesh%ntri))
    allocate(mesh%edge_node(2, mesh%nedge))
    allocate(mesh%edge_tri(2, mesh%nedge))
    allocate(mesh%tri_edge(3, mesh%ntri))
    mesh%tri_node = 0
    mesh%tri_node_F = 0
    mesh%tri_theta2pi = .false.
    mesh%tri_theta_extent = 0d0
    mesh%edge_node = 0
    mesh%edge_tri = 0
    mesh%tri_edge = 0
    kedge = 0
    ktri = mesh%kt_low(1)
    ! define triangles on innermost flux surface
    kf = 1
    do kp = 1, mesh%kp_max(kf)
      ! node indices for triangle with poloidal edge on outer flux surface
      ktri = ktri + 1
      mesh%orient(ktri) = .true.
      mesh%tri_node(1, ktri) = mesh%kp_low(kf) + kp
      mesh%tri_node(2, ktri) = mesh%kp_low(kf) + mod(kp, mesh%kp_max(kf)) + 1
      mesh%tri_node(3, ktri) = mesh%kp_low(kf)
      mesh%tri_node_F(ktri) = 3
      mesh%tri_theta_extent(1, ktri) = mesh%node_theta_geom(mesh%tri_node(1, ktri))
      mesh%tri_theta_extent(2, ktri) = upper_branch(mesh%node_theta_geom(mesh%tri_node(2, ktri)))
      ! triangle and nodes indices for poloidal edge
      kedge = mesh%kp_low(kf) + kp - 1
      mesh%edge_node(:, kedge) = [mesh%tri_node(1, ktri), mesh%tri_node(2, ktri)]
      mesh%edge_tri(1, kedge) = ktri
      mesh%tri_edge(1, ktri) = kedge
      ! triangle indices for radial edge
      kt = ktri - mesh%kt_low(kf)
      ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
      kedge = mesh%npoint + mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf))
      mesh%edge_tri(:, kedge) = [ktri, ktri_adj]
      mesh%tri_edge(2, ktri) = kedge
      mesh%tri_edge(3, ktri_adj) = kedge
    end do
    mesh%tri_theta2pi(mesh%kt_low(1) + mesh%kp_max(1)) = .true.
    ! define triangles on outer flux surfaces
    do kf = 2, mesh%nflux
      kp_lo = 1
      kp_hi = 1
      theta2pi = .false.
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        ! edge i is fixed by nodes(1) and nodes(2)
        ! nodes(3) and nodes(4) are candidates for node I
        nodes(1) = mesh%kp_low(kf - 1) + kp_lo
        nodes(2) = mesh%kp_low(kf) + kp_hi
        nodes(3) = mesh%kp_low(kf) + mod(kp_hi, mesh%kp_max(kf)) + 1
        nodes(4) = mesh%kp_low(kf - 1) + mod(kp_lo, mesh%kp_max(kf - 1)) + 1
        if (kt < mesh%kt_max(kf)) then
          ! test Delaunay condition for nodes(3) as node I
          mat(:, 1) = mesh%node_R(nodes(1:3)) - mesh%node_R(nodes(4))
          mat(:, 2) = mesh%node_Z(nodes(1:3)) - mesh%node_Z(nodes(4))
          mat(:, 3) = mat(:, 1) * mat(:, 1) + mat(:, 2) * mat(:, 2)
          mesh%orient(ktri) = 0d0 >= &
            mat(1, 1) * mat(2, 2) * mat(3, 3) + mat(1, 2) * mat(2, 3) * mat(3, 1) + mat(1, 3) * mat(2, 1) * mat(3, 2) - &
            mat(3, 1) * mat(2, 2) * mat(1, 3) - mat(3, 2) * mat(2, 3) * mat(1, 1) - mat(3, 3) * mat(2, 1) * mat(1, 2)
        else
          ! the last edge is already fixed, so we ignore the Delaunay condition
          if (kp_lo == 1 .and. kp_hi == mesh%kp_max(kf)) then
            mesh%orient(ktri) = .true.
          elseif (kp_lo == mesh%kp_max(kf - 1) .and. kp_hi == 1) then
            mesh%orient(ktri) = .false.
          else
            write (logger%msg, '("Cannot close triangle loop correctly: ' // &
              'kf = ", i0, ", kp_lo = ", i0, ", kp_hi = ", i0)') kf, kp_lo, kp_hi
            if (logger%err) call logger%write_msg
            error stop
          end if
        end if
        if (mesh%orient(ktri)) then
          ! node indices for triangle with poloidal edge on outer flux surface
          mesh%tri_node(:, ktri) = nodes([2, 3, 1])
          mesh%tri_node_F(ktri) = 3
          ! triangle and nodes indices for poloidal edge
          kedge = mesh%kp_low(kf) + kp_hi - 1
          mesh%edge_node(:, kedge) = [mesh%tri_node(1, ktri), mesh%tri_node(2, ktri)]
          mesh%edge_tri(1, kedge) = ktri
          mesh%tri_edge(1, ktri) = kedge
          kp_hi = mod(kp_hi, mesh%kp_max(kf)) + 1
          theta2pi = theta2pi .or. kp_hi == 1
        else
          ! node indices for triangle with poloidal edge on inner flux surface
          mesh%tri_node(:, ktri) = nodes([2, 4, 1])
          mesh%tri_node_F(ktri) = 1
          ! triangle index for poloidal edge
          kedge = mesh%kp_low(kf - 1) + kp_lo - 1
          mesh%edge_tri(2, kedge) = ktri
          mesh%tri_edge(1, ktri) = kedge
          kp_lo = mod(kp_lo, mesh%kp_max(kf - 1)) + 1
          theta2pi = theta2pi .or. kp_lo == 1
        end if
        mesh%tri_theta2pi(ktri) = theta2pi
        ! triangle indices for radial edge
        ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
        kedge = mesh%npoint + mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf))
        mesh%edge_tri(:, kedge) = [ktri, ktri_adj]
        mesh%tri_edge(2, ktri) = kedge
        mesh%tri_edge(3, ktri_adj) = kedge
        ! cache poloidal extent of triangles for point_location_check
        if (theta2pi) then
          mesh%tri_theta_extent(1, ktri) = minval(upper_branch(mesh%node_theta_geom(mesh%tri_node(:, ktri))))
          mesh%tri_theta_extent(2, ktri) = maxval(upper_branch(mesh%node_theta_geom(mesh%tri_node(:, ktri))))
        else
          mesh%tri_theta_extent(1, ktri) = minval(mesh%node_theta_geom(mesh%tri_node(:, ktri)))
          mesh%tri_theta_extent(2, ktri) = maxval(mesh%node_theta_geom(mesh%tri_node(:, ktri)))
        end if
      end do
    end do
    ! set theta = 2 pi exactly
    where (mesh%tri_theta_extent(2, :) < epsilon(1d0))
      mesh%tri_theta_extent(2, :) = 2d0 * pi
    end where
    ! nodes for radial edges
    mesh%edge_node(:, mesh%npoint:) = mesh%tri_node([3, 1], :)
    ! cache bounding boxes
    allocate(mesh%tri_RZ_extent(2, 2, mesh%ntri))
    do ktri = 1, mesh%ntri
      nodes(:3) = mesh%tri_node(:, ktri)
      mesh%tri_RZ_extent(:, 1, ktri) = [minval(mesh%node_R(nodes(:3))), minval(mesh%node_Z(nodes(:3)))]
      mesh%tri_RZ_extent(:, 2, ktri) = [maxval(mesh%node_R(nodes(:3))), maxval(mesh%node_Z(nodes(:3)))]
    end do
    call rtree_init(mesh%ntri, mesh%tri_RZ_extent)
    ! cache areas and 'centroids'
    allocate(mesh%cntr_R(mesh%ntri))
    allocate(mesh%cntr_Z(mesh%ntri))
    allocate(mesh%area(mesh%ntri))
    do ktri = 1, mesh%ntri
      call ring_centered_avg_coord(ktri, mesh%cntr_R(ktri), mesh%cntr_Z(ktri))
      mesh%area(ktri) = tri_area(ktri)
    end do
    ! cache edge midpoints and Gauss-Legendre quadrature evaluation points
    allocate(mesh%mid_R(mesh%nedge))
    allocate(mesh%mid_Z(mesh%nedge))
    allocate(mesh%edge_R(mesh%nedge))
    allocate(mesh%edge_Z(mesh%nedge))
    allocate(mesh%GL_weights(mesh%GL_order))
    allocate(mesh%GL_R(mesh%GL_order, mesh%nedge), mesh%GL_Z(mesh%GL_order, mesh%nedge))
    call gauss_legendre_unit_interval(mesh%GL_order, points, mesh%GL_weights)
    do kedge = 1, mesh%nedge
      mesh%mid_R(kedge) = sum(mesh%node_R(mesh%edge_node(:, kedge))) * 0.5d0
      mesh%mid_Z(kedge) = sum(mesh%node_Z(mesh%edge_node(:, kedge))) * 0.5d0
      mesh%edge_R(kedge) = mesh%node_R(mesh%edge_node(2, kedge)) - mesh%node_R(mesh%edge_node(1, kedge))
      mesh%edge_Z(kedge) = mesh%node_Z(mesh%edge_node(2, kedge)) - mesh%node_Z(mesh%edge_node(1, kedge))
      do k = 1, mesh%GL_order
        mesh%GL_R(k, kedge) = mesh%node_R(mesh%edge_node(2, kedge)) * points(k) + &
          mesh%node_R(mesh%edge_node(1, kedge)) * points(mesh%GL_order - k + 1)
        mesh%GL_Z(k, kedge) = mesh%node_Z(mesh%edge_node(2, kedge)) * points(k) + &
          mesh%node_Z(mesh%edge_node(1, kedge)) * points(mesh%GL_order - k + 1)
      end do
    end do
    allocate(mesh%GL2_weights(mesh%GL2_order))
    allocate(mesh%GL2_R(mesh%GL2_order, mesh%ntri), mesh%GL2_Z(mesh%GL2_order, mesh%ntri))
    mesh%GL2_weights(:) = [0.109951743655322d0, 0.109951743655322d0, 0.109951743655322d0, &
      0.223381589678011d0, 0.223381589678011d0, 0.223381589678011d0]
    points2(:, :) = reshape([ &
      0.816847572980459d0, 0.091576213509771d0, 0.091576213509771d0, &
      0.091576213509771d0, 0.816847572980459d0, 0.091576213509771d0, &
      0.091576213509771d0, 0.091576213509771d0, 0.816847572980459d0, &
      0.108103018168070d0, 0.445948490915965d0, 0.445948490915965d0, &
      0.445948490915965d0, 0.108103018168070d0, 0.445948490915965d0, &
      0.445948490915965d0, 0.445948490915965d0, 0.108103018168070d0 &
      ], [3, mesh%GL2_order])
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        mesh%GL2_R(k, ktri) = sum(mesh%node_R(mesh%tri_node(:, ktri)) * points2(:, k))
        mesh%GL2_Z(k, ktri) = sum(mesh%node_Z(mesh%tri_node(:, ktri)) * points2(:, k))
      end do
    end do
  end subroutine connect_mesh_points

  pure elemental function upper_branch(angle)
    use mephit_util, only: pi
    real(dp), intent(in) :: angle
    real(dp) :: upper_branch
    if (angle <= 0d0) then
      upper_branch = angle + 2d0 * pi
    else
      upper_branch = angle
    end if
  end function upper_branch

  subroutine extend_over_branch(angles, ext, angles_ext)
    use mephit_conf, only: logger
    use mephit_util, only: pi
    real(dp), dimension(:), intent(in) :: angles
    integer, intent(in) :: ext
    real(dp), dimension(-ext:), intent(out) :: angles_ext
    integer :: N

    if (size(angles) + 2 * ext /= size(angles_ext)) then
      call logger%msg_arg_size('extend_over_branch', &
        'size(angles) + 2 * ext', 'size(angles_ext)', &
        size(angles) + 2 * ext, size(angles_ext))
      if (logger%err) call logger%write_msg
      error stop
    end if
    N = size(angles) - 1
    angles_ext(0:N) = angles
    angles_ext(-ext:-1) = angles_ext(N - ext + 1:N) - 2d0 * pi
    angles_ext(N + 1:N + ext) = angles_ext(0:ext - 1) + 2d0 * pi
  end subroutine extend_over_branch

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

  subroutine check_mesh
    use mephit_conf, only: logger
    integer :: kedge, ktri, ktri_check, k_min, k_max, discrepancies, &
      edge_count(mesh%nedge), tri_count(0:mesh%ntri), common_tri(2)

    ! check node indices
    k_min = minval(mesh%tri_node)
    k_max = maxval(mesh%tri_node)
    if (k_min /= 1 .or. k_max /= mesh%npoint) then
      write (logger%msg, '("mesh%tri_node values are out of range [1, ", i0, "]: [", i0, ", ", i0, "]")') &
        mesh%npoint, k_min, k_max
      if (logger%debug) call logger%write_msg
    end if
    k_min = minval(mesh%edge_node)
    k_max = maxval(mesh%edge_node)
    if (k_min /= 1 .or. k_max /= mesh%npoint) then
      write (logger%msg, '("mesh%edge_node values are out of range [1, ", i0, "]: [", i0, ", ", i0, "]")') &
        mesh%npoint, k_min, k_max
      if (logger%debug) call logger%write_msg
    end if
    ! check edge indices
    k_min = minval(mesh%tri_edge)
    k_max = maxval(mesh%tri_edge)
    if (k_min /= 1 .or. k_max /= mesh%nedge) then
      write (logger%msg, '("mesh%tri_edge values are out of range ' // &
        '[1, ", i0, "]: [", i0, ", ", i0, "]")') mesh%nedge, k_min, k_max
      if (logger%debug) call logger%write_msg
    end if
    edge_count = 0
    do ktri = 1, mesh%ntri
      edge_count(mesh%tri_edge(:, ktri)) = edge_count(mesh%tri_edge(:, ktri)) + 1
    end do
    discrepancies = count(2 /= edge_count(:(mesh%kp_low(mesh%nflux) - 1)))
    if (discrepancies > 0) then
      write (logger%msg, '(i0, " discrepancies in edge index counts ' // &
        'from mesh%tri_edge for non-boundary poloidal edges")') discrepancies
      if (logger%debug) call logger%write_msg
    end if
    discrepancies = count(1 /= edge_count(mesh%kp_low(mesh%nflux):(mesh%npoint - 1)))
    if (discrepancies > 0) then
      write (logger%msg, '(i0, " discrepancies in edge index counts ' // &
        'from mesh%tri_edge for boundary poloidal edges")') discrepancies
      if (logger%debug) call logger%write_msg
    end if
    discrepancies = count(2 /= edge_count(mesh%npoint:))
    if (discrepancies > 0) then
      write (logger%msg, '(i0, " discrepancies in edge index counts ' // &
        'from mesh%tri_edge for radial edges")') discrepancies
      if (logger%debug) call logger%write_msg
    end if
    ! check triangle indices
    k_min = minval(mesh%edge_tri)
    k_max = maxval(mesh%edge_tri)
    if (k_min /= 0 .or. k_max /= mesh%ntri) then
      write (logger%msg, '("mesh%edge_tri values are out of range ' // &
        '[0, ", i0, "]: [", i0, ", ", i0, "]")') mesh%ntri, k_min, k_max
      if (logger%debug) call logger%write_msg
    end if
    tri_count = 0
    do kedge = 1, mesh%nedge
      tri_count(mesh%edge_tri(:, kedge)) = tri_count(mesh%edge_tri(:, kedge)) + 1
    end do
    if (tri_count(0) /= mesh%kp_max(mesh%nflux)) then
      write (logger%msg, '("expected ", i0, " boundary triangles, but counted ", i0, " in mesh%edge_tri")') &
        mesh%kp_max(mesh%nflux), tri_count(0)
      if (logger%debug) call logger%write_msg
    end if
    discrepancies = count(3 /= tri_count(1:))
    if (discrepancies > 0) then
      write (logger%msg, '(i0, " discrepancies in triangle index counts ' // &
        'from mesh%edge_tri")') discrepancies
      if (logger%debug) call logger%write_msg
    end if
    ! check consistency of connections
    do kedge = 1, mesh%nedge
      call common_triangles(mesh%edge_node(1, kedge), mesh%edge_node(2, kedge), common_tri)
      if (any(mesh%edge_tri(:, kedge) /= common_tri)) then
        write (logger%msg, '("mesh%edge_tri(:, ", i0, ") = [", i0, ", ", i0, "], ' // &
          'but common_triangle yields [", i0, ", ", i0, "]")') kedge, mesh%edge_tri(:, kedge), common_tri
        if (logger%debug) call logger%write_msg
      end if
    end do
    ! check point locations
    do ktri = 1, mesh%ntri
      ktri_check = point_location(mesh%cntr_R(ktri), mesh%cntr_Z(ktri))
      if (ktri /= ktri_check) then
        write (logger%msg, '("point_location returns ", i0, " for centroid of triangle ", i0)') ktri_check, ktri
        if (logger%debug) call logger%write_msg
      end if
      ktri_check = point_location_check(mesh%cntr_R(ktri), mesh%cntr_Z(ktri))
      if (ktri /= ktri_check) then
        write (logger%msg, '("point_location_check returns ", i0, " for centroid of triangle ", i0)') ktri_check, ktri
        if (logger%debug) call logger%write_msg
      end if
    end do
  end subroutine check_mesh

  subroutine compute_shielding_auxiliaries(s, m)
    use mephit_util, only: pi, clight, interp1d
    type(shielding_t), intent(inout) :: s
    integer, intent(in) :: m
    integer :: kf
    real(dp) :: dq_dpsi, normalized_distance

    kf = mesh%res_ind(m)
    dq_dpsi = interp1d(equil%psi_eqd, equil%qpsi, fs_half%psi(kf), 3, .true.)
    s%coeff = clight * mesh%n / (4d0 * pi * mesh%R_O) * &
      abs(dq_dpsi / (fs_half%q(kf) * fs_half%dp_dpsi(kf))) / &
      (mesh%n * abs(fs%q(kf) - fs%q(kf - 1)))
    s%cross_fade(:) = 0d0
    do kf = 0, mesh%nflux
      normalized_distance = abs(fs%psi(kf) - mesh%psi_res(m)) / mesh%delta_psi_mn(m)
      if (normalized_distance <= 0d0) then
        s%cross_fade(kf) = 1d0
      elseif (normalized_distance >= 1d0) then
        s%cross_fade(kf) = 0d0
      else
        s%cross_fade(kf) = exp(-2d0 * pi / (1 - normalized_distance) * &
          exp(-sqrt(2d0) / normalized_distance))
      end if
    end do
  end subroutine compute_shielding_auxiliaries

  subroutine compute_gpec_jacfac
    integer :: kf, kt, ktri

    allocate(mesh%gpec_jacfac(mesh%nflux))
    mesh%gpec_jacfac(:) = 0d0
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        associate (s => cache%sample_polmodes_half(ktri))
          mesh%gpec_jacfac(kf) = mesh%gpec_jacfac(kf) + &
            s%sqrt_g * s%R * hypot(s%B0_Z, -s%B0_R)
        end associate
      end do
      mesh%gpec_jacfac(kf) = mesh%gpec_jacfac(kf) / dble(mesh%kt_max(kf))
    end do
  end subroutine compute_gpec_jacfac

  !> Compute coarse grid for poloidal mode sampling points
  subroutine compute_sample_polmodes(sample, half_grid)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_conf, only: conf, logger
    use mephit_util, only: pi
    use field_sub, only : field
    type(coord_cache_t), dimension(:), intent(out) :: sample
    logical, intent(in) :: half_grid
    integer :: kf, kpol, k
    integer, dimension(mesh%nflux) :: npol, k_low
    real(dp) :: dum, q
    real(dp), allocatable, dimension(:) :: psi, rad

    if (half_grid) then
      if (size(sample) /= mesh%ntri) then
        call logger%msg_arg_size('compute_sample_polmodes', &
          'size(sample)', 'mesh%ntri', size(sample), mesh%ntri)
        if (logger%err) call logger%write_msg
        error stop
      end if
      allocate(psi, source = fs_half%psi)
      allocate(rad, source = fs_half%rad)
      npol = mesh%kt_max
      k_low = mesh%kt_low
    else
      if (size(sample) /= mesh%npoint) then
        call logger%msg_arg_size('compute_sample_polmodes', &
          'size(sample)', 'mesh%npoint', size(sample), mesh%npoint)
        if (logger%err) call logger%write_msg
        error stop
      end if
      allocate(psi, source = fs%psi)
      allocate(rad, source = fs%rad)
      npol = mesh%kp_max
      k_low = mesh%kp_low
    end if
    do kf = 1, mesh%nflux
      do kpol = 1, npol(kf)
        k = k_low(kf) + kpol
        associate (s => sample(k))
          s%psi = psi(kf)
          s%theta = 2d0 * pi * dble(kpol - 1) / dble(npol(kf))
          if (conf%kilca_scale_factor /= 0 .and. conf%debug_kilca_geom_theta) then
            s%R = mesh%R_O + rad(kf) * cos(s%theta)
            s%Z = mesh%Z_O + rad(kf) * sin(s%theta)
            s%dR_dtheta = -rad(kf) * sin(s%theta)
            s%dZ_dtheta =  rad(kf) * cos(s%theta)
          else
            ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
            call magdata_in_symfluxcoord_ext(2, dum, s%psi - fs%psi(0), s%theta, &
              q, dum, s%sqrt_g, dum, dum, &
              s%R, dum, s%dR_dtheta, s%Z, dum, s%dZ_dtheta)
          end if
          s%ktri = point_location(s%R, s%Z)
          call field(s%R, 0d0, s%Z, s%B0_R, s%B0_phi, s%B0_Z, &
            dum, dum, dum, dum, dum, dum, dum, dum, dum)
          if (conf%kilca_scale_factor /= 0 .and. conf%debug_kilca_geom_theta) then
            s%sqrt_g = equil%cocos%sgn_dpsi * rad(kf) / &
              (-s%B0_R * sin(s%theta) + s%B0_Z * cos(s%theta))
          else
            ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
            ! taken together, these three always yield a positive sign in COCOS 3
            s%sqrt_g = s%sqrt_g * abs(q)
          end if
        end associate
      end do
    end do
    deallocate(psi, rad)
  end subroutine compute_sample_polmodes

  subroutine compute_sample_jnperp(sample_jnperp)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_util, only: interp1d
    type(coord_cache_t), dimension(:), intent(out) :: sample_jnperp
    integer :: kf, kp, kpoi_lo, kpoi_hi, kedge
    real(dp) :: dum, theta_geom_mid, q
    real(dp), dimension(:), allocatable :: theta_flux_ext, theta_geom_ext

    do kf = 1, mesh%nflux
      kpoi_lo = mesh%kp_low(kf) + 1
      kpoi_hi = mesh%kp_low(kf) + mesh%kp_max(kf)
      allocate(theta_flux_ext(mesh%kp_max(kf) + 4), theta_geom_ext(mesh%kp_max(kf) + 4))
      call extend_over_branch(mesh%node_theta_flux(kpoi_lo:kpoi_hi), 2, theta_flux_ext)
      call extend_over_branch(mesh%node_theta_geom(kpoi_lo:kpoi_hi), 2, theta_geom_ext)
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        associate (s => sample_jnperp(kedge), f => cache%mid_fields(kedge))
          theta_geom_mid = upper_branch(atan2(mesh%mid_Z(kedge) - mesh%Z_O, &
                                              mesh%mid_R(kedge) - mesh%R_O))
          s%theta = interp1d(theta_geom_ext, theta_flux_ext, theta_geom_mid, 3)
          ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
          call magdata_in_symfluxcoord_ext(2, dum, fs%psi(kf) - fs%psi(0), s%theta, &
            q, dum, s%sqrt_g, dum, dum, s%R, dum, s%dR_dtheta, s%Z, dum, s%dZ_dtheta)
          s%psi = f%psi
          s%B0_R = f%B0(1)
          s%B0_phi = f%B0(2)
          s%B0_Z = f%B0(3)
          ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
          ! taken together, these three always yield a positive sign in COCOS 3
          s%sqrt_g = s%sqrt_g * abs(q)
        end associate
      end do
      deallocate(theta_flux_ext, theta_geom_ext)
    end do
  end subroutine compute_sample_jnperp

  !> Compute fine grid for parallel current sampling points.
  subroutine compute_sample_Ires(sample_Ires, GL_weights, GL_order, m)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_util, only: pi, binsearch, interp_psi_pol
    use field_sub, only : field
    type(coord_cache_t), dimension(:, :), intent(inout), allocatable :: sample_Ires
    real(dp), dimension(:), intent(inout), allocatable :: GL_weights
    integer, intent(in) :: GL_order
    integer, intent(in) :: m
    integer :: nrad, krad, npol, kpol, kf_min, kf_max, kf, kGL
    real(dp) :: weights(GL_order), points(GL_order), q, dum
    real(dp), allocatable :: theta(:), psi(:)

    npol = mesh%kp_max(mesh%res_ind(m))
    allocate(theta(npol))
    theta(:) = 2d0 * pi * [(dble(kpol - 1), kpol = 1, npol)] / dble(npol)
    call gauss_legendre_unit_interval(GL_order, points, weights)
    call binsearch(fs%rad, lbound(fs%rad, 1), fs%rad(mesh%nflux) * &
      mesh%rad_norm_res(m) - mesh%Delta_rad_mn(m), kf_min)
    call binsearch(fs%rad, lbound(fs%rad, 1), fs%rad(mesh%nflux) * &
      mesh%rad_norm_res(m) + mesh%Delta_rad_mn(m), kf_max)
    nrad = (kf_max - kf_min + 1) * GL_order
    allocate(psi(nrad))
    if (allocated(GL_weights)) deallocate(GL_weights)
    allocate(GL_weights(nrad))
    do kf = kf_min, kf_max
      do kGL = 1, GL_order
        krad = (kf - kf_min) * GL_order + kGL
        psi(krad) = fs%psi(kf) * points(kGL) + &
          fs%psi(kf - 1) * points(GL_order - kGL + 1)
        GL_weights(krad) = weights(kGL) * abs(fs%psi(kf) - fs%psi(kf - 1))
      end do
    end do
    if (allocated(sample_Ires)) deallocate(sample_Ires)
    allocate(sample_Ires(npol, nrad))
    do krad = 1, nrad
      do kpol = 1, npol
        associate (s => sample_Ires(kpol, krad))
          s%psi = psi(krad)
          s%theta = theta(kpol)
          ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad) - fs%psi(0), theta(kpol), &
            q, dum, s%sqrt_g, dum, dum, s%R, dum, s%dR_dtheta, s%Z, dum, s%dZ_dtheta)
          s%ktri = point_location(s%R, s%Z)
          call field(s%R, 0d0, s%Z, s%B0_R, s%B0_phi, s%B0_Z, &
            dum, dum, dum, dum, dum, dum, dum, dum, dum)
          ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
          ! taken together, these three always yield a positive sign in COCOS 3
          s%sqrt_g = s%sqrt_g * abs(q)
        end associate
      end do
    end do
    deallocate(theta, psi)
  end subroutine compute_sample_Ires

  subroutine compute_kilca_auxiliaries
    use mephit_util, only: resample1d
    integer :: kf, kp, kpoi, m
    real(dp), dimension(0:mesh%nflux) :: q_prime

    allocate(mesh%avg_R2gradpsi2(0:mesh%nflux))
    mesh%avg_R2gradpsi2(:) = 0d0
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kpoi = mesh%kp_low(kf) + kp
        associate (s => cache%sample_polmodes(kpoi))
          mesh%avg_R2gradpsi2(kf) = mesh%avg_R2gradpsi2(kf) + &
            s%R ** 4 * (s%B0_Z ** 2 + s%B0_R ** 2)
        end associate
      end do
      mesh%avg_R2gradpsi2(kf) = mesh%avg_R2gradpsi2(kf) / dble(mesh%kp_max(kf))
    end do
    call resample1d(fs%psi, fs%q, fs%psi, q_prime, 3, .true.)
    if (allocated(mesh%damping)) deallocate(mesh%damping)
    allocate(mesh%damping(0:mesh%nflux))
    mesh%damping(:) = 0d0
    do m = mesh%m_res_min, mesh%m_res_max
      mesh%damping(:) = mesh%damping + q_prime / fs%q * mesh%delta_psi_mn(m) * &
        exp(-(fs%psi - mesh%psi_res(m)) ** 2 / mesh%delta_psi_mn(m) ** 2) * mesh%n
    end do
  end subroutine compute_kilca_auxiliaries

  function point_location(R, Z) result(location)
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer
    real(dp), intent(in) :: R, Z
    integer :: location
    type(c_ptr) :: results_p
    integer(c_int) :: result_size
    integer(c_int), dimension(:), pointer :: results
    integer :: k
    real(dp) :: thickness

    location = -1
    thickness = fs%rad(mesh%nflux) * sqrt(epsilon(1d0)) * 8d0
    call Rtree_query(R, Z, result_size, results_p)
    call c_f_pointer(results_p, results, [result_size])
    do k = 1, result_size
      if (point_in_triangle(results(k), R, Z, thickness)) then
        location = results(k)
        return
      end if
    end do
  end function point_location

  function point_location_geom(R, Z, hint_psi) result(location)
    use mephit_util, only: interp_psi_pol, binsearch, pos_angle
    real(dp), intent(in) :: R, Z
    real(dp), intent(in), optional :: hint_psi
    integer :: location
    integer :: kf, kp, kedge, ktri_min, ktri_max, ktri
    real(dp) :: thickness, psi, theta, theta_interval(0:maxval(mesh%kp_max))

    ! return -1 if point is outside computational box
    location = -1
    if (R < mesh%R_min .or. R > mesh%R_max .or. Z < mesh%Z_min .or. Z > mesh%Z_max) return
    ! return -2 if point is outside separatrix
    location = -2
    if (present(hint_psi)) then
      psi = hint_psi
    else
      psi = interp_psi_pol(R, Z)
    end if
    if (equil%cocos%sgn_dpsi == +1) then
      if (psi > fs%psi(mesh%nflux)) return
    else
      if (psi < fs%psi(mesh%nflux)) return
    end if
    ! estimate radial and poloidal position
    call binsearch(fs%psi, lbound(fs%psi, 1), psi, kf)
    theta = pos_angle(atan2(Z - mesh%Z_O, R - mesh%R_O))

    theta_interval = 0d0
    ktri_min = mesh%ntri
    ktri_max = 1
    thickness = fs%rad(mesh%nflux) * sqrt(epsilon(1d0)) * 8d0
    location = -3
    ! if inner flux surface is not the magnetic axis, scan its edge nodes
    if (kf > 1) then
      theta_interval(0:(mesh%kp_max(kf-1) - 1)) = &
        mesh%node_theta_geom((mesh%kp_low(kf-1) + 1):(mesh%kp_low(kf-1) + mesh%kp_max(kf-1)))
      theta_interval(mesh%kp_max(kf-1)) = upper_branch(mesh%node_theta_geom(mesh%kp_low(kf-1) + 1))
      call binsearch(theta_interval(0:mesh%kp_max(kf-1)), lbound(theta_interval, 1), theta, kp)
      kedge = mesh%kp_low(kf-1) + kp - 1
      ktri_min = mesh%edge_tri(2, kedge)
      ktri_max = mesh%edge_tri(2, kedge)
    end if
    ! scan edge nodes of outer flux surface
    theta_interval(0:(mesh%kp_max(kf) - 1)) = &
      mesh%node_theta_geom((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))
    theta_interval(mesh%kp_max(kf)) = upper_branch(mesh%node_theta_geom(mesh%kp_low(kf) + 1))
    call binsearch(theta_interval(0:mesh%kp_max(kf)), lbound(theta_interval, 1), theta, kp)
    kedge = mesh%kp_low(kf) + kp - 1
    ktri_min = min(ktri_min, mesh%edge_tri(1, kedge))
    ktri_max = max(ktri_max, mesh%edge_tri(1, kedge))
    ! scan triangle interval
    do ktri = ktri_min, ktri_max
      if (point_in_triangle(ktri, R, Z, thickness)) then
        location = ktri
        return
      end if
    end do
    ! flux surfaces circumscribe poloidal edges, so parts of the outer triangles may be included
    if (kf < mesh%nflux) then
      kf = kf + 1
      ! re-use previous results for inner flux surface
      ktri_min = mesh%edge_tri(2, kedge)
      ktri_max = mesh%edge_tri(2, kedge)
      ! scan edge nodes of outer flux surface
      theta_interval(0:(mesh%kp_max(kf) - 1)) = &
        mesh%node_theta_geom((mesh%kp_low(kf) + 1):(mesh%kp_low(kf) + mesh%kp_max(kf)))
      theta_interval(mesh%kp_max(kf)) = upper_branch(mesh%node_theta_geom(mesh%kp_low(kf) + 1))
      call binsearch(theta_interval(0:mesh%kp_max(kf)), lbound(theta_interval, 1), theta, kp)
      kedge = mesh%kp_low(kf) + kp - 1
      ktri_min = min(ktri_min, mesh%edge_tri(1, kedge))
      ktri_max = max(ktri_max, mesh%edge_tri(1, kedge))
      ! scan triangle interval
      do ktri = ktri_min, ktri_max
        if (point_in_triangle(ktri, R, Z, thickness)) then
          location = ktri
          return
        end if
      end do
    end if
  end function point_location_geom

  function point_location_check(R, Z, hint_psi) result(location)
    use mephit_util, only: interp_psi_pol, binsearch, pos_angle
    real(dp), intent(in) :: R, Z
    real(dp), intent(in), optional :: hint_psi
    integer :: location
    integer :: kf, ktri_min, ktri_max, ktri
    real(dp) :: thickness, psi, theta
    logical :: tri_mask(mesh%ntri)

    ! return -1 if point is outside computational box
    location = -1
    if (R < mesh%R_min .or. R > mesh%R_max .or. Z < mesh%Z_min .or. Z > mesh%Z_max) return
    ! return -2 if point is outside separatrix
    location = -2
    if (present(hint_psi)) then
      psi = hint_psi
    else
      psi = interp_psi_pol(R, Z)
    end if
    if (equil%cocos%sgn_dpsi == +1) then
      if (psi > fs%psi(mesh%nflux)) return
    else
      if (psi < fs%psi(mesh%nflux)) return
    end if
    ! estimate radial and poloidal position
    call binsearch(fs%psi, lbound(fs%psi, 1), psi, kf)
    theta = pos_angle(atan2(Z - mesh%Z_O, R - mesh%R_O))

    ktri_min = mesh%ntri
    ktri_max = 1
    thickness = fs%rad(mesh%nflux) * sqrt(epsilon(1d0)) * 8d0
    location = -3
    ! scan flux surface
    tri_mask = .false.
    tri_mask((mesh%kt_low(kf) + 1):(mesh%kt_low(kf) + mesh%kt_max(kf))) = .true.
    ktri_min = findloc(mesh%tri_theta_extent(1, :) <= theta .and. &
      theta <= mesh%tri_theta_extent(2, :), .true., 1, tri_mask)
    ktri_max = findloc(mesh%tri_theta_extent(1, :) <= theta .and. &
      theta <= mesh%tri_theta_extent(2, :), .true., 1, tri_mask, back = .true.)
    ! scan triangle interval
    do ktri = ktri_min, ktri_max
      if (point_in_triangle(ktri, R, Z, thickness)) then
        location = ktri
        return
      end if
    end do
    ! flux surfaces circumscribe poloidal edges, so parts of the outer triangles may be included
    if (kf < mesh%nflux) then
      kf = kf + 1
      tri_mask = .false.
      tri_mask((mesh%kt_low(kf) + 1):(mesh%kt_low(kf) + mesh%kt_max(kf))) = .true.
      ktri_min = findloc(mesh%tri_theta_extent(1, :) <= theta .and. &
        theta <= mesh%tri_theta_extent(2, :), .true., 1, tri_mask)
      ktri_max = findloc(mesh%tri_theta_extent(1, :) <= theta .and. &
        theta <= mesh%tri_theta_extent(2, :), .true., 1, tri_mask, back = .true.)
      ! scan triangle interval
      do ktri = ktri_min, ktri_max
        if (point_in_triangle(ktri, R, Z, thickness)) then
          location = ktri
          return
        end if
      end do
    end if
  end function point_location_check

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

  function mesh_interp_theta_flux(R, Z, hint_ktri) result(theta)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mephit_util, only: pi
    real(dp), intent(in) :: R, Z
    integer, intent(in), optional :: hint_ktri
    real(dp) :: theta
    integer :: ktri, nodes(3)
    real(dp) :: DOF(3), Delta_R(3), Delta_Z(3)

    if (present(hint_ktri)) then
      ktri = hint_ktri
    else
      ktri = point_location(R, Z)
    end if
    if (ktri <= 0 .or. ktri > mesh%ntri) then
      theta = ieee_value(1d0, ieee_quiet_nan)
      return
    end if
    nodes = mesh%tri_node(:, ktri)
    DOF = mesh%node_theta_flux(nodes)
    if (mesh%tri_theta2pi(ktri)) then
      where (abs(DOF) <= 0d0)  ! suppress compiler warning about exact comparison
        DOF = 2d0 * pi
      end where
    end if
    if (nodes(3) == mesh%kp_low(1)) then
      ! avoid singular point at axis
      DOF(3) = 0.5d0 * (DOF(1) + DOF(2))
    end if
    Delta_R = R - mesh%node_R(nodes)
    Delta_Z = Z - mesh%node_Z(nodes)
    theta = sum(DOF * (Delta_R([2, 3, 1]) * Delta_Z([3, 1, 2]) - &
      Delta_R([3, 1, 2]) * Delta_Z([2, 3, 1]))) * 0.5d0 / mesh%area(ktri)
  end function mesh_interp_theta_flux

  subroutine mesh_write(mesh, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(mesh_t), intent(in) :: mesh
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer :: orient(mesh%ntri), tri_theta2pi(mesh%ntri)

    where (mesh%orient)
      orient = 1
    elsewhere
      orient = 0
    end where
    where (mesh%tri_theta2pi)
      tri_theta2pi = 1
    elsewhere
      tri_theta2pi = 0
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
      comment = 'normalized small radius (outboard from O point, or towards X point)' // &
      ' in resonance with given poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rsmall_res', mesh%rsmall_res, &
      lbound(mesh%rsmall_res), ubound(mesh%rsmall_res), unit = 'cm', &
      comment = 'normalized small radius (of equivalent-area circle) in resonance with given poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/delta_mn', mesh%delta_mn, &
      lbound(mesh%delta_mn), ubound(mesh%delta_mn), &
      comment = 'resonant layer width', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/delta_psi_mn', mesh%delta_psi_mn, &
      lbound(mesh%delta_psi_mn), ubound(mesh%delta_psi_mn), &
      comment = 'resonant layer width (in units of psi)', unit = 'Mx')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/delta_rad_mn', mesh%delta_rad_mn, &
      lbound(mesh%delta_rad_mn), ubound(mesh%delta_rad_mn), &
      comment = 'resonant layer width (in units of rad)', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/damping', mesh%damping, &
      lbound(mesh%damping), ubound(mesh%damping), &
      comment = 'damping factors', unit = '1')
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_theta2pi', tri_theta2pi, &
      lbound(tri_theta2pi), ubound(tri_theta2pi), &
      comment = 'true (1) if theta = 2 pi instead of theta = 0')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_theta_extent', mesh%tri_theta_extent, &
      lbound(mesh%tri_theta_extent), ubound(mesh%tri_theta_extent), &
      comment = 'range of geometrical poloidal angle covered by triangle')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_RZ_extent', mesh%tri_RZ_extent, &
      lbound(mesh%tri_RZ_extent), ubound(mesh%tri_RZ_extent), &
      comment = 'bounding box (R_min, Z_min, R_max, Z_max) of triangle')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/orient', orient, &
      lbound(orient), ubound(orient), &
      comment = 'triangle orientation: true (1) if edge f lies on outer flux surface')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_node', mesh%edge_node, &
      lbound(mesh%edge_node), ubound(mesh%edge_node), &
      comment = 'nodes connected by edge with global index')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_tri', mesh%edge_tri, &
      lbound(mesh%edge_tri), ubound(mesh%edge_tri), &
      comment = 'triangles connected by edge with global index')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_edge', mesh%tri_edge, &
      lbound(mesh%tri_edge), ubound(mesh%tri_edge), &
      comment = 'global edge indices for triangles')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/mid_R', mesh%mid_R, &
      lbound(mesh%mid_R), ubound(mesh%mid_R), &
      comment = 'R coordinate of edge midpoint', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/mid_Z', mesh%mid_Z, &
      lbound(mesh%mid_Z), ubound(mesh%mid_Z), &
      comment = 'Z coordinate of edge midpoint', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_R', mesh%edge_R, &
      lbound(mesh%edge_R), ubound(mesh%edge_R), &
      comment = 'R coordinate of edge vector', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/edge_Z', mesh%edge_Z, &
      lbound(mesh%edge_Z), ubound(mesh%edge_Z), &
      comment = 'Z coordinate of edge vector', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/area', mesh%area, &
      lbound(mesh%area), ubound(mesh%area), &
      comment = 'triangle area', unit = 'cm^2')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cntr_R', mesh%cntr_R, &
      lbound(mesh%cntr_R), ubound(mesh%cntr_R), &
      comment = 'R coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/cntr_Z', mesh%cntr_Z, &
      lbound(mesh%cntr_Z), ubound(mesh%cntr_Z), &
      comment = 'Z coordinate of triangle ''centroid''', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL_order', mesh%GL_order, &
      comment = 'order of Gauss-Legendre quadrature')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL_weights', mesh%GL_weights, &
      lbound(mesh%GL_weights), ubound(mesh%GL_weights), &
      comment = 'weights of Gauss-Legendre quadrature', unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL_R', mesh%GL_R, &
      lbound(mesh%GL_R), ubound(mesh%GL_R), &
      comment = 'R coordinate of Gauss-Legendre points on edge', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL_Z', mesh%GL_Z, &
      lbound(mesh%GL_Z), ubound(mesh%GL_Z), &
      comment = 'Z coordinate of Gauss-Legendre points on edge', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL2_order', mesh%GL2_order, &
      comment = 'number of Gauss-Legendre quadrature points on triangle')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL2_weights', mesh%GL2_weights, &
      lbound(mesh%GL2_weights), ubound(mesh%GL2_weights), &
      comment = 'weights of 2D Gauss-Legendre quadrature', unit = '1')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL2_R', mesh%GL2_R, &
      lbound(mesh%GL2_R), ubound(mesh%GL2_R), &
      comment = 'R coordinate of Gauss-Legendre points on triangle', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/GL2_Z', mesh%GL2_Z, &
      lbound(mesh%GL2_Z), ubound(mesh%GL2_Z), &
      comment = 'Z coordinate of Gauss-Legendre points on triangle', unit = 'cm')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac, &
      lbound(mesh%gpec_jacfac), ubound(mesh%gpec_jacfac), unit = 'cm^2', &
      comment = 'Jacobian surface factor between flux surfaces')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/avg_R2gradpsi2', mesh%avg_R2gradpsi2, &
      lbound(mesh%avg_R2gradpsi2), ubound(mesh%avg_R2gradpsi2), unit = 'Mx^2', &
      comment = 'Flux surface average (on flux surfaces) for KiLCA coupling')
    call h5_close(h5id_root)
  end subroutine mesh_write

#ifdef USE_MFEM
  subroutine mesh_write_MFEM
    use mephit_conf, only: basename_suffix, decorate_filename
    integer :: fid, ktri, kp, kpoi

    open(newunit = fid, file = decorate_filename('core_plasma.mesh', '', basename_suffix), &
      status = 'replace', form = 'formatted', action = 'write')
    write (fid, '("MFEM mesh v1.0", /)')
    write (fid, '("dimension", /, "2", /)')
    write (fid, '("elements", /, i0)') mesh%ntri
    do ktri = 1, mesh%ntri
      ! <element attribute> <geometry type> <vertex indices ...>
      ! attribute: 1 for plasma core, geometry type: 2 for triangle
      write (fid, '("1 2", 3(1x, i0))') mesh%tri_node(:, ktri) - 1
    end do
    write (fid, '(/, "boundary", /, i0)') mesh%kp_max(mesh%nflux)
    do kp = 1, mesh%kp_max(mesh%nflux)
      ! <boundary element attribute> <geometry type> <vertex indices ...>
      ! attribute: 1 for plasma core, geometry type: 1 for segment
      write (fid, '("1 1", 2(1x, i0))') mesh%kp_low(mesh%nflux) + kp - 1, &
        mesh%kp_low(mesh%nflux) + mod(kp, mesh%kp_max(mesh%nflux))
    end do
    write (fid, '(/, "vertices", /, i0, /, "2")') mesh%npoint
    do kpoi = 1, mesh%npoint
      write (fid, '(es24.16e3, 1x, es24.16e3)') mesh%node_R(kpoi), mesh%node_Z(kpoi)
    end do
    close(fid)
  end subroutine mesh_write_MFEM
#endif

  subroutine write_FreeFem_mesh
    use iso_c_binding, only: c_null_char
    use mephit_util, only: pi, linspace
    use mephit_conf, only: basename_suffix, decorate_filename
    integer :: fid, kpoi, ktri, kp, kedge, npt_inner, npt_outer
    real(dp) :: R_min, R_max, R_mid, R_rad, Z_min, Z_max, Z_mid, Z_rad
    real(dp), parameter :: outer_border_refinement = 0.125d0, outer_box_scale = 2d0
    real(dp), allocatable :: bdry_R(:), bdry_Z(:), theta(:)

    open(newunit = fid, file = decorate_filename('core_plasma.msh', '', basename_suffix), &
      status = 'replace', form = 'formatted', action = 'write')
    write (fid, '(3(1x, i0))') mesh%npoint, mesh%ntri, mesh%kp_max(mesh%nflux) - 1
    do kpoi = 1, mesh%npoint
      write (fid, '(2(1x, es23.15e3), 1x, i0)') &
        mesh%node_R(kpoi), mesh%node_Z(kpoi), 0
    end do
    do ktri = 1, mesh%ntri
      write (fid, '(4(1x, i0))') mesh%tri_node(:, ktri), 0
    end do
    do kp = 1, mesh%kp_max(mesh%nflux)
      write (fid, '(4(1x, i0))') mesh%kp_low(mesh%nflux) + kp, &
        mesh%kp_low(mesh%nflux) + mod(kp, mesh%kp_max(mesh%nflux)) + 1, 1
    end do
    close(fid)
    ! edge numbering as defined in mephit_mesh::connect_mesh_points
    open(newunit = fid, file = decorate_filename('edgemap.dat', '', basename_suffix), &
      status = 'replace', form = 'formatted', action = 'write')
    ! poloidal edges: node indices are sorted in ascending order except for the last triangle
    do kedge = 1, mesh%npoint - 1
      write (fid, '(2(1x, i0))') mesh%edge_tri(1, kedge), sign(1, &
        mesh%edge_node(2, kedge) - mesh%edge_node(1, kedge))
    end do
    ! radial edges: use edge i for natural ordering
    do kedge = mesh%npoint, mesh%nedge
      write (fid, '(2(1x, i0))') mesh%edge_tri(2, kedge), 3
    end do
    close(fid)

    ! extend mesh
    R_min = minval(mesh%node_R)
    R_max = maxval(mesh%node_R)
    R_mid = 0.5d0 * (R_max + R_min)
    R_rad = 0.5d0 * (R_max - R_min) * outer_box_scale
    Z_min = minval(mesh%node_Z)
    Z_max = maxval(mesh%node_Z)
    Z_mid = 0.5d0 * (Z_max + Z_min)
    Z_rad = 0.5d0 * (Z_max - Z_min) * outer_box_scale
    npt_inner = mesh%kp_max(mesh%nflux)
    npt_outer = nint(outer_border_refinement * npt_inner)
    allocate(bdry_R(npt_inner + npt_outer), bdry_Z(npt_inner + npt_outer))
    kpoi = mesh%kp_low(mesh%nflux) + 1
    bdry_R(:npt_inner) = mesh%node_R(kpoi:)
    bdry_Z(:npt_inner) = mesh%node_Z(kpoi:)
    allocate(theta(npt_outer))
    theta(:) = linspace(0d0, 2d0 * pi, npt_outer, 0, 1)
    bdry_R(npt_inner+1:) = R_mid + R_rad * cos(theta)
    bdry_Z(npt_inner+1:) = Z_mid + Z_rad * sin(theta)
    call FEM_triangulate_external(npt_inner, npt_outer, bdry_R, bdry_Z, R_mid, Z_mid, &
      decorate_filename('outer.msh', '', basename_suffix) // c_null_char)
    deallocate(bdry_R, bdry_Z, theta)
  end subroutine write_FreeFem_mesh

  subroutine mesh_read(mesh, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use mephit_conf, only: conf
    type(mesh_t), intent(inout) :: mesh
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer, allocatable :: orient(:), tri_theta2pi(:)

    call mesh_deinit(mesh)
    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_O', mesh%R_O)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_O', mesh%Z_O)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/R_X', mesh%R_X)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/Z_X', mesh%Z_X)
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL_order', mesh%GL_order)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL2_order', mesh%GL2_order)
    ! TODO: allocate deferred-shape arrays in hdf5_tools and skip allocation here
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%psi_res(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%rad_norm_res(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%rsmall_res(mesh%m_res_min:mesh%m_res_max))
    if (conf%kilca_scale_factor /= 0 .and. conf%kilca_pol_mode /= 0) then
      allocate(mesh%res_modes(1))
    else
      allocate(mesh%res_modes(mesh%m_res_max - mesh%m_res_min + 1))
    end if
    allocate(mesh%delta_mn(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%delta_psi_mn(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%delta_rad_mn(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%damping(0:mesh%nflux))
    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux))
    allocate(mesh%kt_low(mesh%nflux))
    allocate(mesh%node_R(mesh%npoint))
    allocate(mesh%node_Z(mesh%npoint))
    allocate(mesh%node_theta_flux(mesh%npoint))
    allocate(mesh%node_theta_geom(mesh%npoint))
    allocate(mesh%tri_node(3, mesh%ntri))
    allocate(mesh%tri_node_F(mesh%ntri))
    allocate(tri_theta2pi(mesh%ntri))
    allocate(mesh%tri_theta2pi(mesh%ntri))
    allocate(mesh%tri_theta_extent(2, mesh%ntri))
    allocate(mesh%tri_RZ_extent(2, 2, mesh%ntri))
    allocate(mesh%orient(mesh%ntri))
    allocate(orient(mesh%ntri))
    allocate(mesh%edge_node(2, mesh%nedge))
    allocate(mesh%edge_tri(2, mesh%nedge))
    allocate(mesh%tri_edge(3, mesh%ntri))
    allocate(mesh%mid_R(mesh%nedge))
    allocate(mesh%mid_Z(mesh%nedge))
    allocate(mesh%edge_R(mesh%nedge))
    allocate(mesh%edge_Z(mesh%nedge))
    allocate(mesh%GL_weights(mesh%GL_order))
    allocate(mesh%GL_R(mesh%GL_order, mesh%nedge))
    allocate(mesh%GL_Z(mesh%GL_order, mesh%nedge))
    allocate(mesh%GL2_weights(mesh%GL2_order))
    allocate(mesh%GL2_R(mesh%GL2_order, mesh%nedge))
    allocate(mesh%GL2_Z(mesh%GL2_order, mesh%nedge))
    allocate(mesh%gpec_jacfac(mesh%nflux))
    allocate(mesh%avg_R2gradpsi2(0:mesh%nflux))
    allocate(mesh%area(mesh%ntri))
    allocate(mesh%cntr_R(mesh%ntri))
    allocate(mesh%cntr_Z(mesh%ntri))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_ind', mesh%res_ind)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi_res', mesh%psi_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rad_norm_res', mesh%rad_norm_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rsmall_res', mesh%rsmall_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_modes', mesh%res_modes)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/delta_mn', mesh%delta_mn)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/delta_psi_mn', mesh%delta_psi_mn)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/delta_rad_mn', mesh%delta_rad_mn)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/damping', mesh%damping)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kp_max', mesh%kp_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kp_low', mesh%kp_low)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kt_max', mesh%kt_max)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/kt_low', mesh%kt_low)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_R', mesh%node_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_Z', mesh%node_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_theta_flux', mesh%node_theta_flux)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/node_theta_geom', mesh%node_theta_geom)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_node', mesh%tri_node)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_node_F', mesh%tri_node_F)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_theta2pi', tri_theta2pi)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_theta_extent', mesh%tri_theta_extent)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_RZ_extent', mesh%tri_RZ_extent)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/orient', orient)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_node', mesh%edge_node)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_tri', mesh%edge_tri)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_edge', mesh%tri_edge)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/mid_R', mesh%mid_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/mid_Z', mesh%mid_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_R', mesh%edge_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/edge_Z', mesh%edge_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/area', mesh%area)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cntr_R', mesh%cntr_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/cntr_Z', mesh%cntr_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL_weights', mesh%GL_weights)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL_R', mesh%GL_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL_Z', mesh%GL_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL2_weights', mesh%GL2_weights)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL2_R', mesh%GL2_R)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/GL2_Z', mesh%GL2_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/avg_R2gradpsi2', mesh%avg_R2gradpsi2)
    call h5_close(h5id_root)
    where (orient == 1)
      mesh%orient = .true.
    elsewhere
      mesh%orient = .false.
    end where
    deallocate(orient)
    where (tri_theta2pi == 1)
      mesh%tri_theta2pi = .true.
    elsewhere
      mesh%tri_theta2pi = .false.
    end where
    deallocate(tri_theta2pi)
    call rtree_init(mesh%ntri, mesh%tri_RZ_extent)
  end subroutine mesh_read

  subroutine mesh_deinit(this)
    class(mesh_t), intent(inout) :: this

    if (allocated(this%res_ind)) deallocate(this%res_ind)
    if (allocated(this%psi_res)) deallocate(this%psi_res)
    if (allocated(this%rad_norm_res)) deallocate(this%rad_norm_res)
    if (allocated(this%rsmall_res)) deallocate(this%rsmall_res)
    if (allocated(this%res_modes)) deallocate(this%res_modes)
    if (allocated(this%delta_mn)) deallocate(this%delta_mn)
    if (allocated(this%delta_psi_mn)) deallocate(this%delta_psi_mn)
    if (allocated(this%delta_rad_mn)) deallocate(this%delta_rad_mn)
    if (allocated(this%damping)) deallocate(this%damping)
    if (allocated(this%kp_max)) deallocate(this%kp_max)
    if (allocated(this%kt_max)) deallocate(this%kt_max)
    if (allocated(this%kp_low)) deallocate(this%kp_low)
    if (allocated(this%kt_low)) deallocate(this%kt_low)
    if (allocated(this%node_R)) deallocate(this%node_R)
    if (allocated(this%node_Z)) deallocate(this%node_Z)
    if (allocated(this%node_theta_flux)) deallocate(this%node_theta_flux)
    if (allocated(this%node_theta_geom)) deallocate(this%node_theta_geom)
    if (allocated(this%tri_node)) deallocate(this%tri_node)
    if (allocated(this%tri_node_F)) deallocate(this%tri_node_F)
    if (allocated(this%tri_theta2pi)) deallocate(this%tri_theta2pi)
    if (allocated(this%tri_theta_extent)) deallocate(this%tri_theta_extent)
    if (allocated(this%tri_RZ_extent)) deallocate(this%tri_RZ_extent)
    if (allocated(this%orient)) deallocate(this%orient)
    if (allocated(this%edge_node)) deallocate(this%edge_node)
    if (allocated(this%edge_tri)) deallocate(this%edge_tri)
    if (allocated(this%tri_edge)) deallocate(this%tri_edge)
    if (allocated(this%mid_R)) deallocate(this%mid_R)
    if (allocated(this%mid_Z)) deallocate(this%mid_Z)
    if (allocated(this%edge_R)) deallocate(this%edge_R)
    if (allocated(this%edge_Z)) deallocate(this%edge_Z)
    if (allocated(this%area)) deallocate(this%area)
    if (allocated(this%cntr_R)) deallocate(this%cntr_R)
    if (allocated(this%cntr_Z)) deallocate(this%cntr_Z)
    if (allocated(this%GL_weights)) deallocate(this%GL_weights)
    if (allocated(this%GL_R)) deallocate(this%GL_R)
    if (allocated(this%GL_Z)) deallocate(this%GL_Z)
    if (allocated(this%GL2_weights)) deallocate(this%GL2_weights)
    if (allocated(this%GL2_R)) deallocate(this%GL2_R)
    if (allocated(this%GL2_Z)) deallocate(this%GL2_Z)
    if (allocated(this%gpec_jacfac)) deallocate(this%gpec_jacfac)
    if (allocated(this%avg_R2gradpsi2)) deallocate(this%avg_R2gradpsi2)
  end subroutine mesh_deinit

  subroutine compare_gpec_coordinates
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, nf90_inq_varid, &
      nf90_inquire_dimension, nf90_get_var, nf90_close, nf90_global, nf90_get_att, nf90_strerror
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext, psipol_max
    use mephit_conf, only: conf, logger, datafile
    use mephit_util, only: pi
    use field_sub, only : field
    character(len = *), parameter :: dataset = 'debug_GPEC'
    character(len = 1024) :: filename
    logical :: file_exists
    integer(HID_T) :: h5id_root
    integer :: ncid_file, ncid, nrad, npol, krad, kpol, idum, fid
    integer, parameter :: stepsize = 10
    real(dp) :: dum, B0_R, B0_Z, unit_normal(0:1), q, chi1
    real(dp), allocatable :: psi(:), theta(:), R(:, :), Z(:, :), xi_n(:, :, :), &
      jac(:, :), sqrt_g(:, :), delpsi(:, :), grad_psi(:, :), contradenspsi(:, :)
    complex(dp), allocatable :: xi_n_R(:), xi_n_Z(:), jacfac(:, :)

    write (filename, '("gpec_control_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (logger%msg, '("File ", a, " found, performing GPEC comparison.")') &
      trim(filename)
    if (logger%info) call logger%write_msg
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_get_att", nf90_get_att(ncid_file, nf90_global, 'chi1', chi1))
    call check_error("nf90_close", nf90_close(ncid_file))
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/chi1', chi1, &
      unit = '?', comment = 'GPEC poloidal flux normalization')
    call h5_add(h5id_root, dataset // '/psipol_max', psipol_max, &
      unit = 'Mx', comment = 'MEPHIT poloidal flux normalization')
    call h5_close(h5id_root)

    write (filename, '("gpec_profile_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (logger%msg, '("File ", a, " found, performing GPEC coordinate comparison.")') &
      trim(filename)
    if (logger%info) call logger%write_msg
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "psi_n", ncid))
    call check_error("nf90_inquire_dimension", &
      nf90_inquire_dimension(ncid_file, ncid, len = nrad))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "theta_dcon", ncid))
    call check_error("nf90_inquire_dimension", &
      nf90_inquire_dimension(ncid_file, ncid, len = npol))
    allocate(psi(nrad), theta(npol), R(npol, nrad), Z(npol, nrad), xi_n(0:1, npol, nrad))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "psi_n", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, psi))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "theta_dcon", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, theta))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "R", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, R, &
      count = flip1d(shape(R)), map = flip1d(stride1d(shape(R)))))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "z", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Z, &
      count = flip1d(shape(Z)), map = flip1d(stride1d(shape(Z)))))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "xi_n_fun", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, xi_n, &
      count = flip1d(shape(xi_n)), map = flip1d(stride1d(shape(xi_n)))))
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
      unit = 'cm', comment = 'R(theta, psi) from GPEC')
    call h5_add(h5id_root, dataset // '/Z_GPEC', Z, lbound(Z), ubound(Z), &
      unit = 'cm', comment = 'Z(theta, psi) from GPEC')
    do krad = 1, nrad
      do kpol = 1, npol
        call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
          dum, dum, dum, dum, dum, R(kpol, krad), dum, dum, Z(kpol, krad), dum, dum)
      end do
    end do
    allocate(xi_n_R(npol), xi_n_Z(npol))
    krad = nrad
    do kpol = 1, npol
      call field(R(kpol, krad), 0d0, Z(kpol, krad), B0_R, dum, B0_Z, &
        dum, dum, dum, dum, dum, dum, dum, dum, dum)
      unit_normal = [R(kpol, krad) * B0_Z, -R(kpol, krad) * B0_R] * equil%cocos%sgn_dpsi
      unit_normal = unit_normal / norm2(unit_normal)
      xi_n_R(kpol) = unit_normal(0) * 1d2 * cmplx(xi_n(0, kpol, krad), xi_n(1, kpol, krad), dp)
      xi_n_Z(kpol) = unit_normal(1) * 1d2 * cmplx(xi_n(0, kpol, krad), xi_n(1, kpol, krad), dp)
    end do
    call h5_add(h5id_root, dataset // '/R', R, lbound(R), ubound(R), &
      unit = 'cm', comment = 'R(theta, psi)')
    call h5_add(h5id_root, dataset // '/Z', Z, lbound(Z), ubound(Z), &
      unit = 'cm', comment = 'Z(theta, psi)')
    call h5_add(h5id_root, dataset // '/xi_n_R', xi_n_R, lbound(xi_n_R), ubound(xi_n_R), &
      unit = 'cm', comment = 'Radial component of normal displacement xi_n(theta) at last flux surface')
    call h5_add(h5id_root, dataset // '/xi_n_Z', xi_n_Z, lbound(xi_n_Z), ubound(xi_n_Z), &
      unit = 'cm', comment = 'Axial component of normal displacement xi_n(theta) at last flux surface')
    call h5_close(h5id_root)
    deallocate(psi, theta, R, Z, xi_n, xi_n_R, xi_n_Z)

    nrad = (nrad - 1 + stepsize) / stepsize
    allocate(psi(nrad), theta(npol), R(npol, nrad), Z(npol, nrad), sqrt_g(npol, nrad))
    filename = 'gpec_diagnostics_jacfac_1_fun.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    logger%msg = 'File ' // trim(filename) // ' found, performing GPEC jacfac comparison.'
    if (logger%info) call logger%write_msg
    allocate(jacfac(npol, nrad), contradenspsi(npol, nrad))
    open(newunit = fid, file = filename, status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
      do kpol = 1, npol
        read (fid, '(1x, 4(1x, es16.8))') psi(krad), theta(kpol), jacfac(kpol, krad)
      end do
    end do
    close(fid)
    psi(:) = psipol_max * psi
    theta(:) = 2d0 * pi * theta
    filename = 'gpec_diagnostics_delpsi_fun.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    logger%msg = 'File ' // trim(filename) // ' found, performing GPEC grad psi comparison.'
    if (logger%info) call logger%write_msg
    allocate(delpsi(npol, nrad), grad_psi(npol, nrad))
    open(newunit = fid, file = filename, status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
      do kpol = 1, npol
        read (fid, '(1x, 3(1x, es16.8))') psi(krad), theta(kpol), delpsi(kpol, krad)
      end do
    end do
    close(fid)
    psi(:) = psipol_max * psi
    theta(:) = 2d0 * pi * theta
    delpsi(:, :) = 1d-2 * psipol_max * delpsi
    do krad = 1, nrad
      do kpol = 1, npol
        call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
          q, dum, sqrt_g(kpol, krad), dum, dum, R(kpol, krad), dum, dum, Z(kpol, krad), dum, dum)
        call field(R(kpol, krad), 0d0, Z(kpol, krad), B0_R, dum, B0_Z, &
          dum, dum, dum, dum, dum, dum, dum, dum, dum)
        grad_psi(kpol, krad) = R(kpol, krad) * hypot(B0_Z, -B0_R)
        contradenspsi(kpol, krad) = grad_psi(kpol, krad) * sqrt_g(kpol, krad) * abs(q)
      end do
      ! normalize by flux surface average - dimensionless factor
      ! conversion factors for delpsi and jac cancel in GPEC jacfac
      contradenspsi(:, krad) = contradenspsi(:, krad) / sum(contradenspsi(:, krad)) * dble(npol)
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/jacfac', jacfac, lbound(jacfac), ubound(jacfac), &
      unit = '1', comment = 'GPEC jacfac at (theta, psi)')
    call h5_add(h5id_root, dataset // '/contradenspsi', contradenspsi, lbound(contradenspsi), ubound(contradenspsi), &
      unit = '1', comment = 'MEPHIT jacfac at (theta, psi)')
    call h5_add(h5id_root, dataset // '/delpsi', delpsi, lbound(delpsi), ubound(delpsi), &
      unit = 'G cm', comment = 'GPEC grad psi at (theta, psi)')
    call h5_add(h5id_root, dataset // '/grad_psi', grad_psi, lbound(grad_psi), ubound(grad_psi), &
      unit = 'G cm', comment = 'MEPHIT grad psi at (theta, psi)')
    call h5_close(h5id_root)
    deallocate(psi, theta, R, Z, delpsi, grad_psi, jacfac, sqrt_g)

    filename = '2d.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (filename, '("dcon_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    write (logger%msg, '("Files ", a, " and 2d.out found, performing GPEC Jacobian comparison.")') &
      trim(filename)
    if (logger%info) call logger%write_msg
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_get_att", nf90_get_att(ncid_file, nf90_global, 'mpsi', nrad))
    call check_error("nf90_get_att", nf90_get_att(ncid_file, nf90_global, 'mtheta', npol))
    call check_error("nf90_close", nf90_close(ncid_file))
    nrad = nrad + 1
    npol = npol + 1
    allocate(psi(nrad), theta(npol), jac(npol, nrad), sqrt_g(npol, nrad), R(npol, nrad), Z(npol, nrad))
    open(newunit = fid, file = '2d.out', status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
      read (fid, '(6x, i3, 6x, e11.3)') idum, psi(krad)
      psi(krad) = psi(krad) * psipol_max
      read (fid, *)
      read (fid, *)
      read (fid, *)
      do kpol = 1, npol
        read (fid, '(i6, 1p, 8e11.3)') idum, theta(kpol), dum, dum, dum, dum, &
          R(kpol, krad), Z(kpol, krad), jac(kpol, krad)
        theta(kpol) = theta(kpol) * 2d0 * pi
        call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
          q, dum, sqrt_g(kpol, krad), dum, dum, dum, dum, dum, dum, dum, dum)
        sqrt_g(kpol, krad) = sqrt_g(kpol, krad) * abs(q)
      end do
      read (fid, *)
      read (fid, *)
      read (fid, *)
    end do
    R(:, :) = R * 1d2  ! m to cm
    Z(:, :) = Z * 1d2  ! m to cm
    jac(:, :) = jac / (2d2 * pi * chi1)  ! m per T to cm per G, normalization in GPEC
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/jac/')
    call h5_add(h5id_root, dataset // '/jac/psi', psi, lbound(psi), ubound(psi), &
      unit = 'Mx', comment = 'poloidal flux (shifted to zero at axis)')
    call h5_add(h5id_root, dataset // '/jac/theta', theta, lbound(theta), ubound(theta), &
      unit = 'rad', comment = 'flux poloidal angle')
    call h5_add(h5id_root, dataset // '/jac/R', R, lbound(R), ubound(R), &
      unit = 'cm', comment = 'R(theta, psi)')
    call h5_add(h5id_root, dataset // '/jac/Z', Z, lbound(Z), ubound(Z), &
      unit = 'cm', comment = 'Z(theta, psi)')
    call h5_add(h5id_root, dataset // '/jac/sqrt_g', sqrt_g, lbound(sqrt_g), ubound(sqrt_g), &
      unit = 'cm G^-1', comment = 'MEPHIT Jacobian at (theta, psi)')
    call h5_add(h5id_root, dataset // '/jac/jac', jac, lbound(jac), ubound(jac), &
      unit = 'cm G^-1', comment = 'GPEC Jacobian at (theta, psi)')
    call h5_close(h5id_root)
    deallocate(psi, theta, jac, sqrt_g, R, Z)

  contains

    subroutine check_error(funcname, status)
      character(len = *), intent(in) :: funcname
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
        write (logger%msg, '(a, " returned error ", i0, ": ", a)') funcname, status, nf90_strerror(status)
        if (logger%err) call logger%write_msg
        error stop
      end if
    end subroutine check_error

    pure function flip1d(arr)
      integer, dimension(:), intent(in) :: arr
      integer, dimension(size(arr, 1)) :: flip1d
      flip1d = arr(ubound(arr, 1):lbound(arr, 1):-1)
    end function flip1d

    pure function stride1d(arr)
      integer, dimension(:), intent(in) :: arr
      integer, dimension(size(arr, 1)) :: stride1d
      integer :: k
      stride1d = [1, (product(arr(:k)), k = 1, size(arr, 1) - 1)]
    end function stride1d
  end subroutine compare_gpec_coordinates

  subroutine write_illustration_data(npsi, ntheta, nrad, npol)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext, psipol_max
    use mephit_util, only: pi, linspace
    integer, intent(in) :: npsi, ntheta, nrad, npol
    integer :: fid, krad, kpol
    real(dp) :: dum
    real(dp), allocatable :: psi(:), theta(:), R(:, :), Z(:, :)

    open(newunit = fid, file = 'illustration.asc', status = 'replace', form = 'formatted', action = 'write')
    write (fid, '(6(1x, i0))') npsi, ntheta, nrad, npol, equil%nbbbs, equil%limitr
    allocate(psi(npsi), theta(npol), R(npol, npsi), Z(npol, npsi))
    psi(:) = linspace(0d0, psipol_max, npsi, 1, 1)
    theta(:) = linspace(0d0, 2d0 * pi, npol, 0, 1)
    do krad = 1, npsi
      do kpol = 1, npol
        call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
          dum, dum, dum, dum, dum, R(kpol, krad), dum, dum, Z(kpol, krad), dum, dum)
      end do
    end do
    do krad = 1, npsi
      do kpol = 1, npol
        write (fid, '(es24.16e3, 1x, es24.16e3)') R(kpol, krad), Z(kpol, krad)
      end do
    end do
    deallocate(psi, theta, R, Z)
    allocate(psi(nrad), theta(ntheta), R(nrad, ntheta), Z(nrad, ntheta))
    psi(:) = linspace(0d0, psipol_max, nrad, 0, 0)
    theta(:) = linspace(0d0, 2d0 * pi, ntheta, 0, 1)
    do kpol = 1, ntheta
      do krad = 2, nrad
        call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
          dum, dum, dum, dum, dum, R(krad, kpol), dum, dum, Z(krad, kpol), dum, dum)
      end do
    end do
    do kpol = 1, ntheta
      write (fid, '(es24.16e3, 1x, es24.16e3)') mesh%R_O, mesh%Z_O
      do krad = 2, nrad
        write (fid, '(es24.16e3, 1x, es24.16e3)') R(krad, kpol), Z(krad, kpol)
      end do
    end do
    deallocate(psi, theta, R, Z)
    do kpol = 1, equil%nbbbs
      write (fid, '(es24.16e3, 1x, es24.16e3)') equil%rbbbs(kpol), equil%zbbbs(kpol)
    end do
    do kpol = 1, equil%limitr
      write (fid, '(es24.16e3, 1x, es24.16e3)') equil%rlim(kpol), equil%zlim(kpol)
    end do
    close(fid)
  end subroutine write_illustration_data

  subroutine init_flux_variables
    use mephit_conf, only: conf, logger, pres_prof_eps, pres_prof_par, pres_prof_geqdsk, &
      q_prof_flux, q_prof_rot, q_prof_geqdsk
    use mephit_util, only: resample1d

    select case (conf%pres_prof)
    case (pres_prof_eps)
      call compute_pres_prof_eps
    case (pres_prof_par)
      call compute_pres_prof_par
    case (pres_prof_geqdsk)
      call compute_pres_prof_geqdsk
    case default
      write (logger%msg, '("unknown pressure profile selection", i0)') conf%pres_prof
      if (logger%err) call logger%write_msg
      error stop
    end select

    select case (conf%q_prof)
    case (q_prof_flux)
      write (logger%msg, '("q profile selection ", i0, ' // &
        '" currently not supported.")') conf%q_prof
      if (logger%err) call logger%write_msg
      error stop
    case (q_prof_rot)
      call compute_safety_factor_rot
    case (q_prof_geqdsk)
      call compute_safety_factor_geqdsk
    case default
      write (logger%msg, '("unknown q profile selection: ", i0)') conf%q_prof
      if (logger%err) call logger%write_msg
      error stop
    end select

    call resample1d(equil%psi_eqd, equil%fpol, fs%psi, fs%F, 3)
    call resample1d(equil%psi_eqd, equil%ffprim, fs%psi, fs%FdF_dpsi, 3)
    call resample1d(equil%psi_eqd, equil%fpol, fs_half%psi, fs_half%F, 3)
    call resample1d(equil%psi_eqd, equil%ffprim, fs_half%psi, fs_half%FdF_dpsi, 3)
  end subroutine init_flux_variables

  subroutine compute_pres_prof_eps
    use mephit_conf, only: conf, logger
    use mephit_util, only: ev2erg
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
    write (logger%msg, '("temp@axis: ", es24.16e3, ", dens@axis: ", es24.16e3)') &
      temp(0), dens(0)
    if (logger%info) call logger%write_msg
    fs%p(:) = dens * temp * ev2erg
    fs%dp_dpsi(:) = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
    dens(1:) = (fs_half%psi - psi_ext) / psi_int * conf%dens_max + conf%dens_min
    temp(1:) = (fs_half%psi - psi_ext) / psi_int * conf%temp_max + conf%temp_min
    fs_half%p(:) = dens(1:) * temp(1:) * ev2erg
    fs_half%dp_dpsi(:) = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
  end subroutine compute_pres_prof_eps

  subroutine compute_pres_prof_par
    use mephit_conf, only: conf
    use mephit_util, only: ev2erg
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
    fs%p(:) = dens * temp * ev2erg
    fs%dp_dpsi(:) = (dens * dtemp_dpsi + ddens_dpsi * temp) * ev2erg
    dens(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
      (conf%dens_max - conf%dens_min) + conf%dens_min
    temp(1:) = (fs_half%psi - psi_ext) / (psi_int - psi_ext) * &
      (conf%temp_max - conf%temp_min) + conf%temp_min
    fs_half%p(:) = dens(1:) * temp(1:) * ev2erg
    fs_half%dp_dpsi(:) = (dens(1:) * dtemp_dpsi + ddens_dpsi * temp(1:)) * ev2erg
  end subroutine compute_pres_prof_par

  subroutine compute_pres_prof_geqdsk
    use mephit_util, only: resample1d

    call resample1d(equil%psi_eqd, equil%pres, fs%psi, fs%p, 3)
    call resample1d(equil%psi_eqd, equil%pprime, fs%psi, fs%dp_dpsi, 3)
    call resample1d(equil%psi_eqd, equil%pres, fs_half%psi, fs_half%p, 3)
    call resample1d(equil%psi_eqd, equil%pprime, fs_half%psi, fs_half%dp_dpsi, 3)
  end subroutine compute_pres_prof_geqdsk

  subroutine compute_safety_factor_flux
    use mephit_util, only: pi, resample1d
    integer :: kf, kt, ktri

    fs_half%q = 0d0
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        fs_half%q(kf) = fs_half%q(kf) + cache%cntr_fields(ktri)%B0(2) * mesh%area(ktri)
      end do
      fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    ! Lagrange polynomial extrapolation for values at separatrix and magnetic axis
    call resample1d(fs_half%psi, fs_half%q, fs%psi, fs%q, 3)
  end subroutine compute_safety_factor_flux

  subroutine compute_safety_factor_rot
    use magdata_in_symfluxcoor_mod, only: qsaf
    use mephit_util, only: resample1d

    ! Lagrange polynomial extrapolation for value at magnetic axis
    call resample1d(psi_fine, qsaf, fs%psi, fs%q, 3)
    call resample1d(psi_fine, qsaf, fs_half%psi, fs_half%q, 3)
  end subroutine compute_safety_factor_rot

  subroutine compute_safety_factor_geqdsk
    use mephit_util, only: resample1d

    call resample1d(equil%psi_eqd, equil%qpsi, fs%psi, fs%q, 3)
    call resample1d(equil%psi_eqd, equil%qpsi, fs_half%psi, fs_half%q, 3)
  end subroutine compute_safety_factor_geqdsk

  subroutine check_resonance_positions
    use mephit_conf, only: logger
    integer :: m, kf
    real(dp), dimension(mesh%nflux) :: abs_err

    do m = mesh%m_res_max, mesh%m_res_min, -1
      abs_err = [(abs(abs(fs_half%q(kf)) - dble(m) / dble(mesh%n)), kf = 1, mesh%nflux)]
      kf = minloc(abs_err, 1)
      if (kf /= mesh%res_ind(m)) then
        write (logger%msg, '("m = ", i0, ": q is resonant at index ", i0, ' // &
          '", but psi is resonant at index ", i0)') m, kf, mesh%res_ind(m)
        if (logger%warn) call logger%write_msg
      end if
    end do
  end subroutine check_resonance_positions

  subroutine check_safety_factor
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdata_in_symfluxcoor_mod, only: psisurf, qsaf
    use mephit_conf, only: datafile
    use mephit_util, only: pi
    character(len = *), parameter :: grp = 'debug_q'
    integer(HID_T) :: h5id_root

    integer :: kf, kt, ktri
    real(dp) :: step_q(mesh%nflux), step_psi_norm(mesh%nflux), eqd_psi_norm(equil%nw)

    step_q(:) = 0d0
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        step_q(kf) = step_q(kf) + cache%cntr_fields(ktri)%B0(2) * mesh%area(ktri)
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

  subroutine equilibrium_field(R, Z, B0, dB0_dR, dB0_dZ, psi, Bmod, dBmod_dR, dBmod_dZ)
    use field_eq_mod, only: psib
    use field_sub, only : field_eq, psif
    real(dp), intent(in) :: R, Z
    real(dp), intent(out), dimension(3) :: B0, dB0_dR, dB0_dZ
    real(dp), intent(out) :: psi, Bmod, dBmod_dR, dBmod_dZ
    real(dp), dimension(3) :: dB0_dphi

    call field_eq(R, 0d0, Z, B0(1), B0(2), B0(3), dB0_dR(1), dB0_dphi(1), dB0_dZ(1), &
      dB0_dR(2), db0_dphi(2), dB0_dZ(2), dB0_dR(3), dB0_dphi(3), dB0_dZ(3))
    psi = psif - psib  ! see intperp_psi_pol in mephit_util
    Bmod = sqrt(sum(B0 ** 2))
    dBmod_dR = (sum(B0 * dB0_dR)) / Bmod
    dBmod_dZ = (sum(B0 * dB0_dZ)) / Bmod
  end subroutine equilibrium_field

  subroutine cache_equilibrium_field
    integer :: ktri, kedge, k

    ! edge midpoints
    do kedge = 1, mesh%nedge
      associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge), Z => mesh%mid_Z(kedge))
        call equilibrium_field(R, Z, f%B0, f%dB0_dR, f%dB0_dZ, f%psi, f%Bmod, f%dBmod_dR, f%dBmod_dZ)
        f%theta = mesh_interp_theta_flux(R, Z, mesh%edge_tri(1, kedge))
      end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
      associate (f => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri), Z => mesh%cntr_Z(ktri))
        call equilibrium_field(R, Z, f%B0, f%dB0_dR, f%dB0_dZ, f%psi, f%Bmod, f%dBmod_dR, f%dBmod_dZ)
        f%theta = mesh_interp_theta_flux(R, Z, ktri)
      end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
      do k = 1, mesh%GL_order
        associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
          call equilibrium_field(R, Z, f%B0, f%dB0_dR, f%dB0_dZ, f%psi, f%Bmod, f%dBmod_dR, f%dBmod_dZ)
          f%theta = mesh_interp_theta_flux(R, Z, mesh%edge_tri(1, kedge))
        end associate
      end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
          call equilibrium_field(R, Z, f%B0, f%dB0_dR, f%dB0_dZ, f%psi, f%Bmod, f%dBmod_dR, f%dBmod_dZ)
          f%theta = mesh_interp_theta_flux(R, Z, ktri)
        end associate
      end do
    end do
  end subroutine cache_equilibrium_field

  subroutine compute_curr0
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mephit_conf, only: conf, logger, curr_prof_ps, curr_prof_rot, curr_prof_geqdsk
    real(dp) :: nan
    integer :: kedge, ktri, k

    nan = ieee_value(1d0, ieee_quiet_nan)
    do kedge = 1, mesh%nedge
      cache%mid_fields(kedge)%j0(:) = nan
      cache%mid_fields(kedge)%dj0_dR(:) = nan
      cache%mid_fields(kedge)%dj0_dZ(:) = nan
    end do
    do ktri = 1, mesh%ntri
      cache%cntr_fields(ktri)%j0(:) = nan
      cache%cntr_fields(ktri)%dj0_dR(:) = nan
      cache%cntr_fields(ktri)%dj0_dZ(:) = nan
    end do
    do kedge = 1, mesh%nedge
      do k = 1, mesh%GL_order
        cache%edge_fields(k, kedge)%j0(:) = nan
        cache%edge_fields(k, kedge)%dj0_dR(:) = nan
        cache%edge_fields(k, kedge)%dj0_dZ(:) = nan
      end do
    end do
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        cache%area_fields(k, ktri)%j0(:) = nan
        cache%area_fields(k, ktri)%dj0_dR(:) = nan
        cache%area_fields(k, ktri)%dj0_dZ(:) = nan
      end do
    end do

    select case (conf%curr_prof)
    case (curr_prof_ps)
      write (logger%msg, '("current profile selection ", i0, " currently not supported.")') conf%curr_prof
      if (logger%err) call logger%write_msg
      error stop
    case (curr_prof_rot)
      call compute_curr0_rot
    case (curr_prof_geqdsk)
      call compute_curr0_geqdsk
    case default
      write (logger%msg, '("unknown current profile selection: ", i0)') conf%curr_prof
      if (logger%err) call logger%write_msg
      error stop
    end select
  end subroutine compute_curr0

  ! TODO: debugging routine to calculate PS current
  subroutine compute_curr0_ps
    use mephit_util, only: clight
    integer :: kf, kt, ktri, kp, kedge
    real(dp), dimension(mesh%nflux) :: B2avg, B2avg_half

    B2avg = 0d0
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        associate (f => cache%mid_fields(kedge))
          B2avg(kf) = B2avg(kf) + f%Bmod ** 2
        end associate
      end do
      B2avg(kf) = B2avg(kf) / mesh%kp_max(kf)
    end do
    B2avg_half = 0d0
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
        associate (f => cache%mid_fields(kedge))
          B2avg_half(kf) = B2avg_half(kf) + f%Bmod ** 2
        end associate
      end do
      B2avg_half(kf) = B2avg_half(kf) / mesh%kt_max(kf)
    end do
    ! edges in poloidal direction
    do kf = 1, mesh%nflux
      do kp = 1, mesh%kp_max(kf)
        kedge = mesh%kp_low(kf) + kp - 1
        associate (f => cache%mid_fields(kedge))
          f%j0(2) = clight * mesh%mid_R(kedge) * fs%dp_dpsi(kf) * &
            (1d0 - f%B0(2) ** 2 / B2avg(kf))
        end associate
      end do
    end do
    ! edges in radial direction
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
        associate (f => cache%mid_fields(kedge))
          f%j0(2) = clight * mesh%mid_R(kedge) * fs_half%dp_dpsi(kf) * &
            (1d0 - f%B0(2) ** 2 / B2avg_half(kf))
        end associate
      end do
    end do
    ! weighted triangle centroids
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        associate (f => cache%cntr_fields(ktri))
          f%j0(2) = clight * mesh%cntr_R(ktri) * fs_half%dp_dpsi(kf) * &
            (1d0 - f%B0(2) ** 2 / B2avg_half(kf))
        end associate
      end do
    end do
  end subroutine compute_curr0_ps

  subroutine compute_curr0_rot
    integer :: ktri, kedge, k

    ! edge midpoints
    do kedge = 1, mesh%nedge
      associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge))
        f%j0(:) = curr0_rot(f, R)
      end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
      associate (f => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri))
        f%j0(:) = curr0_rot(f, R)
      end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
      do k = 1, mesh%GL_order
        associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge))
          f%j0(:) = curr0_rot(f, R)
        end associate
      end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri))
          f%j0(:) = curr0_rot(f, R)
        end associate
      end do
    end do

  contains
    function curr0_rot(f, R)
      use mephit_util, only: clight, pi
      type(field_cache_t), intent(in) :: f
      real(dp), intent(in) :: R
      real(dp), dimension(3) :: curr0_rot
      curr0_rot = 0.25d0 / pi * clight * [-f%dB0_dZ(2), f%dB0_dZ(1) - f%dB0_dR(3), f%dB0_dR(2) + f%B0(2) / R]
    end function curr0_rot
  end subroutine compute_curr0_rot

  subroutine curr0_geqdsk(R, psi, B0, dB0_dR, dB0_dZ, j0, dj0_dR, dj0_dZ)
    use mephit_util, only: clight, pi
    real(dp), intent(in) :: R, psi
    real(dp), intent(in), dimension(3) :: B0, dB0_dR, dB0_dZ
    real(dp), intent(out), dimension(3) :: j0, dj0_dR, dj0_dZ
    real(dp) :: dp0_dpsi, d2p0_dpsi2, F, dF_dpsi, FdF_dpsi, d2F_dpsi2

    dp0_dpsi = interp1d(equil%psi_eqd, equil%pprime, psi, 3)
    d2p0_dpsi2 = interp1d(equil%psi_eqd, equil%pprime, psi, 3, .true.)
    F = interp1d(equil%psi_eqd, equil%fpol, psi, 3)
    FdF_dpsi = interp1d(equil%psi_eqd, equil%ffprim, psi, 3)
    dF_dpsi = interp1d(equil%psi_eqd, equil%fprime, psi, 3)
    d2F_dpsi2 = interp1d(equil%psi_eqd, equil%fprime, psi, 3, .true.)
    j0(1) = 0.25d0 / pi * clight * dF_dpsi * B0(1)
    j0(3) = 0.25d0 / pi * clight * dF_dpsi * B0(3)
    j0(2) = clight * (dp0_dpsi * R + 0.25d0 / (pi * R) * FdF_dpsi)
    dj0_dR(1) = 0.25d0 / pi * clight * (dF_dpsi * dB0_dR(1) &
      + R * B0(3) * B0(1) * d2F_dpsi2)
    dj0_dZ(1) = 0.25d0 / pi * clight * (dF_dpsi * dB0_dZ(1) &
      - R * B0(1) * B0(1) * d2F_dpsi2)
    dj0_dR(3) = 0.25d0 / pi * clight * (dF_dpsi * dB0_dR(3) &
      + R * B0(3) * B0(3) * d2F_dpsi2)
    dj0_dZ(3) = 0.25d0 / pi * clight * (dF_dpsi * dB0_dZ(3) &
      - R * B0(1) * B0(3) * d2F_dpsi2)
    dj0_dR(2) = clight * (dp0_dpsi + d2p0_dpsi2 * R * R * B0(3) + &
      0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2 - FdF_dpsi / R))
    dj0_dZ(2) = clight * (-d2p0_dpsi2 * R * R * B0(1) + &
      0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2))
  end subroutine curr0_geqdsk

  subroutine compute_curr0_geqdsk
    integer :: ktri, kedge, k

    ! edge midpoints
    do kedge = 1, mesh%nedge
      associate (c => cache%mid_fields(kedge), R => mesh%mid_R(kedge))
        call curr0_geqdsk(R, c%psi, c%B0, c%dB0_dR, c%dB0_dZ, c%j0, c%dj0_dR, c%dj0_dZ)
      end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
      associate (c => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri))
        call curr0_geqdsk(R, c%psi, c%B0, c%dB0_dR, c%dB0_dZ, c%j0, c%dj0_dR, c%dj0_dZ)
      end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
      do k = 1, mesh%GL_order
        associate (c => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge))
          call curr0_geqdsk(R, c%psi, c%B0, c%dB0_dR, c%dB0_dZ, c%j0, c%dj0_dR, c%dj0_dZ)
        end associate
      end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        associate (c => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri))
          call curr0_geqdsk(R, c%psi, c%B0, c%dB0_dR, c%dB0_dZ, c%j0, c%dj0_dR, c%dj0_dZ)
        end associate
      end do
    end do
  end subroutine compute_curr0_geqdsk

  subroutine check_curr0
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_conf, only: datafile
    use mephit_util, only: clight, pi, linspace
    use field_sub, only : field
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

  subroutine flux_func_cache_check
    use mephit_conf, only: logger
    logger%msg = 'checking flux_func_cache...'
    if (logger%debug) call logger%write_msg
    write (logger%msg, '("array bounds: fs%psi(", i0, ":", i0, "), ' // &
      'fs%rad(", i0, ":", i0, "), fs_half%psi(", i0, ":", i0, "), ' // &
      'fs_half%rad(", i0, ":", i0, ")")') lbound(fs%psi, 1), ubound(fs%psi, 1), &
      lbound(fs%rad, 1), ubound(fs%rad, 1), lbound(fs_half%psi, 1), &
      ubound(fs_half%psi, 1), lbound(fs_half%rad, 1), ubound(fs_half%rad, 1)
    if (logger%debug) call logger%write_msg
    write (logger%msg, '("expected sign of psi''(r): ", sp, i0, ss)') equil%cocos%sgn_dpsi
    if (logger%debug) call logger%write_msg
    write (logger%msg, '(i0, " ordering violations for psi")') &
      count((fs%psi(1:) - fs_half%psi) * equil%cocos%sgn_dpsi <= 0d0) + &
      count([(fs_half%psi(1) - fs%psi(0)) * equil%cocos%sgn_dpsi] <= 0d0)
    if (logger%debug) call logger%write_msg
    write (logger%msg, '(i0, " ordering violations for rad")') &
      count(fs%rad(1:) <= fs_half%rad) + count([fs_half%rad(1)] <= [fs%rad(0)])
    if (logger%debug) call logger%write_msg
  end subroutine flux_func_cache_check

end module mephit_mesh
