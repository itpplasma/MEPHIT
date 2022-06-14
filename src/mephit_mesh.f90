module mephit_mesh

  use iso_fortran_env, only: dp => real64
  use mephit_util, only: g_eqdsk, interp1d

  implicit none

  private

  ! types and associated procedures
  public :: flux_func_cache, flux_func_cache_init, flux_func_cache_deinit
  public :: mesh_t, mesh_write, mesh_read, mesh_deinit, &
       generate_mesh, write_cache, read_cache, point_location
  public :: coord_cache_t, coord_cache_ext_t, cache_t, cache_init, cache_deinit, &
       cache_write, cache_read, compute_sample_Ipar

  ! testing and debugging procedures
  public :: check_mesh, write_illustration_data, flux_func_cache_check, &
       check_safety_factor, check_curr0

  ! module variables
  public :: equil, psi_interpolator, psi_fine_interpolator, fs, fs_half, mesh, &
       cache

  type(g_eqdsk) :: equil
  type(interp1d) :: psi_interpolator
  type(interp1d) :: psi_fine_interpolator

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

     !> Area of poloidal cross-section, in cm^2
     real(dp), dimension(:), allocatable, public :: area

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

     !> Number of unrefined flux surfaces to be replaced by refined ones.
     integer, allocatable :: deletions(:)

     !> Width ratio of neighbouring refined flux surfaces.
     real(dp), allocatable :: refinement(:)

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
     real(dp), allocatable :: tri_theta_extent(:, :)

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

     !> Triangle indexing for resonant flux surfaces only.
     integer, allocatable :: shielding_kt_low(:)
     !> Local node index for shielding pressure evaluation points.
     integer, allocatable :: shielding_L1_kp(:, :, :)
     !> Weighting factor for shielding pressure evaluation points.
     real(dp), allocatable :: shielding_L1_weight(:, :, :)
     !> R coordinate of shielding pressure evaluation points.
     real(dp), allocatable :: shielding_L1_R(:, :, :)
     !> Z coordinate of shielding pressure evaluation points.
     real(dp), allocatable :: shielding_L1_Z(:, :, :)
     !> Free parameter in the compensated scheme for shielding
     real(dp), allocatable :: shielding_coeff(:)

     !> Surface integral Jacobian used for normalization of poloidal modes in GPEC
     complex(dp), allocatable :: gpec_jacfac(:, :)

  end type mesh_t

  type(mesh_t) :: mesh

  type :: coord_cache_t
     integer :: ktri
     real(dp) :: R, Z, psi, theta, sqrt_g, B0_R, B0_Z, dR_dtheta, dZ_dtheta
  end type coord_cache_t

  type, extends(coord_cache_t) :: coord_cache_ext_t
     real(dp) :: rad, q, dq_dpsi, dR_dpsi, dZ_dpsi, d2R_dpsi_dtheta, d2Z_dpsi_dtheta, &
          B0_phi, dB0R_dR, dB0R_dZ, dB0phi_dR, dB0phi_dZ, dB0Z_dR, dB0Z_dZ, B0_2
  end type coord_cache_ext_t

  type :: field_cache_t
     real(dp) :: psi, B0(3), j0(3), Bmod, dBmod_dR, dBmod_dZ, &
          dB0_dR(3), dB0_dZ(3), dj0_dR(3), dj0_dZ(3)
  end type field_cache_t

  type :: cache_t
     integer :: nrad, npol
     type(coord_cache_t), allocatable :: sample_polmodes_half(:, :), sample_polmodes(:, :)
     type(field_cache_t), allocatable :: edge_fields(:, :), area_fields(:, :), mid_fields(:), cntr_fields(:)
     real(dp), allocatable :: B0_flux(:)
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
       allocate(this%area(nflux))
       allocate(this%perimeter(nflux))
    else
       allocate(this%psi(0:nflux))
       allocate(this%rad(0:nflux))
       allocate(this%F(0:nflux))
       allocate(this%p(0:nflux))
       allocate(this%FdF_dpsi(0:nflux))
       allocate(this%dp_dpsi(0:nflux))
       allocate(this%q(0:nflux))
       allocate(this%area(0:nflux))
       allocate(this%perimeter(0:nflux))
    end if
    this%psi = 0d0
    this%rad = 0d0
    this%F = 0d0
    this%p = 0d0
    this%FdF_dpsi = 0d0
    this%dp_dpsi = 0d0
    this%q = 0d0
    this%area = 0d0
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
    if (allocated(this%area)) deallocate(this%area)
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/area', &
         cache%area, lbound(cache%area), ubound(cache%area), unit = 'cm^2', &
         comment = 'poloidal cross-section area ' // trim(adjustl(comment)))
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/area', cache%area)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/perimeter', cache%perimeter)
    call h5_close(h5id_root)
  end subroutine flux_func_cache_read

  subroutine coord_cache_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(coord_cache_t), dimension(:, :), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/ktri', cache%ktri, &
         lbound(cache), ubound(cache), &
         comment = 'triangle index of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/R', cache%R, &
         lbound(cache), ubound(cache), unit = 'cm', &
         comment = 'R coordinate of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/Z', cache%Z, &
         lbound(cache), ubound(cache), unit = 'cm', &
         comment = 'Z coordinate of ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', cache%theta, &
         lbound(cache), ubound(cache), unit = 'Mx', &
         comment = 'poloidal flux at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/theta', cache%theta, &
         lbound(cache), ubound(cache), unit = 'rad', &
         comment = 'flux poloidal angle at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/sqrt_g', cache%sqrt_g, &
         lbound(cache), ubound(cache), unit = 'cm G^-1', &
         comment = 'Jacobian at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache%B0_R, &
         lbound(cache), ubound(cache), unit = 'G', &
         comment = 'R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache%B0_Z, &
         lbound(cache), ubound(cache), unit = 'G', &
         comment = 'Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dR_dtheta', cache%dR_dtheta, &
         lbound(cache), ubound(cache), unit = 'cm rad^-1', &
         comment = 'Jacobian element (R, theta) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dZ_dtheta', cache%dZ_dtheta, &
         lbound(cache), ubound(cache), unit = 'cm rad^-1', &
         comment = 'Jacobian element (Z, theta) at ' // trim(adjustl(comment)))
    call h5_close(h5id_root)
  end subroutine coord_cache_write

  subroutine coord_cache_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(coord_cache_t), dimension(:, :), intent(inout) :: cache
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

  subroutine coord_cache_ext_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    type(coord_cache_ext_t), dimension(:, :), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call coord_cache_write(cache, file, dataset, comment)
    call h5_open_rw(file, h5id_root)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/rad', cache%rad, &
         lbound(cache), ubound(cache), unit = 'cm', &
         comment = 'minor radius (along theta = 0) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/q', cache%q, &
         lbound(cache), ubound(cache), unit = '1', &
         comment = 'safety factor q at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dq_dpsi', cache%dq_dpsi, &
         lbound(cache), ubound(cache), unit = 'Mx^-1', &
         comment = 'q''(psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dR_dpsi', cache%dR_dpsi, &
         lbound(cache), ubound(cache), unit = 'cm Mx^-1', &
         comment = 'Jacobian element (R, psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dZ_dpsi', cache%dZ_dpsi, &
         lbound(cache), ubound(cache), unit = 'cm Mx^-1', &
         comment = 'Jacobian element (Z, psi) at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2R_dpsi_dtheta', cache%d2R_dpsi_dtheta, &
         lbound(cache), ubound(cache), unit = 'cm Mx^-1 rad^-1', &
         comment = 'Hessian element (psi, theta) of R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/d2Z_dpsi_dtheta', cache%d2Z_dpsi_dtheta, &
         lbound(cache), ubound(cache), unit = 'cm Mx^-1 rad^-1', &
         comment = 'Hessian element (psi, theta) of Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache%B0_phi, &
         lbound(cache), ubound(cache), unit = 'G', &
         comment = 'physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache%dB0R_dR, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'R derivative of B0_R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache%dB0R_dZ, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_R at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache%dB0phi_dR, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'R derivative of B0_phi at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache%dB0phi_dZ, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_phi at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache%dB0Z_dR, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'R derivative of B0_Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache%dB0Z_dZ, &
         lbound(cache), ubound(cache), unit = 'G cm^-1', &
         comment = 'Z derivative of B0_Z at ' // trim(adjustl(comment)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_2', cache%B0_2, &
         lbound(cache), ubound(cache), unit = 'G^2', &
         comment = 'square of equilibrium magnetic field  at ' // trim(adjustl(comment)))
    call h5_close(h5id_root)
  end subroutine coord_cache_ext_write

  subroutine coord_cache_ext_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(coord_cache_ext_t), dimension(:, :), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call coord_cache_read(cache, file, dataset)
    call h5_open(file, h5id_root)
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

  subroutine field_cache_write(cache, file, dataset, comment)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    class(field_cache_t), dimension(..), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    select rank (cache)
    rank (1)
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', cache(:)%psi, &
            lbound(cache), ubound(cache), unit = 'Mx', &
            comment = 'poloidal flux at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache(:)%B0(1), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache(:)%B0(3), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache(:)%B0(2), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_R', cache(:)%j0(1), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_Z', cache(:)%j0(3), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'Z component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_phi', cache(:)%j0(2), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0', cache(:)%Bmod, &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0_dR', cache(:)%dBmod_dR, &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0_dZ', cache(:)%dBmod_dZ, &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache(:)%dB0_dR(1), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache(:)%dB0_dZ(1), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache(:)%dB0_dR(2), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache(:)%dB0_dZ(2), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache(:)%dB0_dR(3), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache(:)%dB0_dZ(3), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0R_dR', cache(:)%dj0_dR(1), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0R_dZ', cache(:)%dj0_dZ(1), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dR', cache(:)%dj0_dR(2), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dZ', cache(:)%dj0_dZ(2), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dR', cache(:)%dj0_dR(3), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of Z component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dZ', cache(:)%dj0_dZ(3), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of Z component of equilibrium current density at ' // trim(adjustl(comment)))
    rank (2)
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/psi', cache(:, :)%psi, &
            lbound(cache), ubound(cache), unit = 'Mx', &
            comment = 'poloidal flux at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache(:, :)%B0(1), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache(:, :)%B0(3), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache(:, :)%B0(2), &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_R', cache(:, :)%j0(1), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_Z', cache(:, :)%j0(3), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'Z component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/j0_phi', cache(:, :)%j0(2), &
            lbound(cache), ubound(cache), unit = 'statA', &
            comment = 'physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0', cache(:, :)%Bmod, &
            lbound(cache), ubound(cache), unit = 'G', &
            comment = 'magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0_dR', cache(:, :)%dBmod_dR, &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0_dZ', cache(:, :)%dBmod_dZ, &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of magnitude of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache(:, :)%dB0_dR(1), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache(:, :)%dB0_dZ(1), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of R component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache(:, :)%dB0_dR(2), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache(:, :)%dB0_dZ(2), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of physical phi component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache(:, :)%dB0_dR(3), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'R derivative of Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache(:, :)%dB0_dZ(3), &
            lbound(cache), ubound(cache), unit = 'G cm^-1', &
            comment = 'Z derivative of Z component of equilibrium magnetic field at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0R_dR', cache(:, :)%dj0_dR(1), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0R_dZ', cache(:, :)%dj0_dZ(1), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of R component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dR', cache(:, :)%dj0_dR(2), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dZ', cache(:, :)%dj0_dZ(2), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of physical phi component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dR', cache(:, :)%dj0_dR(3), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'R derivative of Z component of equilibrium current density at ' // trim(adjustl(comment)))
       call h5_add(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dZ', cache(:, :)%dj0_dZ(3), &
            lbound(cache), ubound(cache), unit = 'statA cm^-1', &
            comment = 'Z derivative of Z component of equilibrium current density at ' // trim(adjustl(comment)))
    rank default
       error stop 'field_cache_write: rank-1 or rank-2 array expected for argument ''cache'''
    end select
    call h5_close(h5id_root)
  end subroutine field_cache_write

  subroutine field_cache_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    class(field_cache_t), dimension(..), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    select rank (cache)
    rank (1)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi', cache(:)%psi)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache(:)%B0(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache(:)%B0(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache(:)%B0(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_R', cache(:)%j0(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_Z', cache(:)%j0(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_phi', cache(:)%j0(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0', cache(:)%Bmod)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0_dR', cache(:)%dBmod_dR)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0_dZ', cache(:)%dBmod_dZ)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache(:)%dB0_dR(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache(:)%dB0_dZ(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache(:)%dB0_dR(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache(:)%dB0_dZ(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache(:)%dB0_dR(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache(:)%dB0_dZ(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0R_dR', cache(:)%dj0_dR(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0R_dZ', cache(:)%dj0_dZ(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dR', cache(:)%dj0_dR(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dZ', cache(:)%dj0_dZ(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dR', cache(:)%dj0_dR(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dZ', cache(:)%dj0_dZ(3))
    rank (2)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi', cache(:, :)%psi)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_R', cache(:, :)%B0(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_Z', cache(:, :)%B0(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_phi', cache(:, :)%B0(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_R', cache(:, :)%j0(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_Z', cache(:, :)%j0(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/j0_phi', cache(:, :)%j0(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0', cache(:, :)%Bmod)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0_dR', cache(:, :)%dBmod_dR)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0_dZ', cache(:, :)%dBmod_dZ)
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dR', cache(:, :)%dB0_dR(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0R_dZ', cache(:, :)%dB0_dZ(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dR', cache(:, :)%dB0_dR(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0phi_dZ', cache(:, :)%dB0_dZ(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dR', cache(:, :)%dB0_dR(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dB0Z_dZ', cache(:, :)%dB0_dZ(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0R_dR', cache(:, :)%dj0_dR(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0R_dZ', cache(:, :)%dj0_dZ(1))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dR', cache(:, :)%dj0_dR(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0phi_dZ', cache(:, :)%dj0_dZ(2))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dR', cache(:, :)%dj0_dR(3))
       call h5_get(h5id_root, trim(adjustl(dataset)) // '/dj0Z_dZ', cache(:, :)%dj0_dZ(3))
    rank default
       error stop 'field_cache_read: rank-1 or rank-2 array expected for argument ''cache'''
    end select
    call h5_close(h5id_root)
  end subroutine field_cache_read

  subroutine cache_init(cache, nrad, npol)
    type(cache_t), intent(inout) :: cache
    integer, intent(in) :: nrad, npol

    call cache_deinit(cache)
    cache%nrad = nrad
    cache%npol = npol
    allocate(cache%sample_polmodes(npol, nrad), cache%sample_polmodes_half(npol, nrad))
    allocate(cache%edge_fields(mesh%GL_order, mesh%nedge), cache%area_fields(mesh%GL2_order, mesh%ntri))
    allocate(cache%mid_fields(mesh%nedge), cache%cntr_fields(mesh%ntri), cache%B0_flux(mesh%nedge))
  end subroutine cache_init

  subroutine cache_deinit(cache)
    type(cache_t), intent(inout) :: cache

    cache%nrad = 0
    cache%npol = 0
    if (allocated(cache%sample_polmodes)) deallocate(cache%sample_polmodes)
    if (allocated(cache%sample_polmodes_half)) deallocate(cache%sample_polmodes_half)
    if (allocated(cache%edge_fields)) deallocate(cache%edge_fields)
    if (allocated(cache%area_fields)) deallocate(cache%area_fields)
    if (allocated(cache%mid_fields)) deallocate(cache%mid_fields)
    if (allocated(cache%cntr_fields)) deallocate(cache%cntr_fields)
    if (allocated(cache%B0_flux)) deallocate(cache%B0_flux)
  end subroutine cache_deinit

  subroutine cache_write(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    type(cache_t), intent(in) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call coord_cache_write(cache%sample_polmodes_half, file, &
         trim(adjustl(dataset)) // '/sample_polmodes_half', &
         'poloidal mode sampling points between flux surfaces')
    call coord_cache_write(cache%sample_polmodes, file, &
         trim(adjustl(dataset)) // '/sample_polmodes', &
         'poloidal mode sampling points on flux surfaces')
    call field_cache_write(cache%edge_fields, file, &
         trim(adjustl(dataset)) // '/edge_fields', &
         'GL quadrature points on triangle edges')
    call field_cache_write(cache%area_fields, file, &
         trim(adjustl(dataset)) // '/area_fields', &
         'GL quadrature points on triangle areas')
    call field_cache_write(cache%mid_fields, file, &
         trim(adjustl(dataset)) // '/mid_fields', &
         'triangle edge midpoints')
    call field_cache_write(cache%cntr_fields, file, &
         trim(adjustl(dataset)) // '/cntr_fields', &
         'weighted triangle centroids')
    call h5_open_rw(file, h5id_root)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/nrad', cache%nrad, &
         comment = 'number of radial divisions for poloidal mode sampling points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/npol', cache%npol, &
         comment = 'number of poloidal divisions for poloidal mode sampling points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/B0_flux', cache%B0_flux, &
         lbound(cache%B0_flux), ubound(cache%B0_flux), &
         comment = 'equilibrium magnetic flux through triangle edge', unit = 'G cm^2')
    call h5_close(h5id_root)
  end subroutine cache_write

  subroutine cache_read(cache, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(cache_t), intent(inout) :: cache
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer :: nrad, npol

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/nrad', nrad)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/npol', npol)
    call cache_init(cache, nrad, npol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/B0_flux', cache%B0_flux)
    call h5_close(h5id_root)
    call coord_cache_read(cache%sample_polmodes_half, file, &
         trim(adjustl(dataset)) // '/sample_polmodes_half')
    call coord_cache_read(cache%sample_polmodes, file, &
         trim(adjustl(dataset)) // '/sample_polmodes')
    call field_cache_read(cache%edge_fields, file, trim(adjustl(dataset)) // '/edge_fields')
    call field_cache_read(cache%area_fields, file, trim(adjustl(dataset)) // '/area_fields')
    call field_cache_read(cache%mid_fields, file, trim(adjustl(dataset)) // '/mid_fields')
    call field_cache_read(cache%cntr_fields, file, trim(adjustl(dataset)) // '/cntr_fields')
  end subroutine cache_read

  subroutine generate_mesh
    use mephit_conf, only: conf

    if (conf%kilca_scale_factor /= 0) then
       mesh%n = conf%n * conf%kilca_scale_factor
    else
       mesh%n = conf%n
    end if
    call create_mesh_points
    call compare_gpec_coordinates
    call connect_mesh_points
    call write_FreeFem_mesh
    call cache_init(cache, mesh%nflux, 512)
    call compute_sample_polmodes(cache%sample_polmodes_half, cache%npol, .true.)
    call compute_sample_polmodes(cache%sample_polmodes, cache%npol, .false.)
    call compute_gpec_jacfac
    call cache_equilibrium_field
    call init_flux_variables
    call check_resonance_positions
    call compute_shielding_auxiliaries
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

  subroutine compute_resonance_positions(psi_sample, q_sample, psi2rho_norm)
    use mephit_conf, only: conf, logger
    use mephit_util, only: interp1d
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
    type(interp1d) :: psi_sample_interpolator

    if (size(psi_sample) /= size(q_sample)) then
       call logger%msg_arg_size('refine_resonant_surfaces', 'size(psi_sample)', &
            'size(q_sample)', size(psi_sample), size(q_sample))
       if (logger%err) call logger%write_msg
       error stop
    end if
    psi_min = minval(psi_sample)
    psi_max = maxval(psi_sample)
    call psi_sample_interpolator%init(4, psi_sample)
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
    logger%msg = 'resonance positions:'
    if (logger%debug) call logger%write_msg
    do m = mesh%m_res_min, mesh%m_res_max
       mesh%psi_res(m) = zeroin(psi_min, psi_max, q_interp_resonant, 1d-9)
       mesh%rad_norm_res(m) = psi2rho_norm(mesh%psi_res(m))
       write (logger%msg, '("m = ", i2, ", psi_m = ", es24.16e3, ", rho_m = ", f19.16)') &
            m, mesh%psi_res(m), mesh%rad_norm_res(m)
       if (logger%debug) call logger%write_msg
    end do

  contains
    function q_interp_resonant(psi)
      real(dp), intent(in) :: psi
      real(dp) :: q_interp_resonant
      q_interp_resonant = psi_sample_interpolator%eval(abs(q_sample), psi) - dble(m) / dble(mesh%n)
    end function q_interp_resonant
  end subroutine compute_resonance_positions

  subroutine refine_eqd_partition(coarse_sep, nref, deletions, refinement, resonances, diverging_q, partition, ref_ind)
    use mephit_conf, only: logger
    use mephit_util, only: linspace
    real(dp), intent(in) :: coarse_sep
    integer, intent(in) :: nref
    integer, dimension(:), intent(in) :: deletions
    real(dp), dimension(:), intent(in) :: refinement, resonances
    logical, intent(in) :: diverging_q
    real(dp), dimension(:), allocatable, intent(out) :: partition
    integer, dimension(:), intent(out) :: ref_ind
    integer :: kref, k, inter(nref + 1)
    integer, dimension(nref) :: additions, add_lo, add_hi, fine_lo, fine_hi
    real(dp), dimension(nref) :: fine_sep
    real(dp), dimension(:, :), allocatable :: geom_ser, fine_pos_lo, fine_pos_hi

    if (allocated(partition)) deallocate(partition)
    if (nref < 1) then
       mesh%nflux = ceiling(1d0 / coarse_sep)
       allocate(partition(0:mesh%nflux))
       partition(:) = linspace(0d0, 1d0, mesh%nflux + 1, 0, 0)
       return
    end if
    if (nref /= size(deletions)) then
       call logger%msg_arg_size('refine_eqd_partition', 'nref', 'size(deletions)', nref, &
            size(deletions))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (nref /= size(refinement)) then
       call logger%msg_arg_size('refine_eqd_partition', 'nref', 'size(refinement)', nref, &
            size(refinement))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (nref /= size(resonances)) then
       call logger%msg_arg_size('refine_eqd_partition', 'nref', 'size(resonances)', nref, &
            size(resonances))
       if (logger%err) call logger%write_msg
       error stop
    end if
    if (nref /= size(ref_ind)) then
       call logger%msg_arg_size('refine_eqd_partition', 'nref', 'size(ref_ind)', nref, &
            size(ref_ind))
       if (logger%err) call logger%write_msg
       error stop
    end if
    ! construct geometric series
    if (any(refinement >= dble(deletions + 3) / dble(deletions + 1))) then
       logger%msg = 'Requested refinement factor is too high'
       if (logger%err) call logger%write_msg
       error stop
    end if
    additions = ceiling((log(refinement + 1d0) - &
         log(3d0 + dble(2 * deletions) * (1d0 - refinement) - refinement)) / log(refinement)) - 1
    fine_sep = coarse_sep * refinement ** (-additions - 1)
    allocate(geom_ser(maxval(additions) + 1, nref))
    geom_ser(:, :) = reshape([(((refinement(kref) ** k - refinement(kref)) / (refinement(kref) - 1d0), &
         k = 1, maxval(additions) + 1), kref = 1, nref)], [maxval(additions) + 1, nref])
    ! determine positions and numbers of additional points in refined regions
    allocate(fine_pos_lo(maxval(additions) + 1, nref), fine_pos_hi(maxval(additions) + 1, nref))
    fine_pos_lo(:, :) = reshape([((resonances(kref) - fine_sep(kref) * (0.5d0 + geom_ser(k, kref)), &
         k = 1, maxval(additions) + 1), kref = 1, nref)], [maxval(additions) + 1, nref])
    fine_pos_hi(:, :) = reshape([((resonances(kref) + fine_sep(kref) * (0.5d0 + geom_ser(k, kref)), &
         k = 1, maxval(additions) + 1), kref = 1, nref)], [maxval(additions) + 1, nref])
    ! reduce number of additional points if refined regions overlap
    add_lo = additions
    add_hi = additions
    do while (fine_pos_lo(add_lo(1) + 1, 1) < 0d0)
       add_lo(1) = add_lo(1) - 1
    end do
    do kref = 2, nref
       do while (fine_pos_lo(add_lo(kref) + 1, kref) < fine_pos_hi(add_hi(kref-1) + 1, kref - 1))
          if (add_lo(kref) <= add_hi(kref-1)) then
             add_hi(kref-1) = add_hi(kref-1) - 1
          else
             add_lo(kref) = add_lo(kref) - 1
          end if
       end do
    end do
    if (diverging_q) then
       ! continue with fine separation outside last refined region
       add_hi(nref) = 0
    else
       do while (1d0 < fine_pos_hi(add_hi(nref) + 1, nref))
          add_hi(nref) = add_hi(nref) - 1
       end do
    end if
    ! compute number of intervals between refined regions
    inter(1) = ceiling(fine_pos_lo(add_lo(1) + 1, 1) / (fine_sep(1) * refinement(1) ** (add_lo(1) + 1)))
    do kref = 2, nref
       inter(kref) = ceiling((fine_pos_lo(add_lo(kref) + 1, kref) - fine_pos_hi(add_hi(kref-1) + 1, kref-1)) / &
            max(fine_sep(kref) * refinement(kref) ** (add_lo(kref) + 1), &
            fine_sep(kref-1) * refinement(kref-1) ** (add_hi(kref-1) + 1)))
    end do
    if (diverging_q) then
       ! continue with fine separation outside last refined region
       inter(nref + 1) = ceiling((1d0 - fine_pos_hi(add_hi(nref) + 1, nref)) / fine_sep(nref))
    else
       inter(nref + 1) = ceiling((1d0 - fine_pos_hi(add_hi(nref) + 1, nref)) / &
            (fine_sep(nref) * refinement(nref) ** (add_hi(nref) + 1)))
    end if
    ! compute upper and lower array indices of refined regions
    fine_lo(1) = inter(1)
    fine_hi(1) = fine_lo(1) + add_lo(1) + add_hi(1) + 1
    do kref = 2, nref
       fine_lo(kref) = fine_hi(kref-1) + inter(kref)
       fine_hi(kref) = fine_lo(kref) + add_lo(kref) + add_hi(kref) + 1
    end do
    ref_ind = fine_hi - add_hi
    ! compute refined regions around resonant flux surfaces
    mesh%nflux = sum(add_lo) + sum(add_hi) + sum(inter) + nref
    allocate(partition(0:mesh%nflux))
    partition(:) = 0d0
    do kref = 1, nref
       partition(fine_lo(kref) + add_lo(kref):fine_lo(kref):-1) = resonances(kref) - &
            fine_sep(kref) * (0.5d0 + geom_ser(:add_lo(kref) + 1, kref))
       partition(fine_hi(kref) - add_hi(kref):fine_hi(kref)) = resonances(kref) + &
            fine_sep(kref) * (0.5d0 + geom_ser(:add_hi(kref) + 1, kref))
    end do
    ! compute equidistant positions between refined regions
    partition(:fine_lo(1)) = linspace(0d0, partition(fine_lo(1)), inter(1) + 1, 0, 0)
    do kref = 2, nref
       partition(fine_hi(kref-1):fine_lo(kref)) = &
            linspace(partition(fine_hi(kref-1)), partition(fine_lo(kref)), inter(kref) + 1, 0, 0)
    end do
    partition(fine_hi(nref):) = linspace(partition(fine_hi(nref)), 1d0, inter(nref+1) + 1, 0, 0)
  end subroutine refine_eqd_partition

  subroutine refine_resonant_surfaces(coarse_sep, rho_norm_ref)
    use mephit_conf, only: conf, conf_arr, logger
    real(dp), intent(in) :: coarse_sep
    real(dp), dimension(:), allocatable, intent(out) :: rho_norm_ref
    logical :: diverging_q
    integer :: m, m_dense, kref
    integer, dimension(:), allocatable :: ref_ind
    logical, dimension(:), allocatable :: mask

    if (allocated(mesh%refinement)) deallocate(mesh%refinement)
    if (allocated(mesh%deletions)) deallocate(mesh%deletions)
    allocate(mesh%refinement(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%deletions(mesh%m_res_min:mesh%m_res_max))
    mesh%refinement(:) = conf_arr%refinement
    mesh%deletions(:) = conf_arr%deletions
    allocate(mask(mesh%m_res_min:mesh%m_res_max))
    mask(:) = 1d0 < mesh%refinement .and. mesh%refinement < 3d0
    if (conf%kilca_scale_factor /= 0) then
       diverging_q = .false.
    else
       diverging_q = .true.
       ! heuristic: if distance between resonances is less than coarse grid separation,
       ! take inner resonance as last to be refined; outside, only the fine separation is used
       m_dense = mesh%m_res_min + 1
       do while (m_dense <= mesh%m_res_max)
          if (mesh%rad_norm_res(m_dense) - mesh%rad_norm_res(m_dense - 1) < coarse_sep) exit
          m_dense = m_dense + 1
       end do
       mask(m_dense:mesh%m_res_max) = .false.
    end if
    allocate(ref_ind(count(mask)))
    call refine_eqd_partition(coarse_sep, count(mask), pack(mesh%deletions, mask), &
         pack(mesh%refinement, mask), pack(mesh%rad_norm_res, mask), diverging_q, rho_norm_ref, ref_ind)
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
    use mephit_util, only: interp_psi_pol, pi, pos_angle
    use magdata_in_symfluxcoor_mod, only: nlabel, rbeg, psisurf, psipol_max, qsaf, &
         rsmall, circumf, raxis, zaxis, load_magdata_in_symfluxcoord
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
    ! calculates points on a fine grid in the core region by integrating along field lines
    call preload_for_SYNCH
    ! loads points that are calculated in preload_for_SYNCH into module variables
    ! and spline interpolates them for use with magdata_in_symfluxcoord_ext
    call load_magdata_in_symfluxcoord
    mesh%R_O = raxis
    mesh%Z_O = zaxis
    mesh%R_X = X_point(1)
    mesh%Z_X = X_point(2)
    ! field_line_integration_for_SYNCH subtracts psi_axis from psisurf and
    ! load_magdata_in_symfluxcoord_ext divides by psipol_max
    psi_axis = interp_psi_pol(raxis, zaxis)
    call psi_fine_interpolator%init(4, psisurf(1:) * psipol_max + psi_axis)
    ! interpolate between psi and rho
    rad_max = rbeg(nlabel)
    allocate(rho_norm_eqd(nlabel))
    rho_norm_eqd(:) = rbeg / rad_max

    call compute_resonance_positions(psisurf(1:) * psipol_max + psi_axis, qsaf, psi2rho_norm)
    call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)
    call refine_resonant_surfaces(conf%max_Delta_rad / rad_max, rho_norm_ref)
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

    fs%area(:) = [(pi * psi_fine_interpolator%eval(rsmall, fs%psi(kf)) ** 2, &
         kf = 0, mesh%nflux)]
    fs%perimeter(:) = [(psi_fine_interpolator%eval(circumf, fs%psi(kf)), &
         kf = 0, mesh%nflux)]
    fs_half%area(:) = [(pi * psi_fine_interpolator%eval(rsmall, fs_half%psi(kf)) ** 2, &
         kf = 1, mesh%nflux)]
    fs_half%perimeter(:) = [(psi_fine_interpolator%eval(circumf, fs_half%psi(kf)), &
         kf = 1, mesh%nflux)]
    allocate(opt_pol_edge_len(mesh%nflux + 1))
    ! cache averaged radius at half-grid steps
    opt_pol_edge_len(:mesh%nflux) = sqrt(fs_half%area / pi)
    ! extrapolate linearly
    opt_pol_edge_len(mesh%nflux + 1) = 2d0 * opt_pol_edge_len(mesh%nflux) - opt_pol_edge_len(mesh%nflux - 1)
    ! compute successive difference
    opt_pol_edge_len(:mesh%nflux) = opt_pol_edge_len(2:) - opt_pol_edge_len(:mesh%nflux)
    ! averaged radius corresponds to altitude - factor for edge length of equilateral triangle
    opt_pol_edge_len(:mesh%nflux) = 2d0 / sqrt(3d0) * opt_pol_edge_len(:mesh%nflux)
    allocate(mesh%kp_max(mesh%nflux))
    allocate(mesh%kt_max(mesh%nflux))
    allocate(mesh%kp_low(mesh%nflux))
    allocate(mesh%kt_low(mesh%nflux))
    ! round to even numbers
    mesh%kp_max(:) = 2 * nint(0.5d0 * fs%perimeter(1:) / opt_pol_edge_len(:mesh%nflux))
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
    function psi2rho_norm(psi) result(rho_norm)
      real(dp), intent(in) :: psi
      real(dp) :: rho_norm
      rho_norm = psi_fine_interpolator%eval(rho_norm_eqd, psi)
    end function psi2rho_norm
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
    use mephit_util, only: pi, gauss_legendre_unit_interval
    integer :: kf, kp, kp_lo, kp_hi, kt, ktri, ktri_adj, kedge, nodes(4), k
    real(dp) :: mat(3, 3), points(mesh%GL_order), points2(3, mesh%GL2_order)

    allocate(mesh%tri_node(3, mesh%ntri))
    allocate(mesh%tri_node_F(mesh%ntri))
    allocate(mesh%tri_theta_extent(2, mesh%ntri))
    allocate(mesh%orient(mesh%ntri))
    allocate(mesh%edge_node(2, mesh%nedge))
    allocate(mesh%edge_tri(2, mesh%nedge))
    allocate(mesh%tri_edge(3, mesh%ntri))
    mesh%tri_node = 0
    mesh%tri_node_F = 0
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
    ! define triangles on outer flux surfaces
    do kf = 2, mesh%nflux
       kp_lo = 1
       kp_hi = 1
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
          else
             ! node indices for triangle with poloidal edge on inner flux surface
             mesh%tri_node(:, ktri) = nodes([2, 4, 1])
             mesh%tri_node_F(ktri) = 1
             ! triangle index for poloidal edge
             kedge = mesh%kp_low(kf - 1) + kp_lo - 1
             mesh%edge_tri(2, kedge) = ktri
             mesh%tri_edge(1, ktri) = kedge
             kp_lo = mod(kp_lo, mesh%kp_max(kf - 1)) + 1
          end if
          ! triangle indices for radial edge
          ktri_adj = mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf)) + 1
          kedge = mesh%npoint + mesh%kt_low(kf) + mod(kt, mesh%kt_max(kf))
          mesh%edge_tri(:, kedge) = [ktri, ktri_adj]
          mesh%tri_edge(2, ktri) = kedge
          mesh%tri_edge(3, ktri_adj) = kedge
          ! cache poloidal extent of triangles for point_location_check
          if (all(mesh%node_theta_geom(mesh%tri_node(:, ktri)) > epsilon(1d0)) .or. kt < (mesh%kt_max(kf) - kt)) then
             mesh%tri_theta_extent(1, ktri) = minval(mesh%node_theta_geom(mesh%tri_node(:, ktri)))
             mesh%tri_theta_extent(2, ktri) = maxval(mesh%node_theta_geom(mesh%tri_node(:, ktri)))
          else
             mesh%tri_theta_extent(1, ktri) = minval(upper_branch(mesh%node_theta_geom(mesh%tri_node(:, ktri))))
             mesh%tri_theta_extent(2, ktri) = maxval(upper_branch(mesh%node_theta_geom(mesh%tri_node(:, ktri))))
          end if
       end do
    end do
    ! set theta = 2 pi exactly
    where (mesh%tri_theta_extent(2, :) < epsilon(1d0))
       mesh%tri_theta_extent(2, :) = 2d0 * pi
    end where
    ! nodes for radial edges
    mesh%edge_node(:, mesh%npoint:) = mesh%tri_node([3, 1], :)
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
          mesh%GL_R(k, kedge) = mesh%node_R(mesh%edge_node(1, kedge)) * points(k) + &
               mesh%node_R(mesh%edge_node(2, kedge)) * points(mesh%GL_order - k + 1)
          mesh%GL_Z(k, kedge) = mesh%node_Z(mesh%edge_node(1, kedge)) * points(k) + &
               mesh%node_Z(mesh%edge_node(2, kedge)) * points(mesh%GL_order - k + 1)
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

  subroutine compute_shielding_auxiliaries
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use points_2d, only: theta_geom2theta_flux
    use mephit_util, only: pi, clight, interp_psi_pol, binsearch, pos_angle
    integer :: m, kf, kt, kedge, kgl, kp, k
    real(dp) :: s, psi, dum, R, Z, theta, dq_dpsi
    real(dp), dimension(:), allocatable :: theta_geom, theta_flux

    allocate(mesh%shielding_kt_low(mesh%m_res_min:mesh%m_res_max))
    mesh%shielding_kt_low(mesh%m_res_min) = 0
    do m = mesh%m_res_min + 1, mesh%m_res_max
       mesh%shielding_kt_low(m) = mesh%shielding_kt_low(m-1) + mesh%kt_max(mesh%res_ind(m-1))
    end do
    allocate(mesh%shielding_L1_kp(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    allocate(mesh%shielding_L1_weight(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    allocate(mesh%shielding_L1_R(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    allocate(mesh%shielding_L1_Z(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    allocate(theta_geom(maxval(mesh%kt_max(mesh%res_ind)) + 1), &
         theta_flux(maxval(mesh%kt_max(mesh%res_ind)) + 1))
    do m = mesh%m_res_min, mesh%m_res_max
       kf = mesh%res_ind(m)
       do kgl = 1, mesh%GL_order
          theta_geom(:) = 0d0
          do kt = 2, mesh%kt_max(kf)
             kedge = mesh%npoint + mesh%kt_low(kf) + kt - 1
             theta_geom(kt) = pos_angle(atan2(mesh%GL_Z(kgl, kedge) - mesh%Z_O, &
                  mesh%GL_R(kgl, kedge) - mesh%R_O))
          end do
          theta_geom(mesh%kt_max(kf) + 1) = 2d0 * pi
          kedge = mesh%npoint + mesh%kt_low(kf)
          psi = interp_psi_pol(mesh%GL_R(kgl, kedge), mesh%GL_Z(kgl, kedge)) - fs%psi(0)
          ! inp_label = 2 to use poloidal flux psi
          call theta_geom2theta_flux(2, s, psi, theta_geom(:mesh%kt_max(kf) + 1), &
               theta_flux(:mesh%kt_max(kf) + 1))
          ! loop over inner and outer flux surface
          do k = -1, 0
             psi = fs%psi(kf+k) - fs%psi(0)
             theta_geom(:mesh%kp_max(kf+k)) = mesh%node_theta_geom((mesh%kp_low(kf+k) + 1):&
                  (mesh%kp_low(kf+k) + mesh%kp_max(kf+k)))
             theta_geom(mesh%kp_max(kf+k) + 1) = 2d0 * pi
             do kt = 1, mesh%kt_max(kf)
                ! inp_label = 2 to use poloidal flux psi
                call magdata_in_symfluxcoord_ext(2, dum, psi, theta_flux(kt), &
                     dum, dum, dum, dum, dum, R, dum, dum, Z, dum, dum)
                mesh%shielding_L1_R(k, kgl, mesh%shielding_kt_low(m) + kt) = R
                mesh%shielding_L1_Z(k, kgl, mesh%shielding_kt_low(m) + kt) = Z
                theta = pos_angle(atan2(Z - mesh%Z_O, R - mesh%R_O))
                call binsearch(theta_geom(:mesh%kp_max(kf+k) + 1), 0, theta, kp)
                mesh%shielding_L1_kp(k, kgl, mesh%shielding_kt_low(m) + kt) = kp
                mesh%shielding_L1_weight(k, kgl, mesh%shielding_kt_low(m) + kt) = &
                     (theta - theta_geom(kp)) / (theta_geom(kp + 1) - theta_geom(kp))
             end do
          end do
       end do
    end do
    deallocate(theta_geom, theta_flux)
    allocate(mesh%shielding_coeff(mesh%m_res_min:mesh%m_res_max))
    do m = mesh%m_res_min, mesh%m_res_max
       kf = mesh%res_ind(m)
       dq_dpsi = psi_interpolator%eval(equil%qpsi, fs_half%psi(kf), .true.)
       mesh%shielding_coeff(m) = clight * mesh%n / (4d0 * pi * mesh%R_O) * &
            abs(dq_dpsi / (fs_half%q(kf) * fs_half%dp_dpsi(kf))) / &
            (mesh%n * abs(fs%q(kf) - fs%q(kf-1)))
    end do
  end subroutine compute_shielding_auxiliaries

  subroutine compute_gpec_jacfac
    use mephit_util, only: imun
    integer, parameter :: m_max = 16
    integer :: kf, kpol, m
    complex(dp) :: fourier_basis(-m_max:m_max)

    allocate(mesh%gpec_jacfac(-m_max:m_max, mesh%nflux))
    mesh%gpec_jacfac(:, :) = (0d0, 0d0)
    do kf = 1, mesh%nflux
       do kpol = 1, cache%npol
          associate (s => cache%sample_polmodes_half(kpol, kf))
            fourier_basis = [(exp(-imun * m * s%theta), m = -m_max, m_max)]
            mesh%gpec_jacfac(:, kf) = mesh%gpec_jacfac(:, kf) + s%sqrt_g * &
                 s%R * hypot(s%B0_Z, -s%B0_R) * fourier_basis
          end associate
       end do
    end do
    mesh%gpec_jacfac(:, :) = mesh%gpec_jacfac / cache%npol
  end subroutine compute_gpec_jacfac

  !> Compute coarse grid for poloidal mode sampling points
  subroutine compute_sample_polmodes(sample, npol, half_grid)
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_conf, only: conf
    use mephit_util, only: pi
    type(coord_cache_t), dimension(:, :), allocatable, intent(inout) :: sample
    integer, intent(in) :: npol
    logical, intent(in) :: half_grid
    integer :: kf, kpol
    real(dp) :: dum, q
    real(dp), allocatable, dimension(:) :: psi, rad

    if (allocated(sample)) deallocate(sample)
    allocate(sample(npol, mesh%nflux))
    if (half_grid) then
       allocate(psi, source = fs_half%psi)
       allocate(rad, source = fs_half%rad)
    else
       allocate(psi, source = fs%psi)
       allocate(rad, source = fs%rad)
    end if
    do kf = 1, mesh%nflux
       do kpol = 1, npol
          associate (s => sample(kpol, kf))
            s%psi = psi(kf)
            s%theta = 2d0 * pi * dble(kpol - 1) / dble(npol)
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
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
            call field(s%R, 0d0, s%Z, s%B0_R, dum, s%B0_Z, &
                 dum, dum, dum, dum, dum, dum, dum, dum, dum)
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
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

  !> Compute fine grid for parallel current sampling points.
  subroutine compute_sample_Ipar(sample_Ipar, m)
    use magdata_in_symfluxcoor_mod, only: psipol_max, magdata_in_symfluxcoord_ext
    use mephit_conf, only: conf
    use mephit_util, only: pi, linspace, interp_psi_pol
    type(coord_cache_ext_t), dimension(:, :), intent(inout) :: sample_Ipar
    integer, intent(in) :: m
    integer :: nrad, npol, krad, kpol
    real(dp) :: rad_min, rad_max, rad_eqd(equil%nw), drad_dpsi, dum
    real(dp), dimension(:), allocatable :: rad, psi, theta

    npol = size(sample_Ipar, 1)
    nrad = size(sample_Ipar, 2)
    rad_eqd(:) = linspace(fs%rad(0), fs%rad(mesh%nflux), equil%nw, 0, 0)
    rad_min = mesh%rad_norm_res(m) * fs%rad(mesh%nflux) - 2d0 * conf%max_Delta_rad
    rad_max = mesh%rad_norm_res(m) * fs%rad(mesh%nflux) + 2d0 * conf%max_Delta_rad
    if (rad_max > fs%rad(mesh%nflux)) then
       rad_max = fs%rad(mesh%nflux)
       rad_min = (2d0 * mesh%rad_norm_res(m) - 1d0) * fs%rad(mesh%nflux)
    end if
    allocate(rad(nrad), psi(nrad), theta(npol))
    rad(:) = linspace(rad_min, rad_max, nrad, 1, 1)
    ! TODO: replace by dedicated interpolation function
    if (conf%kilca_scale_factor /= 0) then
       psi(:) = [(interp_psi_pol(mesh%R_O, mesh%Z_O + rad(krad)), krad = 1, nrad)]
    else
       psi(:) = [(interp_psi_pol(mesh%R_O + rad(krad), mesh%Z_O), krad = 1, nrad)]
    end if
    theta(:) = 2d0 * pi * [(dble(kpol - 1), kpol = 1, npol)] / dble(npol)
    do krad = 1, nrad
       do kpol = 1, npol
          associate (s => sample_Ipar(kpol, krad))
            s%rad = rad(krad)
            s%psi = psi(krad)
            s%theta = theta(kpol)
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%R = mesh%R_O + rad(krad) * cos(theta(kpol))
               s%Z = mesh%Z_O + rad(krad) * sin(theta(kpol))
               s%dR_dtheta = -rad(krad) * sin(theta(kpol))
               s%dZ_dtheta =  rad(krad) * cos(theta(kpol))
               s%dq_dpsi = psi_interpolator%eval(equil%qpsi, psi(krad))
               drad_dpsi = 1d0 / psi_interpolator%eval(rad_eqd, psi(krad))
               s%dR_dpsi = drad_dpsi * cos(theta(kpol))
               s%dZ_dpsi = drad_dpsi * sin(theta(kpol))
               s%d2R_dpsi_dtheta = -drad_dpsi * sin(theta(kpol))
               s%d2Z_dpsi_dtheta =  drad_dpsi * cos(theta(kpol))
            else
               ! psi is shifted by -psi_axis in magdata_in_symfluxcoor_mod
               call magdata_in_symfluxcoord_ext(2, dum, psi(krad) - fs%psi(0), theta(kpol), &
                    s%q, s%dq_dpsi, s%sqrt_g, dum, dum, &
                    s%R, s%dR_dpsi, s%dR_dtheta, &
                    s%Z, s%dZ_dpsi, s%dZ_dtheta, &
                    s%d2R_dpsi_dtheta, s%d2Z_dpsi_dtheta)
               ! psi is normalized in derivatives - rescale
               s%dq_dpsi = s%dq_dpsi / psipol_max
               s%dR_dpsi = s%dR_dpsi / psipol_max
               s%dZ_dpsi = s%dZ_dpsi / psipol_max
               s%d2R_dpsi_dtheta = s%d2R_dpsi_dtheta / psipol_max
               s%d2Z_dpsi_dtheta = s%d2Z_dpsi_dtheta / psipol_max
            end if
            s%ktri = point_location(s%R, s%Z)
            call field(s%R, 0d0, s%Z, s%B0_R, s%B0_phi, s%B0_Z, &
                 s%dB0R_dR, dum, s%dB0R_dZ, s%dB0phi_dR, dum, s%dB0phi_dZ, &
                 s%dB0Z_dR, dum, s%dB0Z_dZ)
            s%B0_2 = s%B0_R ** 2 + s%B0_Z ** 2 + s%B0_phi ** 2
            if (conf%kilca_pol_mode /= 0 .and. conf%debug_kilca_geom_theta) then
               s%sqrt_g = equil%cocos%sgn_dpsi * s%rad / &
                    (-s%B0_R * sin(s%theta) + s%B0_Z * cos(s%theta))
            else
               ! sqrt_g misses a factor of q and the signs of dpsi_drad and B0_phi
               ! taken together, these three always yield a positive sign in COCOS 3
               s%sqrt_g = s%sqrt_g * abs(s%q)
            end if
          end associate
       end do
    end do
    deallocate(rad, psi, theta)
  end subroutine compute_sample_Ipar

  function point_location(R, Z, hint_psi) result(location)
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
  end function point_location

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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/refinement', mesh%refinement, &
         lbound(mesh%refinement), ubound(mesh%refinement), &
         comment = 'width ratio of neighbouring refined flux surfaces')
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/tri_theta_extent', mesh%tri_theta_extent, &
         lbound(mesh%tri_theta_extent), ubound(mesh%tri_theta_extent), &
         comment = 'range of geometrical poloidal angle covered by triangle')
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
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_kt_low', mesh%shielding_kt_low, &
         lbound(mesh%shielding_kt_low), ubound(mesh%shielding_kt_low), &
         comment = 'triangle indexing for resonant flux surfaces only')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_kp', mesh%shielding_L1_kp, &
         lbound(mesh%shielding_L1_kp), ubound(mesh%shielding_L1_kp), &
         comment = 'local node index for shielding pressure evaluation points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_weight', mesh%shielding_L1_weight, &
         lbound(mesh%shielding_L1_weight), ubound(mesh%shielding_L1_weight), unit = '1', &
         comment = 'weighting factor for shielding pressure evaluation points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_R', mesh%shielding_L1_R, &
         lbound(mesh%shielding_L1_R), ubound(mesh%shielding_L1_R), unit = 'cm', &
         comment = 'R coordinate of shielding pressure evaluation points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_Z', mesh%shielding_L1_Z, &
         lbound(mesh%shielding_L1_Z), ubound(mesh%shielding_L1_Z), unit = 'cm', &
         comment = 'Z coordinate of shielding pressure evaluation points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/shielding_coeff', mesh%shielding_coeff, &
         lbound(mesh%shielding_coeff), ubound(mesh%shielding_coeff), unit = 'g^-1 cm s', &
         comment = 'Coefficient in the compensated scheme for shielding')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac, &
         lbound(mesh%gpec_jacfac), ubound(mesh%gpec_jacfac), unit = 'cm^2', &
         comment = 'Jacobian surface factor between flux surfaces')
    call h5_close(h5id_root)
  end subroutine mesh_write

  subroutine write_FreeFem_mesh
    integer :: fid, kpoi, ktri, kp, kedge

    open(newunit = fid, file = 'inputformaxwell.msh', status = 'replace', form = 'formatted', action = 'write')
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
    open(newunit = fid, file = 'edgemap.dat', status = 'replace', form = 'formatted', action = 'write')
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
  end subroutine write_FreeFem_mesh

  subroutine mesh_read(mesh, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use mephit_conf, only: conf
    type(mesh_t), intent(inout) :: mesh
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root
    integer, allocatable :: orient(:)

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
    allocate(mesh%deletions(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%refinement(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%res_ind(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%psi_res(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%rad_norm_res(mesh%m_res_min:mesh%m_res_max))
    if (conf%kilca_scale_factor /= 0) then
       allocate(mesh%res_modes(1))
    else
       allocate(mesh%res_modes(mesh%m_res_max - mesh%m_res_min + 1))
    end if
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
    allocate(mesh%tri_theta_extent(2, mesh%ntri))
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
    allocate(mesh%shielding_kt_low(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%shielding_coeff(mesh%m_res_min:mesh%m_res_max))
    allocate(mesh%gpec_jacfac(-16:16, mesh%nflux))
    allocate(mesh%area(mesh%ntri))
    allocate(mesh%cntr_R(mesh%ntri))
    allocate(mesh%cntr_Z(mesh%ntri))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/deletions', mesh%deletions)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/refinement', mesh%refinement)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_ind', mesh%res_ind)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/psi_res', mesh%psi_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/rad_norm_res', mesh%rad_norm_res)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/res_modes', mesh%res_modes)
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/tri_theta_extent', mesh%tri_theta_extent)
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
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_kt_low', mesh%shielding_kt_low)
    allocate(mesh%shielding_L1_kp(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_kp', mesh%shielding_L1_kp)
    allocate(mesh%shielding_L1_weight(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_weight', mesh%shielding_L1_weight)
    allocate(mesh%shielding_L1_R(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_R', mesh%shielding_L1_R)
    allocate(mesh%shielding_L1_Z(-1:0, mesh%GL_order, sum(mesh%kt_max(mesh%res_ind))))
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_L1_Z', mesh%shielding_L1_Z)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/shielding_coeff', mesh%shielding_coeff)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/gpec_jacfac', mesh%gpec_jacfac)
    call h5_close(h5id_root)
    where (orient == 1)
       mesh%orient = .true.
    elsewhere
       mesh%orient = .false.
    end where
    deallocate(orient)
  end subroutine mesh_read

  subroutine mesh_deinit(this)
    class(mesh_t), intent(inout) :: this

    if (allocated(this%deletions)) deallocate(this%deletions)
    if (allocated(this%refinement)) deallocate(this%refinement)
    if (allocated(this%res_ind)) deallocate(this%res_ind)
    if (allocated(this%psi_res)) deallocate(this%psi_res)
    if (allocated(this%rad_norm_res)) deallocate(this%rad_norm_res)
    if (allocated(this%res_modes)) deallocate(this%res_modes)
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
    if (allocated(this%tri_theta_extent)) deallocate(this%tri_theta_extent)
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
    if (allocated(this%shielding_kt_low)) deallocate(this%shielding_kt_low)
    if (allocated(this%shielding_L1_kp)) deallocate(this%shielding_L1_kp)
    if (allocated(this%shielding_L1_weight)) deallocate(this%shielding_L1_weight)
    if (allocated(this%shielding_L1_R)) deallocate(this%shielding_L1_R)
    if (allocated(this%shielding_L1_Z)) deallocate(this%shielding_L1_Z)
    if (allocated(this%shielding_coeff)) deallocate(this%shielding_coeff)
    if (allocated(this%gpec_jacfac)) deallocate(this%gpec_jacfac)
  end subroutine mesh_deinit

  subroutine compare_gpec_coordinates
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, nf90_inq_varid, &
         nf90_inquire_dimension, nf90_get_var, nf90_close, nf90_global, nf90_get_att
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext, psipol_max
    use mephit_conf, only: conf, logger, datafile
    use mephit_util, only: pi
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
    do kpol = 1, npol
       do krad = 1, nrad
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
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
    call h5_add(h5id_root, dataset // '/R', R, lbound(R), ubound(R), &
         unit = 'cm', comment = 'R(psi, theta)')
    call h5_add(h5id_root, dataset // '/Z', Z, lbound(Z), ubound(Z), &
         unit = 'cm', comment = 'Z(psi, theta)')
    call h5_add(h5id_root, dataset // '/xi_n_R', xi_n_R, lbound(xi_n_R), ubound(xi_n_R), &
         unit = 'cm', comment = 'Radial component of normal displacement xi_n(theta) at last flux surface')
    call h5_add(h5id_root, dataset // '/xi_n_Z', xi_n_Z, lbound(xi_n_Z), ubound(xi_n_Z), &
         unit = 'cm', comment = 'Axial component of normal displacement xi_n(theta) at last flux surface')
    call h5_close(h5id_root)
    deallocate(psi, theta, R, Z, xi_n, xi_n_R, xi_n_Z)

    nrad = (nrad - 1 + stepsize) / stepsize
    allocate(psi(nrad), theta(npol), R(nrad, npol), Z(nrad, npol), sqrt_g(nrad, npol))
    filename = 'gpec_diagnostics_jacfac_1_fun.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    logger%msg = 'File ' // trim(filename) // ' found, performing GPEC jacfac comparison.'
    if (logger%info) call logger%write_msg
    allocate(jacfac(nrad, npol), contradenspsi(nrad, npol))
    open(newunit = fid, file = filename, status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
       do kpol = 1, npol
          read (fid, '(1x, 4(1x, es16.8))') psi(krad), theta(kpol), jacfac(krad, kpol)
       end do
    end do
    close(fid)
    psi(:) = psipol_max * psi
    theta(:) = 2d0 * pi * theta
    jacfac(:, :) = psipol_max / (2d4 * pi * chi1) * jacfac
    filename = 'gpec_diagnostics_delpsi_fun.out'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) return
    logger%msg = 'File ' // trim(filename) // ' found, performing GPEC grad psi comparison.'
    if (logger%info) call logger%write_msg
    allocate(delpsi(nrad, npol), grad_psi(nrad, npol))
    open(newunit = fid, file = filename, status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
       do kpol = 1, npol
          read (fid, '(1x, 3(1x, es16.8))') psi(krad), theta(kpol), delpsi(krad, kpol)
       end do
    end do
    close(fid)
    psi(:) = psipol_max * psi
    theta(:) = 2d0 * pi * theta
    delpsi(:, :) = 1d-2 * psipol_max * delpsi
    do kpol = 1, npol
       do krad = 1, nrad
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
               q, dum, sqrt_g(krad, kpol), dum, dum, R(krad, kpol), dum, dum, Z(krad, kpol), dum, dum)
          call field(R(krad, kpol), 0d0, Z(krad, kpol), B0_R, dum, B0_Z, &
               dum, dum, dum, dum, dum, dum, dum, dum, dum)
          grad_psi(krad, kpol) = R(krad, kpol) * hypot(B0_Z, -B0_R)
          contradenspsi(krad, kpol) = grad_psi(krad, kpol) * sqrt_g(krad, kpol) * abs(q)
       end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/jacfac', jacfac, lbound(jacfac), ubound(jacfac), &
         unit = '?', comment = 'GPEC jacfac at (psi, theta)')
    call h5_add(h5id_root, dataset // '/contradenspsi', contradenspsi, lbound(contradenspsi), ubound(contradenspsi), &
         unit = 'cm^2', comment = 'MEPHIT jacfac at (psi, theta)')
    call h5_add(h5id_root, dataset // '/delpsi', delpsi, lbound(delpsi), ubound(delpsi), &
         unit = 'G cm', comment = 'GPEC grad psi at (psi, theta)')
    call h5_add(h5id_root, dataset // '/grad_psi', grad_psi, lbound(grad_psi), ubound(grad_psi), &
         unit = 'G cm', comment = 'MEPHIT grad psi at (psi, theta)')
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
    allocate(psi(nrad), theta(npol), jac(nrad, npol), sqrt_g(nrad, npol), R(nrad, npol), Z(nrad, npol))
    open(newunit = fid, file = '2d.out', status = 'old', form = 'formatted', action = 'read')
    read (fid, *)
    read (fid, *)
    do krad = 1, nrad
       read (fid, '(6x, i3, 6x, e24.16)') idum, psi(krad)
       psi(krad) = psi(krad) * psipol_max
       read (fid, *)
       read (fid, *)
       read (fid, *)
       do kpol = 1, npol
          read (fid, '(i6, 1p, 8e24.16)') idum, theta(kpol), dum, dum, dum, dum, &
               R(krad, kpol), Z(krad, kpol), jac(krad, kpol)
          theta(kpol) = theta(kpol) * 2d0 * pi
          call magdata_in_symfluxcoord_ext(2, dum, psi(krad), theta(kpol), &
               q, dum, sqrt_g(krad, kpol), dum, dum, dum, dum, dum, dum, dum, dum)
          sqrt_g(krad, kpol) = sqrt_g(krad, kpol) * abs(q)
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
         unit = 'cm', comment = 'R(psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/Z', Z, lbound(Z), ubound(Z), &
         unit = 'cm', comment = 'Z(psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/sqrt_g', sqrt_g, lbound(sqrt_g), ubound(sqrt_g), &
         unit = 'cm G^-1', comment = 'MEPHIT Jacobian at (psi, theta)')
    call h5_add(h5id_root, dataset // '/jac/jac', jac, lbound(jac), ubound(jac), &
         unit = 'cm G^-1', comment = 'GPEC Jacobian at (psi, theta)')
    call h5_close(h5id_root)
    deallocate(psi, theta, jac, sqrt_g, R, Z)

  contains
    subroutine check_error(funcname, status)
      character(len = *), intent(in) :: funcname
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
         write (logger%msg, '(a, " returned error ", i0)') funcname, status
         if (logger%err) call logger%write_msg
         error stop
      end if
    end subroutine check_error
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
    integer :: kf

    ! initialize fluxvar with equidistant psi values
    call psi_interpolator%init(4, equil%psi_eqd)

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
       call compute_safety_factor_flux
    case (q_prof_rot)
       call compute_safety_factor_rot
    case (q_prof_geqdsk)
       call compute_safety_factor_geqdsk
    case default
       write (logger%msg, '("unknown q profile selection: ", i0)') conf%q_prof
       if (logger%err) call logger%write_msg
       error stop
    end select

    fs%F(:) = [(psi_interpolator%eval(equil%fpol, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs%FdF_dpsi(:) = [(psi_interpolator%eval(equil%ffprim, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%F(:) = [(psi_interpolator%eval(equil%fpol, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    fs_half%FdF_dpsi(:) = [(psi_interpolator%eval(equil%ffprim, fs_half%psi(kf)), kf = 1, mesh%nflux)]
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
    integer :: kf

    fs%p(:) = [(psi_interpolator%eval(equil%pres, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs%dp_dpsi(:) = [(psi_interpolator%eval(equil%pprime, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%p(:) = [(psi_interpolator%eval(equil%pres, fs_half%psi(kf)), kf = 1, mesh%nflux)]
    fs_half%dp_dpsi(:) = [(psi_interpolator%eval(equil%pprime, fs_half%psi(kf)), kf = 1, mesh%nflux)]
  end subroutine compute_pres_prof_geqdsk

  subroutine compute_safety_factor_flux
    use mephit_util, only: pi, interp1d
    integer :: kf, kt, ktri
    type(interp1d) :: psi_half_interpolator

    fs_half%q = 0d0
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          fs_half%q(kf) = fs_half%q(kf) + cache%cntr_fields(ktri)%B0(2) * mesh%area(ktri)
       end do
       fs_half%q(kf) = fs_half%q(kf) * 0.5d0 / pi / (fs%psi(kf) - fs%psi(kf-1))
    end do
    call psi_half_interpolator%init(4, fs_half%psi)
    ! Lagrange polynomial extrapolation for values at separatrix and magnetic axis
    fs%q(:) = [(psi_half_interpolator%eval(fs_half%q, fs%psi(kf)), kf = 0, mesh%nflux)]
    call psi_half_interpolator%deinit
  end subroutine compute_safety_factor_flux

  subroutine compute_safety_factor_rot
    use magdata_in_symfluxcoor_mod, only: qsaf
    integer :: kf

    ! Lagrange polynomial extrapolation for value at magnetic axis
    fs%q(:) = [(psi_fine_interpolator%eval(qsaf, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%q(:) = [(psi_fine_interpolator%eval(qsaf, fs_half%psi(kf)), kf = 1, mesh%nflux)]
  end subroutine compute_safety_factor_rot

  subroutine compute_safety_factor_geqdsk
    integer :: kf

    fs%q(:) = [(psi_interpolator%eval(equil%qpsi, fs%psi(kf)), kf = 0, mesh%nflux)]
    fs_half%q(:) = [(psi_interpolator%eval(equil%qpsi, fs_half%psi(kf)), kf = 1, mesh%nflux)]
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

  subroutine cache_equilibrium_field
    use field_eq_mod, only: psif, psib
    integer :: ktri, kedge, k
    real(dp) :: dum

    ! edge midpoints
    do kedge = 1, mesh%nedge
       associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge), Z => mesh%mid_Z(kedge))
         call field(R, 0d0, Z, f%B0(1), f%B0(2), f%B0(3), &
              f%dB0_dR(1), dum, f%dB0_dZ(1), f%dB0_dR(2), dum, f%dB0_dZ(2), f%dB0_dR(3), dum, f%dB0_dZ(3))
         f%psi = psif - psib  ! see intperp_psi_pol in mephit_util
         f%Bmod = sqrt(sum(f%B0 ** 2))
         f%dBmod_dR = (sum(f%B0 * f%dB0_dR)) / f%Bmod
         f%dBmod_dZ = (sum(f%B0 * f%dB0_dZ)) / f%Bmod
         cache%B0_flux(kedge) = R * (f%B0(1) * mesh%edge_Z(kedge) - f%B0(3) * mesh%edge_R(kedge))
       end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
       associate (f => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri), Z => mesh%cntr_Z(ktri))
         call field(R, 0d0, Z, f%B0(1), f%B0(2), f%B0(3), &
              f%dB0_dR(1), dum, f%dB0_dZ(1), f%dB0_dR(2), dum, f%dB0_dZ(2), f%dB0_dR(3), dum, f%dB0_dZ(3))
         f%psi = psif - psib  ! see intperp_psi_pol in mephit_util
         f%Bmod = sqrt(sum(f%B0 ** 2))
         f%dBmod_dR = (sum(f%B0 * f%dB0_dR)) / f%Bmod
         f%dBmod_dZ = (sum(f%B0 * f%dB0_dZ)) / f%Bmod
       end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
       do k = 1, mesh%GL_order
          associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
            call field(R, 0d0, Z, f%B0(1), f%B0(2), f%B0(3), &
                 f%dB0_dR(1), dum, f%dB0_dZ(1), f%dB0_dR(2), dum, f%dB0_dZ(2), f%dB0_dR(3), dum, f%dB0_dZ(3))
            f%psi = psif - psib  ! see intperp_psi_pol in mephit_util
            f%Bmod = sqrt(sum(f%B0 ** 2))
            f%dBmod_dR = (sum(f%B0 * f%dB0_dR)) / f%Bmod
            f%dBmod_dZ = (sum(f%B0 * f%dB0_dZ)) / f%Bmod
          end associate
       end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
       do k = 1, mesh%GL2_order
          associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
            call field(R, 0d0, Z, f%B0(1), f%B0(2), f%B0(3), &
                 f%dB0_dR(1), dum, f%dB0_dZ(1), f%dB0_dR(2), dum, f%dB0_dZ(2), f%dB0_dR(3), dum, f%dB0_dZ(3))
            f%psi = psif - psib  ! see intperp_psi_pol in mephit_util
            f%Bmod = sqrt(sum(f%B0 ** 2))
            f%dBmod_dR = (sum(f%B0 * f%dB0_dR)) / f%Bmod
            f%dBmod_dZ = (sum(f%B0 * f%dB0_dZ)) / f%Bmod
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
    use mephit_util, only: clight, pi
    integer :: ktri, kedge, k
    real(dp) :: dum, B0_phi, dB0R_dZ, dB0Z_dR, dB0phi_dR, dB0phi_dZ

    ! edge midpoints
    do kedge = 1, mesh%nedge
       associate (f => cache%mid_fields(kedge), R => mesh%mid_R(kedge), Z => mesh%mid_Z(kedge))
         call field(R, 0d0, Z, dum, B0_phi, dum, dum, dum, dB0R_dZ, dB0phi_dR, dum, dB0phi_dZ, dB0Z_dR, dum, dum)
         f%j0(1) = -0.25d0 / pi * clight * dB0phi_dZ
         f%j0(3) = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
         f%j0(2) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
       end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
       associate (f => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri), Z => mesh%cntr_Z(ktri))
         call field(R, 0d0, Z, dum, B0_phi, dum, dum, dum, dB0R_dZ, dB0phi_dR, dum, dB0phi_dZ, dB0Z_dR, dum, dum)
         f%j0(1) = -0.25d0 / pi * clight * dB0phi_dZ
         f%j0(3) = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
         f%j0(2) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
       end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
       do k = 1, mesh%GL_order
          associate (f => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge), Z => mesh%GL_Z(k, kedge))
            call field(R, 0d0, Z, dum, B0_phi, dum, dum, dum, dB0R_dZ, dB0phi_dR, dum, dB0phi_dZ, dB0Z_dR, dum, dum)
            f%j0(1) = -0.25d0 / pi * clight * dB0phi_dZ
            f%j0(3) = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
            f%j0(2) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          end associate
       end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
       do k = 1, mesh%GL2_order
          associate (f => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri), Z => mesh%GL2_Z(k, ktri))
            call field(R, 0d0, Z, dum, B0_phi, dum, dum, dum, dB0R_dZ, dB0phi_dR, dum, dB0phi_dZ, dB0Z_dR, dum, dum)
            f%j0(1) = -0.25d0 / pi * clight * dB0phi_dZ
            f%j0(3) = 0.25d0 / pi * clight * (dB0phi_dR + B0_phi / R)
            f%j0(2) = 0.25d0 / pi * clight * (dB0R_dZ - dB0Z_dR)
          end associate
       end do
    end do
  end subroutine compute_curr0_rot

  subroutine compute_curr0_geqdsk
    use mephit_util, only: clight, pi
    integer :: ktri, kedge, k
    real(dp) :: dp0_dpsi, d2p0_dpsi2, F, dF_dpsi, FdF_dpsi, d2F_dpsi2, fprime(equil%nw)

    fprime = equil%ffprim / equil%fpol
    ! edge midpoints
    do kedge = 1, mesh%nedge
       associate (c => cache%mid_fields(kedge), R => mesh%mid_R(kedge))
         dp0_dpsi = psi_interpolator%eval(equil%pprime, c%psi)
         d2p0_dpsi2 = psi_interpolator%eval(equil%pprime, c%psi, .true.)
         F = psi_interpolator%eval(equil%fpol, c%psi)
         FdF_dpsi = psi_interpolator%eval(equil%ffprim, c%psi)
         dF_dpsi = psi_interpolator%eval(fprime, c%psi)
         d2F_dpsi2 = psi_interpolator%eval(fprime, c%psi, .true.)
         c%j0(1) = 0.25d0 / pi * clight * dF_dpsi * c%B0(1)
         c%j0(3) = 0.25d0 / pi * clight * dF_dpsi * c%B0(3)
         c%j0(2) = clight * (dp0_dpsi * R + 0.25d0 / (pi * R) * FdF_dpsi)
         c%dj0_dR(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(1) &
              + R * c%B0(3) * c%B0(1) * d2F_dpsi2)
         c%dj0_dZ(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(1) &
              - R * c%B0(1) * c%B0(1) * d2F_dpsi2)
         c%dj0_dR(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(3) &
              + R * c%B0(3) * c%B0(3) * d2F_dpsi2)
         c%dj0_dZ(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(3) &
              - R * c%B0(1) * c%B0(3) * d2F_dpsi2)
         c%dj0_dR(2) = clight * (dp0_dpsi + d2p0_dpsi2 * R * R * c%B0(3) + &
              0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2 - FdF_dpsi / R))
         c%dj0_dZ(2) = clight * (-d2p0_dpsi2 * R * R * c%B0(1) + &
              0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2))
       end associate
    end do
    ! weighted triangle centroids
    do ktri = 1, mesh%ntri
       associate (c => cache%cntr_fields(ktri), R => mesh%cntr_R(ktri))
         dp0_dpsi = psi_interpolator%eval(equil%pprime, c%psi)
         d2p0_dpsi2 = psi_interpolator%eval(equil%pprime, c%psi, .true.)
         F = psi_interpolator%eval(equil%fpol, c%psi)
         FdF_dpsi = psi_interpolator%eval(equil%ffprim, c%psi)
         dF_dpsi = psi_interpolator%eval(fprime, c%psi)
         d2F_dpsi2 = psi_interpolator%eval(fprime, c%psi, .true.)
         c%j0(1) = 0.25d0 / pi * clight * dF_dpsi * c%B0(1)
         c%j0(3) = 0.25d0 / pi * clight * dF_dpsi * c%B0(3)
         c%j0(2) = clight * (dp0_dpsi * R + 0.25d0 / (pi * R) * FdF_dpsi)
         c%dj0_dR(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(1) &
              + R * c%B0(3) * c%B0(1) * d2F_dpsi2)
         c%dj0_dZ(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(1) &
              - R * c%B0(1) * c%B0(1) * d2F_dpsi2)
         c%dj0_dR(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(3) &
              + R * c%B0(3) * c%B0(3) * d2F_dpsi2)
         c%dj0_dZ(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(3) &
              - R * c%B0(1) * c%B0(3) * d2F_dpsi2)
         c%dj0_dR(2) = clight * (dp0_dpsi + d2p0_dpsi2 * R * R * c%B0(3) + &
              0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2 - FdF_dpsi / R))
         c%dj0_dZ(2) = clight * (-d2p0_dpsi2 * R * R * c%B0(1) + &
              0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2))
       end associate
    end do
    ! Gauss-Legendre evaluation points on triangle edges
    do kedge = 1, mesh%nedge
       do k = 1, mesh%GL_order
          associate (c => cache%edge_fields(k, kedge), R => mesh%GL_R(k, kedge))
            dp0_dpsi = psi_interpolator%eval(equil%pprime, c%psi)
            d2p0_dpsi2 = psi_interpolator%eval(equil%pprime, c%psi, .true.)
            F = psi_interpolator%eval(equil%fpol, c%psi)
            FdF_dpsi = psi_interpolator%eval(equil%ffprim, c%psi)
            dF_dpsi = psi_interpolator%eval(fprime, c%psi)
            d2F_dpsi2 = psi_interpolator%eval(fprime, c%psi, .true.)
            c%j0(1) = 0.25d0 / pi * clight * dF_dpsi * c%B0(1)
            c%j0(3) = 0.25d0 / pi * clight * dF_dpsi * c%B0(3)
            c%j0(2) = clight * (dp0_dpsi * R + 0.25d0 / (pi * R) * FdF_dpsi)
            c%dj0_dR(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(1) &
                 + R * c%B0(3) * c%B0(1) * d2F_dpsi2)
            c%dj0_dZ(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(1) &
                 - R * c%B0(1) * c%B0(1) * d2F_dpsi2)
            c%dj0_dR(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(3) &
                 + R * c%B0(3) * c%B0(3) * d2F_dpsi2)
            c%dj0_dZ(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(3) &
                 - R * c%B0(1) * c%B0(3) * d2F_dpsi2)
            c%dj0_dR(2) = clight * (dp0_dpsi + d2p0_dpsi2 * R * R * c%B0(3) + &
                 0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2 - FdF_dpsi / R))
            c%dj0_dZ(2) = clight * (-d2p0_dpsi2 * R * R * c%B0(1) + &
                 0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2))
          end associate
       end do
    end do
    ! Gauss-Legendre evaluation points on triangle areas
    do ktri = 1, mesh%ntri
       do k = 1, mesh%GL2_order
          associate (c => cache%area_fields(k, ktri), R => mesh%GL2_R(k, ktri))
            dp0_dpsi = psi_interpolator%eval(equil%pprime, c%psi)
            d2p0_dpsi2 = psi_interpolator%eval(equil%pprime, c%psi, .true.)
            F = psi_interpolator%eval(equil%fpol, c%psi)
            FdF_dpsi = psi_interpolator%eval(equil%ffprim, c%psi)
            dF_dpsi = psi_interpolator%eval(fprime, c%psi)
            d2F_dpsi2 = psi_interpolator%eval(fprime, c%psi, .true.)
            c%j0(1) = 0.25d0 / pi * clight * dF_dpsi * c%B0(1)
            c%j0(3) = 0.25d0 / pi * clight * dF_dpsi * c%B0(3)
            c%j0(2) = clight * (dp0_dpsi * R + 0.25d0 / (pi * R) * FdF_dpsi)
            c%dj0_dR(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(1) &
                 + R * c%B0(3) * c%B0(1) * d2F_dpsi2)
            c%dj0_dZ(1) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(1) &
                 - R * c%B0(1) * c%B0(1) * d2F_dpsi2)
            c%dj0_dR(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dR(3) &
                 + R * c%B0(3) * c%B0(3) * d2F_dpsi2)
            c%dj0_dZ(3) = 0.25d0 / pi * clight * (dF_dpsi * c%dB0_dZ(3) &
                 - R * c%B0(1) * c%B0(3) * d2F_dpsi2)
            c%dj0_dR(2) = clight * (dp0_dpsi + d2p0_dpsi2 * R * R * c%B0(3) + &
                 0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2 - FdF_dpsi / R))
            c%dj0_dZ(2) = clight * (-d2p0_dpsi2 * R * R * c%B0(1) + &
                 0.25d0 / (pi * R) * (dF_dpsi ** 2 + F * d2F_dpsi2))
          end associate
       end do
    end do
  end subroutine compute_curr0_geqdsk

  subroutine check_curr0
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdata_in_symfluxcoor_mod, only: magdata_in_symfluxcoord_ext
    use mephit_conf, only: datafile
    use mephit_util, only: clight, pi, linspace
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
