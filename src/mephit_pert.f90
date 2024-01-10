module mephit_pert

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  ! types and associated procedures
  public :: L1_t, L1_init, L1_deinit, L1_write, L1_read, L1_interp
  public :: RT0_t, RT0_init, RT0_deinit, RT0_write, RT0_read, RT0_interp, &
    RT0_project_pol_comp, RT0_project_tor_comp, RT0_tor_comp_from_zero_div, &
    RT0_L2int_num, RT0_L2int, RT0_triplot, RT0_rectplot
  public :: polmodes_t, polmodes_init, polmodes_deinit, &
    polmodes_write, polmodes_read, L1_poloidal_modes
  public :: vec_polmodes_t, vec_polmodes_init, vec_polmodes_deinit, &
    vec_polmodes_write, vec_polmodes_read, RT0_poloidal_modes
  public :: vac_t, vac_init, vac_deinit, vac_write, vac_read, generate_vacfield

  ! module variables
  public :: vac

  ! interfaces
  abstract interface
    function vector_element_projection(ktri, weight, R, Z, n, f)
      use mephit_mesh, only: field_cache_t
      import :: dp
      integer, intent(in) :: ktri
      real(dp), intent(in) :: weight
      real(dp), intent(in) :: R
      real(dp), intent(in) :: Z
      real(dp), intent(in) :: n(:)
      type(field_cache_t), intent(in) :: f
      complex(dp) :: vector_element_projection
    end function vector_element_projection
  end interface

  type :: L1_t
    !> Number of points on which the L1 are defined
    integer :: npoint

    !> Degrees of freedom, given in base units.
    complex(dp), allocatable :: DOF(:)
  end type L1_t

  type :: RT0_t
    !> Number of triangles on which the RT0 elements are defined
    integer :: ntri

    !> Edge fluxes \f$ R \vec{v} \cdot \vec{n} \f$, given in base units of
    !> \f$ \vec{v} \f$ times cm^2.
    !>
    !> The normal vector points to the right of the edge vector. Thus on poloidal edges,
    !> it points radially outward, and on radial edges, it points in poloidally clockwise direction.
    complex(dp), allocatable :: DOF(:)

    !> Physical toroidal component \f$ v_{n (\phi)} \f$ in base units of \f$ \vec{v} \f$.
    !>
    !> Values are taken at each triangle.
    complex(dp), allocatable :: comp_phi(:)
  end type RT0_t

  type :: polmodes_t
    !> Highest absolute poloidal mode number.
    integer :: m_max

    !> Number of flux surfaces, i.e., radial divisions.
    integer :: nflux

    !> Fourier coefficient indexed by poloidal mode number and flux surface.
    complex(dp), allocatable :: coeff(:, :)
  end type polmodes_t

  type :: vec_polmodes_t
    !> Highest absolute poloidal mode number.
    integer :: m_max

    !> Number of flux surfaces, i.e., radial divisions.
    integer :: nflux

    !> Component in radial direction, indexed by poloidal mode number and flux surface.
    !>
    !> For tokamak geometry, this is the contravariant density psi component; for KiLCA
    !> geometry, this is the r component.
    complex(dp), allocatable :: coeff_rad(:, :)

    !> Component normal to flux surface, indexed by poloidal mode number and flux surface.
    !>
    !> For tokamak and KiLCA geometry, this is the physical component perpendicular to the
    !> flux surface.
    complex(dp), allocatable :: coeff_n(:, :)

    !> Component in poloidal direction, indexed by poloidal mode number and flux surface.
    !>
    !> For tokamak geometry, this is the covariant theta component; for KiLCA geometry,
    !> this is the physical theta component.
    complex(dp), allocatable :: coeff_pol(:, :)

    !> Component in radial direction, indexed by poloidal mode number and flux surface.
    !>
    !> For tokamak and KiLCA geometry, this is the physical phi component.
    complex(dp), allocatable :: coeff_tor(:, :)
  end type vec_polmodes_t

  type :: vac_t
    !> Vacuum perturbation field in units of Gauss.
    type(RT0_t) :: Bn

    !> Lower and upper bound of of vac_t::kilca_vac_coeff and vac_t::kilca_pol_modes.
    integer :: m_min, m_max

    !> Single poloidal mode number used with KiLCA. If zero, use whole available range.
    integer :: kilca_pol_mode

    !> poloidal modes excited in KiLCA
    integer, allocatable :: kilca_pol_modes(:)

    !> Integration constant for resonant vacuum perturbation in KiLCA comparison.
    complex(dp), allocatable :: kilca_vac_coeff(:)
  end type vac_t

  type(vac_t) :: vac

contains

  subroutine L1_init(this, npoint)
    type(L1_t), intent(inout) :: this
    integer, intent(in) :: npoint

    call L1_deinit(this)
    this%npoint = npoint
    allocate(this%DOF(npoint))
    this%DOF = (0d0, 0d0)
  end subroutine L1_init

  subroutine L1_deinit(this)
    type(L1_t), intent(inout) :: this

    this%npoint = 0
    if (allocated(this%DOF)) deallocate(this%DOF)
  end subroutine L1_deinit

  !> Real and imaginary part of \p scalar_dof (e.g. #presn) are written, in that order,
  !> to \p outfile (e.g. #mephit_conf::presn_file), where line number corresponds to the
  !> knot index in #mesh_mod::mesh_point.
  subroutine L1_write(elem, file, dataset, comment, unit)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close

    type(L1_t), intent(in) :: elem
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    character(len = *), intent(in) :: unit
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/L1_DOF', &
      elem%DOF, lbound(elem%DOF), ubound(elem%DOF), &
      comment = 'degrees of freedom of ' // trim(adjustl(comment)), &
      unit = trim(adjustl(unit)))
    call h5_close(h5id_root)
  end subroutine L1_write

  subroutine L1_read(elem, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(L1_t), intent(inout) :: elem
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/L1_DOF', elem%DOF)
    call h5_close(h5id_root)
  end subroutine L1_read

  subroutine L1_interp(elem, ktri, R, Z, val, grad_val)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mephit_conf, only: logger
    use mephit_util, only: imun
    use mephit_mesh, only: mesh
    integer, intent(in) :: ktri
    type(L1_t), intent(in) :: elem
    real(dp), intent(in) :: R, Z
    complex(dp), intent(out) :: val
    complex(dp), intent(out), dimension(:), optional :: grad_val
    integer :: nodes(3)
    real(dp) :: nan, Delta_R(3), Delta_Z(3)

    if (present(grad_val)) then
      if (3 /= size(grad_val)) then
        call logger%msg_arg_size('L1_interp', '3', 'size(grad_val)', 3, size(grad_val))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end if
    if (ktri <= 0 .or. ktri > mesh%ntri) then
      nan = ieee_value(1d0, ieee_quiet_nan)
      val = cmplx(nan, nan, dp)
      if (present(grad_val)) then
        grad_val(:) = cmplx(nan, nan, dp)
      end if
      return
    end if
    nodes = mesh%tri_node(:, ktri)
    Delta_R = R - mesh%node_R(nodes)
    Delta_Z = Z - mesh%node_Z(nodes)
    val = sum(elem%DOF(nodes) * (Delta_R([2, 3, 1]) * Delta_Z([3, 1, 2]) - &
      Delta_R([3, 1, 2]) * Delta_Z([2, 3, 1]))) * 0.5d0 / mesh%area(ktri)
    if (present(grad_val)) then
      grad_val(1) = sum(elem%DOF(nodes) * (Delta_Z([3, 1, 2]) - Delta_Z([2, 3, 1]))) * &
        0.5d0 / mesh%area(ktri)
      grad_val(2) = imun * mesh%n / R * val
      grad_val(3) = sum(elem%DOF(nodes) * (Delta_R([2, 3, 1]) - Delta_R([3, 1, 2]))) * &
        0.5d0 / mesh%area(ktri)
    end if
  end subroutine L1_interp

  subroutine RT0_init(this, nedge, ntri)
    type(RT0_t), intent(inout) :: this
    integer, intent(in) :: nedge, ntri

    call RT0_deinit(this)
    this%ntri = ntri
    allocate(this%DOF(nedge))
    allocate(this%comp_phi(ntri))
    this%DOF = (0d0, 0d0)
    this%comp_phi = (0d0, 0d0)
  end subroutine RT0_init

  subroutine RT0_deinit(this)
    type(RT0_t), intent(inout) :: this

    this%ntri = 0
    if (allocated(this%DOF)) deallocate(this%DOF)
    if (allocated(this%comp_phi)) deallocate(this%comp_phi)
  end subroutine RT0_deinit

  subroutine RT0_interp(elem, ktri, R, Z, vec, dvec_dR, dvec_dphi, dvec_dZ)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mephit_conf, only: logger
    use mephit_util, only: imun
    use mephit_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
    integer, intent(in) :: ktri
    real(dp), intent(in) :: R, Z
    complex(dp), dimension(:), intent(out) :: vec
    complex(dp), dimension(:), intent(out), optional :: dvec_dR, dvec_dphi, dvec_dZ
    integer :: nodes(3)
    real(dp) :: nan
    complex(dp) :: DOF(3)

    if (3 /= size(vec)) then
      call logger%msg_arg_size('RT0_interp', '3', 'size(vec)', 3, size(vec))
      if (logger%err) call logger%write_msg
      error stop
    end if
    if (present(dvec_dR)) then
      if (3 /= size(dvec_dR)) then
        call logger%msg_arg_size('RT0_interp', '3', 'size(dvec_dR)', 3, size(dvec_dR))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end if
    if (present(dvec_dphi)) then
      if (3 /= size(dvec_dphi)) then
        call logger%msg_arg_size('RT0_interp', '3', 'size(dvec_dphi)', 3, size(dvec_dphi))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end if
    if (present(dvec_dZ)) then
      if (3 /= size(dvec_dZ)) then
        call logger%msg_arg_size('RT0_interp', '3', 'size(dvec_dZ)', 3, size(dvec_dZ))
        if (logger%err) call logger%write_msg
        error stop
      end if
    end if
    if (ktri <= 0 .or. ktri > mesh%ntri) then
      nan = ieee_value(1d0, ieee_quiet_nan)
      vec(:) = cmplx(nan, nan, dp)
      if (present(dvec_dR)) then
        dvec_dR(:) = cmplx(nan, nan, dp)
      end if
      if (present(dvec_dphi)) then
        dvec_dphi(:) = cmplx(nan, nan, dp)
      end if
      if (present(dvec_dZ)) then
        dvec_dZ(:) = cmplx(nan, nan, dp)
      end if
      return
    end if
    ! indices of nodes and opposite edges as defined in mephit_mesh::connect_mesh_points
    nodes = mesh%tri_node(:, ktri)
    DOF = elem%DOF(mesh%tri_edge(:, ktri))
    if (mesh%orient(ktri)) then
      nodes = nodes([3, 1, 2])
      DOF(2) = -DOF(2)
    else
      nodes = nodes([1, 3, 2])
      DOF(1:2) = -DOF(1:2)
    end if
    vec(1) = 0.5d0 / mesh%area(ktri) / R * sum(DOF * (R - mesh%node_R(nodes)))
    vec(2) = elem%comp_phi(ktri)
    vec(3) = 0.5d0 / mesh%area(ktri) / R * sum(DOF * (Z - mesh%node_Z(nodes)))
    if (present(dvec_dR)) then
      dvec_dR(1) = 0.5d0 / mesh%area(ktri) / R ** 2 * sum(DOF * mesh%node_R(nodes))
      dvec_dR(2) = (0d0, 0d0)
      dvec_dR(3) = -vec(3) / R
    end if
    if (present(dvec_dphi)) then
      dvec_dphi(:) = imun * mesh%n / R * vec
    end if
    if (present(dvec_dZ)) then
      dvec_dZ(1) = (0d0, 0d0)
      dvec_dZ(2) = (0d0, 0d0)
      dvec_dZ(3) = 0.5d0 / mesh%area(ktri) / R * sum(DOF)
    end if
  end subroutine RT0_interp

  subroutine RT0_project_pol_comp(elem, proj)
    use mephit_mesh, only: mesh, cache
    type(RT0_t), intent(inout) :: elem
    procedure(vector_element_projection) :: proj
    integer :: kedge, k, ktri
    real(dp) :: n_f(3)

    do kedge = 1, mesh%nedge
      ktri = mesh%edge_tri(1, kedge)
      n_f = [mesh%edge_Z(kedge), 0d0, -mesh%edge_R(kedge)]
      do k = 1, mesh%GL_order
        elem%DOF(kedge) = elem%DOF(kedge) + proj(ktri, mesh%GL_weights(k) * mesh%GL_R(k, kedge), &
          mesh%GL_R(k, kedge), mesh%GL_Z(k, kedge), n_f, cache%edge_fields(k, kedge))
      end do
    end do
  end subroutine RT0_project_pol_comp

  subroutine RT0_project_tor_comp(elem, proj)
    use mephit_mesh, only: mesh, cache
    type(RT0_t), intent(inout) :: elem
    procedure(vector_element_projection) :: proj
    integer :: ktri, k
    real(dp), parameter :: n_f(3) = [0d0, 1d0, 0d0]

    do ktri = 1, mesh%ntri
      do k = 1, mesh%GL2_order
        elem%comp_phi(ktri) = elem%comp_phi(ktri) + proj(ktri, mesh%GL2_weights(k), &
          mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), n_f, cache%area_fields(k, ktri))
      end do
    end do
  end subroutine RT0_project_tor_comp

  pure subroutine RT0_tor_comp_from_zero_div(elem)
    use mephit_util, only: imun
    use mephit_mesh, only: mesh
    type(RT0_t), intent(inout) :: elem
    integer :: ktri

    forall (ktri = 1:mesh%ntri, mesh%orient(ktri))
      elem%comp_phi(ktri) = imun / mesh%n / mesh%area(ktri) * (elem%DOF(mesh%tri_edge(1, ktri)) &
        - elem%DOF(mesh%tri_edge(2, ktri)) + elem%DOF(mesh%tri_edge(3, ktri)))
    end forall
    forall (ktri = 1:mesh%ntri, .not. mesh%orient(ktri))
      elem%comp_phi(ktri) = imun / mesh%n / mesh%area(ktri) * (-elem%DOF(mesh%tri_edge(1, ktri)) &
        - elem%DOF(mesh%tri_edge(2, ktri)) + elem%DOF(mesh%tri_edge(3, ktri)))
    end forall
  end subroutine RT0_tor_comp_from_zero_div

  function RT0_L2int_num(elem)
    use mephit_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
    real(dp) :: RT0_L2int_num
    integer :: ktri, k
    complex(dp) :: vec(3)
    real(dp) :: series

    RT0_L2int_num = 0d0
    do ktri = 1, mesh%ntri
      series = 0d0
      do k = 1, mesh%GL2_order
        call RT0_interp(elem, ktri, mesh%GL2_R(k, ktri), mesh%GL2_Z(k, ktri), vec)
        series = series + mesh%GL2_weights(k) * mesh%GL2_R(k, ktri) ** 2 * &
          (vec(1)%re ** 2 + vec(1)%im ** 2 + vec(3)%re ** 2 + vec(3)%im ** 2)
      end do
      RT0_L2int_num = RT0_L2int_num + series * mesh%area(ktri)
    end do
    RT0_L2int_num = sqrt(RT0_L2int_num)
  end function RT0_L2int_num

  function RT0_L2int(elem)
    use mephit_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
    real(dp) :: RT0_L2int
    integer :: ktri
    real(dp) :: prods(6)

    RT0_L2int = 0d0
    do ktri = 1, mesh%ntri
      associate (R => mesh%edge_R, Z => mesh%edge_Z, f => mesh%tri_edge(1, ktri), &
        o => mesh%tri_edge(2, ktri), i => mesh%tri_edge(3, ktri))
        prods = [elem%DOF(o)%re ** 2 + elem%DOF(o)%im ** 2, &
          elem%DOF(f)%re ** 2 + elem%DOF(f)%im ** 2, &
          elem%DOF(i)%re ** 2 + elem%DOF(i)%im ** 2, &
          2d0 * (elem%DOF(f)%re * elem%DOF(i)%re + elem%DOF(f)%im * elem%DOF(i)%im), &
          2d0 * (elem%DOF(i)%re * elem%DOF(o)%re + elem%DOF(i)%im * elem%DOF(o)%im), &
          2d0 * (elem%DOF(o)%re * elem%DOF(f)%re + elem%DOF(o)%im * elem%DOF(f)%im)]
        if (mesh%orient(ktri)) then
          prods(5:6) = -prods(5:6)  ! negate DOF(o)
        else
          prods(4:5) = -prods(4:5)  ! negate DOF(o) and DOF(f)
        end if
        ! negate R(o) and Z(o) because edge o points in clockwise direction locally
        RT0_L2int = RT0_L2int + ( &
          (sum([3, 1, 1, 1, -1, -1] * prods)) * (R(i) ** 2 + Z(i) ** 2) - &
          (sum([3, -1, 3, 1, -3, 1] * prods)) * (R(i) * R(o) + Z(i) * Z(o)) + &
          (sum([1, 1, 3, -1, -1, 1] * prods)) * (R(o) ** 2 + Z(o) ** 2) &
          ) / (24 * mesh%area(ktri))
      end associate
    end do
    RT0_L2int = sqrt(RT0_L2int)
  end function RT0_L2int

  subroutine RT0_read(elem, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(RT0_t), intent(inout) :: elem
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/RT0_DOF', elem%DOF)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/RT0_comp_phi', elem%comp_phi)
    call h5_close(h5id_root)
  end subroutine RT0_read

  pure function jacobian(kf, kt, R) result(metric_det)
    use mephit_mesh, only: equil, fs_half, mesh, cache
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: R
    real(dp) :: metric_det
    metric_det = equil%cocos%sgn_dpsi * fs_half%q(kf) * R / &
      cache%cntr_fields(mesh%kt_low(kf) + kt)%B0(2)
  end function jacobian

  subroutine RT0_write(elem, file, dataset, comment, unit, plots)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(RT0_t), intent(in) :: elem
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    character(len = *), intent(in) :: unit
    integer, intent(in) :: plots
    integer(HID_T) :: h5id_root
    real(dp), allocatable :: rect_R(:), rect_Z(:)
    complex(dp), allocatable :: comp_R(:), comp_Z(:), comp_psi_contravar_dens(:), &
      comp_theta_covar(:), rect_comp_R(:, :), rect_comp_phi(:, :), rect_comp_Z(:, :)

    if (plots >= 1) then
      call RT0_triplot(elem, comp_R, comp_Z, comp_psi_contravar_dens, comp_theta_covar)
    end if
    if (plots >= 2) then
      call RT0_rectplot(elem, rect_R, rect_Z, rect_comp_R, rect_comp_phi, rect_comp_Z)
    end if
    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/RT0_DOF', &
      elem%DOF, lbound(elem%DOF), ubound(elem%DOF), &
      comment = 'degrees of freedom of ' // trim(adjustl(comment)), &
      unit = trim(adjustl(unit)) // ' cm^2')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/RT0_comp_phi', &
      elem%comp_phi, lbound(elem%comp_phi), ubound(elem%comp_phi), &
      comment = 'physical phi component of ' // trim(adjustl(comment)), &
      unit = trim(adjustl(unit)))
    if (plots >= 1) then
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/comp_R', &
        comp_R, lbound(comp_R), ubound(comp_R), &
        comment = 'R component of ' // trim(adjustl(comment)) // ' at centroid', &
        unit = trim(adjustl(unit)))
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/comp_Z', &
        comp_Z, lbound(comp_Z), ubound(comp_Z), &
        comment = 'Z component of ' // trim(adjustl(comment)) // ' at centroid', &
        unit = trim(adjustl(unit)))
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/comp_psi_contravar_dens', &
        comp_psi_contravar_dens, lbound(comp_psi_contravar_dens), ubound(comp_psi_contravar_dens), &
        comment = 'contravariant density psi component of ' // trim(adjustl(comment)) // ' at centroid', &
        unit = trim(adjustl(unit)) // ' cm^2')
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/comp_theta_covar', &
        comp_theta_covar, lbound(comp_theta_covar), ubound(comp_theta_covar), &
        comment = 'covariant theta component of ' // trim(adjustl(comment)) // ' at centroid', &
        unit = trim(adjustl(unit)) // ' cm')
      deallocate(comp_R, comp_Z, comp_psi_contravar_dens, comp_theta_covar)
    end if
    if (plots >= 2) then
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/rect_R', rect_R, &
        lbound(rect_R), ubound(rect_R), comment = 'R coordinate of rectangular grid', unit = 'cm')
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/rect_Z', rect_Z, &
        lbound(rect_Z), ubound(rect_Z), comment = 'Z coordinate of rectangular grid', unit = 'cm')
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/rect_comp_R', &
        rect_comp_R, lbound(rect_comp_R), ubound(rect_comp_R), &
        comment = 'R component of ' // trim(adjustl(comment)) // ' on GEQSDK grid', &
        unit = trim(adjustl(unit)))
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/rect_comp_phi', &
        rect_comp_phi, lbound(rect_comp_phi), ubound(rect_comp_phi), &
        comment = 'physical phi component of ' // trim(adjustl(comment)) // ' on GEQSDK grid', &
        unit = trim(adjustl(unit)))
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/rect_comp_Z', &
        rect_comp_Z, lbound(rect_comp_Z), ubound(rect_comp_Z), &
        comment = 'Z component of ' // trim(adjustl(comment)) // ' on GEQDSK grid', &
        unit = trim(adjustl(unit)))
      deallocate(rect_R, rect_Z, rect_comp_R, rect_comp_phi, rect_comp_Z)
    end if
    call h5_close(h5id_root)
  end subroutine RT0_write

  subroutine RT0_triplot(elem, comp_R, comp_Z, comp_psi_contravar_dens, comp_theta_covar)
    use mephit_mesh, only: equil, mesh, cache
    type(RT0_t), intent(in) :: elem
    complex(dp), intent(out), dimension(:), allocatable :: comp_R, comp_Z, &
      comp_psi_contravar_dens, comp_theta_covar
    integer :: kf, kt, ktri
    complex(dp) :: vec(3)

    allocate(comp_R(mesh%ntri), comp_Z(mesh%ntri), &
      comp_psi_contravar_dens(mesh%ntri), comp_theta_covar(mesh%ntri))
    do kf = 1, mesh%nflux
      do kt = 1, mesh%kt_max(kf)
        ktri = mesh%kt_low(kf) + kt
        associate (R => mesh%cntr_R(ktri), Z => mesh%cntr_Z(ktri), f => cache%cntr_fields(ktri))
          call RT0_interp(elem, ktri, R, Z, vec)
          comp_R(ktri) = vec(1)
          comp_Z(ktri) = vec(3)
          ! projection to contravariant psi component
          comp_psi_contravar_dens(ktri) = (comp_R(ktri) * f%B0(3) - &
            comp_Z(ktri) * f%B0(1)) * R * jacobian(kf, kt, R)
          ! projection to covariant theta component
          comp_theta_covar(ktri) = equil%cocos%sgn_dpsi * jacobian(kf, kt, R) * &
            (comp_R(ktri) * f%B0(1) + comp_Z(ktri) * f%B0(3))
        end associate
      end do
    end do
  end subroutine RT0_triplot

  subroutine RT0_rectplot(elem, rect_R, rect_Z, comp_R, comp_phi, comp_Z)
    use mephit_mesh, only: equil, mesh, point_location
    type(RT0_t), intent(in) :: elem
    real(dp), intent(out), dimension(:), allocatable :: rect_R, rect_Z
    complex(dp), intent(out), dimension(:, :), allocatable :: comp_R, comp_phi, comp_Z
    integer :: kw, kh, ktri
    complex(dp) :: vec(3)

    allocate(rect_R(equil%nw), rect_Z(equil%nh))
    rect_R(:) = equil%R_eqd
    rect_Z(:) = equil%Z_eqd
    allocate(comp_R(equil%nw, equil%nh), comp_phi(equil%nw, equil%nh), comp_Z(equil%nw, equil%nh))
    do kw = 1, equil%nw
      do kh = 1, equil%nh
        ktri = point_location(rect_R(kw), rect_Z(kh))
        if (ktri > mesh%kt_low(1) .and. ktri <= mesh%ntri) then
          call RT0_interp(elem, ktri, rect_R(kw), rect_Z(kh), vec)
          comp_R(kw, kh) = vec(1)
          comp_phi(kw, kh) = vec(2)
          comp_Z(kw, kh) = vec(3)
        else
          comp_R(kw, kh) = (0d0, 0d0)
          comp_phi(kw, kh) = (0d0, 0d0)
          comp_Z(kw, kh) = (0d0, 0d0)
        end if
      end do
    end do
  end subroutine RT0_rectplot

  subroutine polmodes_init(this, m_max, nflux)
    type(polmodes_t), intent(inout) :: this
    integer, intent(in) :: m_max
    integer, intent(in) :: nflux

    call polmodes_deinit(this)
    this%m_max = abs(m_max)
    this%nflux = nflux
    allocate(this%coeff(-this%m_max:this%m_max, 0:nflux))
  end subroutine polmodes_init

  subroutine polmodes_deinit(this)
    type(polmodes_t), intent(inout) :: this

    this%m_max = 0
    this%nflux = 0
    if (allocated(this%coeff)) deallocate(this%coeff)
  end subroutine polmodes_deinit

  subroutine polmodes_read(polmodes, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(polmodes_t), intent(inout) :: polmodes
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/coeff', polmodes%coeff)
    call h5_close(h5id_root)
  end subroutine polmodes_read

  subroutine polmodes_write(polmodes, file, dataset, comment, unit)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(polmodes_t), intent(in) :: polmodes
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    character(len = *), intent(in) :: unit
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_max', polmodes%m_max, &
      'maximal absolute poloidal mode number')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff', &
      polmodes%coeff, lbound(polmodes%coeff), ubound(polmodes%coeff), &
      comment = trim(adjustl(comment)), unit = trim(adjustl(unit)))
    call h5_close(h5id_root)
  end subroutine polmodes_write

  subroutine vec_polmodes_init(this, m_max, nflux)
    type(vec_polmodes_t), intent(inout) :: this
    integer, intent(in) :: m_max
    integer, intent(in) :: nflux

    call vec_polmodes_deinit(this)
    this%m_max = abs(m_max)
    this%nflux = nflux
    allocate(this%coeff_rad(-this%m_max:this%m_max, 1:nflux))
    allocate(this%coeff_n(-this%m_max:this%m_max, 1:nflux))
    allocate(this%coeff_pol(-this%m_max:this%m_max, 1:nflux))
    allocate(this%coeff_tor(-this%m_max:this%m_max, 1:nflux))
  end subroutine vec_polmodes_init

  subroutine vec_polmodes_deinit(this)
    type(vec_polmodes_t), intent(inout) :: this

    this%m_max = 0
    this%nflux = 0
    if (allocated(this%coeff_rad)) deallocate(this%coeff_rad)
    if (allocated(this%coeff_n)) deallocate(this%coeff_n)
    if (allocated(this%coeff_pol)) deallocate(this%coeff_pol)
    if (allocated(this%coeff_tor)) deallocate(this%coeff_tor)
  end subroutine vec_polmodes_deinit

  subroutine vec_polmodes_read(vec_polmodes, file, dataset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(vec_polmodes_t), intent(inout) :: vec_polmodes
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/coeff_rad', vec_polmodes%coeff_rad)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/coeff_n', vec_polmodes%coeff_n)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/coeff_pol', vec_polmodes%coeff_pol)
    call h5_get(h5id_root, trim(adjustl(dataset)) // '/coeff_tor', vec_polmodes%coeff_tor)
    call h5_close(h5id_root)
  end subroutine vec_polmodes_read

  subroutine vec_polmodes_write(vec_polmodes, file, dataset, comment, unit)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf
    type(vec_polmodes_t), intent(in) :: vec_polmodes
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: dataset
    character(len = *), intent(in) :: comment
    character(len = *), intent(in) :: unit
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/m_max', vec_polmodes%m_max, &
      'maximal absolute poloidal mode number')
    if (conf%kilca_scale_factor /= 0) then
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_rad', &
        vec_polmodes%coeff_rad, lbound(vec_polmodes%coeff_rad), ubound(vec_polmodes%coeff_rad), &
        comment = 'r component of ' // trim(adjustl(comment)), &
        unit = trim(adjustl(unit)))
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_pol', &
        vec_polmodes%coeff_pol, lbound(vec_polmodes%coeff_pol), ubound(vec_polmodes%coeff_pol), &
        comment = 'physical theta component of  ' // trim(adjustl(comment)), &
        unit = trim(adjustl(unit)))
    else
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_rad', &
        vec_polmodes%coeff_rad, lbound(vec_polmodes%coeff_rad), ubound(vec_polmodes%coeff_rad), &
        comment = 'contravariant density psi component of ' // trim(adjustl(comment)), &
        unit = trim(adjustl(unit)) // ' cm^2')
      call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_pol', &
        vec_polmodes%coeff_pol, lbound(vec_polmodes%coeff_pol), ubound(vec_polmodes%coeff_pol), &
        comment = 'covariant theta component of ' // trim(adjustl(comment)), &
        unit = trim(adjustl(unit)) // ' cm')
    end if
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_n', &
      vec_polmodes%coeff_n, lbound(vec_polmodes%coeff_n), ubound(vec_polmodes%coeff_n), &
      comment = 'physical component perpendicular to flux surface of ' // trim(adjustl(comment)), &
      unit = trim(adjustl(unit)))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/coeff_tor', &
      vec_polmodes%coeff_tor, lbound(vec_polmodes%coeff_tor), ubound(vec_polmodes%coeff_tor), &
      comment = 'physical phi component of ' // trim(adjustl(comment)), &
      unit = trim(adjustl(unit)))
    call h5_close(h5id_root)
  end subroutine vec_polmodes_write

  subroutine L1_poloidal_modes(elem, polmodes)
    use mephit_mesh, only: mesh, cache
    type(L1_t), intent(in) :: elem
    type(polmodes_t), intent(inout) :: polmodes
    integer :: kf, log2, kpol, k

    polmodes%coeff(:, :) = (0d0, 0d0)
    polmodes%coeff(0, 0) = elem%DOF(1)
    do kf = 1, mesh%nflux
      log2 = cache%log2_kp_max(kf)
      do kpol = 1, shiftl(1, log2)
        k = cache%kp_low(kf) + kpol
        associate (s => cache%sample_polmodes(k))
          call L1_interp(elem, s%ktri, s%R, s%Z, cache%fft(log2)%samples(kpol))
        end associate
      end do
      call cache%fft(log2)%apply(-polmodes%m_max, polmodes%m_max, polmodes%coeff(:, kf))
    end do
  end subroutine L1_poloidal_modes

  subroutine RT0_poloidal_modes(elem, vec_polmodes)
    use mephit_conf, only: conf
    use mephit_util, only: bent_cyl2straight_cyl
    use mephit_mesh, only: equil, mesh, cache
    type(RT0_t), intent(in) :: elem
    type(vec_polmodes_t), intent(inout) :: vec_polmodes
    integer :: npol, log2, kf, kpol, k
    complex(dp) :: cyl_vec(3)
    complex(dp), dimension(:), allocatable :: comp_rad, comp_n, comp_pol, comp_tor

    npol = shiftl(1, cache%max_log2)
    allocate(comp_rad(npol), comp_n(npol), comp_pol(npol), comp_tor(npol))
    vec_polmodes%coeff_rad(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_n(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_pol(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_tor(:, :) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      log2 = cache%log2_kt_max(kf)
      npol = shiftl(1, log2)
      do kpol = 1, npol
        k = cache%kt_low(kf) + kpol
        associate (s => cache%sample_polmodes_half(k))
          call RT0_interp(elem, s%ktri, s%R, s%Z, cyl_vec)
          if (conf%kilca_scale_factor /= 0) then
            call bent_cyl2straight_cyl(cyl_vec(1), cyl_vec(2), cyl_vec(3), &
              s%theta, comp_rad(kpol), comp_pol(kpol), comp_tor(kpol))
            comp_n(kpol) = comp_rad(kpol)
          else
            comp_rad(kpol) = s%sqrt_g * (cyl_vec(1) * s%B0_Z - cyl_vec(3) * s%B0_R) * s%R
            comp_n(kpol) = (cyl_vec(1) * s%B0_Z - cyl_vec(3) * s%B0_R) / hypot(s%B0_R, s%B0_Z)
            comp_pol(kpol) = cyl_vec(1) * s%dR_dtheta + cyl_vec(3) * s%dZ_dtheta
            comp_tor(kpol) = cyl_vec(2)
          end if
        end associate
      end do
      cache%fft(log2)%samples(:) = comp_rad(:npol)
      call cache%fft(log2)%apply(-vec_polmodes%m_max, vec_polmodes%m_max, &
        vec_polmodes%coeff_rad(:, kf))
      cache%fft(log2)%samples(:) = comp_n(:npol)
      call cache%fft(log2)%apply(-vec_polmodes%m_max, vec_polmodes%m_max, &
        vec_polmodes%coeff_n(:, kf))
      vec_polmodes%coeff_n(:, kf) = vec_polmodes%coeff_n(:, kf) * equil%cocos%sgn_dpsi
      cache%fft(log2)%samples(:) = comp_pol(:npol)
      call cache%fft(log2)%apply(-vec_polmodes%m_max, vec_polmodes%m_max, &
        vec_polmodes%coeff_pol(:, kf))
      cache%fft(log2)%samples(:) = comp_tor(:npol)
      call cache%fft(log2)%apply(-vec_polmodes%m_max, vec_polmodes%m_max, &
        vec_polmodes%coeff_tor(:, kf))
    end do
    deallocate(comp_rad, comp_n, comp_pol, comp_tor)
  end subroutine RT0_poloidal_modes

  subroutine vac_init(this, nedge, ntri, m_min, m_max)
    use mephit_conf, only: conf
    class(vac_t), intent(inout) :: this
    integer, intent(in) :: nedge, ntri, m_min, m_max
    integer :: m

    call vac_deinit(this)
    call RT0_init(this%Bn, nedge, ntri)
    if (conf%kilca_scale_factor /= 0) then
      this%m_min = m_min
      this%m_max = m_max
      if (conf%kilca_pol_mode /= 0) then
        allocate(vac%kilca_pol_modes(2))
        vac%kilca_pol_modes(:) = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
      else
        allocate(vac%kilca_pol_modes(2 * (m_max - m_min + 1)))
        vac%kilca_pol_modes(:) = [([m, -m], m = m_min, m_max)]
      end if
      allocate(this%kilca_vac_coeff(m_min:m_max))
      this%kilca_vac_coeff = (0d0, 0d0)
    end if
  end subroutine vac_init

  subroutine vac_deinit(this)
    class(vac_t), intent(inout) :: this

    call RT0_deinit(this%Bn)
    this%m_min = 0
    this%m_max = 0
    this%kilca_pol_mode = 0
    if (allocated(this%kilca_pol_modes)) deallocate(this%kilca_pol_modes)
    if (allocated(this%kilca_vac_coeff)) deallocate(this%kilca_vac_coeff)
  end subroutine vac_deinit

  subroutine vac_read(this, file, grp_name)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use mephit_conf, only: conf
    class(vac_t), intent(inout) :: this
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: grp_name
    integer(HID_T) :: h5id_root

    call RT0_read(vac%Bn, file, trim(grp_name) // '/Bn')
    if (conf%kilca_scale_factor /= 0) then
      call h5_open(file, h5id_root)
      call h5_get(h5id_root, trim(adjustl(grp_name)) // '/m_min', this%m_max)
      call h5_get(h5id_root, trim(adjustl(grp_name)) // '/m_max', this%m_max)
      call h5_get(h5id_root, trim(adjustl(grp_name)) // '/kilca_pol_mode', this%kilca_pol_mode)
      call h5_get(h5id_root, trim(adjustl(grp_name)) // '/kilca_pol_modes', this%kilca_pol_modes)
      call h5_get(h5id_root, trim(adjustl(grp_name)) // '/kilca_vac_coeff', this%kilca_vac_coeff)
      call h5_close(h5id_root)
    end if
  end subroutine vac_read

  subroutine vac_write(this, file, grp_name)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_add, h5_close
    use mephit_conf, only: conf
    class(vac_t), intent(in) :: this
    character(len = *), intent(in) :: file
    character(len = *), intent(in) :: grp_name
    integer(HID_T) :: h5id_root

    call RT0_write(vac%Bn, file, trim(adjustl(grp_name)) // '/Bn', &
      'magnetic field (vacuum)', 'G', 1)
    if (conf%kilca_scale_factor /= 0) then
      ! already called h5_create_parent_groups
      call h5_open_rw(file, h5id_root)
      call h5_add(h5id_root, trim(adjustl(grp_name)) // '/m_min', this%m_max, &
        'minimum poloidal mode number')
      call h5_add(h5id_root, trim(adjustl(grp_name)) // '/m_max', this%m_max, &
        'maximal poloidal mode number')
      call h5_add(h5id_root, trim(adjustl(grp_name)) // '/kilca_pol_mode', this%kilca_pol_mode, &
        'poloidal mode number of resonance to be considered in pre- and post-processing')
      call h5_add(h5id_root, trim(adjustl(grp_name)) // '/kilca_pol_modes', &
        this%kilca_pol_modes, lbound(this%kilca_pol_modes), ubound(this%kilca_pol_modes), &
        'full range of poloidal mode numbers')
      call h5_add(h5id_root, trim(adjustl(grp_name)) // '/kilca_vac_coeff', &
        this%kilca_vac_coeff, lbound(this%kilca_vac_coeff), ubound(this%kilca_vac_coeff), &
        'coefficient of modified Bessel function of first kind', 'G')
      call h5_close(h5id_root)
    end if
  end subroutine vac_write

  subroutine read_Bnvac_Nemov(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bnvac_R, Bnvac_Z)
    use input_files, only: pfile
    use mephit_conf, only: conf
    use mephit_util, only: imun, pi, linspace
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    integer :: fid, nphi, idum, kR, kphi, kZ
    real(dp) :: B_R, B_phi, B_Z
    complex(dp), allocatable :: fourier_basis(:)

    call read_field_input
    open(newunit = fid, file = pfile, status = 'old', action = 'read', form = 'formatted')
    read (fid, *) nR, nphi, nZ, idum
    read (fid, *) Rmin, Rmax
    read (fid, *) ! phimin, phimax
    read (fid, *) Zmin, Zmax
    ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
    allocate(fourier_basis(nphi - 1))
    fourier_basis(:) = exp(-imun * conf%n * linspace(0d0, 2d0 * pi, nphi - 1, 0, 1)) &
      / dble(nphi - 1)
    allocate(Bnvac_R(nR, nZ), Bnvac_Z(nR, nZ))
    Bnvac_R = (0d0, 0d0)
    Bnvac_Z = (0d0, 0d0)
    do kR = 1, nR
      do kphi = 1, nphi
        do kZ = 1, nZ
          read (fid, *) B_R, B_phi, B_Z
          ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
          if (kphi /= nphi) then
            Bnvac_R(kR, kZ) = Bnvac_R(kR, kZ) + fourier_basis(kphi) * B_R
            Bnvac_Z(kR, kZ) = Bnvac_Z(kR, kZ) + fourier_basis(kphi) * B_Z
          end if
        end do
      end do
    end do
    close(fid)
    deallocate(fourier_basis)
  end subroutine read_Bnvac_Nemov

  subroutine read_Bnvac_GPEC(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bnvac_R, Bnvac_Z)
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, nf90_inq_varid, &
      nf90_inquire_dimension, nf90_get_var, nf90_close
    use mephit_conf, only: conf, logger
    use mephit_mesh, only: equil
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    character(len = 1024) :: filename
    logical :: file_exists
    integer :: ncid_file, ncid
    real(dp), dimension(:), allocatable :: R, Z
    real(dp), dimension(:, :, :), allocatable :: Bn_R, Bn_Z, Bnplas_R, Bnplas_Z

    write (filename, '("gpec_cylindrical_output_n", i0, ".nc")') conf%n
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) then
      write (logger%msg, '("File ", a, " not found, cannot read vacuum perturbation field.")') &
        trim(filename)
      if (logger%err) call logger%write_msg
      error stop
    end if
    call check_error("nf90_open", nf90_open(filename, nf90_nowrite, ncid_file))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "R", ncid))
    call check_error("nf90_inquire_dimension", &
      nf90_inquire_dimension(ncid_file, ncid, len = nR))
    call check_error("nf90_inq_dimid", nf90_inq_dimid(ncid_file, "z", ncid))
    call check_error("nf90_inquire_dimension", &
      nf90_inquire_dimension(ncid_file, ncid, len = nZ))
    allocate(R(nR), Z(nZ), &
      Bn_R(nR, nZ, 0:1), Bn_Z(nR, nZ, 0:1), Bnplas_R(nR, nZ, 0:1), Bnplas_Z(nR, nZ, 0:1))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "R", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, R))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "z", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Z))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "b_r", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Bn_R))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "b_z", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Bn_Z))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "b_r_plasma", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Bnplas_R))
    call check_error("nf90_inq_varid", nf90_inq_varid(ncid_file, "b_z_plasma", ncid))
    call check_error("nf90_get_var", nf90_get_var(ncid_file, ncid, Bnplas_Z))
    call check_error("nf90_close", nf90_close(ncid_file))
    ! m to cm
    Rmin = 1d2 * R(1)
    Rmax = 1d2 * R(nR)
    ! GPEC omits offset
    Zmin = 1d2 * Z(1) + equil%zmid
    Zmax = 1d2 * Z(nZ) + equil%zmid
    allocate(Bnvac_R(nR, nZ), Bnvac_Z(nR, nZ))
    ! T to G, factor 1/2 from Fourier series, complex conjugate
    Bnvac_R(:, :) = 0.5d4 * cmplx(Bn_R(:, :, 0) - Bnplas_R(:, :, 0), &
      -(Bn_R(:, :, 1) - Bnplas_R(:, :, 1)), dp)
    Bnvac_Z(:, :) = 0.5d4 * cmplx(Bn_Z(:, :, 0) - Bnplas_Z(:, :, 0), &
      -(Bn_Z(:, :, 1) - Bnplas_Z(:, :, 1)), dp)
    deallocate(R, Z, Bn_R, Bn_Z, Bnplas_R, Bnplas_Z)

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
  end subroutine read_Bnvac_GPEC

  subroutine compute_Bnvac(Bn)
    use coil_tools, only: read_currents, read_Bnvac_Fourier
    use mephit_conf, only: conf, logger, vac_src_nemov, vac_src_gpec, vac_src_fourier
    type(RT0_t), intent(inout) :: Bn
    integer :: nR, nZ
    real(dp) :: Rmin, Rmax, Zmin, Zmax
    real(dp), dimension(:), allocatable :: Ic
    complex(dp), dimension(:, :), allocatable :: Bn_R, Bn_Z

    ! initialize vacuum field
    select case (conf%vac_src)
    case (vac_src_nemov)
      call read_Bnvac_Nemov(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    case (vac_src_gpec)
      call read_Bnvac_GPEC(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    case (vac_src_fourier)
      call read_currents(conf%currents_file, Ic)
      call read_Bnvac_Fourier(conf%coil_file, conf%n, conf%Biot_Savart_prefactor * Ic, &
        nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
      deallocate(Ic)
    case default
      write (logger%msg, '("unknown vacuum field source selection", i0)') conf%vac_src
      if (logger%err) call logger%write_msg
      error stop
    end select
    call vector_potential_single_mode(conf%n, nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    deallocate(Bn_R, Bn_Z)
    ! project to finite elements
    Bn%DOF = (0d0, 0d0)
    call RT0_project_pol_comp(Bn, project_splined)
    call RT0_tor_comp_from_zero_div(Bn)

  contains
    function project_splined(ktri, weight, R, Z, n_f, f)
      use mephit_mesh, only: field_cache_t
      complex(dp) :: project_splined
      integer, intent(in) :: ktri
      real(dp), intent(in) :: weight, R, Z, n_f(:)
      type(field_cache_t), intent(in) :: f
      complex(dp) :: B(3)

      call spline_bn(conf%n, R, Z, B(1), B(2), B(3))
      project_splined = weight * sum(B * n_f)
    end function project_splined
  end subroutine compute_Bnvac

  subroutine generate_vacfield(vac)
    use bdivfree_mod, only: Rpoi, Zpoi, ipoint, AZnRe, AZnIm, ARnRe, ARnIm
    use mephit_conf, only: conf
    type(vac_t), intent(inout) :: vac

    if (conf%kilca_scale_factor /= 0) then
      call compute_kilca_vac_coeff(vac)
      call compute_kilca_vacuum(vac%Bn)
      call debug_Bmnvac_kilca
      call debug_RT0(vac%Bn)
    else
      if (conf%nonres) then
        call compute_Bn_nonres(vac%Bn)
      else
        call compute_Bnvac(vac%Bn)
        call debug_Bnvac_rectplot
        call debug_B0_rectplot
        call debug_Bmnvac
        call debug_RT0(vac%Bn)
        call debug_fouriermodes
        ! deallocate arrays allocated in vector_potential_single_mode
        if (allocated(Rpoi)) deallocate(Rpoi)
        if (allocated(Zpoi)) deallocate(Zpoi)
        if (allocated(ipoint)) deallocate(ipoint)
        if (allocated(AZnRe)) deallocate(AZnRe)
        if (allocated(AZnIm)) deallocate(AZnIm)
        if (allocated(ARnRe)) deallocate(ARnRe)
        if (allocated(ARnIm)) deallocate(ARnIm)
      end if
    end if
  end subroutine generate_vacfield

  subroutine debug_Bnvac_rectplot
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use bdivfree_mod, only: nR, nZ, Rpoi, Zpoi
    use mephit_conf, only: conf, datafile
    character(len=*), parameter :: dataset = 'vac/Bn'
    integer :: kR, kZ
    integer(HID_T) :: h5id_root
    complex(dp), dimension(nR, nZ) :: Bn_R, Bn_phi, Bn_Z

    do kZ = 1, nZ
      do kR = 1, nR
        call spline_bn(conf%n, Rpoi(kR), Zpoi(kZ), Bn_R(kR, kZ), Bn_phi(kR, kZ), Bn_Z(kR, kZ))
      end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/rect_R', Rpoi, lbound(Rpoi), ubound(Rpoi), &
      comment = 'R coordinate of rectangular grid', unit = 'cm')
    call h5_add(h5id_root, dataset // '/rect_Z', Zpoi, lbound(Zpoi), ubound(Zpoi), &
      comment = 'Z coordinate of rectangular grid', unit = 'cm')
    call h5_add(h5id_root, dataset // '/rect_comp_R', Bn_R, &
      lbound(Bn_R), ubound(Bn_R), unit = 'Mx', &
      comment = 'R component of magnetic field (vacuum) on rectangular grid')
    call h5_add(h5id_root, dataset // '/rect_comp_phi', Bn_phi, &
      lbound(Bn_phi), ubound(Bn_phi), unit = 'Mx', &
      comment = 'physical phi component of magnetic field (vacuum) on rectangular grid')
    call h5_add(h5id_root, dataset // '/rect_comp_Z', Bn_Z, &
      lbound(Bn_Z), ubound(Bn_Z), unit = 'Mx', &
      comment = 'Z component of magnetic field (vacuum) on rectangular grid')
    call h5_close(h5id_root)
  end subroutine debug_Bnvac_rectplot

  subroutine debug_B0_rectplot
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use bdivfree_mod, only: nR, nZ, Rpoi, Zpoi
    use mephit_conf, only: datafile
    character(len = *), parameter :: dataset = 'equil'
    integer :: kR, kZ
    integer(HID_T) :: h5id_root
    real(dp) :: dum
    real(dp), dimension(nR, nZ) :: B0_R, B0_phi, B0_Z

    do kZ = 1, nZ
      do kR = 1, nR
        call field(Rpoi(kR), 0d0, Zpoi(kZ), B0_R(kR, kZ), B0_phi(kR, kZ), B0_Z(kR, kZ), &
          dum, dum, dum, dum, dum, dum, dum, dum, dum)
      end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/rect_R', Rpoi, lbound(Rpoi), ubound(Rpoi), &
      comment = 'R coordinate of rectangular grid', unit = 'cm')
    call h5_add(h5id_root, dataset // '/rect_Z', Zpoi, lbound(Zpoi), ubound(Zpoi), &
      comment = 'Z coordinate of rectangular grid', unit = 'cm')
    call h5_add(h5id_root, dataset // '/B0_R', B0_R, &
      lbound(B0_R), ubound(B0_R), unit = 'Mx', &
      comment = 'R component of magnetic field (equilibrium) on rectangular grid')
    call h5_add(h5id_root, dataset // '/B0_phi', B0_phi, &
      lbound(B0_phi), ubound(B0_phi), unit = 'Mx', &
      comment = 'physical phi component of magnetic field (equilibrium) on rectangular grid')
    call h5_add(h5id_root, dataset // '/B0_Z', B0_Z, &
      lbound(B0_Z), ubound(B0_Z), unit = 'Mx', &
      comment = 'Z component of magnetic field (equilibrium) on rectangular grid')
    call h5_close(h5id_root)
  end subroutine debug_B0_rectplot

  subroutine debug_Bmnvac
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_mesh, only: equil, mesh, cache
    character(len = *), parameter :: dataset = 'vac/Bn'
    integer :: npol, kf, log2, kpol, k
    integer(HID_T) :: h5id_root
    complex(dp) :: Bn_R, Bn_Z, Bn_phi
    complex(dp), dimension(:), allocatable :: Bn_contradenspsi, Bn_n
    complex(dp), dimension(-conf%m_max:conf%m_max, mesh%nflux) :: Bmn_contradenspsi, Bmn_n

    npol = shiftl(1, cache%max_log2)
    allocate(Bn_contradenspsi(npol), Bn_n(npol))
    Bmn_contradenspsi(:, :) = (0d0, 0d0)
    Bmn_n(:, :) = (0d0, 0d0)
    do kf = 1, mesh%nflux
      log2 = cache%log2_kt_max(kf)
      npol = shiftl(1, log2)
      do kpol = 1, npol
        k = cache%kt_low(kf) + kpol
        associate (s => cache%sample_polmodes_half(k))
          call spline_bn(conf%n, s%R, s%Z, Bn_R, Bn_phi, Bn_Z)
          Bn_contradenspsi(kpol) = s%sqrt_g * (Bn_R * s%B0_Z - Bn_Z * s%B0_R) * s%R
          Bn_n(kpol) = (Bn_R * s%B0_Z - Bn_Z * s%B0_R) / hypot(s%B0_R, s%B0_Z)
        end associate
      end do
      cache%fft(log2)%samples(:) = Bn_contradenspsi(:npol)
      call cache%fft(log2)%apply(-conf%m_max, conf%m_max, Bmn_contradenspsi(:, kf))
      cache%fft(log2)%samples(:) = Bn_n(:npol)
      call cache%fft(log2)%apply(-conf%m_max, conf%m_max, Bmn_n(:, kf))
      Bmn_n(:, kf) = Bmn_n(:, kf) * equil%cocos%sgn_dpsi
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/comp_psi_contravar_dens', Bmn_contradenspsi, &
      lbound(Bmn_contradenspsi), ubound(Bmn_contradenspsi), unit = 'Mx', &
      comment = 'poloidal modes of contravariant density component of vacuum perturbation')
    call h5_add(h5id_root, dataset // '/comp_n', Bmn_n, &
      lbound(Bmn_n), ubound(Bmn_n), unit = 'G', &
      comment = 'poloidal modes of perpendicular component of vacuum perturbation')
    call h5_close(h5id_root)
  end subroutine debug_Bmnvac

  subroutine compute_Bn_nonres(Bn)
    use mephit_conf, only: conf
    use mephit_mesh, only: mesh, cache
    type(RT0_t), intent(inout) :: Bn
    integer :: kedge, k
    complex(dp) :: Bnpsi

    ! fluxes in poloidal direction are set to zero
    Bn%DOF = (0d0, 0d0)
    do kedge = 1, mesh%npoint - 1
      do k = 1, mesh%GL_order
        associate (f => cache%edge_fields(k, kedge))
          Bnpsi = -mesh%R_O * f%B0(2) / mesh%GL_R(k, kedge)
          Bn%DOF(kedge) = Bn%DOF(kedge) + mesh%GL_weights(k) * &
            Bnpsi * (mesh%edge_R(kedge) ** 2 + mesh%edge_Z(kedge) ** 2) / &
            (f%B0(1) * mesh%edge_R(kedge) + f%B0(3) * mesh%edge_Z(kedge))
        end associate
      end do
    end do
    call RT0_tor_comp_from_zero_div(Bn)
    if (conf%quad_avg) call avg_flux_on_quad(Bn)
  end subroutine compute_Bn_nonres

  subroutine avg_flux_on_quad(elem)
    use mephit_util, only: imun
    use mephit_mesh, only: mesh
    type(RT0_t), intent(inout) :: elem

    integer :: kf, kt, kedge, ktri1, ktri2
    complex(dp) :: tor_flux_avg, tor_flux_diff

    do kf = 2, mesh%nflux
      do kt = 1, mesh%kt_max(kf), 2
        kedge = mesh%npoint + mesh%kt_low(kf) + kt  ! every odd edge
        ktri1 = mesh%kt_low(kf) + kt
        ktri2 = mesh%kt_low(kf) + kt + 1
        tor_flux_avg = 0.5d0 * (elem%comp_phi(ktri2) * mesh%area(ktri2) + &
          elem%comp_phi(ktri1) * mesh%area(ktri1))
        tor_flux_diff = 0.5d0 * (elem%comp_phi(ktri2) * mesh%area(ktri2) - &
          elem%comp_phi(ktri1) * mesh%area(ktri1))
        ! DoF points in the direction opposite to the difference above
        elem%DOF(kedge) = elem%DOF(kedge) + imun * mesh%n * tor_flux_diff
        elem%comp_phi(ktri1) = tor_flux_avg / mesh%area(ktri1)
        elem%comp_phi(ktri2) = tor_flux_avg / mesh%area(ktri2)
      end do
    end do
  end subroutine avg_flux_on_quad

  ! calculate resonant vacuum perturbation
  subroutine compute_kilca_vacuum(Bn)
    use mephit_mesh, only: mesh
    type(RT0_t), intent(inout) :: Bn
    integer :: kedge, k
    real(dp) :: rho, theta
    complex(dp) :: B_R, B_phi, B_Z

    Bn%DOF = (0d0, 0d0)
    do kedge = 1, mesh%nedge
      do k = 1, mesh%GL_order
        rho = hypot(mesh%GL_R(k, kedge) - mesh%R_O, mesh%GL_Z(k, kedge) - mesh%Z_O)
        theta = atan2(mesh%GL_Z(k, kedge) - mesh%Z_O, mesh%GL_R(k, kedge) - mesh%R_O)
        call kilca_vacuum(mesh%n, vac%kilca_pol_modes, mesh%R_O, rho, theta, B_R, B_phi, B_Z)
        Bn%DOF(kedge) = Bn%DOF(kedge) + mesh%GL_weights(k) * &
          (B_R * mesh%edge_Z(kedge) - B_Z * mesh%edge_R(kedge)) * mesh%GL_R(k, kedge)
      end do
    end do
    ! toroidal flux via zero divergence
    call RT0_tor_comp_from_zero_div(Bn)
  end subroutine compute_kilca_vacuum

  !> Calculate the vacuum perturbation field in cylindrical coordinates from the Fourier
  !> series of all given modes.
  !>
  !> @param tor_mode toroidal mode number, scaled by mephit_conf::kilca_scale_factor
  !> (usually mephit_conf::n)
  !> @param pol_modes array of poloidal mode numbers
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> mephit_util::g_eqdsk::rcentr)
  !> @param r radial distance \f$ r \f$ from magnetic axis
  !> @param theta geometrical poloidal angle \f$ theta \f$ (coinciding with symmetry flux
  !> coordinates' poloidal angle in this geometry)
  !> @param B_R physical component \f$ B_{R} (r, theta, n) \f$ of the vacuum perturbation
  !> field
  !> @param B_phi physical component \f$ B_{(\varphi)} (r, theta, n) \f$ of the vacuum
  !> perturbation field
  !> @param B_Z physical component \f$ B_{Z} (r, theta, n) \f$ of the vacuum perturbation
  !> field
  subroutine kilca_vacuum(tor_mode, pol_modes, R_0, r, theta, B_R, B_phi, B_Z)
    use mephit_util, only: imun, straight_cyl2bent_cyl
    integer, intent(in) :: tor_mode, pol_modes(1:)
    real(dp), intent(in) :: R_0, r, theta
    complex(dp), intent(out) :: B_R, B_phi, B_Z
    complex(dp) :: B_rad, B_pol, B_tor, temp_B_rad, temp_B_pol, temp_B_tor
    complex(dp), dimension(size(pol_modes)) :: fourier_basis
    integer :: k

    B_rad = (0d0, 0d0)
    B_pol = (0d0, 0d0)
    B_tor = (0d0, 0d0)
    fourier_basis = exp(imun * pol_modes * theta)
    do k = 1, ubound(pol_modes, 1)
      call kilca_vacuum_fourier(tor_mode, pol_modes(k), R_0, r, &
        vac%kilca_vac_coeff(abs(pol_modes(k))), temp_B_rad, temp_B_pol, temp_B_tor)
      B_rad = B_rad + temp_B_rad * fourier_basis(k)
      B_pol = B_pol + temp_B_pol * fourier_basis(k)
      B_tor = B_tor + temp_B_tor * fourier_basis(k)
    end do
    call straight_cyl2bent_cyl(B_rad, B_pol, B_tor, theta, B_R, B_phi, B_Z)
  end subroutine kilca_vacuum

  subroutine compute_kilca_vac_coeff(vac)
    use h5lt, only: h5ltget_dataset_info_f
    use hdf5_tools, only: hid_t, hsize_t, size_t, h5_open, h5_get, h5_exists, h5_close
    use mephit_conf, only: conf, logger, cmplx_fmt
    use mephit_mesh, only: mesh
    type(vac_t), intent(inout) :: vac
    character(len = *), parameter :: grp_ptrn = '("output/postprocessor", i0)'
    character(len = 32) :: grp_name  ! /output/postprocessor1234567890
    integer(hid_t) :: h5id_root
    integer(hsize_t) :: dims(2)
    integer(size_t) :: type_size
    integer :: k, m, n, type_id, h5err
    real(dp) :: r2dum(1, 2), r
    real(dp), allocatable :: kilca_r(:, :), kilca_Bz(:, :)
    complex(dp) :: B_rad, B_pol, B_tor, Bz

    call h5_open(trim(adjustl(conf%kilca_vac_output)), h5id_root)
    k = 0
    do while (.true.)
      k = k + 1
      write (grp_name, grp_ptrn) k
      if (.not. h5_exists(h5id_root, grp_name)) exit
      call h5_get(h5id_root, trim(grp_name) // '/mode', r2dum)
      m = int(r2dum(1, 1))
      n = int(r2dum(1, 2))
      if (n /= conf%n .or. m > vac%m_max .or. m < vac%m_min) cycle
      call h5ltget_dataset_info_f(h5id_root, trim(grp_name) // '/Bz', &
        dims, type_id, type_size, h5err)
      allocate(kilca_r(dims(1), 1), kilca_Bz(dims(1), dims(2)))
      call h5_get(h5id_root, trim(grp_name) // '/r', kilca_r)
      call h5_get(h5id_root, trim(grp_name) // '/Bz', kilca_Bz)
      r = kilca_r(dims(1), 1)
      if (dims(2) > 1) then
        Bz = cmplx(kilca_Bz(dims(1), 1), kilca_Bz(dims(1), 2), dp)
      else
        Bz = cmplx(kilca_Bz(dims(1), 1), 0d0, dp)
      end if
      call kilca_vacuum_fourier(mesh%n, m, mesh%R_O, r, (1d0, 0d0), B_rad, B_pol, B_tor)
      vac%kilca_vac_coeff(m) = Bz / B_tor
      deallocate(kilca_r, kilca_Bz)
    end do
    call h5_close(h5id_root)
    logger%msg = 'effective vacuum perturbation field coefficients:'
    if (logger%info) call logger%write_msg
    do m = mesh%m_res_min, mesh%m_res_max
      write (logger%msg, '("vac%kilca_vac_coeff(", i0, ") = ", ' // cmplx_fmt // ')') m, &
        vac%kilca_vac_coeff(m)
      if (logger%info) call logger%write_msg
    end do
  end subroutine compute_kilca_vac_coeff

  !> Calculate the Fourier coefficient of the vacuum perturbation field for a given
  !> toroidal-poloidal mode.
  !>
  !> @param tor_mode toroidal mode number, scaled by mephit_conf::kilca_scale_factor
  !> (usually mephit_conf::n)
  !> @param pol_mode poloidal mode number
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> mephit_util::g_eqdsk::rcentr)
  !> @param r radial distance \f$ r \f$ from magnetic axis
  !> @param vac_coeff coefficient of modified Bessel functions (integration constant)
  !> @param B_rad physical component \f$ B_{r} (r, m, n) \f$ of the vacuum perturbation
  !> field
  !> @param B_pol physical component \f$ B_{(\theta)} (r, m, n) \f$ of the vacuum
  !> perturbation field
  !> @param B_tor physical component \f$ B_{z} (r, m, n) \f$ of the vacuum perturbation
  !> field
  subroutine kilca_vacuum_fourier(tor_mode, pol_mode, R_0, r, vac_coeff, &
    B_rad, B_pol, B_tor)
    use fgsl, only: fgsl_double, fgsl_int, fgsl_success, fgsl_sf_bessel_icn_array
    use mephit_conf, only: logger
    use mephit_util, only: imun
    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r
    complex(dp), intent(in) :: vac_coeff
    complex(dp), intent(out) :: B_rad, B_pol, B_tor
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(abs(pol_mode)-1, abs(pol_mode)+1, k_z_r, I_m)
    if (status /= fgsl_success .and. logger%err) then
      write (logger%msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
      call logger%write_msg
    end if
    B_rad = 0.5d0 * (I_m(-1) + I_m(1)) * vac_coeff
    B_pol = imun * pol_mode / k_z_r * I_m(0) * vac_coeff
    B_tor = imun * I_m(0) * vac_coeff
  end subroutine kilca_vacuum_fourier

  subroutine debug_Bmnvac_kilca
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: datafile
    use mephit_mesh, only: fs_half, mesh
    character(len = *), parameter :: dataset = 'Bmnvac'
    integer(HID_T) :: h5id_root
    complex(dp), dimension(mesh%nflux, vac%m_min:vac%m_max) :: &
      B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    integer :: kf, m

    do m = vac%m_min, vac%m_max
      do kf = 1, mesh%nflux
        call kilca_vacuum_fourier(mesh%n, -m, mesh%R_O, fs_half%rad(kf), &
          vac%kilca_vac_coeff(m), B_rad_neg(kf, m), B_pol_neg(kf, m), B_tor_neg(kf, m))
        call kilca_vacuum_fourier(mesh%n, m, mesh%R_O, fs_half%rad(kf), &
          vac%kilca_vac_coeff(m), B_rad_pos(kf, m), B_pol_pos(kf, m), B_tor_pos(kf, m))
      end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/comp_rad_neg', B_rad_neg, lbound(B_rad_neg), ubound(B_rad_neg), &
      comment = 'radial component of negative poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_add(h5id_root, dataset // '/comp_pol_neg', B_pol_neg, lbound(B_pol_neg), ubound(B_pol_neg), &
      comment = 'poloidal component of negative poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_add(h5id_root, dataset // '/comp_tor_neg', B_tor_neg, lbound(B_tor_neg), ubound(B_tor_neg), &
      comment = 'toroidal component of negative poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_add(h5id_root, dataset // '/comp_rad_pos', B_rad_pos, lbound(B_rad_pos), ubound(B_rad_pos), &
      comment = 'radial component of positive poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_add(h5id_root, dataset // '/comp_pol_pos', B_pol_pos, lbound(B_pol_pos), ubound(B_pol_pos), &
      comment = 'poloidal component of positive poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_add(h5id_root, dataset // '/comp_tor_pos', B_tor_pos, lbound(B_tor_pos), ubound(B_tor_pos), &
      comment = 'toroidal component of positive poloidal mode of KiLCA vacuum at half-grid steps', unit = 'G')
    call h5_close(h5id_root)
  end subroutine debug_Bmnvac_kilca

  subroutine debug_RT0(Bn)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use mephit_conf, only: conf, datafile
    use mephit_mesh, only: mesh, fs_half, cache
    type(RT0_t), intent(in) :: Bn
    character(len = *), parameter :: dataset = 'debug_RT0'
    integer(HID_T) :: h5id_root
    integer :: kf, log2, kpol, k
    complex(dp) :: B_interp(3)
    complex(dp), dimension(size(cache%sample_polmodes_half)) :: &
      B_R, B_Z, B_phi, B_R_interp, B_Z_interp, B_phi_interp

    do kf = 1, mesh%nflux
      log2 = cache%log2_kt_max(kf)
      do kpol = 1, shiftl(1, log2)
        k = cache%kt_low(kf) + kpol
        associate (s => cache%sample_polmodes_half(k))
          if (conf%kilca_scale_factor /= 0) then
            call kilca_vacuum(mesh%n, vac%kilca_pol_modes, mesh%R_O, &
              fs_half%rad(kf), s%theta, B_R(k), B_phi(k), B_Z(k))
          else
            call spline_bn(conf%n, s%R, s%Z, B_R(k), B_phi(k), B_Z(k))
          end if
          call RT0_interp(Bn, s%ktri, s%R, s%Z, B_interp)
          B_R_interp(k) = B_interp(1)
          B_phi_interp(k) = B_interp(2)
          B_Z_interp(k) = B_interp(3)
        end associate
      end do
    end do
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/Bn_R', B_R, lbound(B_R), ubound(B_R), &
      comment = 'R component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_add(h5id_root, dataset // '/Bn_phi', B_phi, lbound(B_phi), ubound(B_phi), &
      comment = 'phi component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_add(h5id_root, dataset // '/Bn_Z', B_Z, lbound(B_Z), ubound(B_Z), &
      comment = 'Z component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_add(h5id_root, dataset // '/Bn_R_RT0', B_R_interp, lbound(B_R_interp), ubound(B_R_interp), &
      comment = 'RT0-interpolated R component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_add(h5id_root, dataset // '/Bn_phi_RT0', B_phi_interp, lbound(B_phi_interp), ubound(B_phi_interp), &
      comment = 'RT0-interpolated phi component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_add(h5id_root, dataset // '/Bn_Z_RT0', B_Z_interp, lbound(B_Z_interp), ubound(B_Z_interp), &
      comment = 'RT0-interpolated Z component of magnetic field (vacuum perturbation)', unit = 'G')
    call h5_close(h5id_root)
  end subroutine debug_RT0

  subroutine debug_fouriermodes
    use mephit_conf, only: conf, datafile, logger
    use mephit_util, only: resample1d, linspace, imun
    use mephit_mesh, only: equil
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    logical :: file_exists
    integer :: fid, ntor, mpol, nlabel, nsqpsi, k
    integer(HID_T) :: h5id_root
    real(dp) :: flabel_min, flabel_max, ddum, psimax, phimax, sgn_dpsi
    real(dp), allocatable :: qsaf(:), rsmall(:), psisurf(:), rq_eqd(:), psi_n(:)
    complex(dp), allocatable :: z3dum(:, :, :), Amn_theta(:, :, :), Bmn_contradenspsi(:, :)
    character(len = *), parameter :: ptrn_nsqpsi = 'nrad = ', ptrn_psimax = 'psimax', &
      ptrn_phimax = 'phimax', dataset = 'debug_fouriermodes'
    character(len = 1024) :: line

    inquire(file = 'amn.dat', exist = file_exists)
    if (.not. file_exists) return
    inquire(file = 'equil_r_q_psi.dat', exist = file_exists)
    if (.not. file_exists) return
    logger%msg = 'Files amn.dat and equil_r_q_psi.dat found, performing Fourier mode comparison.'
    if (logger%info) call logger%write_msg
    open(newunit = fid, form = 'unformatted', file = 'amn.dat', action = 'read')
    read (fid) ntor, mpol, nlabel, flabel_min, flabel_max
    allocate(z3dum(-mpol:mpol, ntor, nlabel), Amn_theta(-mpol:mpol, ntor, nlabel))
    read (fid) z3dum, Amn_theta
    close(fid)
    deallocate(z3dum)
    open(newunit = fid, form = 'formatted', file = 'equil_r_q_psi.dat', action = 'read')
    read (fid, '(a)') line
    call extract(line, ptrn_nsqpsi)
    read (line, *) nsqpsi
    read (fid, '(a)') line
    call extract(line, ptrn_psimax)
    read (line, *) psimax
    call extract(line, ptrn_phimax)
    read (line, *) phimax
    read (fid, *)
    allocate(qsaf(nsqpsi), psisurf(nsqpsi), rsmall(nsqpsi))
    do k = 1, nsqpsi
      read (fid, *) ddum, qsaf(k), psisurf(k), ddum, ddum, rsmall(k)
    end do
    close(fid)
    sgn_dpsi = sign(1d0, psimax)
    allocate(rq_eqd(nlabel), psi_n(nlabel), Bmn_contradenspsi(-mpol:mpol, nlabel))
    rq_eqd(:) = linspace(abs(flabel_min), abs(flabel_max), nlabel, 0, 0)
    call resample1d(rsmall * abs(qsaf), psisurf / psimax, rq_eqd, psi_n, 3)
    ! if qsaf does not have the expected sign, theta points in the wrong direction,
    ! and we have to reverse the index m and the overall sign
    if (equil%cocos%sgn_q * qsaf(nsqpsi) < 0d0) then
      Bmn_contradenspsi(:, :) = -imun * conf%n * sgn_dpsi * Amn_theta(mpol:-mpol:-1, conf%n, :)
    else
      Bmn_contradenspsi(:, :) = imun * conf%n * sgn_dpsi * Amn_theta(:, conf%n, :)
    end if
    call h5_open_rw(datafile, h5id_root)
    call h5_create_parent_groups(h5id_root, dataset // '/')
    call h5_add(h5id_root, dataset // '/m_max', mpol, 'maximal absolute poloidal mode number')
    call h5_add(h5id_root, dataset // '/psi_n', psi_n, lbound(psi_n), ubound(psi_n), &
      comment = 'normalized poloidal flux of interpolation points', unit = '1')
    call h5_add(h5id_root, dataset // '/comp_psi_contravar_dens', Bmn_contradenspsi, &
      lbound(Bmn_contradenspsi), ubound(Bmn_contradenspsi), unit = 'Mx', &
      comment = 'normalized poloidal flux of interpolation points')
    call h5_close(h5id_root)
    deallocate(Amn_theta, qsaf, psisurf, rsmall, rq_eqd, psi_n, Bmn_contradenspsi)

  contains
    subroutine extract(string, preceding)
      character(len = *), intent(inout) :: string
      character(len = *), intent(in) :: preceding
      integer :: i
      i = index(string, preceding) + len(preceding) - 1
      string(:i) = ' '
      string = adjustl(string)
    end subroutine extract
  end subroutine debug_fouriermodes

end module mephit_pert
