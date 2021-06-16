module magdif_pert

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: L1_init, L1_deinit, L1_read, L1_write, RT0_init, RT0_deinit, RT0_read, &
       RT0_write, RT0_interp, RT0_check_div_free, RT0_check_redundant_edges, RT0_triplot, &
       RT0_rectplot, RT0_poloidal_modes, vec_polmodes_init, vec_polmodes_deinit, &
       vec_polmodes_read, vec_polmodes_write, AUG_coils_read, AUG_coils_write_Nemov, &
       AUG_coils_read_Nemov, AUG_coils_write_GPEC, AUG_coils_read_GPEC, AUG_coils_write_Fourier, &
       read_currents_Nemov, Biot_Savart_sum_coils, write_Bvac_Nemov, generate_vacfield, compute_Bnvac, &
       compute_Bn_nonres, avg_flux_on_quad, &
       compute_kilca_vacuum, kilca_vacuum, compute_kilca_vac_coeff, kilca_vacuum_fourier, &
       check_kilca_vacuum, check_RT0, debug_fouriermodes, debug_Bnvac_rectplot, debug_Bmnvac

  type, public :: L1_t
     !> Number of points on which the L1 are defined
     integer :: npoint

     !> Degrees of freedom, given in base units.
     complex(dp), allocatable :: DOF(:)
  end type L1_t

  type, public :: RT0_t
     !> Number of triangles on which the RT0 elements are defined
     integer :: ntri

     !> Edge fluxes \f$ R \vec{v} \cdot \vec{n} \f$, given in base units of
     !> \f$ \vec{v} \f$ times cm^2.
     !>
     !> Values are stored seprately for each triangle, i.e. twice per edge. The first index
     !> refers to the edge and the second index refers to the triangle.
     complex(dp), allocatable :: DOF(:, :)

     !> Physical toroidal component \f$ v_{n (\phi)} \f$ in base units of \f$ \vec{v} \f$.
     !>
     !> Values are taken at each triangle.
     complex(dp), allocatable :: comp_phi(:)
  end type RT0_t

  type, public :: vec_polmodes_t
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
  !> to \p outfile (e.g. #magdif_conf::presn_file), where line number corresponds to the
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

  subroutine RT0_init(this, ntri)
    type(RT0_t), intent(inout) :: this
    integer, intent(in) :: ntri

    call RT0_deinit(this)
    this%ntri = ntri
    allocate(this%DOF(3, ntri))
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

  subroutine RT0_interp(ktri, elem, R, Z, comp_R, comp_Z, comp_phi, &
       comp_R_dR, comp_R_dZ, comp_Z_dR, comp_Z_dZ)
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use magdif_mesh, only: mesh
    integer, intent(in) :: ktri
    type(RT0_t), intent(in) :: elem
    real(dp), intent(in) :: r, z
    complex(dp), intent(out) :: comp_R, comp_Z
    complex(dp), intent(out), optional :: comp_phi, comp_R_dR, comp_R_dZ, &
         comp_Z_dR, comp_Z_dZ
    integer :: nodes(3)
    real(dp) :: nan

    if (ktri <= 0) then
       nan = ieee_value(1d0, ieee_quiet_nan)
       comp_R = cmplx(nan, nan, dp)
       comp_Z = cmplx(nan, nan, dp)
       if (present(comp_phi)) then
          comp_phi = cmplx(nan, nan, dp)
       end if
       if (present(comp_R_dZ)) then
          comp_R_dZ = (0d0, 0d0)
       end if
       if (present(comp_Z_dR)) then
          comp_Z_dR = cmplx(nan, nan, dp)
       end if
       if (present(comp_Z_dZ)) then
          comp_Z_dZ = cmplx(nan, nan, dp)
       end if
       return
    end if
    nodes = mesh%tri_node(:, ktri)
    ! edge 1 lies opposite to knot 3, etc.
    comp_R = 0.5d0 / mesh%area(ktri) / R * ( &
         elem%DOF(1, ktri) * (R - mesh%node_R(nodes(3))) + &
         elem%DOF(2, ktri) * (R - mesh%node_R(nodes(1))) + &
         elem%DOF(3, ktri) * (R - mesh%node_R(nodes(2))))
    comp_Z = 0.5d0 / mesh%area(ktri) / R * ( &
         elem%DOF(1, ktri) * (Z - mesh%node_Z(nodes(3))) + &
         elem%DOF(2, ktri) * (Z - mesh%node_Z(nodes(1))) + &
         elem%DOF(3, ktri) * (Z - mesh%node_Z(nodes(2))))
    if (present(comp_phi)) then
       comp_phi = elem%comp_phi(ktri)
    end if
    if (present(comp_R_dR)) then
       comp_R_dR = 0.5d0 / mesh%area(ktri) / R ** 2 * ( &
            elem%DOF(1, ktri) * mesh%node_R(nodes(3)) + &
            elem%DOF(2, ktri) * mesh%node_R(nodes(1)) + &
            elem%DOF(3, ktri) * mesh%node_R(nodes(2)))
    end if
    if (present(comp_R_dZ)) then
       comp_R_dZ = (0d0, 0d0)
    end if
    if (present(comp_Z_dR)) then
       comp_Z_dR = -comp_Z / R
    end if
    if (present(comp_Z_dZ)) then
       comp_Z_dZ = sum(elem%DOF(:, ktri)) * 0.5d0 / mesh%area(ktri) / R
    end if
  end subroutine RT0_interp

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
  subroutine RT0_check_div_free(elem, n, rel_err, field_name)
    use magdif_conf, only: log
    use magdif_util, only: imun
    use magdif_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
    integer, intent(in) :: n
    real(dp), intent(in) :: rel_err
    character(len = *), intent(in) :: field_name

    integer :: ktri
    real(dp) :: div, abs_flux

    do ktri = 1, mesh%ntri
       abs_flux = sum(abs(elem%DOF(:, ktri))) + abs(imun * n * elem%comp_phi(ktri) * &
            mesh%area(ktri))
       if (abs_flux > 0d0) then
          div = abs((sum(elem%DOF(:, ktri)) + imun * n * elem%comp_phi(ktri) * &
               mesh%area(ktri))) / abs_flux
          if (div > rel_err) then
             write (log%msg, '("divergence of ", a, ' // &
                  '" above threshold in triangle ", i0, ": ", es24.16e3)') &
                  trim(field_name), ktri, div
              if (log%err) call log%write
              error stop
          end if
       end if
    end do
  end subroutine RT0_check_div_free

  subroutine RT0_check_redundant_edges(elem, name)
    use magdif_conf, only: log, cmplx_fmt
    use magdif_mesh, only: mesh
    type(RT0_t), intent(in) :: elem
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
       if (abs(real(elem%DOF(ke, ktri))) < small) then
          inconsistent = inconsistent .or. abs(real(elem%DOF(ke_adj, ktri_adj))) >= small
       else
          inconsistent = inconsistent .or. eps < abs(1d0 + &
               real(elem%DOF(ke_adj, ktri_adj)) / real(elem%DOF(ke, ktri)))
       end if
       if (abs(aimag(elem%DOF(ke, ktri))) < small) then
          inconsistent = inconsistent .or. abs(aimag(elem%DOF(ke_adj, ktri_adj))) >= small
       else
          inconsistent = inconsistent .or. eps < abs(1d0 + &
               aimag(elem%DOF(ke_adj, ktri_adj)) / aimag(elem%DOF(ke, ktri)))
       end if
       if (inconsistent) then
          write (log%msg, '("inconsistent redundant edges: ", ' // &
               'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ', ", ", ' // &
               'a, "(", i0, ", ", i0, ") = ", ' // cmplx_fmt // ')') &
               trim(name), ke, ktri, elem%DOF(ke, ktri), &
               trim(name), ke_adj, ktri_adj, elem%DOF(ke_adj, ktri_adj)
          if (log%err) call log%write
          error stop
       end if
    end do
  end subroutine RT0_check_redundant_edges

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

  pure function jacobian(kf, kt, r) result(metric_det)
    use magdif_mesh, only: equil, fs_half, mesh, B0phi_Omega
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: r
    real(dp) :: metric_det
    metric_det = equil%cocos%sgn_dpsi * fs_half%q(kf) * r / B0phi_Omega(mesh%kt_low(kf) + kt)
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
    end if
    call h5_close(h5id_root)
    if (allocated(comp_R)) deallocate(comp_R)
    if (allocated(comp_Z)) deallocate(comp_Z)
    if (allocated(comp_psi_contravar_dens)) deallocate(comp_psi_contravar_dens)
    if (allocated(comp_theta_covar)) deallocate(comp_theta_covar)
    if (allocated(rect_R)) deallocate(rect_R)
    if (allocated(rect_Z)) deallocate(rect_Z)
    if (allocated(rect_comp_R)) deallocate(rect_comp_R)
    if (allocated(rect_comp_phi)) deallocate(rect_comp_phi)
    if (allocated(rect_comp_Z)) deallocate(rect_comp_Z)
  end subroutine RT0_write

  subroutine RT0_triplot(elem, comp_R, comp_Z, comp_psi_contravar_dens, comp_theta_covar)
    use magdif_mesh, only: equil, mesh, B0R_Omega, B0Z_Omega
    type(RT0_t), intent(in) :: elem
    complex(dp), intent(out), dimension(:), allocatable :: comp_R, comp_Z, &
         comp_psi_contravar_dens, comp_theta_covar
    integer :: kf, kt, ktri
    real(dp) :: R, Z

    allocate(comp_R(mesh%ntri))
    allocate(comp_Z(mesh%ntri))
    allocate(comp_psi_contravar_dens(mesh%ntri))
    allocate(comp_theta_covar(mesh%ntri))
    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          R = mesh%R_Omega(ktri)
          Z = mesh%Z_Omega(ktri)
          call RT0_interp(ktri, elem, R, Z, comp_R(ktri), comp_Z(ktri))
          ! projection to contravariant psi component
          comp_psi_contravar_dens(ktri) = (comp_R(ktri) * B0Z_Omega(ktri) - &
               comp_Z(ktri) * B0R_Omega(ktri)) * R * jacobian(kf, kt, R)
          ! projection to covariant theta component
          comp_theta_covar(ktri) = equil%cocos%sgn_dpsi * (comp_R(ktri) * B0R_Omega(ktri) + &
               comp_Z(ktri) * B0Z_Omega(ktri)) * jacobian(kf, kt, R)
       end do
    end do
  end subroutine RT0_triplot

  subroutine RT0_rectplot(elem, rect_R, rect_Z, comp_R, comp_phi, comp_Z)
    use magdif_mesh, only: equil, mesh, point_location
    type(RT0_t), intent(in) :: elem
    real(dp), intent(out), dimension(:), allocatable :: rect_R, rect_Z
    complex(dp), intent(out), dimension(:, :), allocatable :: comp_R, comp_phi, comp_Z
    integer :: kw, kh, ktri

    allocate(rect_R(equil%nw))
    allocate(rect_Z(equil%nh))
    rect_R(:) = equil%R_eqd
    rect_Z(:) = equil%Z_eqd
    allocate(comp_R(equil%nw, equil%nh))
    allocate(comp_phi(equil%nw, equil%nh))
    allocate(comp_Z(equil%nw, equil%nh))
    do kw = 1, equil%nw
       do kh = 1, equil%nh
          ktri = point_location(rect_R(kw), rect_Z(kh))
          if (ktri > mesh%kt_low(1) .and. ktri <= mesh%ntri) then
             call RT0_interp(ktri, elem, rect_R(kw), rect_Z(kh), &
                  comp_R(kw, kh), comp_Z(kw, kh), comp_phi(kw, kh))
          else
             comp_R(kw, kh) = (0d0, 0d0)
             comp_phi(kw, kh) = (0d0, 0d0)
             comp_Z(kw, kh) = (0d0, 0d0)
          end if
       end do
    end do
  end subroutine RT0_rectplot

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
    use magdif_conf, only: conf
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
    if (conf%kilca_pol_mode /= 0) then
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

  subroutine RT0_poloidal_modes(elem, vec_polmodes)
    use magdif_conf, only: conf
    use magdif_util, only: imun, bent_cyl2straight_cyl
    use magdif_mesh, only: equil, mesh, sample_polmodes
    type(RT0_t), intent(in) :: elem
    type(vec_polmodes_t), intent(inout) :: vec_polmodes
    integer :: kf, kt, ktri, m
    complex(dp) :: fourier_basis(-vec_polmodes%m_max:vec_polmodes%m_max)
    complex(dp) :: comp_R, comp_Z, comp_phi, comp_rad, comp_n, comp_pol, comp_tor

    vec_polmodes%coeff_rad(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_n(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_pol(:, :) = (0d0, 0d0)
    vec_polmodes%coeff_tor(:, :) = (0d0, 0d0)
    associate (s => sample_polmodes, v => vec_polmodes)
      do kf = 1, mesh%nflux
         do kt = 1, mesh%kt_max(kf)
            ktri = mesh%kt_low(kf) + kt
            fourier_basis = [(exp(-imun * m * s%theta(ktri)), m = -v%m_max, v%m_max)]
            call RT0_interp(s%ktri(ktri), elem, s%R(ktri), s%Z(ktri), comp_R, comp_Z, comp_phi)
            if (conf%kilca_scale_factor /= 0) then
               call bent_cyl2straight_cyl(comp_R, comp_phi, comp_Z, &
                    s%theta(ktri), comp_rad, comp_pol, comp_tor)
               comp_n = comp_rad
            else
               comp_rad = s%sqrt_g(ktri) * (comp_R * s%B0_Z(ktri) - &
                    comp_Z * s%B0_R(ktri)) * s%R(ktri)
               comp_n = (comp_R * s%B0_Z(ktri) - comp_Z * s%B0_R(ktri)) &
                    / hypot(s%B0_R(ktri), s%B0_Z(ktri))
               comp_pol = comp_R * s%dR_dtheta(ktri) + comp_Z * s%dZ_dtheta(ktri)
               comp_tor = comp_phi
            end if
            v%coeff_rad(:, kf) = v%coeff_rad(:, kf) + comp_rad * fourier_basis
            v%coeff_n(:, kf) = v%coeff_n(:, kf) + comp_n * fourier_basis
            v%coeff_pol(:, kf) = v%coeff_pol(:, kf) + comp_pol * fourier_basis
            v%coeff_tor(:, kf) = v%coeff_tor(:, kf) + comp_tor * fourier_basis
         end do
         v%coeff_rad(:, kf) = v%coeff_rad(:, kf) / mesh%kt_max(kf)
         v%coeff_n(:, kf) = equil%cocos%sgn_dpsi * v%coeff_n(:, kf) / mesh%kt_max(kf)
         v%coeff_pol(:, kf) = v%coeff_pol(:, kf) / mesh%kt_max(kf)
         v%coeff_tor(:, kf) = v%coeff_tor(:, kf) / mesh%kt_max(kf)
      end do
    end associate
  end subroutine RT0_poloidal_modes

  subroutine AUG_coils_read(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), dimension(:, :, :), allocatable :: XYZ
    character(len = 8) :: filename
    integer :: fid, status, kc, ks
    real(dp), dimension(:, :), allocatable :: R_Z_phi

    ncoil = 8
    nwind = 5
    ! count number of lines = number of coil segments
    filename = 'Bu1n.asc'
    open(newunit = fid, file = directory // '/' // filename, status = 'old', &
         action = 'read', form = 'formatted')
    nseg = 0
    do
       read(fid, *, iostat = status)
       if (status /= 0) exit
       nseg = nseg + 1
    end do
    close(fid)
    ! allocate coordinates and read all coil data, starting with upper B coil set
    allocate(R_Z_phi(nseg, 3))
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, ncoil
       write (filename, '("Bu", i1, "n.asc")') kc
       open(newunit = fid, file = directory // '/' // filename, status = 'old', &
            action = 'read', form = 'formatted')
       do ks = 1, nseg
          read (fid, '(f8.6, 1x, f8.5, 1x, f8.5)') R_Z_phi(ks, 1), R_Z_phi(ks, 2), R_Z_phi(ks, 3)
       end do
       close(fid)
       XYZ(1, :, kc) = 1d2 * R_Z_phi(:, 1) * cos(R_Z_phi(:, 3))
       XYZ(2, :, kc) = 1d2 * R_Z_phi(:, 1) * sin(R_Z_phi(:, 3))
       XYZ(3, :, kc) = 1d2 * R_Z_phi(:, 2)
    end do
    do kc = 1, ncoil
       write (filename, '("Bl", i1, "n.asc")') kc
       open(newunit = fid, file = directory // '/' // filename, status = 'old', &
            action = 'read', form = 'formatted')
       do ks = 1, nseg
          read (fid, '(f8.6, 1x, f8.5, 1x, f8.5)') R_Z_phi(ks, 1), R_Z_phi(ks, 2), R_Z_phi(ks, 3)
       end do
       close(fid)
       XYZ(1, :, kc + ncoil) = 1d2 * R_Z_phi(:, 1) * cos(R_Z_phi(:, 3))
       XYZ(2, :, kc + ncoil) = 1d2 * R_Z_phi(:, 1) * sin(R_Z_phi(:, 3))
       XYZ(3, :, kc + ncoil) = 1d2 * R_Z_phi(:, 2)
    end do
    if (allocated(R_Z_phi)) deallocate(R_Z_phi)
  end subroutine AUG_coils_read

  subroutine AUG_coils_write_Nemov(directory, ncoil, nseg, XYZ)
    use magdif_conf, only: log
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks

    if (3 /= size(XYZ, 1)) then
       call log%msg_arg_size('AUG_coils_write_Nemov', '3', 'size(XYZ, 1)', 3, size(XYZ, 1))
       if (log%err) call log%write
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       call log%msg_arg_size('AUG_coils_write_Nemov', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       call log%msg_arg_size('AUG_coils_write_Nemov', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3))
       if (log%err) call log%write
       error stop
    end if
    open(newunit = fid, file = directory // '/co_asd.dd', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(1x, i6)') 2 * ncoil * (nseg + 1)
    do kc = 1, 2 * ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es17.9e2), 1x, es12.4e2, 1x, i3)') &
               XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc), 1.0, kc
       end do
       write (fid, '(3(1x, es17.9e2), 1x, es12.4e2, 1x, i3)') &
            XYZ(1, 1, kc), XYZ(2, 1, kc), XYZ(3, 1, kc), 0.0, kc
    end do
    close(fid)
  end subroutine AUG_coils_write_Nemov

  subroutine AUG_coils_read_Nemov(directory, ncoil, nseg, nwind, XYZ)
    use magdif_conf, only: log
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), allocatable, dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks, temp
    real(dp) :: rdum

    ncoil = 8
    nwind = 5
    open(newunit = fid, file = directory // '/co_asd.dd', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, *) temp
    nseg = temp / (2 * ncoil) - 1
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, 2 * ncoil
       do ks = 1, nseg
          read (fid, *) XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc), rdum, temp
          if (temp /= kc) then
             write (log%msg, '("Expected coil index ", i0, " in co_asd.dd at line ", ' // &
                  'i0, ", but got ", i0)') kc, (kc - 1) * (nseg + 1) + ks + 1, temp
             if (log%err) call log%write
             error stop
          end if
       end do
       read (fid, *)
    end do
    close(fid)
  end subroutine AUG_coils_read_Nemov

  subroutine AUG_coils_write_GPEC(directory, ncoil, nseg, nwind, XYZ)
    use magdif_conf, only: log
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks

    if (3 /= size(XYZ, 1)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', '3', 'size(XYZ, 1)', 3, size(XYZ, 1))
       if (log%err) call log%write
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3))
       if (log%err) call log%write
       error stop
    end if
    open(newunit = fid, file = directory // '/aug_bu.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(3(1x, i4), 1x, f7.2)') ncoil, 1, nseg + 1, dble(nwind)
    do kc = 1, ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es12.4e2))') &
               1d-2 * XYZ(1, ks, kc), 1d-2 * XYZ(2, ks, kc), 1d-2 * XYZ(3, ks, kc)
       end do
       write (fid, '(3(1x, es12.4e2))') &
            1d-2 * XYZ(1, 1, kc), 1d-2 * XYZ(2, 1, kc), 1d-2 * XYZ(3, 1, kc)
    end do
    close(fid)
    open(newunit = fid, file = directory // '/aug_bl.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    write (fid, '(3(1x, i4), 1x, f7.2)') ncoil, 1, nseg + 1, dble(nwind)
    do kc = ncoil + 1, 2 * ncoil
       do ks = 1, nseg
          write (fid, '(3(1x, es12.4e2))') &
               1d-2 * XYZ(1, ks, kc), 1d-2 * XYZ(2, ks, kc), 1d-2 * XYZ(3, ks, kc)
       end do
       write (fid, '(3(1x, es12.4e2))') &
            1d-2 * XYZ(1, 1, kc), 1d-2 * XYZ(2, 1, kc), 1d-2 * XYZ(3, 1, kc)
    end do
    close(fid)
  end subroutine AUG_coils_write_GPEC

  subroutine AUG_coils_read_GPEC(directory, ncoil, nseg, nwind, XYZ)
    character(len = *), intent(in) :: directory
    integer, intent(out) :: ncoil, nseg, nwind
    real(dp), intent(out), allocatable, dimension(:, :, :) :: XYZ
    integer :: fid, kc, ks, idum
    real(dp) :: ddum

    open(newunit = fid, file = directory // '/aug_bu.dat', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, '(3(1x, i4), 1x, f7.2)') ncoil, idum, nseg, ddum
    nseg = nseg - 1
    nwind = int(ddum)
    allocate(XYZ(3, nseg, 2 * ncoil))
    do kc = 1, ncoil
       do ks = 1, nseg
          read (fid, '(3(1x, es12.4e2))') XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc)
       end do
       read (fid, *)
    end do
    close(fid)
    open(newunit = fid, file = directory // '/aug_bl.dat', status = 'old', &
         action = 'read', form = 'formatted')
    read (fid, *)
    do kc = ncoil + 1, 2 * ncoil
       do ks = 1, nseg
          read (fid, '(3(1x, es12.4e2))') XYZ(1, ks, kc), XYZ(2, ks, kc), XYZ(3, ks, kc)
       end do
       read (fid, *)
    end do
    close(fid)
    XYZ(:, :, :) = 1d2 * XYZ
  end subroutine AUG_coils_read_GPEC

  subroutine AUG_coils_write_Fourier(directory, ncoil, nseg, nwind, XYZ, nmax, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use magdif_conf, only: log
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), dimension(:, :, :, :, :), allocatable :: Bn
    character(len = 1024) :: filename
    integer :: ntor, kc
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname

    if (3 /= size(XYZ, 1)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', '3', 'size(XYZ, 1)', 3, size(XYZ, 1))
       if (log%err) call log%write
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       call log%msg_arg_size('AUG_coils_write_GPEC', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3))
       if (log%err) call log%write
       error stop
    end if
    if (nmax > nphi / 4) then
       write (log%msg, '("Requested nmax = ", i0, ", but only ", i0, " modes available.")') &
            nmax, nphi / 4
       if (log%err) call log%write
       error stop
    end if
    call Biot_Savart_Fourier(ncoil, nseg, nwind, XYZ, nmax, &
         Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    filename = directory // '/AUG_B_coils.h5'
    call h5_open_rw(trim(filename), h5id_root)
    do ntor = 0, nmax
       write (modename, '("ntor_", i2.2)') ntor
       call h5_create_parent_groups(h5id_root, modename // '/')
       call h5_add(h5id_root, modename // '/R_min', Rmin, &
            unit = 'cm', comment = 'minimal R coordinate of computational grid')
       call h5_add(h5id_root, modename // '/R_max', Rmax, &
            unit = 'cm', comment = 'maximal R coordinate of computational grid')
       call h5_add(h5id_root, modename // '/Z_min', Zmin, &
            unit = 'cm', comment = 'minimal Z coordinate of computational grid')
       call h5_add(h5id_root, modename // '/Z_max', Zmax, &
            unit = 'cm', comment = 'maximal Z coordinate of computational grid')
       call h5_add(h5id_root, modename // '/nR', nR, 'number of grid points in R direction')
       call h5_add(h5id_root, modename // '/nZ', nZ, 'number of grid points in Z direction')
       call h5_add(h5id_root, modename // '/ncoil', ncoil, 'number of upper and lower B coils')
       do kc = 1, 2 * ncoil
          write (coilname, '("coil_", i2.2)') kc
          call h5_create_parent_groups(h5id_root, modename // '/' // coilname // '/')
          call h5_add(h5id_root, modename // '/' // coilname // '/Bn_R', &
               Bn(ntor, 1, :, :, kc), [1, 1], [nR, nZ], &
               unit = 'G', comment = 'R component of coil field for I_c = 0.1 A')
          call h5_add(h5id_root, modename // '/' // coilname // '/Bn_phi', &
               Bn(ntor, 2, :, :, kc), [1, 1], [nR, nZ], &
               unit = 'G', comment = 'physical phi component of coil field for I_c = 0.1 A')
          call h5_add(h5id_root, modename // '/' // coilname // '/Bn_Z', &
               Bn(ntor, 3, :, :, kc), [1, 1], [nR, nZ], &
               unit = 'G', comment = 'Z component of coil field for I_c = 0.1 A')
       end do
    end do
    call h5_close(h5id_root)
    if (allocated(Bn)) deallocate(Bn)
  end subroutine AUG_coils_write_Fourier

  subroutine read_currents_Nemov(directory, Ic)
    character(len = *), intent(in) :: directory
    real(dp), intent(out), allocatable :: Ic(:)
    integer, parameter :: nwind = 5, ncoil = 8
    integer :: fid

    allocate(Ic(2 * ncoil))
    open(newunit = fid, file = directory // '/cur_asd.dd', status = 'old', action = 'read', form = 'formatted')
    read (fid, *) Ic(:)
    close(fid)
    Ic(:) = Ic / real(nwind)
  end subroutine read_currents_Nemov

  subroutine Biot_Savart_sum_coils(ncoil, nseg, nwind, XYZ, Ic, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: log
    use magdif_util, only: linspace
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    real(dp), intent(in), dimension(:) :: Ic
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    real(dp), intent(out), dimension(:, :, :, :), allocatable :: Bvac
    integer :: kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), rXYZ(3), crXYZ(3), BXYZ(3), BRfZ(3)
    real(dp), dimension(size(XYZ, 1), size(XYZ, 2), size(XYZ, 3)) :: cXYZ, dXYZ

    if (3 /= size(XYZ, 1)) then
       call log%msg_arg_size('Biot_Savart_sum_coils', '3', 'size(XYZ, 1)', 3, size(XYZ, 1))
       if (log%err) call log%write
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       call log%msg_arg_size('Biot_Savart_sum_coils', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       call log%msg_arg_size('Biot_Savart_sum_coils', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(Ic)) then
       call log%msg_arg_size('Biot_Savart_sum_coils', '2 * ncoil', 'size(Ic)', &
            2 * ncoil, size(Ic))
       if (log%err) call log%write
       error stop
    end if
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    cXYZ(:, :nseg-1, :) = 0.5d0 * (XYZ(:, 2:, :) + XYZ(:, :nseg-1, :))
    cXYZ(:, nseg, :) = 0.5d0 * (XYZ(:, 1, :) + XYZ(:, nseg, :))
    dXYZ(:, :nseg-1, :) = XYZ(:, 2:, :) - XYZ(:, :nseg-1, :)
    dXYZ(:, nseg, :) = XYZ(:, 1, :) - XYZ(:, nseg, :)
    allocate(Bvac(3, nZ, nphi, nR))
    Bvac(:, :, :, :) = 0d0
    do kR = 1, nR
       do kphi = 1, nphi
          rXYZ(1:2) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi)]
          do kZ = 1, nZ
             rXYZ(3) = Z(kZ)
             ! Biot-Savart integral over coil segments
             do kc = 1, 2 * ncoil
                BRfZ(:) = 0d0
                do ks = 1, nseg
                   crXYZ(:) = rXYZ(:) - cXYZ(:, ks, kc)
                   BXYZ(:) = &
                        [dXYZ(2, ks, kc) * crXYZ(3) - dXYZ(3, ks, kc) * crXYZ(2), &
                        dXYZ(3, ks, kc) * crXYZ(1) - dXYZ(1, ks, kc) * crXYZ(3), &
                        dXYZ(1, ks, kc) * crXYZ(2) - dXYZ(2, ks, kc) * crXYZ(1)] &
                        / sqrt(sum(crXYZ ** 2)) ** 3
                   BRfZ(:) = BRfZ + &
                        [BXYZ(2) * sinphi(kphi) + BXYZ(1) * cosphi(kphi), &
                        BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi), &
                        BXYZ(3)]
                end do
                Bvac(:, kZ, kphi, kR) = Bvac(:, kZ, kphi, kR) + nwind * Ic(kc) * BRfZ
             end do
          end do
       end do
    end do
  end subroutine Biot_Savart_sum_coils

  subroutine Biot_Savart_Fourier(ncoil, nseg, nwind, XYZ, nmax, &
       Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)
    use iso_c_binding, only: c_ptr, c_double, c_double_complex, c_size_t, c_f_pointer
    use FFTW3, only: fftw_alloc_real, fftw_alloc_complex, fftw_plan_dft_r2c_1d, FFTW_PATIENT, &
         FFTW_DESTROY_INPUT, fftw_execute_dft_r2c, fftw_destroy_plan, fftw_free
    use constants, only: pi  ! orbit_mod.f90
    use magdif_conf, only: log
    use magdif_util, only: linspace
    integer, intent(in) :: ncoil, nseg, nwind
    real(dp), intent(in), dimension(:, :, :) :: XYZ
    integer, intent(in) :: nmax
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in) :: nR, nphi, nZ
    complex(dp), intent(out), dimension(:, :, :, :, :), allocatable :: Bn
    integer :: nfft, kc, ks, kR, kphi, kZ
    real(dp), dimension(nphi) :: phi, cosphi, sinphi
    real(dp) :: R(nR), Z(nZ), rXYZ(3), crXYZ(3), BXYZ(3)
    real(dp), dimension(size(XYZ, 1), size(XYZ, 2), size(XYZ, 3)) :: cXYZ, dXYZ
    type(c_ptr) :: plan_R, plan_phi, plan_Z, p_BR, p_Bphi, p_BZ, p_BnR, p_Bnphi, p_BnZ
    real(c_double), dimension(:), pointer :: BR, Bphi, BZ
    complex(c_double_complex), dimension(:), pointer :: BnR, Bnphi, BnZ

    if (3 /= size(XYZ, 1)) then
       call log%msg_arg_size('Biot_Savart_Fourier', '3', 'size(XYZ, 1)', 3, size(XYZ, 1))
       if (log%err) call log%write
       error stop
    end if
    if (nseg /= size(XYZ, 2)) then
       call log%msg_arg_size('Biot_Savart_Fourier', 'nseg', 'size(XYZ, 2)', nseg, size(XYZ, 2))
       if (log%err) call log%write
       error stop
    end if
    if (2 * ncoil /= size(XYZ, 3)) then
       call log%msg_arg_size('Biot_Savart_Fourier', '2 * ncoil', 'size(XYZ, 3)', &
            2 * ncoil, size(XYZ, 3))
       if (log%err) call log%write
       error stop
    end if
    if (nmax > nphi / 4) then
       write (log%msg, '("Requested nmax = ", i0, ", but only ", i0, " modes available.")') &
            nmax, nphi / 4
       if (log%err) call log%write
       error stop
    end if
    nfft = nphi / 2 + 1
    R(:) = linspace(Rmin, Rmax, nR, 0, 0)
    phi(:) = linspace(0d0, 2d0 * pi, nphi, 0, 1)
    cosphi(:) = cos(phi)
    sinphi(:) = sin(phi)
    Z(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    cXYZ(:, :nseg-1, :) = 0.5d0 * (XYZ(:, 2:, :) + XYZ(:, :nseg-1, :))
    cXYZ(:, nseg, :) = 0.5d0 * (XYZ(:, 1, :) + XYZ(:, nseg, :))
    dXYZ(:, :nseg-1, :) = XYZ(:, 2:, :) - XYZ(:, :nseg-1, :)
    dXYZ(:, nseg, :) = XYZ(:, 1, :) - XYZ(:, nseg, :)
    allocate(Bn(0:nmax, 3, nR, nZ, 2 * ncoil))
    ! prepare FFTW
    p_BR = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_BR, BR, [nphi])
    p_Bphi = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_Bphi, Bphi, [nphi])
    p_BZ = fftw_alloc_real(int(nphi, c_size_t))
    call c_f_pointer(p_BZ, BZ, [nphi])
    p_BnR = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_BnR, BnR, [nfft])
    p_Bnphi = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_Bnphi, Bnphi, [nfft])
    p_BnZ = fftw_alloc_complex(int(nfft, c_size_t))
    call c_f_pointer(p_BnZ, BnZ, [nfft])
    plan_R = fftw_plan_dft_r2c_1d(nphi, BR, BnR, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))
    plan_phi = fftw_plan_dft_r2c_1d(nphi, Bphi, Bnphi, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))
    plan_Z = fftw_plan_dft_r2c_1d(nphi, BZ, BnZ, ior(FFTW_PATIENT, FFTW_DESTROY_INPUT))
    do kc = 1, 2 * ncoil
       do kZ = 1, nZ
          rXYZ(3) = Z(kZ)
          do kR = 1, nR
             ! sample points toroidally
             do kphi = 1, nphi
                rXYZ(1:2) = [R(kR) * cosphi(kphi), R(kR) * sinphi(kphi)]
                BR(kphi) = 0d0
                Bphi(kphi) = 0d0
                BZ(kphi) = 0d0
                ! Biot-Savart integral over coil segments
                do ks = 1, nseg
                   crXYZ(:) = rXYZ(:) - cXYZ(:, ks, kc)
                   BXYZ(:) = &
                        [dXYZ(2, ks, kc) * crXYZ(3) - dXYZ(3, ks, kc) * crXYZ(2), &
                        dXYZ(3, ks, kc) * crXYZ(1) - dXYZ(1, ks, kc) * crXYZ(3), &
                        dXYZ(1, ks, kc) * crXYZ(2) - dXYZ(2, ks, kc) * crXYZ(1)] &
                        / sqrt(sum(crXYZ ** 2)) ** 3
                   BR(kphi) = BR(kphi) + BXYZ(2) * sinphi(kphi) + BXYZ(1) * cosphi(kphi)
                   Bphi(kphi) = Bphi(kphi) + BXYZ(2) * cosphi(kphi) - BXYZ(1) * sinphi(kphi)
                   BZ(kphi) = BZ(kphi) + BXYZ(3)
                end do
             end do
             call fftw_execute_dft_r2c(plan_R, BR, BnR)
             call fftw_execute_dft_r2c(plan_phi, Bphi, Bnphi)
             call fftw_execute_dft_r2c(plan_Z, BZ, BnZ)
             Bn(0:nmax, 1, kR, kZ, kc) = nwind * BnR(1:nmax+1) / dble(nphi)
             Bn(0:nmax, 2, kR, kZ, kc) = nwind * Bnphi(1:nmax+1) / dble(nphi)
             Bn(0:nmax, 3, kR, kZ, kc) = nwind * BnZ(1:nmax+1) / dble(nphi)
          end do
       end do
    end do
    call fftw_destroy_plan(plan_R)
    call fftw_destroy_plan(plan_phi)
    call fftw_destroy_plan(plan_Z)
    call fftw_free(p_BR)
    call fftw_free(p_Bphi)
    call fftw_free(p_BZ)
    call fftw_free(p_BnR)
    call fftw_free(p_Bnphi)
    call fftw_free(p_BnZ)
    ! nullify pointers past this point
  end subroutine Biot_Savart_Fourier

  subroutine vector_potential_single_mode(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    use bdivfree_mod, only: nR_mod => nr, nZ_mod => nz, Rmin_mod => rmin, Zmin_mod => zmin, &
         ntor, icp, ipoint, hr, hz, rpoi, zpoi, aznre, aznim, arnre, arnim
    use magdif_util, only: imun, linspace
    use magdif_conf, only: conf
    integer, intent(in) :: nR, nZ
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(in) :: Bn_R(:, :), Bn_Z(:, :)
    integer :: iz, imi(nz), ima(nz), jmi(nr), jma(nr)
    complex(dp) :: An_R(size(Bn_Z, 1), size(Bn_Z, 2)), An_Z(size(Bn_R, 1), size(Bn_R, 2))

    nR_mod = nR
    nZ_mod = nZ
    ntor = conf%n
    icp = nr * nz
    Rmin_mod = Rmin
    Zmin_mod = Zmin
    hr = (Rmax - Rmin) / dble(nR - 1)
    hz = (Zmax - Zmin) / dble(nZ - 1)
    allocate(rpoi(nR), zpoi(nZ), ipoint(nR, nZ))
    rpoi(:) = linspace(Rmin, Rmax, nR, 0, 0)
    zpoi(:) = linspace(Zmin, Zmax, nZ, 0, 0)
    do iZ = 1, nZ
       An_Z(:, iZ) = -imun * Bn_R(:, iZ) * rpoi(:) / dble(ntor)
       An_R(:, iZ) = imun * Bn_Z(:, iZ) * rpoi(:) / dble(ntor)
    end do
    allocate(aznre(6, 6, icp, ntor), aznim(6, 6, icp, ntor))
    allocate(arnre(6, 6, icp, ntor), arnim(6, 6, icp, ntor))
    aznre(:, :, :, :ntor-1) = (0d0, 0d0)
    aznim(:, :, :, :ntor-1) = (0d0, 0d0)
    arnre(:, :, :, :ntor-1) = (0d0, 0d0)
    arnim(:, :, :, :ntor-1) = (0d0, 0d0)
    imi(:) = 1
    ima(:) = nR
    jmi(:) = 1
    jma(:) = nZ
    call s2dcut(nr, nz, hr, hz, An_Z%re, imi, ima, jmi, jma, icp, aznre(:, :, :, ntor), ipoint)
    call s2dcut(nr, nz, hr, hz, An_Z%im, imi, ima, jmi, jma, icp, aznim(:, :, :, ntor), ipoint)
    call s2dcut(nr, nz, hr, hz, An_R%re, imi, ima, jmi, jma, icp, arnre(:, :, :, ntor), ipoint)
    call s2dcut(nr, nz, hr, hz, An_R%im, imi, ima, jmi, jma, icp, arnim(:, :, :, ntor), ipoint)
  end subroutine vector_potential_single_mode

  subroutine write_Bvac_Nemov(directory, Rmin, Rmax, Zmin, Zmax, Bvac)
    use magdif_conf, only: log
    use constants, only: pi  ! orbit_mod.f90
    character(len = *), intent(in) :: directory
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    real(dp), intent(in), dimension(:, :, :, :) :: Bvac
    integer :: fid, nR, nphi, nZ, kR, kphi, kZ

    if (3 /= size(Bvac, 1)) then
       call log%msg_arg_size('write_Bnvac_Nemov', '3', 'size(Bvac, 1)', 3, size(Bvac, 1))
       if (log%err) call log%write
       error stop
    end if
    nZ = size(Bvac, 2)
    nphi = size(Bvac, 3)
    nR = size(Bvac, 4)
    open(newunit = fid, file = directory // '/field.dat', status = 'replace', &
         action = 'write', form = 'formatted')
    ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
    write (fid, '(i0, 1x, i0, 1x, i0, 1x, i0)') nR, nphi + 1, nZ, 1
    write (fid, '(es24.16e3, 1x, es24.16e3)') Rmin, Rmax
    write (fid, '(es24.16e3, 1x, es24.16e3)') 0d0, 2d0 * pi
    write (fid, '(es24.16e3, 1x, es24.16e3)') Zmin, Zmax
    do kR = 1, nR
       do kphi = 1, nphi
          do kZ = 1, nZ
             write (fid, '(es24.16e3, 1x, es24.16e3, 1x, es24.16e3)') Bvac(:, kZ, kphi, kR)
          end do
       end do
       ! kisslinger_asdex.f90 writes points at phi = 0 and at phi = 2 pi
       kphi = 1
       do kZ = 1, nZ
          write (fid, '(es24.16e3, 1x, es24.16e3, 1x, es24.16e3)') Bvac(:, kZ, kphi, kR)
       end do
    end do
    close(fid)
  end subroutine write_Bvac_Nemov

  subroutine read_Bnvac_Nemov(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bnvac_R, Bnvac_Z)
    use input_files, only: pfile
    use magdif_conf, only: conf
    use magdif_util, only: imun, linspace, get_field_filenames
    use constants, only: pi  ! orbit_mod.f90
    character(len = 1024) :: gfile, convexfile
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    integer :: fid, nphi, idum, kR, kphi, kZ
    real(dp) :: B_R, B_phi, B_Z
    complex(dp), allocatable :: fourier_basis(:)

    call get_field_filenames(gfile, pfile, convexfile)
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
    use magdif_conf, only: conf, log
    use magdif_mesh, only: equil
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
       write (log%msg, '("File ", a, " not found, cannot read vacuum perturbation field.")') &
         trim(filename)
       if (log%err) call log%write
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
         write (log%msg, '(a, " returned error ", i0)') funcname, status
         if (log%err) call log%write
         error stop
      end if
    end subroutine check_error
  end subroutine read_Bnvac_GPEC

  subroutine read_Bnvac_Fourier(directory, ntor, Ic, nR, nZ, Rmin, Rmax, Zmin, Zmax, &
       Bnvac_R, Bnvac_Z)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    use magdif_conf, only: log
    character(len = *), intent(in) :: directory
    integer, intent(in) :: ntor
    real(dp), intent(in), dimension(:) :: Ic
    integer, intent(out) :: nR, nZ
    real(dp), intent(out) :: Rmin, Rmax, Zmin, Zmax
    complex(dp), intent(out), dimension(:, :), allocatable :: Bnvac_R, Bnvac_Z
    character(len = 1024) :: filename
    logical :: file_exists
    integer :: kc, ncoil
    integer(HID_T) :: h5id_root
    character(len = 7) :: modename, coilname
    complex(dp), dimension(:, :), allocatable :: Bn

    filename = trim(adjustl(directory)) // '/AUG_B_coils.h5'
    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) then
       write (log%msg, '("File ", a, " not found, cannot read vacuum perturbation field.")') &
         trim(filename)
       if (log%err) call log%write
       error stop
    end if
    call h5_open(filename, h5id_root)
    write (modename, '("ntor_", i2.2)') ntor
    call h5_get(h5id_root, modename // '/R_min', Rmin)
    call h5_get(h5id_root, modename // '/R_max', Rmax)
    call h5_get(h5id_root, modename // '/Z_min', Zmin)
    call h5_get(h5id_root, modename // '/Z_max', Zmax)
    call h5_get(h5id_root, modename // '/nR', nR)
    call h5_get(h5id_root, modename // '/nZ', nZ)
    call h5_get(h5id_root, modename // '/ncoil', ncoil)
    if (2 * ncoil /= size(Ic)) then
       call log%msg_arg_size('read_Bnvac_Fourier', '2 * ncoil', 'size(Ic)', &
            2 * ncoil, size(Ic))
       if (log%err) call log%write
       error stop
    end if
    allocate(Bnvac_R(nR, nZ), Bnvac_Z(nR, nZ), Bn(nR, nZ))
    Bnvac_R(:, :) = (0d0, 0d0)
    Bnvac_Z(:, :) = (0d0, 0d0)
    do kc = 1, 2 * ncoil
       write (coilname, '("coil_", i2.2)') kc
       call h5_get(h5id_root, modename // '/' // coilname // '/Bn_R', Bn)
       Bnvac_R(:, :) = Bnvac_R + Ic(kc) * Bn
       call h5_get(h5id_root, modename // '/' // coilname // '/Bn_Z', Bn)
       Bnvac_Z(:, :) = Bnvac_Z + Ic(kc) * Bn
    end do
    call h5_close(h5id_root)
    if (allocated(Bn)) deallocate(Bn)
  end subroutine read_Bnvac_Fourier

  subroutine compute_Bnvac(Bn)
    use magdif_conf, only: conf, log, vac_src_nemov, vac_src_gpec, vac_src_fourier
    use magdif_util, only: gauss_legendre_unit_interval, imun
    use magdif_mesh, only: mesh
    type(RT0_t), intent(inout) :: Bn
    integer, parameter :: order = 2
    integer :: nR, nZ, ktri, ke, k
    real(dp) :: Rmin, Rmax, Zmin, Zmax, R, Z, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_phi, B_Z
    complex(dp), dimension(:, :), allocatable :: Bn_R, Bn_Z

    ! initialize vacuum field
    select case (conf%vac_src)
    case (vac_src_nemov)
       call read_Bnvac_Nemov(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    case (vac_src_gpec)
       call read_Bnvac_GPEC(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    case (vac_src_fourier)
       call read_Bnvac_Fourier('.', conf%n, conf%Ic, &
            nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    case default
       write (log%msg, '("unknown vacuum field source selection", i0)') conf%vac_src
       if (log%err) call log%write
       error stop
    end select
    call vector_potential_single_mode(nR, nZ, Rmin, Rmax, Zmin, Zmax, Bn_R, Bn_Z)
    deallocate(Bn_R, Bn_Z)
    ! project to finite elements
    call gauss_legendre_unit_interval(order, points, weights)
    do ktri = 1, mesh%ntri
       node_R = mesh%node_R([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       node_Z = mesh%node_Z([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       do ke = 1, 3
          edge_R = node_R(ke + 1) - node_R(ke)
          edge_Z = node_Z(ke + 1) - node_Z(ke)
          do k = 1, order
             R = node_R(ke) * points(k) + node_R(ke + 1) * points(order - k + 1)
             Z = node_Z(ke) * points(k) + node_Z(ke + 1) * points(order - k + 1)
             call spline_bn(conf%n, R, Z, B_R, B_phi, B_Z)
             Bn%DOF(ke, ktri) = Bn%DOF(ke, ktri) + &
                  weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
          end do
       end do
       ! toroidal flux via zero divergence
       Bn%comp_phi(ktri) = imun / mesh%n * sum(Bn%DOF(:, ktri)) / mesh%area(ktri)
    end do
  end subroutine compute_Bnvac

  subroutine generate_vacfield
    use magdif_conf, only: conf, datafile
    use magdif_mesh, only: mesh
    type(RT0_t) :: Bn

    call RT0_init(Bn, mesh%ntri)
    if (conf%kilca_scale_factor /= 0) then
       call compute_kilca_vac_coeff
       call compute_kilca_vacuum(Bn)
       call check_kilca_vacuum
       call check_RT0(Bn)
    else
       if (conf%nonres) then
          call compute_Bn_nonres(Bn)
       else
          call compute_Bnvac(Bn)
          call debug_Bnvac_rectplot
          call debug_B0_rectplot
          call debug_Bmnvac
          call debug_fouriermodes
       end if
    end if
    call RT0_check_redundant_edges(Bn, 'Bnvac')
    call RT0_check_div_free(Bn, mesh%n, 1d-9, 'Bnvac')
    call RT0_write(Bn, datafile, 'Bnvac', 'magnetic field (vacuum)', 'G', 1)
    call RT0_deinit(Bn)
  end subroutine generate_vacfield

  subroutine debug_Bnvac_rectplot
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    use bdivfree_mod, only: nR, nZ, Rpoi, Zpoi
    use magdif_conf, only: conf, datafile
    character(len=*), parameter :: dataset = 'Bnvac'
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
    use magdif_conf, only: datafile
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
    use magdif_conf, only: conf, datafile
    use magdif_util, only: imun
    use magdif_mesh, only: equil, mesh, sample_polmodes
    character(len=*), parameter :: dataset = 'Bmnvac'
    integer, parameter :: m_max = 24
    integer :: kf, kt, ktri, m
    integer(HID_T) :: h5id_root
    complex(dp) :: Bn_R, Bn_Z, Bn_phi, Bn_contradenspsi, Bn_n
    complex(dp) :: fourier_basis(-m_max:m_max), Bmn_contradenspsi(-m_max:m_max, mesh%nflux), &
         Bmn_n(-m_max:m_max, mesh%nflux)

    Bmn_contradenspsi(:, :) = (0d0, 0d0)
    Bmn_n(:, :) = (0d0, 0d0)
    associate (s => sample_polmodes)
      do kf = 1, mesh%nflux
         do kt = 1, mesh%kt_max(kf)
            ktri = mesh%kt_low(kf) + kt
            fourier_basis = [(exp(-imun * m * s%theta(ktri)), m = -m_max, m_max)]
            call spline_bn(conf%n, s%R(ktri), s%Z(ktri), Bn_R, Bn_phi, Bn_Z)
            Bn_contradenspsi = s%sqrt_g(ktri) * (Bn_R * s%B0_Z(ktri) - &
                 Bn_Z * s%B0_R(ktri)) * s%R(ktri)
            Bn_n = (Bn_R * s%B0_Z(ktri) - Bn_Z * s%B0_R(ktri)) / hypot(s%B0_R(ktri), s%B0_Z(ktri))
            Bmn_contradenspsi(:, kf) = Bmn_contradenspsi(:, kf) + Bn_contradenspsi * fourier_basis
            Bmn_n(:, kf) = Bmn_n(:, kf) + Bn_n * fourier_basis
         end do
         Bmn_contradenspsi(:, kf) = Bmn_contradenspsi(:, kf) / mesh%kt_max(kf)
         Bmn_n(:, kf) = equil%cocos%sgn_dpsi * Bmn_n(:, kf) / mesh%kt_max(kf)
      end do
    end associate
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
    use magdif_conf, only: conf
    use magdif_util, only: imun
    use magdif_mesh, only: mesh, B0R, B0phi, B0Z
    type(RT0_t), intent(inout) :: Bn
    integer :: kf, kt, ktri, base, tip
    real(dp) :: r, lr, lz
    complex(dp) :: Bnpsi

    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          ! use midpoint of edge f
          base = mesh%lf(1, ktri)
          tip = mesh%lf(2, ktri)
          R = (mesh%node_R(base) + mesh%node_R(tip)) * 0.5d0
          lR = mesh%node_R(tip) - mesh%node_R(base)
          lZ = mesh%node_Z(tip) - mesh%node_Z(base)
          Bnpsi = -mesh%R_O * B0phi(mesh%ef(ktri), ktri) / r
          Bn%DOF(mesh%ef(ktri), ktri) = Bnpsi * (lR ** 2 + lZ ** 2) / &
               (B0R(mesh%ef(ktri), ktri) * lR + B0Z(mesh%ef(ktri), ktri) * lZ)
          Bn%comp_phi(ktri) = imun / mesh%n * Bn%DOF(mesh%ef(ktri), ktri) / mesh%area(ktri)
          Bn%DOF(mesh%ei(ktri), ktri) = (0d0, 0d0)
          Bn%DOF(mesh%eo(ktri), ktri) = (0d0, 0d0)
       end do
    end do
    if (conf%quad_avg) call avg_flux_on_quad(Bn)
  end subroutine compute_Bn_nonres

  subroutine avg_flux_on_quad(elem)
    use magdif_util, only: imun
    use magdif_mesh, only: mesh
    type(RT0_t), intent(inout) :: elem

    integer :: kf, kt, ktri1, ktri2
    complex(dp) :: tor_flux_avg, tor_flux_diff

    do kf = 2, mesh%nflux
       do kt = 1, mesh%kt_max(kf), 2
          ktri1 = mesh%kt_low(kf) + kt
          ktri2 = mesh%kt_low(kf) + kt + 1
          tor_flux_avg = 0.5d0 * (elem%comp_phi(ktri2) * mesh%area(ktri2) + &
               elem%comp_phi(ktri1) * mesh%area(ktri1))
          tor_flux_diff = 0.5d0 * (elem%comp_phi(ktri2) * mesh%area(ktri2) - &
               elem%comp_phi(ktri1) * mesh%area(ktri1))
          elem%comp_phi(ktri1) = tor_flux_avg / mesh%area(ktri1)
          elem%DOF(mesh%eo(ktri1), ktri1) = elem%DOF(mesh%eo(ktri1), ktri1) - &
               imun * mesh%n * tor_flux_diff
          elem%comp_phi(ktri2) = tor_flux_avg / mesh%area(ktri2)
          elem%DOF(mesh%ei(ktri2), ktri2) = elem%DOF(mesh%ei(ktri2), ktri2) + &
               imun * mesh%n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  ! calculate resonant vacuum perturbation
  subroutine compute_kilca_vacuum(Bn)
    use magdif_conf, only: conf
    use magdif_util, only: imun, gauss_legendre_unit_interval
    use magdif_mesh, only: mesh
    type(RT0_t), intent(inout) :: Bn
    integer, parameter :: order = 2
    integer :: ktri, k, ke, pol_modes(2)
    real(dp) :: R, Z, rho, theta, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_phi, B_Z

    pol_modes = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
    call gauss_legendre_unit_interval(order, points, weights)
    do ktri = 1, mesh%ntri
       node_R = mesh%node_R([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       node_Z = mesh%node_Z([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       do ke = 1, 3
          edge_R = node_R(ke + 1) - node_R(ke)
          edge_Z = node_Z(ke + 1) - node_Z(ke)
          do k = 1, order
             R = node_R(ke) * points(k) + node_R(ke + 1) * points(order - k + 1)
             Z = node_Z(ke) * points(k) + node_Z(ke + 1) * points(order - k + 1)
             rho = hypot(R - mesh%R_O, Z - mesh%Z_O)
             theta = atan2(Z - mesh%Z_O, R - mesh%R_O)
             call kilca_vacuum(mesh%n, pol_modes, mesh%R_O, rho, theta, B_R, B_phi, B_Z)
             Bn%DOF(ke, ktri) = Bn%DOF(ke, ktri) + &
                  weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
          end do
       end do
       ! toroidal flux via zero divergence
       Bn%comp_phi(ktri) = imun / mesh%n * sum(Bn%DOF(:, ktri)) / mesh%area(ktri)
    end do
  end subroutine compute_kilca_vacuum

  !> Calculate the vacuum perturbation field in cylindrical coordinates from the Fourier
  !> series of all given modes.
  !>
  !> @param tor_mode toroidal mode number, scaled by magdif_config::kilca_scale_factor
  !> (usually magdif_config::n)
  !> @param pol_modes array of poloidal mode numbers
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> magdif_util::g_eqdsk::rcentr)
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
    use magdif_conf, only: conf_arr
    use magdif_util, only: imun, straight_cyl2bent_cyl
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
            conf_arr%kilca_vac_coeff(abs(pol_modes(k))), temp_B_rad, temp_B_pol, temp_B_tor)
       B_rad = B_rad + temp_B_rad * fourier_basis(k)
       B_pol = B_pol + temp_B_pol * fourier_basis(k)
       B_tor = B_tor + temp_B_tor * fourier_basis(k)
    end do
    call straight_cyl2bent_cyl(B_rad, B_pol, B_tor, theta, B_R, B_phi, B_Z)
  end subroutine kilca_vacuum

  subroutine compute_kilca_vac_coeff
    use magdif_conf, only: conf_arr, log, cmplx_fmt
    use magdif_mesh, only: mesh
    integer :: m
    complex(dp) :: B_rad, B_pol, B_tor

    do m = mesh%m_res_min, mesh%m_res_max
       if (abs(conf_arr%kilca_vac_r(m)) <= 0d0) then
          write (log%msg, '("ignoring kilca_vac_r(", i0, "), ' // &
               'resorting to kilca_vac_coeff(", i0, ")")') m, m
          if (log%info) call log%write
          cycle
       end if
       call kilca_vacuum_fourier(mesh%n, m, mesh%R_O, conf_arr%kilca_vac_r(m), (1d0, 0d0), &
            B_rad, B_pol, B_tor)
       conf_arr%kilca_vac_coeff(m) = conf_arr%kilca_vac_Bz(m) / B_tor
    end do
    log%msg = 'effective vacuum perturbation field coefficients:'
    if (log%info) call log%write
    do m = mesh%m_res_min, mesh%m_res_max
       write (log%msg, '("kilca_vac_coeff(", i0, ") = ", ' // cmplx_fmt // ')') m, &
            conf_arr%kilca_vac_coeff(m)
       if (log%info) call log%write
    end do
  end subroutine compute_kilca_vac_coeff

  !> Calculate the Fourier coefficient of the vacuum perturbation field for a given
  !> toroidal-poloidal mode.
  !>
  !> @param tor_mode toroidal mode number, scaled by magdif_config::kilca_scale_factor
  !> (usually magdif_config::n)
  !> @param pol_mode poloidal mode number
  !> @param R_0 distance of straight cylinder axis to torus axis (usually
  !> magdif_util::g_eqdsk::rcentr)
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
    use magdif_conf, only: log
    use magdif_util, only: imun
    integer, intent(in) :: tor_mode, pol_mode
    real(dp), intent(in) :: R_0, r
    complex(dp), intent(in) :: vac_coeff
    complex(dp), intent(out) :: B_rad, B_pol, B_tor
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status

    k_z_r = tor_mode / R_0 * r
    status = fgsl_sf_bessel_icn_array(abs(pol_mode)-1, abs(pol_mode)+1, k_z_r, I_m)
    if (status /= fgsl_success .and. log%err) then
       write (log%msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
       call log%write
    end if
    B_rad = 0.5d0 * (I_m(-1) + I_m(1)) * vac_coeff
    B_pol = imun * pol_mode / k_z_r * I_m(0) * vac_coeff
    B_tor = imun * I_m(0) * vac_coeff
  end subroutine kilca_vacuum_fourier

  subroutine check_kilca_vacuum
    use magdif_conf, only: conf, conf_arr, longlines
    use magdif_mesh, only: fs_half, mesh
    complex(dp) :: B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    real(dp) :: rad
    integer :: kf, fid, abs_pol_mode

    abs_pol_mode = abs(conf%kilca_pol_mode)
    open(newunit = fid, file = 'cmp_vac.dat', recl = 3 * longlines)
    do kf = 1, mesh%nflux
       rad = fs_half%rad(kf)
       call kilca_vacuum_fourier(mesh%n, -abs_pol_mode, mesh%R_O, rad, &
            conf_arr%kilca_vac_coeff(abs_pol_mode), B_rad_neg, B_pol_neg, B_tor_neg)
       call kilca_vacuum_fourier(mesh%n, abs_pol_mode, mesh%R_O, rad, &
            conf_arr%kilca_vac_coeff(abs_pol_mode), B_rad_pos, B_pol_pos, B_tor_pos)
       write (fid, '(13(1x, es24.16e3))') rad, &
           B_rad_neg, B_pol_neg, B_tor_neg, B_rad_pos, B_pol_pos, B_tor_pos
    end do
    close(fid)
  end subroutine check_kilca_vacuum

  subroutine check_RT0(Bn)
    use magdif_conf, only: conf, longlines
    use magdif_mesh, only: fs_half, mesh, point_location
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    type(RT0_t), intent(in) :: Bn
    integer :: fid, kf, kpol, ktri
    real(dp) :: rad, theta, R, Z
    complex(dp) :: B_R, B_Z, B_phi, B_R_interp, B_Z_interp, B_phi_interp
    integer :: pol_modes(2)

    pol_modes = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
    open(newunit = fid, file = 'cmp_RT0.dat', recl = longlines)
    do kf = mesh%nflux / 3, mesh%nflux / 3 ! 1, nflux
       rad = fs_half%rad(kf)
       do kpol = 1, 2 * conf%nkpol
          theta = (kpol - 0.5d0) / dble(2 * conf%nkpol) * 2d0 * pi
          call kilca_vacuum(mesh%n, pol_modes, mesh%R_O, rad, theta, B_R, B_phi, B_Z)
          R = mesh%R_O + rad * cos(theta)
          Z = mesh%Z_O + rad * sin(theta)
          ktri = point_location(R, Z)
          call RT0_interp(ktri, Bn, R, Z, B_R_interp, B_Z_interp, B_phi_interp)
          write (fid, '(14(1x, es23.15e3))') rad, theta, B_R, B_phi, B_Z, &
               B_R_interp, B_phi_interp, B_Z_interp
       end do
    end do
    close(fid)
  end subroutine check_RT0

  subroutine debug_fouriermodes
    use magdif_conf, only: conf, datafile, log
    use magdif_util, only: flux_func, linspace, imun
    use magdif_mesh, only: equil
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
    type(flux_func) :: rq_interpolator

    inquire(file = 'amn.dat', exist = file_exists)
    if (.not. file_exists) return
    inquire(file = 'equil_r_q_psi.dat', exist = file_exists)
    if (.not. file_exists) return
    log%msg = 'Files amn.dat and equil_r_q_psi.dat found, performing Fourier mode comparison.'
    if (log%info) call log%write
    open(newunit = fid, form = 'unformatted', file = 'amn.dat')
    read (fid) ntor, mpol, nlabel, flabel_min, flabel_max
    allocate(z3dum(-mpol:mpol, ntor, nlabel))
    allocate(Amn_theta(-mpol:mpol, ntor, nlabel))
    read (fid) z3dum, Amn_theta
    close(fid)
    open(newunit = fid, form = 'formatted', file = 'equil_r_q_psi.dat')
    read (fid, '(a)') line
    call extract(line, ptrn_nsqpsi)
    read (line, *) nsqpsi
    read (fid, '(a)') line
    call extract(line, ptrn_psimax)
    read (line, *) psimax
    call extract(line, ptrn_phimax)
    read (line, *) phimax
    read (fid, *)
    allocate(qsaf(nsqpsi))
    allocate(psisurf(nsqpsi))
    allocate(rsmall(nsqpsi))
    do k = 1, nsqpsi
       read (fid, *) ddum, qsaf(k), psisurf(k), ddum, ddum, rsmall(k)
    end do
    close(fid)
    sgn_dpsi = sign(1d0, psimax)
    call rq_interpolator%init(4, rsmall * abs(qsaf))
    allocate(rq_eqd(nlabel))
    rq_eqd(:) = linspace(abs(flabel_min), abs(flabel_max), nlabel, 0, 0)
    allocate(psi_n(nlabel))
    psi_n(:) = [(rq_interpolator%interp(psisurf / psimax, rq_eqd(k)), k = 1, nlabel)]
    allocate(Bmn_contradenspsi(-mpol:mpol, nlabel))
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
    if (allocated(z3dum)) deallocate(z3dum)
    if (allocated(Amn_theta)) deallocate(Amn_theta)
    if (allocated(qsaf)) deallocate(qsaf)
    if (allocated(psisurf)) deallocate(psisurf)
    if (allocated(rsmall)) deallocate(rsmall)
    if (allocated(rq_eqd)) deallocate(rq_eqd)
    if (allocated(psi_n)) deallocate(psi_n)
    if (allocated(Bmn_contradenspsi)) deallocate(Bmn_contradenspsi)

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

end module magdif_pert
