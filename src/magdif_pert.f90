module magdif_pert

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: write_scalar_dof, interp_RT0, check_div_free, check_redundant_edges, read_vector_dof, &
       write_vector_dof, write_vector_plot, write_vector_plot_rect, compute_Bn_nonres, &
       avg_flux_on_quad, compute_kilca_vacuum, kilca_vacuum, compute_kilca_vac_coeff, &
       kilca_vacuum_fourier, check_kilca_vacuum, check_RT0

contains

  !> Real and imaginary part of \p scalar_dof (e.g. #presn) are written, in that order,
  !> to \p outfile (e.g. #magdif_conf::presn_file), where line number corresponds to the
  !> knot index in #mesh_mod::mesh_point.
  subroutine write_scalar_dof(scalar_dof, outfile)
    use magdif_conf, only: longlines
    use magdif_mesh, only: mesh
    complex(dp), intent(in) :: scalar_dof(:)
    character(len = *), intent(in) :: outfile

    integer :: kpoint, fid

    open(newunit = fid, file = outfile, recl = longlines, status = 'replace')
    do kpoint = 1, mesh%kp_low(mesh%nflux + 1)
       write (fid, '(2(1x, es24.16e3))') scalar_dof(kpoint)
    end do
    close(fid)
  end subroutine write_scalar_dof

  subroutine interp_RT0(ktri, pol_flux, R, Z, comp_R, comp_Z, &
       comp_R_dR, comp_R_dZ, comp_Z_dR, comp_Z_dZ)
    use mesh_mod, only: knot, triangle, mesh_point, mesh_element
    integer, intent(in) :: ktri
    complex(dp), intent(in) :: pol_flux(:,:)
    real(dp), intent(in) :: r, z
    complex(dp), intent(out) :: comp_R, comp_Z
    complex(dp), intent(out), optional :: comp_R_dR, comp_R_dZ, comp_Z_dR, comp_Z_dZ
    type(triangle) :: elem
    type(knot) :: node(3)

    elem = mesh_element(ktri)
    node = mesh_point(elem%i_knot(:))
    ! edge 1 lies opposite to knot 3, etc.
    comp_R = 1d0 / elem%det_3 / R * ( &
         pol_flux(ktri, 1) * (R - node(3)%Rcoord) + &
         pol_flux(ktri, 2) * (R - node(1)%Rcoord) + &
         pol_flux(ktri, 3) * (R - node(2)%Rcoord))
    comp_Z = 1d0 / elem%det_3 / R * ( &
         pol_flux(ktri, 1) * (Z - node(3)%Zcoord) + &
         pol_flux(ktri, 2) * (Z - node(1)%Zcoord) + &
         pol_flux(ktri, 3) * (Z - node(2)%Zcoord))
    if (present(comp_R_dR)) then
       comp_R_dR = 1d0 / elem%det_3 / R ** 2 * ( &
            pol_flux(ktri, 1) * node(3)%Rcoord + &
            pol_flux(ktri, 2) * node(1)%Rcoord + &
            pol_flux(ktri, 3) * node(2)%Rcoord)
    end if
    if (present(comp_R_dZ)) then
       comp_R_dZ = (0d0, 0d0)
    end if
    if (present(comp_Z_dR)) then
       comp_Z_dR = -comp_Z / R
    end if
    if (present(comp_Z_dZ)) then
       comp_Z_dZ = sum(pol_flux(ktri, :)) / elem%det_3 / R
    end if
  end subroutine interp_RT0

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
    use magdif_mesh, only: mesh
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

  subroutine write_vector_dof(pol_flux, tor_comp, outfile)
    use iso_c_binding, only: c_long
    use magdif_mesh, only: mesh
    use mesh_mod, only: mesh_element_rmp
    complex(dp), intent(in) :: pol_flux(:,:)
    complex(dp), intent(in) :: tor_comp(:)
    character(len = *), intent(in) :: outfile
    integer(c_long) :: length
    integer :: ktri, fid

    length = 8 * mesh%kt_low(mesh%nflux + 1)  ! (Re, Im) of RT0 DoFs + toroidal component
    ! status = 'old' for writing to named pipe
    open(newunit = fid, file = outfile, access = 'stream', status = 'old', &
         action = 'write', form = 'unformatted')
    write (fid) length
    do ktri = 1, mesh%kt_low(mesh%nflux + 1)
       write (fid) pol_flux(ktri, :), tor_comp(ktri) * mesh_element_rmp(ktri)%area
    end do
    close(fid)
  end subroutine write_vector_dof

  subroutine read_vector_dof(pol_flux, tor_comp, infile)
    use iso_c_binding, only: c_long
    use magdif_conf, only: log
    use magdif_mesh, only: mesh
    use mesh_mod, only: mesh_element_rmp
    complex(dp), intent(out) :: pol_flux(:,:)
    complex(dp), intent(out) :: tor_comp(:)
    character(len = *), intent(in) :: infile
    integer(c_long) :: length
    integer :: ktri, fid

    open(newunit = fid, file = infile, access = 'stream', status = 'old', &
         action = 'read', form = 'unformatted')
    read (fid) length
    if (length < 8 * mesh%kt_low(mesh%nflux + 1)) then
       ! (Re, Im) of RT0 DoFs + toroidal component
       write (log%msg, '("File ", a, " only contains ", i0, " real values, ' // &
            'expected ", i0, ".")') infile, length, 8 * mesh%kt_low(mesh%nflux + 1)
       if (log%err) call log%write
       error stop
    end if
    do ktri = 1, mesh%kt_low(mesh%nflux + 1)
       read (fid) pol_flux(ktri, :), tor_comp(ktri)
    end do
    tor_comp = tor_comp / mesh_element_rmp%area
    close(fid)
  end subroutine read_vector_dof

  pure function jacobian(kf, kt, r) result(metric_det)
    use magdif_mesh, only: equil, fs_half, mesh, B0phi_Omega
    integer, intent(in) :: kf, kt
    real(dp), intent(in) :: r
    real(dp) :: metric_det
    metric_det = equil%cocos%sgn_dpsi * fs_half%q(kf) * r / B0phi_Omega(mesh%kt_low(kf) + kt)
  end function jacobian

  subroutine write_vector_plot(pol_flux, tor_comp, outfile)
    use magdif_conf, only: conf, longlines
    use magdif_mesh, only: equil, mesh, B0R_Omega, B0Z_Omega
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
       n_cutoff = mesh%nflux - 1
    else
       n_cutoff = mesh%nflux
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
    r = mesh%R_O
    z = mesh%Z_O
    do ktri = mesh%kt_low(n_cutoff+1) + 1, ntri
       write (fid, '(12(1x, es24.16e3))') r, z, (0d0, k = 1, 10)
    end do
    close(fid)
  end subroutine write_vector_plot

  subroutine write_vector_plot_rect(pol_flux, tor_comp, outfile)
    use magdif_conf, only: longlines
    use magdif_mesh, only: equil, mesh, point_location
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
          if (ktri > mesh%kt_low(1) .and. ktri <= mesh%kt_low(mesh%nflux + 1)) then
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

  subroutine compute_Bn_nonres(Bnflux, Bnphi)
    use mesh_mod, only: knot, triangle_rmp, mesh_point, mesh_element_rmp
    use magdif_conf, only: conf
    use magdif_util, only: imun
    use magdif_mesh, only: mesh, B0R, B0phi, B0Z
    complex(dp), intent(inout) :: Bnflux(:,:), Bnphi(:)
    integer :: kf, kt, ktri
    type(triangle_rmp) :: tri
    type(knot) :: base, tip
    real(dp) :: r, lr, lz
    complex(dp) :: Bnpsi

    do kf = 1, mesh%nflux
       do kt = 1, mesh%kt_max(kf)
          ktri = mesh%kt_low(kf) + kt
          tri = mesh_element_rmp(ktri)
          ! use midpoint of edge f
          base = mesh_point(tri%lf(1))
          tip = mesh_point(tri%lf(2))
          r = (base%rcoord + tip%rcoord) * 0.5d0
          lr = tip%rcoord - base%rcoord
          lz = tip%zcoord - base%zcoord
          Bnpsi = -mesh%R_O * B0phi(ktri, tri%ef) / r
          Bnflux(ktri, tri%ef) = Bnpsi * (lr ** 2 + lz ** 2) / &
               (B0r(ktri, tri%ef) * lr + B0z(ktri, tri%ef) * lz)
          Bnphi(ktri) = imun / mesh%n * Bnflux(ktri, tri%ef) / tri%area
          Bnflux(ktri, tri%ei) = (0d0, 0d0)
          Bnflux(ktri, tri%eo) = (0d0, 0d0)
       end do
    end do
    if (conf%quad_avg) call avg_flux_on_quad(Bnflux, Bnphi)
  end subroutine compute_Bn_nonres

  subroutine avg_flux_on_quad(pol_flux, tor_comp)
    use magdif_util, only: imun
    use magdif_mesh, only: mesh
    use mesh_mod, only: triangle_rmp, mesh_element_rmp
    complex(dp), intent(inout) :: pol_flux(:,:)
    complex(dp), intent(inout) :: tor_comp(:)

    integer :: kf, kt, ktri1, ktri2
    complex(dp) :: tor_flux_avg, tor_flux_diff
    type(triangle_rmp) :: tri1, tri2

    do kf = 2, mesh%nflux
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
          pol_flux(ktri1, tri1%eo) = pol_flux(ktri1, tri1%eo) - imun * mesh%n * tor_flux_diff
          tor_comp(ktri2) = tor_flux_avg / tri2%area
          pol_flux(ktri2, tri2%ei) = pol_flux(ktri2, tri2%ei) + imun * mesh%n * tor_flux_diff
       end do
    end do
  end subroutine avg_flux_on_quad

  ! calculate resonant vacuum perturbation
  subroutine compute_kilca_vacuum(Bnflux, Bnphi)
    use mesh_mod, only: ntri, knot, triangle, mesh_point, mesh_element, mesh_element_rmp
    use magdif_conf, only: conf
    use magdif_util, only: imun, gauss_legendre_unit_interval
    use magdif_mesh, only: mesh
    complex(dp), allocatable, intent(out) :: Bnflux(:,:), Bnphi(:)
    integer, parameter :: order = 2
    integer :: ktri, k, ke, pol_modes(2)
    real(dp) :: R, Z, rho, theta, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_phi, B_Z

    call gauss_legendre_unit_interval(order, points, weights)
    allocate(Bnflux(ntri, 3))
    allocate(Bnphi(ntri))
    Bnflux = (0d0, 0d0)
    Bnphi = (0d0, 0d0)
    pol_modes = [conf%kilca_pol_mode, -conf%kilca_pol_mode]
    do ktri = 1, ntri
       associate(tri => mesh_element_rmp(ktri), &
            knots => mesh_point(mesh_element(ktri)%i_knot))
         node_R = [knots(:)%rcoord, knots(1)%rcoord]
         node_Z = [knots(:)%zcoord, knots(1)%zcoord]
         do ke = 1, 3
            edge_R = node_R(ke + 1) - node_R(ke)
            edge_Z = node_Z(ke + 1) - node_Z(ke)
            do k = 1, order
               R = node_R(ke) * points(k) + node_R(ke + 1) * points(order - k + 1)
               Z = node_Z(ke) * points(k) + node_Z(ke + 1) * points(order - k + 1)
               rho = hypot(R - mesh%R_O, Z - mesh%Z_O)
               theta = atan2(Z - mesh%Z_O, R - mesh%R_O)
               call kilca_vacuum(mesh%n, pol_modes, mesh%R_O, rho, theta, B_R, B_phi, B_Z)
               Bnflux(ktri, ke) = Bnflux(ktri, ke) + &
                    weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
            end do
         end do
         ! toroidal flux via zero divergence
         Bnphi(ktri) = imun / mesh%n * sum(Bnflux(ktri, :)) / tri%area
       end associate
    end do
    call check_redundant_edges(Bnflux, .false., 'vacuum B_n')
    call check_div_free(Bnflux, Bnphi, mesh%n, 1d-9, 'vacuum B_n')
    call write_vector_dof(Bnflux, Bnphi, conf%Bn_vac_file)
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

  subroutine check_RT0(Bnflux, Bnphi)
    use magdif_conf, only: conf, longlines
    use magdif_mesh, only: fs_half, mesh, point_location
    use constants, only: pi  ! PRELOAD/SRC/orbit_mod.f90
    complex(dp), intent(in) :: Bnflux(:,:), Bnphi(:)
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
          call interp_RT0(ktri, Bnflux, R, Z, B_R_interp, B_Z_interp)
          B_phi_interp = Bnphi(ktri)
          write (fid, '(14(1x, es23.15e3))') rad, theta, B_R, B_phi, B_Z, &
               B_R_interp, B_phi_interp, B_Z_interp
       end do
    end do
    close(fid)
  end subroutine check_RT0

end module magdif_pert
