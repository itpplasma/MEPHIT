  !> Returns the indices of the two triangles sharing an edge.
  !>
  !> @param knot1 first knot of the edge
  !> @param knot2 second knot of the edge
  !> @param common_tri indices of the triangles sharing the given edge
  !>
  !> The program is halted if the input data is invalid, i.e. if more than two triangles
  !> appear to share the edge.
  subroutine common_triangles(knot1, knot2, common_tri)
    type(knot), intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)

    integer :: k, l, kcom
    kcom = 0
    do k = 1, knot1%n_owners
       do l = 1, knot2%n_owners
          if (knot1%i_owner_tri(k) == knot2%i_owner_tri(l)) then
             kcom = kcom+1
             if (kcom > 2) stop "Error: more than two common triangles for knots"
             common_tri(kcom) = knot1%i_owner_tri(k)
          end if
       end do
    end do
  end subroutine common_triangles

complex(dp) function neumaier_sum(length, array)
  integer, intent(in) :: length
  complex(dp), intent(in), dimension(length) :: array
  real(dp) :: sr, si, cr, ci, tr, ti
  integer :: k

  sr = real(array(1))
  si = aimag(array(1))
  cr = 0.0d0
  ci = 0.0d0
  do k = 2, length
     tr = sr + real(array(k))
     if (abs(sr) >= abs(real(array(k)))) then
        cr = cr + (sr - tr) + real(array(k))
     else
        cr = cr + (real(array(k)) - tr) + sr
     end if
     sr = tr
     ti = si + aimag(array(k))
     if (abs(si) >= abs(aimag(array(k)))) then
        ci = ci + (si - ti) + aimag(array(k))
     else
        ci = ci + (aimag(array(k)) - ti) + si
     end if
     si = ti
  end do
  neumaier_sum = cmplx(sr + cr, si + ci)
end function neumaier_sum


subroutine unshared_knots(common_tri, outer_knot, inner_knot)
  integer, intent(in) :: common_tri(2)
  type(knot), intent(out) :: outer_knot, inner_knot
  integer, dimension(2, 3) :: set
  integer, dimension(2) :: i_unshared
  integer :: tri, k, n_unshared
  character(len = *), parameter :: errmsg = &
       'not exactly one unshared knot between triangles'

  set(1, :) = mesh_element(common_tri(1))%i_knot(:)
  set(2, :) = mesh_element(common_tri(2))%i_knot(:)
  do tri = 1, 2
     n_unshared = 0
     do k = 1, 3
        if (.not. any(set(tri, k) == set(3-tri, :))) then
           n_unshared = n_unshared + 1
           i_unshared(tri) = set(tri, k)
        end if
     end do
     if (n_unshared /= 1) then
        if (log_debug) write(logfile, *) errmsg
        stop errmsg
     end if
  end do
  if (i_unshared(1) > i_unshared(2)) then
     outer_knot = mesh_point(i_unshared(1))
     inner_knot = mesh_point(i_unshared(2))
  else
     outer_knot = mesh_point(i_unshared(2))
     inner_knot = mesh_point(i_unshared(1))
  end if
end subroutine unshared_knots


!> TODO: get rid of this due to redundancy with bnflux
subroutine read_hpsi
  integer :: k
  real(dp) :: dummy_re, dummy_im

  open(1, file = hpsi_file)
  do k = 1, ntri
     read (1, *) dummy_re, dummy_im
     if (.not. nonres) then
        mesh_element_rmp(k)%bnorm_vac = cmplx(dummy_re, dummy_im, dp)
     else
        !Test case: completely non-resonant perturbation
        mesh_element_rmp(k)%bnorm_vac = 3.d0 * R0 * abs(bphicovar) &
             / sum(mesh_point(mesh_element(k)%i_knot(:))%rcoord ** 2 &
             * mesh_point(mesh_element(k)%i_knot(:))%b_mod)
     endif
  end do
  close(1)
end subroutine read_hpsi


    integer, dimension(:,:), allocatable :: connections
    integer :: low, kt, ke, kt_adj, ke_adj
    allocate(connections(ntri, 3))
    connections = 0
    low = 1
    do kt = 1, kt_low(nflux+1)
       do ke = 1, 3
          if (connections(kt, ke) == 0) then
             kt_adj = mesh_element(kt)%neighbour(ke)
             ke_adj = mesh_element(kt)%neighbour_edge(ke)
             if (connections(kt_adj, ke_adj) == 0) then
                connections(kt, ke) = low
                connections(kt_adj, ke_adj) = low
                low = low + 1
             else
                connections(kt, ke) = connections(kt_adj, ke_adj)
             end if
          end if
       end do
    end do
    open(1, file = 'connectivity.dat')
    do kt = 1, ntri
       write (1, *) connections(kt, :)
    end do
    close(1)
    deallocate(connections)

program connectivity
  implicit none
  integer :: fid, kt, ke, low
  integer, parameter :: ntri = 35680
  integer, dimension(ntri, 3) :: connections
  integer, dimension(:), allocatable :: mapping

  open(newunit = fid, file = '../FEM/connectivity.dat')
  do kt = 1, ntri
     read (fid, *) connections(kt, :)
  end do
  close(fid)
  allocate(mapping(minval(connections):maxval(connections)))
  mapping = 0
  low = 1
  do kt = 1, ntri
     do ke = 1, 3
        if (mapping(connections(kt, ke)) == 0) then
           mapping(connections(kt, ke)) = low
           low = low + 1
        end if
        connections(kt, ke) = mapping(connections(kt, ke))
     end do
  end do
  open(newunit = fid, file = '../FEM/connections.dat')
  do kt = 1, ntri
     write (fid, *) connections(kt, :)
  end do
  close(fid)
  deallocate(mapping)
end program connectivity

    do kf = 1, nflux
       if (sheet_current_factor /= 0d0 .and. m_res(kf) > 0) then
          presn(kp_low(kf-1)+1:kp_low(kf-1)+kp_max(kf-1)) = &
               presn(kp_low(kf)+1:kp_low(kf)+kp_max(kf))
       end if
    end do

! get rid of high values at center of non-resonant test field
    C_psi = [(tanh(1d-2 * dble(kf)), kf = 2, nflux+1)]
    ! ...
          if (orient) then
             Bnpsi = -C_psi(kf) * R0 * B0phi(kt_low(kf) + kt, ef) / r
          else
             Bnpsi = -C_psi(kf-1) * R0 * B0phi(kt_low(kf) + kt, ef) / r
          end if
          Bnflux(kt_low(kf) + kt, ef) = Bnpsi * r / Deltapsi * &
               sum(mesh_element(common_tri(:))%det_3)
          Bnphi(kt_low(kf) + kt) = imun / n * Bnflux(kt_low(kf) + kt, ef) &
               / elem%det_3 * 2d0
          Bnflux(kt_low(kf) + kt, ei) = -1d-2 * (psi(kf) - psi(kf-1))
          Bnflux(kt_low(kf) + kt, eo) = 1d-2 * (psi(kf) - psi(kf-1))

  subroutine check_kilca_vacuum
    use fgsl, only: fgsl_double, fgsl_int, fgsl_success, fgsl_sf_bessel_icn_array
    use magdif_config, only: n, nkpol, R0, kilca_pol_mode, log_msg, log_err, log_write, kilca_vac_coeff, longlines
    use magdif_util, only: imun
    use magdif, only: fs_half
    complex(dp) :: B_r, B_theta, B_z
    real(fgsl_double) :: I_m(-1:1), k_z_r
    integer(fgsl_int) :: status
    integer :: kf, fid
    real(dp) :: r

    open(newunit = fid, file = 'cmp_vac.dat', recl = longlines)
    do kf = 1, ubound(fs_half%rad, 1)
       r = fs_half%rad(kf)
       k_z_r = n / R0 * r
       status = fgsl_sf_bessel_icn_array(abs(kilca_pol_mode)-1, abs(kilca_pol_mode)+1, k_z_r, I_m)
       if (status /= fgsl_success .and. log_err) then
          write (log_msg, '("fgsl_sf_bessel_icn_array returned error ", i0)') status
          call log_write
       end if
       B_r = 0.5d0 * (I_m(-1) + I_m(1)) * kilca_vac_coeff
       B_theta = imun * kilca_pol_mode / k_z_r * I_m(0) * kilca_vac_coeff
       B_z = imun * I_m(0) * kilca_vac_coeff
       write (fid, '(7(1x, es23.16))') r, B_r, B_theta, B_z
    end do
    close(fid)
  end subroutine check_kilca_vacuum
