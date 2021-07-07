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

  subroutine interp1d_export_hdf5(this, file, datatset)
    use hdf5_tools, only: HID_T, h5_open_rw, h5_create_parent_groups, h5_add, h5_close
    type(interp1d), intent(in) :: this
    character(len = *) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open_rw(file, h5id_root)
    call h5_create_parent_groups(h5id_root, trim(adjustl(dataset)) // '/')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n_lag', this%n_lag, comment = 'polynomial degree')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n_var', this%n_var, comment = 'number of sample points')
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/indep_var', this%indep_var, &
         lbound(this%indep_var), ubound(this%indep_var), comment = 'sample points')
    call h5_close(h5id_root)
  end subroutine interp1d_export_hdf5

  subroutine interp1d_import_hdf5(this, file, datatset)
    use hdf5_tools, only: HID_T, h5_open, h5_get, h5_close
    type(interp1d), intent(in) :: this
    character(len = *) :: file, dataset
    integer(HID_T) :: h5id_root

    call h5_open(file, h5id_root)
    call interp1d_deinit(this)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n_lag', this%n_lag)
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/n_var', this%n_var)
    allocate(this%indep_var(this%n_var))
    call h5_add(h5id_root, trim(adjustl(dataset)) // '/indep_var', this%indep_var)
    call h5_close(h5id_root)
  end subroutine interp1d_import_hdf5
