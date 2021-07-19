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
