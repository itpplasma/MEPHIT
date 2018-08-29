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
