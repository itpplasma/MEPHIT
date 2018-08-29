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
