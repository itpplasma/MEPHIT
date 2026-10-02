program test_closed_surface_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use closed_surface_boundary_m, only: find_closed_boundary
    implicit none
    real(dp) :: offset, sign_flux, edge, box, target
    integer :: direction, gauge, padding, stat

    do direction = -1, 1, 2
        sign_flux = real(direction, dp)
        do gauge = 0, 1
            offset = real(gauge, dp)*1234.0_dp
            target = offset+sign_flux*0.37_dp**2
            do padding = 1, 3
                box = 2.0_dp+0.5_dp*real(padding, dp)
                call find_closed_boundary(nested_circle, 2.0_dp, box, target, edge, stat)
                if (stat /= 0) error stop 'Nested circle boundary not found'
                if (abs(edge-2.37_dp) > 2.0e-11_dp) &
                    error stop 'Circle radius depends on flux gauge or box padding'
            end do
            call find_closed_boundary(nested_circle, 2.0_dp, 2.1_dp, target, edge, stat)
            if (stat == 0) error stop 'Uncontained boundary was accepted'
            call find_closed_boundary(nested_circle, 2.0_dp, 3.0_dp, offset, edge, stat)
            if (stat == 0) error stop 'Magnetic axis was accepted as boundary'
        end do
    end do
    print *, 'Closed boundary passed exact circle, signed flux and gauge tests'
contains
    function nested_circle(r) result(psi)
        real(dp), intent(in) :: r
        real(dp) :: psi
        psi = offset+sign_flux*(r-2.0_dp)**2
    end function nested_circle
end program test_closed_surface_boundary
