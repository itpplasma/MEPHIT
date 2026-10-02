program test_closed_boundary_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: rad, zet
    use field_sub, only: coefficient, psi_offset, axis_r, axis_z
    use field_line_integration_mod, only: use_eqdsk_boundary, target_boundary_flux, &
        theta0_at_xpoint, x_point
    implicit none
    integer, parameter :: nr = 12, nt = 128
    real(dp), parameter :: a = 37.0_dp, fpol = 620.0_dp*53000.0_dp
    real(dp) :: rmn, rmx, zmn, zmx, raxis, zaxis
    real(dp) :: rbeg(nr), rsmall(nr), q(nr), psi(nr), phitor(nr), perimeter(nr)
    real(dp) :: r(nr, nt), z(nr, nt), b(nr, nt), jac(nr, nt), expected(nr), error
    integer :: direction, gauge, padding, i
    character(len=32) :: probe
    external :: field_line_integration_for_SYNCH

    call get_command_argument(1, probe)
    use_eqdsk_boundary = .true.
    theta0_at_xpoint = .false.
    do direction = -1, 1, 2
        coefficient = real(direction, dp)*53000.0_dp/2.6_dp
        do gauge = 0, 1
            psi_offset = real(gauge, dp)*1.0e8_dp
            axis_z = real(gauge, dp)*17.0_dp
            target_boundary_flux = psi_offset+coefficient*a**2
            do padding = 1, 2
                rad = axis_r+[-a, a]*real(padding+1, dp)
                zet = axis_z+[-a, a]*real(padding+1, dp)
                if (trim(probe) == 'escaped') &
                    zet(2) = axis_z+a*(1.0_dp-1.0e-6_dp)
                if (trim(probe) == 'biased') zet(2) = axis_z+1.5_dp*a
                call field_line_integration_for_SYNCH(3600, 10000, nr, nt, &
                    rmn, rmx, zmn, zmx, raxis, zaxis, rbeg, rsmall, q, psi, phitor, &
                    perimeter, r, z, b, jac)
                expected = [(a*real(i, dp)/real(nr, dp), i = 1, nr)]
                if (maxval(abs(rsmall-expected))/a > 2.0e-7_dp) &
                    error stop 'Nested-circle area radius failed'
                if (maxval(abs(perimeter-2.0_dp*acos(-1.0_dp)*expected))/a > &
                    2.0e-6_dp) error stop 'Nested-circle perimeter failed'
                expected = fpol/(2.0_dp*coefficient*sqrt(axis_r**2-expected**2))
                error = maxval(abs(q-expected))/maxval(abs(expected))
                if (error > 2.0e-7_dp) error stop 'Signed analytic q failed'
                if (maxval(abs(hypot(r(nr, :)-axis_r, z(nr, :)-axis_z)-a))/a > &
                    2.0e-7_dp) error stop 'Edge contour depends on rectangle'
                if (norm2(x_point-[axis_r+a, axis_z])/a > 2.0e-7_dp) &
                    error stop 'Closed boundary anchor is off contour'
                if (abs(psi(nr)-coefficient*a**2)/abs(coefficient*a**2) > &
                    2.0e-7_dp) error stop 'Requested boundary flux failed'
            end do
        end do
    end do
    print *, 'Native closed-boundary geometry passed signed q, area and contour oracles'
end program test_closed_boundary_geometry
