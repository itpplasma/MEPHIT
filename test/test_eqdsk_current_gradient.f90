program test_eqdsk_current_gradient
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mephit_mesh, only: equil, curr0_geqdsk, equilibrium_field
    use mephit_util, only: clight, pi, geqdsk_import_hdf5, init_field, deinit_field
    use geqdsk_tools, only: geqdsk_deinit
    implicit none
    real(dp), parameter :: radii(3) = [590.0_dp, 620.0_dp, 657.0_dp]
    real(dp), parameter :: heights(3) = [-23.0_dp, 0.0_dp, 31.0_dp]
    real(dp) :: curvature, field_axis, slope, R, Z, h, scale, error, maximum
    real(dp) :: j0(3), dR(3), dZ(3), plus1(3), minus1(3), plus2(3), minus2(3)
    real(dp) :: finite_R(3), finite_Z(3)
    integer :: flux_sign, field_sign, profile, i, ir, iz, level, checks

    ! Analytic divergence-free fields and exactly interpolated flux functions.
    ! This tests the derivative identity; force balance is not assumed.
    allocate (equil%psi_eqd(9), equil%pprime(9), equil%fpol(9), &
        equil%ffprim(9), equil%fprime(9))
    do i = 1, 9
        equil%psi_eqd(i) = 25000.0_dp*real(i-5, dp)
    end do
    maximum = 0.0_dp
    checks = 0
    do flux_sign = -1, 1, 2
        curvature = 4.0_dp*real(flux_sign, dp)
        do field_sign = -1, 1, 2
            field_axis = 3.0e6_dp*real(field_sign, dp)
            slope = 50.0_dp*real(field_sign, dp)
            do profile = 1, 3
                equil%fpol = field_axis+slope*equil%psi_eqd
                equil%fprime = slope
                equil%pprime = 1.0e-3_dp+2.0e-9_dp*equil%psi_eqd
                equil%ffprim = slope*equil%fpol
                if (profile == 2) then
                    ! Independently tabulated FF': its own interpolant determines
                    ! the derivative of the actual native toroidal current.
                    equil%ffprim = equil%ffprim+ &
                        700.0_dp*equil%psi_eqd+0.005_dp*equil%psi_eqd**2
                else if (profile == 3) then
                    equil%fprime = 0.0_dp
                    equil%ffprim = 1.0e8_dp
                end if
                do iz = 1, size(heights)
                    Z = heights(iz)
                    do ir = 1, size(radii)
                        R = radii(ir)
                        call evaluate(R, Z, j0, dR, dZ)
                        scale = max(maxval(abs(dR)), maxval(abs(dZ)), &
                            clight*abs(equil%ffprim(5))/(4.0_dp*pi*R**2))
                        do level = 1, 2
                            h = 0.08_dp/real(level, dp)
                            call current(R+h, Z, plus1)
                            call current(R-h, Z, minus1)
                            call current(R+2.0_dp*h, Z, plus2)
                            call current(R-2.0_dp*h, Z, minus2)
                            finite_R = (8.0_dp*(plus1-minus1) &
                                -(plus2-minus2))/(12.0_dp*h)
                            call current(R, Z+h, plus1)
                            call current(R, Z-h, minus1)
                            call current(R, Z+2.0_dp*h, plus2)
                            call current(R, Z-2.0_dp*h, minus2)
                            finite_Z = (8.0_dp*(plus1-minus1) &
                                -(plus2-minus2))/(12.0_dp*h)
                            error = max(maxval(abs(dR-finite_R)), &
                                maxval(abs(dZ-finite_Z)))/scale
                            maximum = max(maximum, error)
                            if (error > 1.0e-8_dp) then
                                print *, 'Failed signs,profile,R,Z,h,error:', &
                                    flux_sign, field_sign, profile, R, Z, h, error
                                error stop 'Native current derivatives fail finite differences'
                            end if
                            checks = checks+1
                        end do
                        if (iz == 2) then
                            if (abs(dZ(2))/scale > 1.0e-12_dp) &
                                error stop 'Toroidal current breaks midplane symmetry'
                        end if
                    end do
                end do
            end do
        end do
    end do
    print *, 'Native current-gradient finite-difference checks:', checks
    print *, 'Maximum normalized derivative discrepancy:', maximum
    if (command_argument_count() == 1) call audit_native_equilibrium()
    if (command_argument_count() > 1) error stop 'Expected at most one native HDF5 path'
contains
    subroutine audit_native_equilibrium()
        character(len=1024) :: filename
        real(dp) :: angle, radius, axis_R, axis_Z, minor, RR, ZZ, step
        real(dp) :: value(3), grad_R(3), grad_Z(3), p1(3), m1(3), p2(3), m2(3)
        real(dp) :: fd_R(3), fd_Z(3), native_scale, discrepancy
        integer :: radial, poloidal, refinement

        call get_command_argument(1, filename)
        call geqdsk_deinit(equil)
        call geqdsk_import_hdf5(equil, trim(filename), 'equil')
        call init_field(equil)
        axis_R = equil%rmaxis
        axis_Z = equil%zmaxis
        minor = min(maxval(equil%rbbbs)-axis_R, axis_R-minval(equil%rbbbs))
        do refinement = 1, 2
            step = 0.04_dp/real(refinement, dp)
            discrepancy = 0.0_dp
            do radial = 1, 3
                radius = minor*real(radial, dp)/4.0_dp
                do poloidal = 0, 7
                    angle = pi*real(poloidal, dp)/4.0_dp
                    RR = axis_R+radius*cos(angle)
                    ZZ = axis_Z+radius*sin(angle)
                    call native_value(RR, ZZ, value, grad_R, grad_Z)
                    call native_value(RR+step, ZZ, p1)
                    call native_value(RR-step, ZZ, m1)
                    call native_value(RR+2.0_dp*step, ZZ, p2)
                    call native_value(RR-2.0_dp*step, ZZ, m2)
                    fd_R = (8.0_dp*(p1-m1)-(p2-m2))/(12.0_dp*step)
                    call native_value(RR, ZZ+step, p1)
                    call native_value(RR, ZZ-step, m1)
                    call native_value(RR, ZZ+2.0_dp*step, p2)
                    call native_value(RR, ZZ-2.0_dp*step, m2)
                    fd_Z = (8.0_dp*(p1-m1)-(p2-m2))/(12.0_dp*step)
                    native_scale = max(maxval(abs(fd_R)), maxval(abs(fd_Z)))
                    if (native_scale <= 0.0_dp) error stop 'Degenerate native audit'
                    discrepancy = max(discrepancy, max(maxval(abs(grad_R-fd_R)), &
                        maxval(abs(grad_Z-fd_Z)))/native_scale)
                end do
            end do
            ! A diagnostic, not a tolerance or equilibrium validity assertion.
            print *, 'Actual EQDSK step [cm], maximum relative gradient error:', &
                step, discrepancy
        end do
        call deinit_field()
    end subroutine audit_native_equilibrium

    subroutine native_value(RR, ZZ, value, grad_R, grad_Z)
        real(dp), intent(in) :: RR, ZZ
        real(dp), intent(out) :: value(3)
        real(dp), intent(out), optional :: grad_R(3), grad_Z(3)
        real(dp) :: B(3), BR(3), BZ(3), psi, Bmod, Bmod_R, Bmod_Z
        real(dp) :: local_R(3), local_Z(3)

        call equilibrium_field(RR, ZZ, B, BR, BZ, psi, Bmod, Bmod_R, Bmod_Z)
        call curr0_geqdsk(RR, psi, B, BR, BZ, value, local_R, local_Z)
        if (.not. all(ieee_is_finite(value))) error stop 'Nonfinite native current'
        if (.not. all(ieee_is_finite(local_R))) error stop 'Nonfinite radial gradient'
        if (.not. all(ieee_is_finite(local_Z))) error stop 'Nonfinite vertical gradient'
        if (present(grad_R)) grad_R = local_R
        if (present(grad_Z)) grad_Z = local_Z
    end subroutine native_value

    subroutine evaluate(R, Z, j0, dR, dZ)
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: j0(3), dR(3), dZ(3)
        real(dp) :: psi, psi_R, psi_Z, F, B(3), BR(3), BZ(3)

        psi = curvature*((R-620.0_dp)**2+Z**2)
        psi_R = 2.0_dp*curvature*(R-620.0_dp)
        psi_Z = 2.0_dp*curvature*Z
        F = field_axis+slope*psi
        B = [-psi_Z/R, F/R, psi_R/R]
        BR = [psi_Z/R**2, slope*psi_R/R-F/R**2, &
            2.0_dp*curvature/R-psi_R/R**2]
        BZ = [-2.0_dp*curvature/R, slope*psi_Z/R, 0.0_dp]
        call curr0_geqdsk(R, psi, B, BR, BZ, j0, dR, dZ)
        if (.not. all(ieee_is_finite(j0))) error stop 'Nonfinite manufactured current'
        if (.not. all(ieee_is_finite(dR))) error stop 'Nonfinite radial gradient'
        if (.not. all(ieee_is_finite(dZ))) error stop 'Nonfinite vertical gradient'
    end subroutine evaluate

    subroutine current(R, Z, j0)
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: j0(3)
        real(dp) :: unused_R(3), unused_Z(3)

        call evaluate(R, Z, j0, unused_R, unused_Z)
    end subroutine current
end program test_eqdsk_current_gradient
