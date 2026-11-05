program test_current_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mephit_iter, only: mde_current_source, perpendicular_current_divergence
    use mephit_util, only: clight, pi
    implicit none
    real(dp), parameter :: radii(3) = [590.0_dp, 620.0_dp, 657.0_dp]
    real(dp), parameter :: heights(3) = [-23.0_dp, 0.0_dp, 31.0_dp]
    real(dp), parameter :: radius_scale = 620.0_dp, height_scale = 31.0_dp
    complex(dp), parameter :: imaginary = (0.0_dp, 1.0_dp)
    real(dp) :: alpha, beta, curvature, axial_curve, R, Z, step, scale, error
    real(dp) :: B0(3), j0(3), half_grad_B2(3), grad_j0B0(3), Bmod
    real(dp) :: maximum_source, maximum_divergence
    complex(dp) :: Bn(3), grad_pn(3), grad_BnB0(3), source, divergence
    complex(dp) :: current(3), plus1(3), minus1(3), plus2(3), minus2(3)
    complex(dp) :: finite_R(3), finite_Z(3), finite_divergence
    integer :: field_sign, axis_sign, toroidal_mode, profile, ir, iz, level, checks

    ! Exact axisymmetric force-balanced pinches, in physical (R,phi,Z) components.
    ! Bphi = alpha*R + curvature*R**3; BZ = beta + axial_curve*R**2.
    ! The independent current below uses B0 x (c*grad(p1) - j0 x Bn)/B0**2.
    ! Its cylindrical divergence is differentiated, not the source expression.
    maximum_source = 0.0_dp
    maximum_divergence = 0.0_dp
    checks = 0
    do field_sign = -1, 1, 2
        alpha = 0.035_dp*real(field_sign, dp)
        do axis_sign = -1, 1, 2
            beta = 7.0_dp*real(axis_sign, dp)
            do toroidal_mode = -3, 3, 6
                do profile = 1, 2
                    curvature = 0.0_dp
                    axial_curve = 0.0_dp
                    if (profile == 2) then
                        curvature = 0.2_dp*alpha/radius_scale**2
                        axial_curve = 0.1_dp*beta/radius_scale**2
                    end if
                    do iz = 1, size(heights)
                        Z = heights(iz)
                        do ir = 1, size(radii)
                            R = radii(ir)
                            call coefficients(R, Z, B0, j0, Bn, grad_pn, &
                                half_grad_B2, grad_j0B0, grad_BnB0)
                            Bmod = sqrt(sum(B0**2))
                            source = mde_current_source(B0, Bmod, j0, &
                                half_grad_B2, Bn, grad_pn, grad_j0B0, grad_BnB0)
                            divergence = perpendicular_current_divergence(B0, &
                                Bmod, j0, half_grad_B2, Bn, grad_pn, &
                                grad_j0B0, grad_BnB0)
                            call physical_current(R, Z, current)
                            do level = 1, 2
                                step = 0.08_dp/real(level, dp)
                                call physical_current(R+step, Z, plus1)
                                call physical_current(R-step, Z, minus1)
                                call physical_current(R+2.0_dp*step, Z, plus2)
                                call physical_current(R-2.0_dp*step, Z, minus2)
                                finite_R = (8.0_dp*(plus1-minus1) &
                                    -(plus2-minus2))/(12.0_dp*step)
                                call physical_current(R, Z+step, plus1)
                                call physical_current(R, Z-step, minus1)
                                call physical_current(R, Z+2.0_dp*step, plus2)
                                call physical_current(R, Z-2.0_dp*step, minus2)
                                finite_Z = (8.0_dp*(plus1-minus1) &
                                    -(plus2-minus2))/(12.0_dp*step)
                                finite_divergence = finite_R(1)+current(1)/R &
                                    +finite_Z(3)+imaginary*real(toroidal_mode, dp) &
                                    *current(2)/R
                                scale = max(abs(finite_R(1)), abs(current(1)/R), &
                                    abs(finite_Z(3)), &
                                    abs(real(toroidal_mode, dp)*current(2)/R))
                                if (.not. finite_complex(source)) &
                                    error stop 'Nonfinite active current source'
                                if (.not. finite_complex(divergence)) &
                                    error stop 'Nonfinite diagnostic divergence'
                                if (.not. finite_complex(finite_divergence)) &
                                    error stop 'Nonfinite finite-difference oracle'
                                if (.not. ieee_is_finite(scale)) &
                                    error stop 'Nonfinite finite-difference scale'
                                if (scale <= tiny(1.0_dp)) &
                                    error stop 'Degenerate finite-difference oracle'
                                error = abs(source+finite_divergence)/scale
                                maximum_source = max(maximum_source, error)
                                if (error > 1.0e-8_dp) then
                                    print *, 'Source signs,profile,R,Z,h,error:', &
                                        field_sign, axis_sign, toroidal_mode, &
                                        profile, R, Z, step, error
                                    error stop 'Current MDE source is not -div(Jperp)'
                                end if
                                error = abs(divergence-finite_divergence)/scale
                                maximum_divergence = max(maximum_divergence, error)
                                if (error > 1.0e-8_dp) &
                                    error stop 'Diagnostic is not div(Jperp)'
                                checks = checks+1
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    print *, 'Native current-source finite-difference checks:', checks
    print *, 'Maximum normalized active-source discrepancy:', maximum_source
    print *, 'Maximum normalized diagnostic discrepancy:', maximum_divergence
contains
    pure logical function finite_complex(value)
        complex(dp), intent(in) :: value

        finite_complex = ieee_is_finite(real(value, dp))
        if (.not. finite_complex) return
        finite_complex = ieee_is_finite(aimag(value))
    end function finite_complex

    subroutine coefficients(R, Z, B0, j0, Bn, grad_pn, half_grad_B2, &
            grad_j0B0, grad_BnB0)
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: B0(3), j0(3), half_grad_B2(3), grad_j0B0(3)
        complex(dp), intent(out) :: Bn(3), grad_pn(3), grad_BnB0(3)
        real(dp) :: x, zeta, dB_dR(3), dj_dR(3)
        complex(dp) :: pressure, dpressure_R, dpressure_Z, f, df_R, ddf_RR, ddf_RZ

        x = R/radius_scale
        zeta = Z/height_scale
        B0 = [0.0_dp, alpha*R+curvature*R**3, beta+axial_curve*R**2]
        dB_dR = [0.0_dp, alpha+3.0_dp*curvature*R**2, 2.0_dp*axial_curve*R]
        j0 = clight/(4.0_dp*pi)*[0.0_dp, -2.0_dp*axial_curve*R, &
            2.0_dp*alpha+4.0_dp*curvature*R**2]
        dj_dR = clight/(4.0_dp*pi)*[0.0_dp, -2.0_dp*axial_curve, &
            8.0_dp*curvature*R]
        half_grad_B2 = [sum(B0*dB_dR), 0.0_dp, 0.0_dp]
        grad_j0B0 = [sum(dj_dR*B0+j0*dB_dR), 0.0_dp, 0.0_dp]
        pressure = 0.02_dp*((1.0_dp, 2.0_dp)*x**3 &
            +(-0.3_dp, 0.7_dp)*zeta+(0.4_dp, -0.6_dp)*x*zeta**2)
        dpressure_R = 0.02_dp/radius_scale*(3.0_dp*(1.0_dp, 2.0_dp)*x**2 &
            +(0.4_dp, -0.6_dp)*zeta**2)
        dpressure_Z = 0.02_dp/height_scale*((-0.3_dp, 0.7_dp) &
            +2.0_dp*(0.4_dp, -0.6_dp)*x*zeta)
        ! Bn = curl(f e_Z exp(in phi)), so div(Bn) vanishes exactly.
        f = 0.5_dp*((0.8_dp, -0.5_dp)*x**3+(0.3_dp, 0.4_dp)*x**2*zeta &
            +(0.2_dp, 0.6_dp)*x**4*zeta**2)
        df_R = 0.5_dp/radius_scale*(3.0_dp*(0.8_dp, -0.5_dp)*x**2 &
            +2.0_dp*(0.3_dp, 0.4_dp)*x*zeta &
            +4.0_dp*(0.2_dp, 0.6_dp)*x**3*zeta**2)
        ddf_RR = 0.5_dp/radius_scale**2*(6.0_dp*(0.8_dp, -0.5_dp)*x &
            +2.0_dp*(0.3_dp, 0.4_dp)*zeta &
            +12.0_dp*(0.2_dp, 0.6_dp)*x**2*zeta**2)
        ddf_RZ = 0.5_dp/(radius_scale*height_scale)* &
            (2.0_dp*(0.3_dp, 0.4_dp)*x &
            +8.0_dp*(0.2_dp, 0.6_dp)*x**3*zeta)
        Bn = [imaginary*real(toroidal_mode, dp)*f/R, -df_R, (0.0_dp, 0.0_dp)]
        grad_pn = [dpressure_R, imaginary*real(toroidal_mode, dp)*pressure/R, &
            dpressure_Z]
        grad_BnB0 = [-ddf_RR*B0(2)-df_R*dB_dR(2), &
            -imaginary*real(toroidal_mode, dp)*df_R*B0(2)/R, -ddf_RZ*B0(2)]
    end subroutine coefficients

    subroutine physical_current(R, Z, current)
        real(dp), intent(in) :: R, Z
        complex(dp), intent(out) :: current(3)
        real(dp) :: B0(3), j0(3), half_grad_B2(3), grad_j0B0(3)
        complex(dp) :: Bn(3), grad_pn(3), grad_BnB0(3), force(3)

        call coefficients(R, Z, B0, j0, Bn, grad_pn, half_grad_B2, &
            grad_j0B0, grad_BnB0)
        force = clight*grad_pn-cross_real_complex(j0, Bn)
        current = cross_real_complex(B0, force)/sum(B0**2)
    end subroutine physical_current

    pure function cross_real_complex(left, right) result(value)
        real(dp), intent(in) :: left(3)
        complex(dp), intent(in) :: right(3)
        complex(dp) :: value(3)

        value(1) = left(2)*right(3)-left(3)*right(2)
        value(2) = left(3)*right(1)-left(1)*right(3)
        value(3) = left(1)*right(2)-left(2)*right(1)
    end function cross_real_complex
end program test_current_source
