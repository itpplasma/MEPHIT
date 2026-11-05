program test_helical_current
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mephit_iter, only: helical_current_vector
    use mephit_pert, only: RT0_t, RT0_init, RT0_interp
    use mephit_mesh, only: mesh
    implicit none
    complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
    real(dp) :: B(3), R, F, BZ, wave, z, step, scale, maximum
    complex(dp) :: a, h, j(3), jp(3), jm(3), jp2(3), jm2(3), derivative
    complex(dp) :: parallel_expected(3)
    integer :: sf, sb, sk, sn, level, checks

    maximum = 0.0_dp
    checks = 0
    do sf = -1, 1, 2
        do sb = -1, 1, 2
            B(1) = 0.7_dp*sb
            B(2) = 3.0_dp*sf
            B(3) = -1.2_dp*sb
            R = 4.3_dp
            a = (0.23_dp, -0.41_dp)
            h = (-0.31_dp, 0.76_dp)
            j = helical_current_vector(R, B, a, h)
            call close_scalar(sum(j*B)/sum(B**2), h, 'parallel projection')
            j = helical_current_vector(R, B, a, (0.0_dp, 0.0_dp))
            call close_scalar(sum(j*B), (0.0_dp, 0.0_dp), &
                'perpendicular correction', maxval(abs(j))*sqrt(sum(B**2)))
            j = helical_current_vector(R, B, (0.0_dp, 0.0_dp), h)
            parallel_expected = h*B
            call close_vector(j, parallel_expected, 'parallel-only current')
            j = helical_current_vector(R, B, (0.0_dp, 0.0_dp), &
                (0.0_dp, 0.0_dp))
            call close_vector(j, [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
                (0.0_dp, 0.0_dp)], 'zero current')
        end do
    end do

    ! Independent cylindrical continuity oracle. BR=0, BZ constant, Bphi=F/R.
    ! Differentiate the native physical current in Z; the Fourier phi derivative
    ! is i*n/R. The coefficient is derived from div J=0, not from the helper.
    R = 3.2_dp
    z = 0.73_dp
    do sf = -1, 1, 2
        F = 12.0_dp*sf
        do sb = -1, 1, 2
            BZ = 2.3_dp*sb
            do sk = -1, 1, 2
                wave = 0.4_dp*sk
                do sn = -1, 1, 2
                    call cylinder(z, j)
                    do level = 1, 2
                        step = 0.005_dp/real(level, dp)
                        call cylinder(z+step, jp)
                        call cylinder(z-step, jm)
                        call cylinder(z+2.0_dp*step, jp2)
                        call cylinder(z-2.0_dp*step, jm2)
                        derivative = (8.0_dp*(jp(3)-jm(3)) &
                            -(jp2(3)-jm2(3)))/(12.0_dp*step)
                        scale = max(abs(derivative), abs(3.0_dp*j(2)/R))
                        maximum = max(maximum, &
                            abs(derivative+ii*3.0_dp*sn*j(2)/R)/scale)
                        if (maximum > 2.0e-10_dp) &
                            error stop 'Cylindrical physical current is not solenoidal'
                        checks = checks+1
                    end do
                end do
            end do
        end do
    end do
    call reconstructed_current()
    print *, 'Helical-current continuity checks:', checks
    print *, 'Maximum normalized divergence:', maximum
contains
    subroutine cylinder(zz, current)
        real(dp), intent(in) :: zz
        complex(dp), intent(out) :: current(3)
        complex(dp) :: hh, aa
        real(dp) :: native_field(3)

        hh = (0.6_dp, -0.7_dp)*exp(ii*wave*zz)
        aa = (wave*BZ+3.0_dp*sn*F/R**2) &
            /(wave*F*BZ-3.0_dp*sn*BZ**2)*hh
        native_field(1) = 0.0_dp
        native_field(2) = F/R
        native_field(3) = BZ
        current = helical_current_vector(R, native_field, aa, hh)
    end subroutine cylinder

    subroutine reconstructed_current()
        type(RT0_t) :: total
        real(dp) :: center_R, center_Z, field(3), normal(2)
        real(dp) :: toroidal_unit(3), inner_cross(3), outer_cross(3)
        complex(dp) :: bg_h, coeff, physical(3)
        complex(dp) :: expected(3), recovered(3), poloidal_flux(2)
        integer :: edge, endpoints(2, 3), first, last

        ! Native RT0 triangle with analytic fluxes of R*Jpol = constant.
        ! The centroid current is exactly the prescribed physical vector.
        mesh%ntri = 1
        allocate (mesh%node_R(3), mesh%node_Z(3), mesh%tri_node(3, 1), &
            mesh%tri_edge(3, 1), mesh%orient(1), mesh%area(1))
        mesh%node_R = [4.0_dp, 3.0_dp, 3.0_dp]
        mesh%node_Z = [0.0_dp, 2.0_dp, 0.0_dp]
        mesh%tri_node(:, 1) = [1, 2, 3]
        mesh%tri_edge(:, 1) = [1, 2, 3]
        mesh%orient = .true.
        mesh%area = 1.0_dp
        endpoints(:, 1) = [1, 2]
        endpoints(:, 2) = [3, 2]
        endpoints(:, 3) = [3, 1]
        center_R = sum(mesh%node_R)/3.0_dp
        center_Z = sum(mesh%node_Z)/3.0_dp
        field = [0.9_dp, -2.3_dp, 1.7_dp]
        bg_h = (0.4_dp, -0.2_dp)
        coeff = (0.17_dp, -0.29_dp)
        physical = helical_current_vector(center_R, field, coeff, bg_h)
        ! The independent vector identity B x (e_phi x B) is perpendicular.
        toroidal_unit = 0.0_dp
        toroidal_unit(2) = 1.0_dp
        inner_cross = cross(toroidal_unit, field)
        outer_cross = cross(field, inner_cross)
        expected = bg_h*field + coeff*center_R*outer_cross
        call RT0_init(total, 3, 1)
        total%comp_phi = physical(2)
        do edge = 1, 3
            first = endpoints(1, edge)
            last = endpoints(2, edge)
            normal(1) = mesh%node_Z(last)-mesh%node_Z(first)
            normal(2) = mesh%node_R(first)-mesh%node_R(last)
            poloidal_flux = center_R*physical([1, 3])
            total%DOF(edge) = sum(poloidal_flux*normal)
        end do
        call RT0_interp(total, 1, center_R, center_Z, recovered)
        call close_vector(recovered, expected, 'reconstructed physical helical current')
        call close_scalar(sum(recovered*field)/sum(field**2), &
            bg_h, 'reconstructed parallel coefficient')
    end subroutine reconstructed_current

    pure function cross(x, y) result(value)
        real(dp), intent(in) :: x(3), y(3)
        real(dp) :: value(3)

        value(1) = x(2)*y(3)-x(3)*y(2)
        value(2) = x(3)*y(1)-x(1)*y(3)
        value(3) = x(1)*y(2)-x(2)*y(1)
    end function cross

    subroutine close_scalar(actual, expected, label, normalization)
        complex(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        real(dp), intent(in), optional :: normalization
        real(dp) :: denominator

        denominator = max(1.0_dp, abs(expected))
        if (present(normalization)) denominator = max(denominator, normalization)
        if (.not. ieee_is_finite(abs(actual))) error stop 'Nonfinite current'
        if (abs(actual-expected) > 5.0e-13_dp*denominator) then
            print *, label, actual, expected
            error stop 'Physical-current invariant failed'
        end if
    end subroutine close_scalar

    subroutine close_vector(actual, expected, label)
        complex(dp), intent(in) :: actual(3), expected(3)
        character(len=*), intent(in) :: label
        integer :: component

        do component = 1, 3
            call close_scalar(actual(component), expected(component), label)
        end do
    end subroutine close_vector
end program test_helical_current
