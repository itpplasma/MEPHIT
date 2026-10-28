program test_resonant_current
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mephit_mesh, only: mesh, fs, equil, coord_cache_t, point_location, &
        point_location_current
    use mephit_pert, only: L1_t, L1_init, compute_Ires
    implicit none
    integer, parameter :: nedge = 16, npol = 256, nrad = 2, m = 11
    real(dp), parameter :: pi = acos(-1.0_dp), center_R = 10.0_dp
    complex(dp), parameter :: a = (0.7_dp, -0.2_dp), b = (-0.03_dp, 0.04_dp)
    complex(dp), parameter :: c = (0.12_dp, -0.07_dp)
    type(L1_t) :: current
    type(coord_cache_t) :: samples(npol, nrad)
    real(dp) :: bounds(2, 2, nedge), weights(nrad), phase, radius, fraction
    real(dp) :: psi_lo, width
    complex(dp) :: actual, expected
    integer :: k, next, nodes(3), krad, kpol, invalid, sign_psi, sign_q, profile
    interface
        subroutine Rtree_init(ntri, tri_bb) bind(C, name='Rtree_init')
            import :: c_int, c_double
            integer(c_int), value, intent(in) :: ntri
            real(c_double), intent(in) :: tri_bb(2, 2, ntri)
        end subroutine Rtree_init
    end interface

    ! A circular LCFS is represented by an inscribed polygon. Quadrature points
    ! inside the physical circle can be outside its straight boundary chords,
    ! as at the TC24 m=11 surface. The current is an exactly representable affine
    ! field; its Fourier integral below is independent of the triangle lookup.
    mesh%nflux = 1
    mesh%ntri = nedge
    allocate (mesh%node_R(nedge+1), mesh%node_Z(nedge+1), &
        mesh%tri_node(3, nedge), mesh%area(nedge), mesh%kp_low(1), &
        mesh%kp_max(1), mesh%edge_node(2, nedge), mesh%edge_tri(2, nedge))
    allocate (fs%rad(0:1), fs%psi(0:1))
    fs%rad = [0.0_dp, 1.0_dp]
    mesh%kp_low = 1
    mesh%kp_max = nedge
    mesh%node_R(1) = center_R
    mesh%node_Z(1) = 0.0_dp
    do k = 1, nedge
        phase = 2.0_dp*pi*real(k-1, dp)/real(nedge, dp)
        mesh%node_R(k+1) = center_R+cos(phase)
        mesh%node_Z(k+1) = sin(phase)
    end do
    do k = 1, nedge
        next = mod(k, nedge)+1
        nodes = [k+1, next+1, 1]
        mesh%tri_node(:, k) = nodes
        mesh%edge_node(:, k) = nodes(1:2)
        mesh%edge_tri(:, k) = [k, 0]
        mesh%area(k) = 0.5_dp*sin(2.0_dp*pi/real(nedge, dp))
        bounds(:, 1, k) = [minval(mesh%node_R(nodes)), &
            minval(mesh%node_Z(nodes))]
        bounds(:, 2, k) = [maxval(mesh%node_R(nodes)), &
            maxval(mesh%node_Z(nodes))]
    end do
    call Rtree_init(nedge, bounds)
    fs%psi = [0.0_dp, 1.0_dp]
    k = point_location(center_R+0.5_dp, 0.1_dp)
    if (point_location_current(center_R+0.5_dp, 0.1_dp, 0.26_dp) /= k) &
        error stop 'Interior element changed'
    if (point_location_current(center_R+1.01_dp, 0.0_dp, 1.02_dp) > 0) &
        error stop 'Point beyond the physical LCFS was accepted'
    if (point_location_current(center_R+3.0_dp, 0.0_dp, 0.99_dp) > 0) &
        error stop 'Unbounded extrapolation was accepted'
    call L1_init(current, nedge+1)
    psi_lo = 0.97_dp**2
    width = 1.0_dp-psi_lo
    weights = 0.5_dp*width
    do sign_psi = -1, 1, 2
        fs%psi = [0.0_dp, real(sign_psi, dp)]
        do sign_q = -1, 1, 2
            equil%cocos%sgn_q = sign_q
            invalid = 0
            do krad = 1, nrad
                fraction = 0.5_dp+(real(2*krad-3, dp)/sqrt(3.0_dp))*0.5_dp
                radius = sqrt(psi_lo+width*fraction)
                do kpol = 1, npol
                    associate (s => samples(kpol, krad))
                        s%psi = real(sign_psi, dp)*radius**2
                        s%theta = 2.0_dp*pi*(real(kpol-1, dp)+0.5_dp)/real(npol, dp)
                        s%R = center_R+radius*cos(s%theta)
                        s%Z = radius*sin(s%theta)
                        ! A positive shaped area measure with a known m=11 term.
                        s%sqrt_g = s%R*(1.0_dp+0.2_dp*cos(m*s%theta))
                        s%B0_R = 0.0_dp
                        s%B0_phi = 2.0_dp
                        s%B0_Z = 0.0_dp
                        s%ktri = point_location(s%R, s%Z)
                        if (s%ktri < 1) invalid = invalid+1
                    end associate
                end do
            end do
            if (invalid == 0) error stop 'Fixture misses the curved boundary gap'
            do profile = 1, 2
                current%DOF = a
                expected = 0.2_dp*width*a
                if (profile == 2) then
                    current%DOF = a+b*mesh%node_R+c*mesh%node_Z
                    expected = 0.2_dp*width*(a+b*center_R)
                end if
                ! Keep the invalid cached indices: old TC24 caches must also work.
                call compute_Ires(samples, weights, current, m, actual)
                if (.not. ieee_is_finite(abs(actual))) &
                    error stop 'NaN in LCFS-clipped resonant current'
                if (abs(actual-expected) > 2.0e-13_dp*abs(expected)) then
                    print *, 'Current integral:', actual, 'expected:', expected
                    error stop 'Analytic resonant-current integral failed'
                end if
            end do
        end do
    end do
    print *, 'Finite m=11 current agrees with constant and affine Fourier oracles'
end program test_resonant_current
