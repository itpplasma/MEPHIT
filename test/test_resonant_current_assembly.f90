program test_resonant_current_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mephit_iter, only: add_resonant_current
    use mephit_pert, only: L1_t, RT0_t, L1_init, RT0_init, RT0_interp
    use mephit_mesh, only: mesh
    implicit none
    call assembled_current()
contains
    subroutine assembled_current()
        type(L1_t) :: parallel, resonant_parallel
        type(RT0_t) :: total, resonant
        real(dp) :: center_R, center_Z, field(3), normal(2)
        complex(dp) :: bg_h, resonant_h, physical(3), resonant_value(3)
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
        physical = cmplx([0.2_dp, -0.5_dp, 0.7_dp], [0.3_dp, 0.4_dp, -0.1_dp], dp)
        resonant_value = cmplx([-0.4_dp, 0.8_dp, 1.2_dp], [0.2_dp, -0.6_dp, 0.5_dp], dp)
        bg_h = sum(physical*field)/sum(field**2)
        resonant_h = sum(resonant_value*field)/sum(field**2)
        expected = physical + resonant_value
        call L1_init(parallel, 3)
        call L1_init(resonant_parallel, 3)
        call RT0_init(total, 3, 1)
        call RT0_init(resonant, 3, 1)
        parallel%DOF = bg_h
        resonant_parallel%DOF = resonant_h
        total%comp_phi = physical(2)
        resonant%comp_phi = resonant_value(2)
        do edge = 1, 3
            first = endpoints(1, edge)
            last = endpoints(2, edge)
            normal(1) = mesh%node_Z(last)-mesh%node_Z(first)
            normal(2) = mesh%node_R(first)-mesh%node_R(last)
            poloidal_flux = center_R*physical([1, 3])
            total%DOF(edge) = sum(poloidal_flux*normal)
            poloidal_flux = center_R*resonant_value([1, 3])
            resonant%DOF(edge) = sum(poloidal_flux*normal)
        end do
        call add_resonant_current(parallel, total, resonant_parallel, resonant)
        call RT0_interp(total, 1, center_R, center_Z, recovered)
        call close_vector(recovered, expected, 'assembled physical current')
        call close_scalar(sum(recovered*field)/sum(field**2), &
            sum(parallel%DOF)/3.0_dp, 'assembled parallel coefficient')
    end subroutine assembled_current
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
end program test_resonant_current_assembly
