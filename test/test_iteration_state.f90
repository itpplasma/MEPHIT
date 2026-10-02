program test_iteration_state
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use mephit_iteration_state, only: select_iteration_state, relative_iteration_norm, &
        iteration_state_flags, iteration_continue, iteration_converged, &
        iteration_limit, iteration_single, iteration_invalid
    implicit none
    real(dp), parameter :: tolerance = 1.0e-7_dp
    complex(dp) :: K(3, 3), P(3, 3), input(3), field(3), drive(3), residual(3)
    complex(dp) :: exact(3), pressure(3), current(3), step(3)
    real(dp) :: relative, forcing, lambda(2), gram_inverse(2, 2), U(3, 2), Lr(2, 2)
    real(dp) :: nan, infinity
    integer :: action, evaluation, updates, checks, i
    logical :: converged, consistent
    character(len=32) :: reason

    checks = 0
    ! An explicitly triangular, nonnormal operator has an independent exact
    ! solution by back substitution. Its preconditioner includes a near-unit
    ! mode as well as a large mode; selected eigenvectors are not orthogonal.
    lambda = [1.0e6_dp, 0.9999_dp]
    K = cmplx(0.0_dp, 0.0_dp, dp)
    K(1, 1) = lambda(1)
    K(2, 2) = lambda(2)
    K(3, 3) = 0.6_dp
    K(1, 2) = 2.0_dp*(lambda(2)-lambda(1))
    drive = [cmplx(1.0_dp, 0.2_dp, dp), cmplx(1.0e-8_dp, -2.0e-8_dp, dp), &
        cmplx(0.1_dp, 0.2_dp, dp)]
    exact(2) = drive(2)/(1.0_dp-lambda(2))
    exact(1) = (drive(1)+K(1, 2)*exact(2))/(1.0_dp-lambda(1))
    exact(3) = drive(3)/(1.0_dp-real(K(3, 3), dp))
    U(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    U(:, 2) = [2.0_dp, 1.0_dp, 0.0_dp]
    gram_inverse(1, :) = [5.0_dp, -2.0_dp]
    gram_inverse(2, :) = [-2.0_dp, 1.0_dp]
    do i = 1, 2
        Lr(i, :) = lambda(i)/(lambda(i)-1.0_dp)*gram_inverse(i, :)
    end do
    P = -matmul(U, matmul(Lr, transpose(U)))
    do i = 1, 3
        P(i, i) = P(i, i)+1.0_dp
    end do
    input = exact+cmplx(1.0e-8_dp, 0.0_dp, dp)*U(:, 1)
    forcing = norm(drive)
    call physical_map(input, field, pressure, current)
    residual = field-input
    step = matmul(P, residual)
    call require(norm(step)/forcing < tolerance, 'Fixture must expose old false pass')
    call require(norm(residual)/forcing > 1.0e-3_dp, &
        'Original equation rejects fixture')
    call require(norm(pressure-1.0e8_dp*(input+step)) > 0.9_dp, &
        'Fixture exposes old saved pressure lag')
    call require(norm(current-(3.0e8_dp+2.0_dp)*(input+step)) > 2.9_dp, &
        'Fixture exposes old saved current lag')
    call select_iteration_state(input, field, norm(residual), forcing, tolerance, &
        .false., .false., action, relative)
    call require(action == iteration_continue, &
        'Original map must reject old false pass')
    updates = 0
    do evaluation = 0, 20
        call physical_map(input, field, pressure, current)
        residual = field-input
        call select_iteration_state(input, field, norm(residual), forcing, tolerance, &
            evaluation == 20, .false., action, relative)
        if (action /= iteration_continue) exit
        input = input+matmul(P, residual)
        updates = updates+1
    end do
    call iteration_state_flags(action, converged, consistent, reason)
    call require(converged, 'Independent linear-map oracle must converge')
    call require(consistent, 'Accepted output must have coherent constitutive state')
    call require(trim(reason) == 'original-map-converged', 'Convergence reason')
    call require(updates > 0, 'False-pass input cannot be accepted without an update')
    call require(norm(field-exact) < 1.0e-12_dp, 'Back-substitution field oracle')
    call require(norm(matmul(K, field)+drive-field)/forcing < tolerance, &
        'Accepted field satisfies independently evaluated original equation')
    call require(norm(pressure-1.0e8_dp*field) < 1.0e-10_dp, &
        'Saved pressure refers to accepted field')
    call require(norm(current-(3.0e8_dp+2.0_dp)*field) < 1.0e-9_dp, &
        'Saved current refers to accepted field')

    ! A fixed iteration budget is a separate oracle: one update of .75*b+v
    ! from zero returns v, with original residual .75*v; no second update.
    K = cmplx(0.0_dp, 0.0_dp, dp)
    do i = 1, 3
        K(i, i) = 0.75_dp
    end do
    drive = cmplx(1.0_dp, 0.0_dp, dp)
    input = cmplx(0.0_dp, 0.0_dp, dp)
    updates = 0
    do evaluation = 0, 1
        call physical_map(input, field, pressure, current)
        residual = field-input
        call select_iteration_state(input, field, norm(residual), norm(drive), &
            tolerance, evaluation == 1, .false., action, relative)
        if (action /= iteration_continue) exit
        input = input+residual
        updates = updates+1
    end do
    call iteration_state_flags(action, converged, consistent, reason)
    call require(action == iteration_limit, 'Exhaustion is unsuccessful')
    call require(.not. converged, 'No convergence flag on exhaustion')
    call require(consistent, 'Failed diagnostic retains coherent input state')
    call require(trim(reason) == 'iteration-limit', 'Exhaustion reason')
    call require(updates == 1, 'Exactly one requested update')
    call require(norm(field-drive) < 1.0e-14_dp, 'No hidden update beyond limit')
    call require(abs(relative-0.75_dp) < 1.0e-14_dp, 'Limit original residual')
    call require(norm(pressure-1.0e8_dp*field) < 1.0e-10_dp, 'Limit pressure coherence')

    ! Convergence on the last allowed input wins over exhaustion.
    input = 4.0_dp*drive
    call physical_map(input, field, pressure, current)
    call select_iteration_state(input, field, norm(field-input), norm(drive), &
        tolerance, .true., .false., action, relative)
    call require(action == iteration_converged, 'Exact solution at update limit')

    ! An accepted nonzero residual still requires pressure at the accepted
    ! input, not the nearby mapped field. Its known constitutive gradient
    ! makes a small field discrepancy observable without another map call.
    input(1) = input(1)+1.0e-8_dp
    call physical_map(input, field, pressure, current)
    call select_iteration_state(input, field, norm(field-input), norm(drive), &
        tolerance, .true., .false., action, relative)
    call require(action == iteration_converged, 'Small original residual accepts input')
    call require(norm(pressure-1.0e8_dp*field) < 1.0e-10_dp, &
        'Accepted nonzero residual preserves constitutive pressure')

    ! Single-update mode returns K*input+v, retains its constitutive input,
    ! and explicitly declares that it is neither accepted nor coherent.
    input = cmplx(2.0_dp, 0.0_dp, dp)
    call physical_map(input, field, pressure, current)
    call select_iteration_state(input, field, norm(field-input), norm(drive), &
        tolerance, .true., .true., action, relative)
    call iteration_state_flags(action, converged, consistent, reason)
    call require(action == iteration_single, 'Explicit single-update action')
    call require(norm(field-cmplx(2.5_dp, 0.0_dp, dp)) < 1.0e-14_dp, &
        'Single-update mapped field oracle')
    call require(.not. converged, 'Single-update is not converged')
    call require(.not. consistent, 'Single-update declares distinct constitutive input')
    call require(trim(reason) == 'single-update', 'Single-update reason')
    call require(norm(pressure-1.0e8_dp*input) < 1.0e-10_dp, &
        'Retained one-shot constitutive input')

    ! No new dimensional tolerance is invented for zero vacuum forcing.
    input = cmplx(0.0_dp, 0.0_dp, dp)
    field = input
    call select_iteration_state(input, field, 0.0_dp, 0.0_dp, tolerance, &
        .true., .false., action, relative)
    call require(action == iteration_converged, 'Exact zero-forcing fixed point')
    call require(relative <= 0.0_dp, 'Defined zero/zero normalized residual')
    field = cmplx(1.0_dp, 0.0_dp, dp)
    call select_iteration_state(input, field, 1.0_dp, 0.0_dp, tolerance, &
        .true., .false., action, relative)
    call require(action == iteration_limit, 'Nonzero residual with zero forcing fails')
    call require(relative > 1.0e100_dp, 'Defined nonzero/zero residual sentinel')
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    call invalid_control(nan, 1.0_dp, tolerance)
    call invalid_control(1.0_dp, infinity, tolerance)
    call invalid_control(1.0_dp, 1.0_dp, nan)
    call invalid_control(1.0_dp, 1.0_dp, 0.0_dp)
    call invalid_control(-1.0_dp, 1.0_dp, tolerance)
    call require(relative_iteration_norm(infinity, 1.0_dp) > 1.0e100_dp, &
        'Nonfinite step has a defined unsuccessful sentinel')
    field = cmplx(nan, 0.0_dp, dp)
    call select_iteration_state(input, field, 0.0_dp, 1.0_dp, tolerance, &
        .true., .false., action, relative)
    call require(action == iteration_invalid, 'Nonfinite mapped field cannot pass')
    print *, 'Independent original-map, state, mode and limit checks:', checks

contains

    subroutine physical_map(value, mapped, p, j)
        complex(dp), intent(in) :: value(:)
        complex(dp), intent(out) :: mapped(:), p(:), j(:)

        mapped = matmul(K, value)+drive
        p = 1.0e8_dp*value
        j = (3.0e8_dp+2.0_dp)*value
    end subroutine physical_map

    function norm(value) result(magnitude)
        complex(dp), intent(in) :: value(:)
        real(dp) :: magnitude

        magnitude = sqrt(sum(abs(value)**2))
    end function norm

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print *, 'FAIL: ', label
            error stop 'Independent iteration oracle failed'
        end if
        checks = checks+1
    end subroutine require

    subroutine invalid_control(residual_norm, forcing_norm, threshold)
        real(dp), intent(in) :: residual_norm, forcing_norm, threshold

        field = cmplx(1.0_dp, 0.0_dp, dp)
        call select_iteration_state(input, field, residual_norm, forcing_norm, &
            threshold, .true., .false., action, relative)
        call iteration_state_flags(action, converged, consistent, reason)
        call require(action == iteration_invalid, 'Invalid control must fail')
        call require(.not. converged, 'Invalid control is not converged')
        call require(trim(reason) == 'invalid-iteration-state', &
            'Invalid control reason')
    end subroutine invalid_control

end program test_iteration_state
