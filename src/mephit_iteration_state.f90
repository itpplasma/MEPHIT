module mephit_iteration_state
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: select_iteration_state, relative_iteration_norm, iteration_state_flags
    public :: finite_iteration_field
    public :: iteration_continue, iteration_converged, iteration_limit
    public :: iteration_single, iteration_invalid

    integer, parameter :: iteration_continue = 0, iteration_converged = 1
    integer, parameter :: iteration_limit = 2, iteration_single = 3
    integer, parameter :: iteration_invalid = 4

contains

    pure function finite_iteration_field(field) result(finite)
        complex(dp), intent(in) :: field(:)
        logical :: finite

        finite = all(ieee_is_finite(real(field, dp)))
        if (.not. finite) return
        finite = all(ieee_is_finite(aimag(field)))
    end function finite_iteration_field

    pure function relative_iteration_norm(value, forcing) result(relative)
        real(dp), intent(in) :: value, forcing
        real(dp) :: relative

        relative = huge(1.0_dp)
        if (.not. ieee_is_finite(value)) return
        if (.not. ieee_is_finite(forcing)) return
        if (value < 0.0_dp) return
        if (forcing < 0.0_dp) return
        if (forcing <= 0.0_dp) then
            ! With zero forcing, only an exactly zero residual is admissible.
            if (value <= 0.0_dp) relative = 0.0_dp
            return
        end if
        if (forcing < 1.0_dp) then
            if (value > huge(1.0_dp)*forcing) return
        end if
        relative = value/forcing
    end function relative_iteration_norm

    pure subroutine select_iteration_state(input, field, residual_norm, forcing_norm, &
            tolerance, exhausted, single, action, map_rel_err)
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(inout) :: field(:)
        real(dp), intent(in) :: residual_norm, forcing_norm, tolerance
        logical, intent(in) :: exhausted, single
        integer, intent(out) :: action
        real(dp), intent(out) :: map_rel_err

        if (size(input) /= size(field)) error stop 'Iteration field size mismatch'
        action = iteration_invalid
        map_rel_err = relative_iteration_norm(residual_norm, forcing_norm)
        if (.not. ieee_is_finite(residual_norm)) then
            field = input
            return
        end if
        if (.not. ieee_is_finite(forcing_norm)) then
            field = input
            return
        end if
        if (.not. ieee_is_finite(tolerance)) then
            field = input
            return
        end if
        if (min(residual_norm, forcing_norm) < 0.0_dp) then
            field = input
            return
        end if
        if (tolerance <= 0.0_dp) then
            field = input
            return
        end if
        if (.not. finite_iteration_field(field)) then
            field = input
            return
        end if
        if (.not. finite_iteration_field(input)) then
            field = input
            return
        end if
        if (single) then
            ! A one-map diagnostic retains its mapped output, without acceptance.
            action = iteration_single
            return
        end if
        field = input
        if (map_rel_err < tolerance) then
            action = iteration_converged
        else if (exhausted) then
            action = iteration_limit
        else
            action = iteration_continue
        end if
    end subroutine select_iteration_state

    pure subroutine iteration_state_flags(action, converged, consistent_state, reason)
        integer, intent(in) :: action
        logical, intent(out) :: converged, consistent_state
        character(len=*), intent(out) :: reason

        converged = .false.
        consistent_state = .true.
        select case (action)
        case (iteration_converged)
            converged = .true.
            reason = 'original-map-converged'
        case (iteration_limit)
            reason = 'iteration-limit'
        case (iteration_single)
            consistent_state = .false.
            reason = 'single-update'
        case (iteration_invalid)
            consistent_state = .false.
            reason = 'invalid-iteration-state'
        case default
            reason = 'continuing'
        end select
    end subroutine iteration_state_flags

end module mephit_iteration_state
