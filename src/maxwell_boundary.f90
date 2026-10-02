module maxwell_boundary_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: validate_outer_polygon, read_mesh_exterior_scale, validate_exterior_scale
contains
    subroutine validate_exterior_scale(scale)
        real(dp), intent(in) :: scale
        if (.not. ieee_is_finite(scale)) error stop 'Maxwell exterior scale must be finite'
        if (scale <= 1.0_dp) error stop 'Maxwell exterior scale must exceed one'
    end subroutine validate_exterior_scale

    subroutine validate_outer_polygon(inner_R, inner_Z, outer_R, outer_Z, scale)
        real(dp), intent(in) :: inner_R(:), inner_Z(:), outer_R(:), outer_Z(:), scale
        integer :: i, j, next, n
        real(dp) :: cross, tolerance, span
        call validate_exterior_scale(scale)
        if (size(inner_R) /= size(inner_Z)) error stop 'Plasma polygon size mismatch'
        if (size(inner_R) < 3) error stop 'Plasma polygon needs three vertices'
        n = size(outer_R)
        if (n /= size(outer_Z)) error stop 'Maxwell polygon size mismatch'
        if (n < 3) error stop 'Maxwell polygon needs three vertices'
        if (.not. all(ieee_is_finite(inner_R))) error stop 'Nonfinite plasma R'
        if (.not. all(ieee_is_finite(inner_Z))) error stop 'Nonfinite plasma Z'
        if (.not. all(ieee_is_finite(outer_R))) error stop 'Nonfinite exterior R'
        if (.not. all(ieee_is_finite(outer_Z))) error stop 'Nonfinite exterior Z'
        if (minval(outer_R) <= 0.0_dp) error stop 'Maxwell exterior crosses R=0'
        span = max(maxval(outer_R)-minval(outer_R), &
            maxval(outer_Z)-minval(outer_Z))
        if (span <= 0.0_dp) error stop 'Degenerate Maxwell exterior'
        tolerance = 64.0_dp*epsilon(1.0_dp)*span**2
        do i = 1, n
            next = mod(i, n)+1
            if (hypot(outer_R(next)-outer_R(i), outer_Z(next)-outer_Z(i)) <= &
                64.0_dp*epsilon(1.0_dp)*span) error stop 'Degenerate Maxwell edge'
            do j = 1, size(inner_R)
                cross = (outer_R(next)-outer_R(i))*(inner_Z(j)-outer_Z(i)) &
                    -(outer_Z(next)-outer_Z(i))*(inner_R(j)-outer_R(i))
                if (cross <= tolerance) &
                    error stop 'Maxwell exterior polygon must strictly contain plasma'
            end do
        end do
    end subroutine validate_outer_polygon

    subroutine read_mesh_exterior_scale(h5id, dataset, requested, effective)
        use hdf5_tools, only: HID_T, h5_obj_exists, h5_get
        integer(HID_T), intent(in) :: h5id
        character(len=*), intent(in) :: dataset
        real(dp), intent(in) :: requested
        real(dp), intent(out) :: effective
        logical :: present_scale
        call validate_exterior_scale(requested)
        ! Reset on every read; historical mesh files used the fixed scale two.
        effective = 2.0_dp
        call h5_obj_exists(h5id, trim(dataset)//'/maxwell_outer_scale', present_scale)
        if (present_scale) &
            call h5_get(h5id, trim(dataset)//'/maxwell_outer_scale', effective)
        if (.not. ieee_is_finite(effective)) error stop 'Invalid cached Maxwell exterior scale'
        if (effective /= requested) &
            error stop 'Maxwell exterior scale mismatch: remesh and rebuild preconditioner'
    end subroutine read_mesh_exterior_scale
end module maxwell_boundary_m
