!> A zero resonant layer width is detected exactly where the layer is used.
program test_resonant_layer_check
    use mephit_mesh, only: zero_width_resonance
    implicit none
    logical :: ok

    ok = .true.
    ! layer used: first zero width (m = 3) is reported, m_res_min = 2
    ok = ok .and. zero_width_resonance(.true., 2, [1d0, 0d0, 0d0]) == 3
    ok = ok .and. zero_width_resonance(.true., 2, [0d0, 1d0]) == 2
    ! all widths positive: none (m_res_min - 1)
    ok = ok .and. zero_width_resonance(.true., 2, [1d0, 2d0]) == 1
    ! layer not used (no damping, currn_model /= 2): zero widths are allowed
    ok = ok .and. zero_width_resonance(.false., 2, [0d0, 0d0]) == 1
    ! a negative or NaN width is not a layer either
    ok = ok .and. zero_width_resonance(.true., 3, [-1d0]) == 3
    ok = ok .and. zero_width_resonance(.true., 3, [ieee_nan()]) == 3
    if (.not. ok) error stop 'resonant layer check test failed'
contains
    function ieee_nan() result(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        double precision :: x
        x = ieee_value(x, ieee_quiet_nan)
    end function ieee_nan
end program test_resonant_layer_check
