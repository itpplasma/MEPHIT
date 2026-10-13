!> Resistive layer width against an independently evaluated reference and its scaling.
program test_layer_width
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mephit_util, only: resistive_layer_width
    implicit none
    ! reference from the formula in resistive_layer_width evaluated in Python for
    ! B = 53 kG, R = 620 cm, n_e = 5e13 cm^-3, T_e = 1 keV, Z = 1, m_i = 2.5, n = 3,
    ! q = 3, dq/dr = 0.1 cm^-1 (v_A = 1.034e7 m/s, eta = 2.489e-8 Ohm m)
    real(dp), parameter :: reference = 0.10589984240932555_dp
    real(dp) :: w
    logical :: ok

    ok = .true.
    w = resistive_layer_width(53000d0, 620d0, 5d13, 1d3, 1d0, 2.5d0, 3, 3d0, 0.1d0)
    ok = ok .and. abs(w/reference - 1d0) < 1d-9
    ! delta ~ (q')^(-1/3): eight times the shear halves the width
    w = resistive_layer_width(53000d0, 620d0, 5d13, 1d3, 1d0, 2.5d0, 3, 3d0, 0.8d0)
    ok = ok .and. abs(w/reference - 0.5d0) < 1d-9
    ! delta ~ n^(-1/3)
    w = resistive_layer_width(53000d0, 620d0, 5d13, 1d3, 1d0, 2.5d0, 24, 3d0, 0.1d0)
    ok = ok .and. abs(w/reference - 0.5d0) < 1d-9
    ! delta ~ R^(1/3) at fixed B, q'
    w = resistive_layer_width(53000d0, 8d0*620d0, 5d13, 1d3, 1d0, 2.5d0, 3, 3d0, 0.1d0)
    ok = ok .and. abs(w/reference - 2d0) < 1d-9
    ! a colder, denser layer is wider: eta ~ T^(-3/2)
    w = resistive_layer_width(53000d0, 620d0, 5d13, 100d0, 1d0, 2.5d0, 3, 3d0, 0.1d0)
    ok = ok .and. w > 3d0*reference
    if (.not. ok) error stop 'resistive layer width test failed'
end program test_layer_width
