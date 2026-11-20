!> Gradient of the equilibrium toroidal current from curr0_geqdsk against central differences
!> of its own j0_phi = c (R p'(psi) + F F'(psi) / (4 pi R)) on a synthetic equilibrium.
program test_curr0_derivatives
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mephit_mesh, only: equil, curr0_geqdsk
    implicit none
    integer, parameter :: nw = 257
    real(dp), parameter :: R0 = 165.0_dp, a = 50.0_dp, F0 = -3.3e6_dp, F1 = 2.0e-3_dp, &
        F2 = -1.0e-9_dp, p1 = 3.0e-4_dp, p2 = -2.0e-10_dp, psi0 = 4.0e5_dp, h = 1.0e-3_dp
    real(dp) :: psi_max, R, Z, B0(3), dB0_dR(3), dB0_dZ(3), j0(3), dj0_dR(3), dj0_dZ(3)
    real(dp) :: jp(3), jm(3), fd_R, fd_Z, err
    integer :: k, kp
    logical :: ok

    psi_max = psi0 * 1.2_dp
    allocate(equil%psi_eqd(nw), equil%pprime(nw), equil%fpol(nw), equil%ffprim(nw), &
        equil%fprime(nw))
    equil%psi_eqd = [(psi_max * (k - 1) / (nw - 1), k = 1, nw)]
    equil%fpol = F0 + F1 * equil%psi_eqd + F2 * equil%psi_eqd ** 2
    equil%fprime = F1 + 2 * F2 * equil%psi_eqd
    equil%ffprim = equil%fpol * equil%fprime
    equil%pprime = p1 + p2 * equil%psi_eqd
    ok = .true.
    do kp = 1, 6
        R = R0 + 0.6_dp * a * cos(1.1_dp * kp)
        Z = 0.6_dp * a * sin(1.1_dp * kp)
        call j0_at(R, Z, j0, dj0_dR, dj0_dZ)
        call j0_at(R + h, Z, jp)
        call j0_at(R - h, Z, jm)
        fd_R = (jp(2) - jm(2)) / (2 * h)
        call j0_at(R, Z + h, jp)
        call j0_at(R, Z - h, jm)
        fd_Z = (jp(2) - jm(2)) / (2 * h)
        err = max(abs(dj0_dR(2) - fd_R), abs(dj0_dZ(2) - fd_Z)) / max(abs(fd_R), abs(fd_Z))
        if (err > 1.0e-5_dp) then
            ok = .false.
            print '("R ", f8.3, " Z ", f8.3, ": dj0phi/dR ", es12.4, " vs ", es12.4, &
                &", dj0phi/dZ ", es12.4, " vs ", es12.4)', R, Z, dj0_dR(2), fd_R, dj0_dZ(2), fd_Z
        end if
    end do
    if (.not. ok) error stop 'equilibrium current gradient does not match j0'

contains

    !> psi = psi0 (1 - ((R - R0)^2 + Z^2) / a^2), B_R = -d_Z psi / R, B_Z = d_R psi / R
    subroutine j0_at(R, Z, j0, dj0_dR, dj0_dZ)
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: j0(3)
        real(dp), intent(out), optional :: dj0_dR(3), dj0_dZ(3)
        real(dp) :: psi, B(3), dB_dR(3), dB_dZ(3), dR(3), dZ(3)

        psi = psi0 * (1 - ((R - R0) ** 2 + Z ** 2) / a ** 2)
        B = [2 * psi0 * Z / (a ** 2 * R), 0.0_dp, -2 * psi0 * (R - R0) / (a ** 2 * R)]
        dB_dR = 0.0_dp
        dB_dZ = 0.0_dp
        call curr0_geqdsk(R, psi, B, dB_dR, dB_dZ, j0, dR, dZ)
        if (present(dj0_dR)) dj0_dR = dR
        if (present(dj0_dZ)) dj0_dZ = dZ
    end subroutine j0_at
end program test_curr0_derivatives
