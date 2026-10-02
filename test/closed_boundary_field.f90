module field_eq_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer :: icall_eq = 0, nrad = 2, nzet = 2
    real(dp) :: rad(2), zet(2), rtf = 620.0_dp, btf = 53000.0_dp
end module field_eq_mod

module magdata_in_symfluxcoor_mod
    implicit none
    double precision :: btor, rbig
end module magdata_in_symfluxcoor_mod

module field_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
    real(dp) :: psi_offset, coefficient, axis_r = 620.0_dp, axis_z
contains
    subroutine field_eq(r, phi, z, br, bp, bz, brr, brp, brz, bpr, bpp, bpz, &
            bzr, bzp, bzz)
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: br, bp, bz, brr, brp, brz, bpr, bpp, bpz, bzr, bzp, bzz
        real(dp), parameter :: fpol = 620.0_dp*53000.0_dp
        psif = psi_offset+coefficient*((r-axis_r)**2+(z-axis_z)**2)
        dpsidr = 2.0_dp*coefficient*(r-axis_r)
        dpsidz = 2.0_dp*coefficient*(z-axis_z)
        d2psidr2 = 2.0_dp*coefficient
        d2psidrdz = 0.0_dp
        d2psidz2 = 2.0_dp*coefficient
        br = -dpsidz/r
        bp = fpol/r
        bz = dpsidr/r
        brr = dpsidz/r**2
        brp = 0.0_dp*phi
        brz = -2.0_dp*coefficient/r
        bpr = -fpol/r**2
        bpp = 0.0_dp
        bpz = 0.0_dp
        bzr = 2.0_dp*coefficient*axis_r/r**2
        bzp = 0.0_dp
        bzz = 0.0_dp
    end subroutine field_eq
end module field_sub
