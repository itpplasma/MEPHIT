module closed_surface_boundary_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: find_closed_boundary
    abstract interface
        function radial_flux_t(r) result(psi)
            import dp
            real(dp), intent(in) :: r
            real(dp) :: psi
        end function radial_flux_t
    end interface
contains
    subroutine find_closed_boundary(flux, axis_r, box_r, target, edge_r, stat)
        procedure(radial_flux_t) :: flux
        real(dp), intent(in) :: axis_r, box_r, target
        real(dp), intent(out) :: edge_r
        integer, intent(out) :: stat
        real(dp) :: left, right, mid, fleft, fright, fmid
        integer :: i

        stat = 1
        edge_r = axis_r
        if (box_r <= axis_r) return
        if (.not. ieee_is_finite(target)) return
        fleft = flux(axis_r)-target
        if (.not. ieee_is_finite(fleft)) return
        if (fleft == 0.0_dp) return
        left = axis_r
        do i = 1, 256
            right = axis_r+(box_r-axis_r)*real(i, dp)/256.0_dp
            fright = flux(right)-target
            if (.not. ieee_is_finite(fright)) return
            if (sign(1.0_dp, fleft) /= sign(1.0_dp, fright)) exit
            if (fright == 0.0_dp) exit
            left = right
            fleft = fright
        end do
        if (i > 256) return
        do i = 1, 80
            mid = 0.5_dp*(left+right)
            fmid = flux(mid)-target
            if (.not. ieee_is_finite(fmid)) return
            if (abs(right-left) <= 1.0e-12_dp*max(1.0_dp, abs(mid))) exit
            if (fmid == 0.0_dp) exit
            if (sign(1.0_dp, fmid) == sign(1.0_dp, fleft)) then
                left = mid
                fleft = fmid
            else
                right = mid
            end if
        end do
        edge_r = mid
        stat = 0
    end subroutine find_closed_boundary
end module closed_surface_boundary_m
