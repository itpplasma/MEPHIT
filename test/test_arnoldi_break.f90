!> Arnoldi eigenvalues of a diagonal operator whose dominant eigenvector is barely in the
!> start vector, so that the first Ritz values all lie below the threshold.
program test_arnoldi_break
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mephit_util, only: arnoldi_break
    implicit none
    integer, parameter :: ndim = 20
    real(dp), parameter :: lambda_max = -1.002_dp
    complex(dp) :: diag(ndim), start_vector(ndim)
    complex(dp), allocatable :: eigvals(:)
    integer :: k, ierr, nritz

    diag(1) = lambda_max
    start_vector(1) = 1d-2
    do k = 2, ndim
        diag(k) = 0.4_dp * real(k - 1, dp) / real(ndim, dp)
        start_vector(k) = 1d0
    end do
    call arnoldi_break(ndim, ndim, 3, 0.5_dp, 1d-9, apply_diag, start_vector, ierr, nritz, eigvals)
    if (ierr /= 0) error stop 'arnoldi_break returned an error'
    if (nritz /= 1) error stop 'arnoldi_break missed the eigenvalue above the threshold'
    if (abs(eigvals(1) - lambda_max) > 1d-9) error stop 'arnoldi_break eigenvalue is wrong'

contains

    subroutine apply_diag(old_val, new_val)
        complex(dp), intent(in) :: old_val(:)
        complex(dp), intent(out) :: new_val(:)

        new_val(:) = diag * old_val
    end subroutine apply_diag
end program test_arnoldi_break
