!> Regression test for refine_unit_partition with a resonance next to the edge.
program mephit_test_refine

  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mephit_mesh, only: refine_unit_partition

  implicit none

  real(dp), parameter :: res(2) = [0.9d0, 0.9996d0]
  real(dp), parameter :: fine_sep(2) = [2d-3, 2d-4]
  integer, parameter :: add_fine(2) = [2, 2]
  real(dp), parameter :: refinement(2) = [1.5d0, 1.5d0]
  real(dp), allocatable :: partition(:)
  integer :: nflux, kf_ref(2), k

  ! the nominal fine band of the second resonance, res + 2.5 * fine_sep = 1.0001,
  ! exceeds the edge at 1 and must be clipped instead of failing
  call refine_unit_partition(2, 2d-2, fine_sep, add_fine, refinement, res, .false., &
    nflux, partition, kf_ref)
  if (size(partition) /= nflux + 1) error stop 'partition size inconsistent with nflux'
  if (abs(partition(0)) > 1d-14) error stop 'partition does not start at 0'
  if (abs(partition(nflux) - 1d0) > 1d-14) error stop 'partition does not end at 1'
  do k = 1, nflux
    if (partition(k) <= partition(k - 1)) error stop 'partition not strictly increasing'
  end do
  do k = 1, 2
    if (.not. (partition(kf_ref(k) - 1) < res(k) .and. res(k) < partition(kf_ref(k)))) &
      error stop 'resonance not bracketed by its interval'
  end do
  ! the clipped region still resolves the layer with the fine separation
  if (partition(kf_ref(2)) - partition(kf_ref(2) - 1) > 1.0001d0 * fine_sep(2)) &
    error stop 'resonant interval coarser than fine separation'
  print '("refine_unit_partition edge clipping ok, nflux = ", i0)', nflux
end program mephit_test_refine
