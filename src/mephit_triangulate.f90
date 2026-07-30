module mephit_triangulate
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: FEM_triangulate_external

contains

  !> Triangulate the region between two closed boundary loops and write the
  !> result in FreeFem msh format. Replaces Shewchuk's Triangle with flags
  !> "BejnpqYz" by fortfem's clean-room reimplementation.
  subroutine FEM_triangulate_external(npt_inner, npt_outer, bdry_R, bdry_Z, &
    R_hole, Z_hole, fname)
    use triangle_compat, only: tc_result_t, triangulate_compat
    integer, intent(in) :: npt_inner, npt_outer
    real(dp), intent(in) :: bdry_R(npt_inner + npt_outer), bdry_Z(npt_inner + npt_outer)
    real(dp), intent(in) :: R_hole, Z_hole
    character(len = *), intent(in) :: fname
    type(tc_result_t) :: res
    real(dp) :: points(2, npt_inner + npt_outer), holes(2, 1)
    integer :: segments(2, npt_inner + npt_outer)
    integer :: fid, k, stat

    points(1, :) = bdry_R
    points(2, :) = bdry_Z
    do k = 1, npt_inner
      segments(:, k) = [k, mod(k, npt_inner) + 1]
    end do
    do k = 1, npt_outer
      segments(:, npt_inner + k) = [npt_inner + k, npt_inner + mod(k, npt_outer) + 1]
    end do
    holes(:, 1) = [R_hole, Z_hole]
    call triangulate_compat(points, segments, holes, res, stat, &
      min_angle = 20d0, quality = .true., nobisect = 1)
    if (stat /= 0) then
      error stop 'FEM_triangulate_external: triangulate_compat failed'
    end if
    open(newunit = fid, file = fname, status = 'replace', form = 'formatted', &
      action = 'write')
    write (fid, '(i0, 2(1x, i0))') res%npoints, res%ntriangles, res%nsegments
    do k = 1, res%npoints
      write (fid, '(a, 1x, a, " 0")') &
        c_double_string(res%points(1, k)), c_double_string(res%points(2, k))
    end do
    do k = 1, res%ntriangles
      write (fid, '(i0, 2(1x, i0), " 1")') res%triangles(:, k)
    end do
    do k = 1, res%nsegments
      write (fid, '(i0, 1x, i0, " 2")') res%segments(:, k)
    end do
    close(fid)
  end subroutine FEM_triangulate_external

  !> Render a double like C's printf("%.16e", ...) to keep the msh file
  !> format byte-identical to the former C writer.
  function c_double_string(val) result(str)
    real(dp), intent(in) :: val
    character(len = :), allocatable :: str
    character(len = 24) :: buf
    integer :: k

    write (buf, '(es23.16e2)') val
    str = trim(adjustl(buf))
    do k = 1, len(str)
      if (str(k:k) == 'E') str(k:k) = 'e'
    end do
  end function c_double_string

end module mephit_triangulate
