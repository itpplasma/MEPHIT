module mephit_triangulation_fortran
    !> Direct Fortran triangulation interface for MEPHIT
    !> Replaces external TRIANGLE library with native Fortran implementation
    
    use iso_fortran_env, only: dp => real64
    use iso_c_binding, only: c_char, c_null_char
    use triangulation_fortran
    implicit none
    
    private
    public :: FEM_triangulate_external_fortran
    public :: write_freefem_mesh_fortran
    
contains

subroutine FEM_triangulate_external_fortran(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    !> Direct Fortran triangulation interface compatible with MEPHIT
    !> Replaces C interface with native Fortran implementation
    
    integer, intent(in) :: npt_inner, npt_outer
    real(dp), intent(in) :: node_R(:), node_Z(:)
    real(dp), intent(in) :: R_O, Z_O  ! Unused in current implementation
    character(len=*), intent(in) :: fname
    
    integer :: total_points, i
    real(dp), allocatable :: points(:,:)
    integer, allocatable :: segments(:,:)
    type(triangulation_result_t) :: result
    
    total_points = npt_inner + npt_outer
    
    ! Allocate input arrays
    allocate(points(2, total_points))
    allocate(segments(2, npt_outer))
    
    ! Pack input points
    do i = 1, total_points
        points(1, i) = node_R(i)
        points(2, i) = node_Z(i)
    end do
    
    ! Create boundary segments (outer boundary only)
    do i = 1, npt_outer
        segments(1, i) = npt_inner + i
        segments(2, i) = npt_inner + mod(i, npt_outer) + 1
    end do
    
    ! Perform triangulation
    call triangulate_fortran(points, segments, result)
    
    ! Write output to FreeFEM format
    call write_freefem_mesh_fortran(result, fname)
    
    ! Clean up
    call cleanup_triangulation(result)
    deallocate(points, segments)
    
    write(*,'(A,I0,A,I0,A,A)') 'Fortran triangulation: ', result%npoints, &
        ' points, ', result%ntriangles, ' triangles → ', trim(fname)
    
end subroutine FEM_triangulate_external_fortran

subroutine write_freefem_mesh_fortran(result, filename)
    !> Write triangulation result to FreeFEM mesh format
    type(triangulation_result_t), intent(in) :: result
    character(len=*), intent(in) :: filename
    
    integer :: unit, i
    
    ! Open file for writing
    open(newunit=unit, file=filename, status='replace', action='write')
    
    ! Write header: npoints ntriangles 0
    write(unit, '(I0,1X,I0,1X,I0)') result%npoints, result%ntriangles, 0
    
    ! Write points: R Z label
    do i = 1, result%npoints
        write(unit, '(ES15.8,1X,ES15.8,1X,I0)') result%points(1, i), result%points(2, i), 0
    end do
    
    ! Write triangles: v1 v2 v3 label
    do i = 1, result%ntriangles
        write(unit, '(I0,1X,I0,1X,I0,1X,I0)') result%triangles(1, i), &
            result%triangles(2, i), result%triangles(3, i), 0
    end do
    
    close(unit)
    
end subroutine write_freefem_mesh_fortran

end module mephit_triangulation_fortran