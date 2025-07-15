module triangulation_fortran
    use iso_fortran_env, only: dp => real64
    use delaunay_types
    use bowyer_watson
    implicit none
    
    private
    public :: triangulation_result_t, triangulate_fortran, triangulate_with_hole_fortran
    public :: triangulate_with_quality_fortran, triangulate_triangle_lib
    public :: cleanup_triangulation
    
    ! Type to hold triangulation results - matches TRIANGLE output format
    type :: triangulation_result_t
        integer :: npoints              ! Number of points
        integer :: ntriangles           ! Number of triangles
        integer :: nsegments            ! Number of segments
        real(dp), allocatable :: points(:,:)      ! Points (2, npoints)
        integer, allocatable :: triangles(:,:)    ! Triangles (3, ntriangles)
        integer, allocatable :: segments(:,:)     ! Segments (2, nsegments)
        integer, allocatable :: neighbors(:,:)    ! Neighbors (3, ntriangles)
    end type triangulation_result_t
    
contains

subroutine triangulate_fortran(points, segments, result, status)
    !> Main triangulation routine using Bowyer-Watson algorithm
    real(dp), intent(in) :: points(:,:)      ! Input points (2, npoints)
    integer, intent(in) :: segments(:,:)     ! Input segments (2, nsegments)
    type(triangulation_result_t), intent(out) :: result
    integer, intent(out), optional :: status
    
    type(mesh_t) :: mesh
    integer :: i, valid_triangles, valid_points
    
    if (present(status)) status = 0
    
    ! Perform Delaunay triangulation
    call delaunay_triangulate(points, mesh)
    
    ! Convert mesh to result format
    call mesh_to_result(mesh, segments, result)
    
    call destroy_mesh(mesh)
    
end subroutine triangulate_fortran

subroutine triangulate_with_hole_fortran(points, segments, hole_point, result, status)
    !> Triangulation with hole - equivalent to TRIANGLE with hole specification
    real(dp), intent(in) :: points(:,:)      ! Input points (2, npoints)
    integer, intent(in) :: segments(:,:)     ! Input segments (2, nsegments)
    real(dp), intent(in) :: hole_point(:)    ! Hole point (2)
    type(triangulation_result_t), intent(out) :: result
    integer, intent(out), optional :: status
    
    if (present(status)) status = 0
    
    ! Placeholder - to be implemented
    call allocate_result(result, size(points, 2), 0, size(segments, 2))
    
    ! Copy input points
    result%points = points
    result%segments = segments
    
    ! TODO: Implement Delaunay triangulation with hole
    if (present(status)) status = -1
    
end subroutine triangulate_with_hole_fortran

subroutine triangulate_with_quality_fortran(points, segments, min_angle, result, status)
    !> Triangulation with quality constraints - equivalent to TRIANGLE's 'q' option
    real(dp), intent(in) :: points(:,:)      ! Input points (2, npoints)
    integer, intent(in) :: segments(:,:)     ! Input segments (2, nsegments)
    real(dp), intent(in) :: min_angle        ! Minimum angle in degrees
    type(triangulation_result_t), intent(out) :: result
    integer, intent(out), optional :: status
    
    if (present(status)) status = 0
    
    ! Placeholder - to be implemented
    call allocate_result(result, size(points, 2), 0, size(segments, 2))
    
    ! Copy input points
    result%points = points
    result%segments = segments
    
    ! TODO: Implement Delaunay triangulation with quality constraints
    if (present(status)) status = -1
    
end subroutine triangulate_with_quality_fortran

subroutine triangulate_triangle_lib(points, segments, result, status)
    !> Wrapper for original TRIANGLE library - for comparison testing
    real(dp), intent(in) :: points(:,:)      ! Input points (2, npoints)
    integer, intent(in) :: segments(:,:)     ! Input segments (2, nsegments)
    type(triangulation_result_t), intent(out) :: result
    integer, intent(out), optional :: status
    
    if (present(status)) status = 0
    
    ! Call the existing TRIANGLE library through C interface
    ! This will be used for comparison in tests
    
    ! TODO: Implement C interface call to TRIANGLE
    ! For now, return error status
    if (present(status)) status = -1
    
end subroutine triangulate_triangle_lib

subroutine allocate_result(result, npoints, ntriangles, nsegments)
    !> Helper to allocate result arrays
    type(triangulation_result_t), intent(out) :: result
    integer, intent(in) :: npoints, ntriangles, nsegments
    
    result%npoints = npoints
    result%ntriangles = ntriangles
    result%nsegments = nsegments
    
    allocate(result%points(2, npoints))
    allocate(result%triangles(3, ntriangles))
    allocate(result%segments(2, nsegments))
    allocate(result%neighbors(3, ntriangles))
    
end subroutine allocate_result

subroutine mesh_to_result(mesh, input_segments, result)
    !> Convert internal mesh format to triangulation result format
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: input_segments(:,:)
    type(triangulation_result_t), intent(out) :: result
    
    integer :: i, valid_points, valid_triangles
    
    ! Count valid points and triangles
    valid_points = 0
    do i = 1, mesh%npoints
        if (mesh%points(i)%valid) valid_points = valid_points + 1
    end do
    
    valid_triangles = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) valid_triangles = valid_triangles + 1
    end do
    
    ! Allocate result arrays
    call allocate_result(result, valid_points, valid_triangles, size(input_segments, 2))
    
    ! Copy valid points
    valid_points = 0
    do i = 1, mesh%npoints
        if (mesh%points(i)%valid) then
            valid_points = valid_points + 1
            result%points(1, valid_points) = mesh%points(i)%x
            result%points(2, valid_points) = mesh%points(i)%y
        end if
    end do
    
    ! Copy valid triangles (need to remap vertex indices)
    valid_triangles = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            valid_triangles = valid_triangles + 1
            ! Map original vertex indices to new indices
            result%triangles(1, valid_triangles) = remap_vertex_index(mesh, mesh%triangles(i)%vertices(1))
            result%triangles(2, valid_triangles) = remap_vertex_index(mesh, mesh%triangles(i)%vertices(2))
            result%triangles(3, valid_triangles) = remap_vertex_index(mesh, mesh%triangles(i)%vertices(3))
        end if
    end do
    
    ! Copy input segments
    result%segments = input_segments
    
end subroutine mesh_to_result

integer function remap_vertex_index(mesh, original_idx)
    !> Map original vertex index to new index (accounting for removed vertices)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: original_idx
    
    integer :: i, valid_count
    
    valid_count = 0
    do i = 1, original_idx
        if (mesh%points(i)%valid) then
            valid_count = valid_count + 1
        end if
    end do
    
    remap_vertex_index = valid_count
end function remap_vertex_index

subroutine cleanup_triangulation(result)
    !> Clean up allocated arrays
    type(triangulation_result_t), intent(inout) :: result
    
    if (allocated(result%points)) deallocate(result%points)
    if (allocated(result%triangles)) deallocate(result%triangles)
    if (allocated(result%segments)) deallocate(result%segments)
    if (allocated(result%neighbors)) deallocate(result%neighbors)
    
    result%npoints = 0
    result%ntriangles = 0
    result%nsegments = 0
    
end subroutine cleanup_triangulation

end module triangulation_fortran