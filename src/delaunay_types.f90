module delaunay_types
    use iso_fortran_env, only: dp => real64
    implicit none
    
    private
    public :: point_t, triangle_t, edge_t, mesh_t
    public :: create_mesh, destroy_mesh, resize_mesh
    public :: add_point, add_triangle, add_edge
    public :: is_valid_triangle, is_valid_edge
    
    ! Point type with coordinates and ID
    type :: point_t
        real(dp) :: x, y
        integer :: id
        logical :: valid = .true.
    end type point_t
    
    ! Triangle type with vertices, neighbors, and validity
    type :: triangle_t
        integer :: vertices(3)   ! Point indices (1-based)
        integer :: neighbors(3)  ! Neighbor triangle indices (1-based, 0 = no neighbor)
        logical :: valid = .true.
    end type triangle_t
    
    ! Edge type with endpoints and constraint flag
    type :: edge_t
        integer :: endpoints(2)  ! Point indices (1-based)
        logical :: constrained = .false.
        logical :: valid = .true.
    end type edge_t
    
    ! Main mesh type containing all geometric data
    type :: mesh_t
        type(point_t), allocatable :: points(:)
        type(triangle_t), allocatable :: triangles(:)
        type(edge_t), allocatable :: edges(:)
        
        integer :: npoints = 0
        integer :: ntriangles = 0
        integer :: nedges = 0
        
        integer :: max_points = 0
        integer :: max_triangles = 0
        integer :: max_edges = 0
        
        ! Super-triangle vertices (for removal after triangulation)
        integer :: super_vertices(3) = 0
    end type mesh_t
    
contains

subroutine create_mesh(mesh, max_points, max_triangles, max_edges)
    !> Initialize mesh with given maximum capacities
    type(mesh_t), intent(out) :: mesh
    integer, intent(in) :: max_points, max_triangles, max_edges
    
    mesh%max_points = max_points
    mesh%max_triangles = max_triangles
    mesh%max_edges = max_edges
    
    mesh%npoints = 0
    mesh%ntriangles = 0
    mesh%nedges = 0
    
    allocate(mesh%points(max_points))
    allocate(mesh%triangles(max_triangles))
    allocate(mesh%edges(max_edges))
    
    mesh%super_vertices = 0
end subroutine create_mesh

subroutine destroy_mesh(mesh)
    !> Clean up mesh memory
    type(mesh_t), intent(inout) :: mesh
    
    if (allocated(mesh%points)) deallocate(mesh%points)
    if (allocated(mesh%triangles)) deallocate(mesh%triangles)
    if (allocated(mesh%edges)) deallocate(mesh%edges)
    
    mesh%npoints = 0
    mesh%ntriangles = 0
    mesh%nedges = 0
    mesh%max_points = 0
    mesh%max_triangles = 0
    mesh%max_edges = 0
    mesh%super_vertices = 0
end subroutine destroy_mesh

subroutine resize_mesh(mesh, new_max_points, new_max_triangles, new_max_edges)
    !> Resize mesh arrays if needed
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: new_max_points, new_max_triangles, new_max_edges
    
    type(point_t), allocatable :: temp_points(:)
    type(triangle_t), allocatable :: temp_triangles(:)
    type(edge_t), allocatable :: temp_edges(:)
    
    ! Resize points if needed
    if (new_max_points > mesh%max_points) then
        allocate(temp_points(new_max_points))
        if (mesh%npoints > 0) then
            temp_points(1:mesh%npoints) = mesh%points(1:mesh%npoints)
        end if
        deallocate(mesh%points)
        mesh%points = temp_points
        mesh%max_points = new_max_points
    end if
    
    ! Resize triangles if needed
    if (new_max_triangles > mesh%max_triangles) then
        allocate(temp_triangles(new_max_triangles))
        if (mesh%ntriangles > 0) then
            temp_triangles(1:mesh%ntriangles) = mesh%triangles(1:mesh%ntriangles)
        end if
        deallocate(mesh%triangles)
        mesh%triangles = temp_triangles
        mesh%max_triangles = new_max_triangles
    end if
    
    ! Resize edges if needed
    if (new_max_edges > mesh%max_edges) then
        allocate(temp_edges(new_max_edges))
        if (mesh%nedges > 0) then
            temp_edges(1:mesh%nedges) = mesh%edges(1:mesh%nedges)
        end if
        deallocate(mesh%edges)
        mesh%edges = temp_edges
        mesh%max_edges = new_max_edges
    end if
end subroutine resize_mesh

function add_point(mesh, x, y, point_id) result(index)
    !> Add a point to the mesh and return its index
    type(mesh_t), intent(inout) :: mesh
    real(dp), intent(in) :: x, y
    integer, intent(in), optional :: point_id
    integer :: index
    
    ! Resize if necessary
    if (mesh%npoints >= mesh%max_points) then
        call resize_mesh(mesh, mesh%max_points * 2, mesh%max_triangles, mesh%max_edges)
    end if
    
    mesh%npoints = mesh%npoints + 1
    index = mesh%npoints
    
    mesh%points(index)%x = x
    mesh%points(index)%y = y
    mesh%points(index)%valid = .true.
    
    if (present(point_id)) then
        mesh%points(index)%id = point_id
    else
        mesh%points(index)%id = index
    end if
end function add_point

function add_triangle(mesh, v1, v2, v3) result(index)
    !> Add a triangle to the mesh and return its index
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: v1, v2, v3
    integer :: index
    
    ! Resize if necessary
    if (mesh%ntriangles >= mesh%max_triangles) then
        call resize_mesh(mesh, mesh%max_points, mesh%max_triangles * 2, mesh%max_edges)
    end if
    
    mesh%ntriangles = mesh%ntriangles + 1
    index = mesh%ntriangles
    
    mesh%triangles(index)%vertices(1) = v1
    mesh%triangles(index)%vertices(2) = v2
    mesh%triangles(index)%vertices(3) = v3
    mesh%triangles(index)%neighbors = 0  ! No neighbors initially
    mesh%triangles(index)%valid = .true.
end function add_triangle

function add_edge(mesh, v1, v2, constrained) result(index)
    !> Add an edge to the mesh and return its index
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: v1, v2
    logical, intent(in), optional :: constrained
    integer :: index
    
    ! Resize if necessary
    if (mesh%nedges >= mesh%max_edges) then
        call resize_mesh(mesh, mesh%max_points, mesh%max_triangles, mesh%max_edges * 2)
    end if
    
    mesh%nedges = mesh%nedges + 1
    index = mesh%nedges
    
    mesh%edges(index)%endpoints(1) = v1
    mesh%edges(index)%endpoints(2) = v2
    mesh%edges(index)%valid = .true.
    
    if (present(constrained)) then
        mesh%edges(index)%constrained = constrained
    else
        mesh%edges(index)%constrained = .false.
    end if
end function add_edge

logical function is_valid_triangle(mesh, tri_idx)
    !> Check if triangle index is valid and triangle is valid
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: tri_idx
    
    is_valid_triangle = .false.
    
    if (tri_idx < 1 .or. tri_idx > mesh%ntriangles) return
    if (.not. mesh%triangles(tri_idx)%valid) return
    
    ! Check vertex indices
    if (any(mesh%triangles(tri_idx)%vertices < 1) .or. &
        any(mesh%triangles(tri_idx)%vertices > mesh%npoints)) return
    
    ! Check vertex validity
    if (.not. all(mesh%points(mesh%triangles(tri_idx)%vertices)%valid)) return
    
    is_valid_triangle = .true.
end function is_valid_triangle

logical function is_valid_edge(mesh, edge_idx)
    !> Check if edge index is valid and edge is valid
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: edge_idx
    
    is_valid_edge = .false.
    
    if (edge_idx < 1 .or. edge_idx > mesh%nedges) return
    if (.not. mesh%edges(edge_idx)%valid) return
    
    ! Check endpoint indices
    if (any(mesh%edges(edge_idx)%endpoints < 1) .or. &
        any(mesh%edges(edge_idx)%endpoints > mesh%npoints)) return
    
    ! Check endpoint validity
    if (.not. all(mesh%points(mesh%edges(edge_idx)%endpoints)%valid)) return
    
    is_valid_edge = .true.
end function is_valid_edge

end module delaunay_types