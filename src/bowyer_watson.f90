module bowyer_watson
    use iso_fortran_env, only: dp => real64
    use delaunay_types
    use geometric_predicates
    implicit none
    
    private
    public :: delaunay_triangulate, insert_point
    public :: create_super_triangle, remove_super_triangle
    public :: find_cavity, fill_cavity
    
contains

subroutine delaunay_triangulate(input_points, mesh)
    !> Main Delaunay triangulation routine using Bowyer-Watson algorithm
    real(dp), intent(in) :: input_points(:,:)  ! (2, npoints)
    type(mesh_t), intent(out) :: mesh
    
    integer :: i, npoints
    integer :: super_idx, point_idx
    
    npoints = size(input_points, 2)
    
    ! Initialize mesh with appropriate size
    call create_mesh(mesh, npoints + 3, 6 * npoints, 3 * npoints)
    
    ! Create super-triangle containing all points
    call create_super_triangle(input_points, mesh, super_idx)
    
    ! Add input points to mesh
    do i = 1, npoints
        point_idx = add_point(mesh, input_points(1, i), input_points(2, i), i)
    end do
    
    ! Insert each point using Bowyer-Watson algorithm
    do i = 4, mesh%npoints  ! Start after super-triangle vertices
        call insert_point(mesh, i)
    end do
    
    ! Remove super-triangle and its associated triangles
    call remove_super_triangle(mesh)
    
end subroutine delaunay_triangulate

subroutine create_super_triangle(input_points, mesh, super_tri_idx)
    !> Create a super-triangle that contains all input points
    real(dp), intent(in) :: input_points(:,:)
    type(mesh_t), intent(inout) :: mesh
    integer, intent(out) :: super_tri_idx
    
    real(dp) :: min_x, max_x, min_y, max_y
    real(dp) :: dx, dy, delta_max, x_mid, y_mid
    integer :: p1, p2, p3
    
    ! Find bounding box of all points
    min_x = minval(input_points(1, :))
    max_x = maxval(input_points(1, :))
    min_y = minval(input_points(2, :))
    max_y = maxval(input_points(2, :))
    
    dx = max_x - min_x
    dy = max_y - min_y
    delta_max = max(dx, dy)
    x_mid = (min_x + max_x) / 2.0_dp
    y_mid = (min_y + max_y) / 2.0_dp
    
    ! Create super-triangle vertices (much larger than bounding box)
    p1 = add_point(mesh, x_mid - 20.0_dp * delta_max, y_mid - delta_max, -1)
    p2 = add_point(mesh, x_mid + 20.0_dp * delta_max, y_mid - delta_max, -2)
    p3 = add_point(mesh, x_mid, y_mid + 20.0_dp * delta_max, -3)
    
    ! Store super-triangle vertex indices
    mesh%super_vertices(1) = p1
    mesh%super_vertices(2) = p2
    mesh%super_vertices(3) = p3
    
    ! Create super-triangle
    super_tri_idx = add_triangle(mesh, p1, p2, p3)
    
end subroutine create_super_triangle

subroutine remove_super_triangle(mesh)
    !> Remove super-triangle vertices and all triangles containing them
    type(mesh_t), intent(inout) :: mesh
    
    integer :: i, j, v
    logical :: contains_super_vertex
    
    ! Mark triangles containing super-triangle vertices as invalid
    do i = 1, mesh%ntriangles
        if (.not. mesh%triangles(i)%valid) cycle
        
        contains_super_vertex = .false.
        do j = 1, 3
            v = mesh%triangles(i)%vertices(j)
            if (any(mesh%super_vertices == v)) then
                contains_super_vertex = .true.
                exit
            end if
        end do
        
        if (contains_super_vertex) then
            mesh%triangles(i)%valid = .false.
        end if
    end do
    
    ! Mark super-triangle vertices as invalid
    do i = 1, 3
        if (mesh%super_vertices(i) > 0 .and. mesh%super_vertices(i) <= mesh%npoints) then
            mesh%points(mesh%super_vertices(i))%valid = .false.
        end if
    end do
    
end subroutine remove_super_triangle

subroutine insert_point(mesh, point_idx)
    !> Insert a point into existing triangulation using Bowyer-Watson algorithm
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: point_idx
    
    integer, allocatable :: cavity_triangles(:)
    integer, allocatable :: cavity_edges(:,:)
    integer :: ncavity_triangles, ncavity_edges
    
    ! Find triangles whose circumcircles contain the new point
    call find_cavity(mesh, point_idx, cavity_triangles, ncavity_triangles)
    
    if (ncavity_triangles == 0) return  ! Point not inside any circumcircle
    
    ! Find the boundary of the cavity
    call find_cavity_boundary(mesh, cavity_triangles, ncavity_triangles, &
                             cavity_edges, ncavity_edges)
    
    ! Remove triangles in the cavity
    call remove_cavity_triangles(mesh, cavity_triangles, ncavity_triangles)
    
    ! Create new triangles connecting the point to the cavity boundary
    call fill_cavity(mesh, point_idx, cavity_edges, ncavity_edges)
    
end subroutine insert_point

subroutine find_cavity(mesh, point_idx, cavity_triangles, ncavity_triangles)
    !> Find all triangles whose circumcircles contain the given point
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: point_idx
    integer, allocatable, intent(out) :: cavity_triangles(:)
    integer, intent(out) :: ncavity_triangles
    
    type(point_t) :: point
    type(point_t) :: pa, pb, pc
    integer :: i
    logical, allocatable :: in_cavity(:)
    
    point = mesh%points(point_idx)
    
    allocate(in_cavity(mesh%ntriangles))
    in_cavity = .false.
    ncavity_triangles = 0
    
    ! Check each valid triangle
    do i = 1, mesh%ntriangles
        if (.not. mesh%triangles(i)%valid) cycle
        
        ! Get triangle vertices
        pa = mesh%points(mesh%triangles(i)%vertices(1))
        pb = mesh%points(mesh%triangles(i)%vertices(2))
        pc = mesh%points(mesh%triangles(i)%vertices(3))
        
        ! Test if point is inside circumcircle
        if (in_circle(pa, pb, pc, point)) then
            in_cavity(i) = .true.
            ncavity_triangles = ncavity_triangles + 1
        end if
    end do
    
    ! Collect cavity triangle indices
    allocate(cavity_triangles(ncavity_triangles))
    ncavity_triangles = 0
    do i = 1, mesh%ntriangles
        if (in_cavity(i)) then
            ncavity_triangles = ncavity_triangles + 1
            cavity_triangles(ncavity_triangles) = i
        end if
    end do
    
end subroutine find_cavity

subroutine find_cavity_boundary(mesh, cavity_triangles, ncavity_triangles, &
                                cavity_edges, ncavity_edges)
    !> Find boundary edges of the cavity (edges that belong to only one cavity triangle)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: cavity_triangles(:), ncavity_triangles
    integer, allocatable, intent(out) :: cavity_edges(:,:)
    integer, intent(out) :: ncavity_edges
    
    integer, allocatable :: all_edges(:,:)
    integer, allocatable :: edge_count(:)
    integer :: nedges_total, i, j, t, v1, v2
    integer :: edge_idx
    logical :: found
    
    ! Collect all edges from cavity triangles
    nedges_total = 3 * ncavity_triangles
    allocate(all_edges(2, nedges_total))
    allocate(edge_count(nedges_total))
    
    nedges_total = 0
    do i = 1, ncavity_triangles
        t = cavity_triangles(i)
        ! Add three edges of triangle t
        do j = 1, 3
            v1 = mesh%triangles(t)%vertices(j)
            v2 = mesh%triangles(t)%vertices(mod(j, 3) + 1)
            
            ! Ensure consistent edge ordering (smaller vertex first)
            if (v1 > v2) then
                call swap_integers(v1, v2)
            end if
            
            ! Check if edge already exists
            found = .false.
            do edge_idx = 1, nedges_total
                if (all_edges(1, edge_idx) == v1 .and. all_edges(2, edge_idx) == v2) then
                    edge_count(edge_idx) = edge_count(edge_idx) + 1
                    found = .true.
                    exit
                end if
            end do
            
            if (.not. found) then
                nedges_total = nedges_total + 1
                all_edges(1, nedges_total) = v1
                all_edges(2, nedges_total) = v2
                edge_count(nedges_total) = 1
            end if
        end do
    end do
    
    ! Count boundary edges (appear only once)
    ncavity_edges = 0
    do i = 1, nedges_total
        if (edge_count(i) == 1) then
            ncavity_edges = ncavity_edges + 1
        end if
    end do
    
    ! Collect boundary edges
    allocate(cavity_edges(2, ncavity_edges))
    ncavity_edges = 0
    do i = 1, nedges_total
        if (edge_count(i) == 1) then
            ncavity_edges = ncavity_edges + 1
            cavity_edges(1, ncavity_edges) = all_edges(1, i)
            cavity_edges(2, ncavity_edges) = all_edges(2, i)
        end if
    end do
    
end subroutine find_cavity_boundary

subroutine remove_cavity_triangles(mesh, cavity_triangles, ncavity_triangles)
    !> Mark cavity triangles as invalid
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: cavity_triangles(:), ncavity_triangles
    
    integer :: i, t
    
    do i = 1, ncavity_triangles
        t = cavity_triangles(i)
        mesh%triangles(t)%valid = .false.
    end do
    
end subroutine remove_cavity_triangles

subroutine fill_cavity(mesh, point_idx, cavity_edges, ncavity_edges)
    !> Create new triangles connecting the point to each cavity boundary edge
    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: point_idx
    integer, intent(in) :: cavity_edges(:,:), ncavity_edges
    
    integer :: i, v1, v2, new_tri
    type(point_t) :: p, p1, p2
    
    p = mesh%points(point_idx)
    
    do i = 1, ncavity_edges
        v1 = cavity_edges(1, i)
        v2 = cavity_edges(2, i)
        p1 = mesh%points(v1)
        p2 = mesh%points(v2)
        
        ! Ensure counter-clockwise orientation
        if (orientation(p1, p2, p) == ORIENTATION_CCW) then
            new_tri = add_triangle(mesh, v1, v2, point_idx)
        else
            new_tri = add_triangle(mesh, v2, v1, point_idx)
        end if
    end do
    
end subroutine fill_cavity

subroutine swap_integers(a, b)
    !> Swap two integers
    integer, intent(inout) :: a, b
    integer :: temp
    
    temp = a
    a = b
    b = temp
end subroutine swap_integers

end module bowyer_watson