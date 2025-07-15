program test_triangulation
    use iso_fortran_env, only: dp => real64
    use delaunay_types
    use geometric_predicates
    use triangulation_fortran
    implicit none
    
    ! Test framework
    integer :: test_count = 0
    integer :: passed_tests = 0
    
    write(*,*) '=== MEPHIT Triangulation Tests ==='
    
    ! Test suite
    call test_data_structures()
    call test_geometric_predicates()
    call test_simple_triangle()
    call test_bowyer_watson()
    call test_square_with_hole()
    call test_complex_boundary()
    call test_quality_constraints()
    call test_edge_cases()
    
    ! Summary
    write(*,*) 
    write(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, '/', test_count
    if (passed_tests == test_count) then
        write(*,*) 'All tests passed!'
    else
        write(*,*) 'Some tests failed!'
        stop 1
    end if
    
contains

subroutine test_data_structures()
    character(len=*), parameter :: test_name = 'Data Structures'
    type(mesh_t) :: mesh
    integer :: p1, p2, p3, t1, e1
    
    call start_test(test_name)
    
    ! Test mesh creation
    call create_mesh(mesh, 100, 200, 150)
    call assert_equal(mesh%max_points, 100, 'Max points set correctly')
    call assert_equal(mesh%max_triangles, 200, 'Max triangles set correctly')
    call assert_equal(mesh%max_edges, 150, 'Max edges set correctly')
    call assert_equal(mesh%npoints, 0, 'Initial point count is zero')
    
    ! Test point addition
    p1 = add_point(mesh, 0.0_dp, 0.0_dp, 1)
    p2 = add_point(mesh, 1.0_dp, 0.0_dp, 2)
    p3 = add_point(mesh, 0.5_dp, 0.866_dp, 3)
    
    call assert_equal(mesh%npoints, 3, 'Three points added')
    call assert_equal(p1, 1, 'First point index')
    call assert_equal(p2, 2, 'Second point index')
    call assert_equal(p3, 3, 'Third point index')
    
    ! Test triangle addition
    t1 = add_triangle(mesh, p1, p2, p3)
    call assert_equal(mesh%ntriangles, 1, 'One triangle added')
    call assert_equal(t1, 1, 'Triangle index')
    call assert_true(is_valid_triangle(mesh, t1), 'Triangle is valid')
    
    ! Test edge addition
    e1 = add_edge(mesh, p1, p2, .true.)
    call assert_equal(mesh%nedges, 1, 'One edge added')
    call assert_equal(e1, 1, 'Edge index')
    call assert_true(is_valid_edge(mesh, e1), 'Edge is valid')
    call assert_true(mesh%edges(e1)%constrained, 'Edge is constrained')
    
    ! Test mesh resizing
    call resize_mesh(mesh, 200, 400, 300)
    call assert_equal(mesh%max_points, 200, 'Points resized')
    call assert_equal(mesh%max_triangles, 400, 'Triangles resized')
    call assert_equal(mesh%max_edges, 300, 'Edges resized')
    call assert_equal(mesh%npoints, 3, 'Point count preserved')
    
    call destroy_mesh(mesh)
    call end_test()
end subroutine

subroutine test_geometric_predicates()
    character(len=*), parameter :: test_name = 'Geometric Predicates'
    type(point_t) :: pa, pb, pc, pd
    type(mesh_t) :: mesh
    integer :: t1, orient_result, p1, p2, p3
    real(dp) :: area
    
    call start_test(test_name)
    
    ! Set up test points
    pa = point_t(0.0_dp, 0.0_dp, 1, .true.)
    pb = point_t(1.0_dp, 0.0_dp, 2, .true.)
    pc = point_t(0.5_dp, 0.866_dp, 3, .true.)
    pd = point_t(0.5_dp, 0.3_dp, 4, .true.)
    
    ! Test orientation
    orient_result = orientation(pa, pb, pc)
    call assert_equal(orient_result, ORIENTATION_CCW, 'CCW orientation')
    
    orient_result = orientation(pa, pc, pb)
    call assert_equal(orient_result, ORIENTATION_CW, 'CW orientation')
    
    orient_result = orientation(pa, pb, point_t(0.5_dp, 0.0_dp, 5, .true.))
    call assert_equal(orient_result, ORIENTATION_COLLINEAR, 'Collinear points')
    
    ! Test in_circle
    call assert_true(in_circle(pa, pb, pc, pd), 'Point inside circumcircle')
    
    pd = point_t(2.0_dp, 2.0_dp, 4, .true.)
    call assert_false(in_circle(pa, pb, pc, pd), 'Point outside circumcircle')
    
    ! Test triangle area
    area = triangle_area(pa, pb, pc)
    call assert_approx_equal(area, 0.433_dp, 0.001_dp, 'Triangle area')
    
    ! Test point in triangle
    call create_mesh(mesh, 10, 10, 10)
    p1 = add_point(mesh, 0.0_dp, 0.0_dp, 1)
    p2 = add_point(mesh, 1.0_dp, 0.0_dp, 2)
    p3 = add_point(mesh, 0.5_dp, 0.866_dp, 3)
    t1 = add_triangle(mesh, p1, p2, p3)
    
    pd = point_t(0.4_dp, 0.2_dp, 4, .true.)  ! Point clearly inside
    call assert_true(point_in_triangle(pd, mesh, t1), 'Point inside triangle')
    
    pd = point_t(0.0_dp, 1.0_dp, 4, .true.)  ! Point clearly outside
    call assert_false(point_in_triangle(pd, mesh, t1), 'Point outside triangle')
    
    call destroy_mesh(mesh)
    call end_test()
end subroutine

subroutine test_simple_triangle()
    character(len=*), parameter :: test_name = 'Simple Triangle'
    type(mesh_t) :: mesh
    type(point_t) :: pa, pb, pc
    integer :: p1, p2, p3, t1
    real(dp) :: area
    
    call start_test(test_name)
    
    ! Create a simple triangle mesh
    call create_mesh(mesh, 10, 10, 10)
    
    ! Add triangle vertices
    p1 = add_point(mesh, 0.0_dp, 0.0_dp, 1)
    p2 = add_point(mesh, 1.0_dp, 0.0_dp, 2)
    p3 = add_point(mesh, 0.5_dp, 0.866_dp, 3)
    
    ! Add triangle
    t1 = add_triangle(mesh, p1, p2, p3)
    
    ! Test triangle properties
    call assert_equal(mesh%npoints, 3, 'Number of points')
    call assert_equal(mesh%ntriangles, 1, 'Number of triangles')
    call assert_true(is_valid_triangle(mesh, t1), 'Triangle is valid')
    
    ! Test triangle area
    pa = mesh%points(p1)
    pb = mesh%points(p2)
    pc = mesh%points(p3)
    area = triangle_area(pa, pb, pc)
    call assert_approx_equal(area, 0.433_dp, 0.001_dp, 'Triangle area')
    
    ! Test orientation (should be CCW)
    call assert_equal(orientation(pa, pb, pc), ORIENTATION_CCW, 'CCW orientation')
    
    call destroy_mesh(mesh)
    call end_test()
end subroutine

subroutine test_bowyer_watson()
    character(len=*), parameter :: test_name = 'Bowyer-Watson Algorithm'
    real(dp), parameter :: points(2,4) = reshape([&
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp], [2, 4])
    integer, parameter :: segments(2,4) = reshape([&
        1, 2, &
        2, 3, &
        3, 4, &
        4, 1], [2, 4])
    
    type(triangulation_result_t) :: result
    
    call start_test(test_name)
    
    ! Test Bowyer-Watson triangulation
    call triangulate_fortran(points, segments, result)
    
    ! Should have 4 points
    call assert_equal(result%npoints, 4, 'Number of points')
    
    ! Should have 2 triangles (for a square)
    call assert_equal(result%ntriangles, 2, 'Number of triangles')
    
    ! Should have 4 segments
    call assert_equal(result%nsegments, 4, 'Number of segments')
    
    ! Check that all triangles are valid (positive area)
    call assert_true(all_triangles_valid(result), 'All triangles have positive area')
    
    call cleanup_triangulation(result)
    call end_test()
end subroutine

logical function all_triangles_valid(result)
    type(triangulation_result_t), intent(in) :: result
    integer :: i
    real(dp) :: area
    
    all_triangles_valid = .true.
    do i = 1, result%ntriangles
        area = compute_triangle_area(result%points, result%triangles(:, i))
        if (area <= 0.0_dp) then
            all_triangles_valid = .false.
            return
        end if
    end do
end function all_triangles_valid

real(dp) function compute_triangle_area(points, triangle)
    real(dp), intent(in) :: points(:,:)
    integer, intent(in) :: triangle(3)
    real(dp) :: x1, y1, x2, y2, x3, y3
    
    x1 = points(1, triangle(1))
    y1 = points(2, triangle(1))
    x2 = points(1, triangle(2))
    y2 = points(2, triangle(2))
    x3 = points(1, triangle(3))
    y3 = points(2, triangle(3))
    
    compute_triangle_area = 0.5_dp * abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)))
end function compute_triangle_area

subroutine test_square_with_hole()
    character(len=*), parameter :: test_name = 'Square with Hole (placeholder)'
    
    call start_test(test_name)
    
    ! Placeholder test - will be implemented with Bowyer-Watson algorithm
    write(*,*) '  TODO: Implement after Bowyer-Watson algorithm'
    
    call end_test()
end subroutine

subroutine test_complex_boundary()
    character(len=*), parameter :: test_name = 'Complex Boundary (placeholder)'
    
    call start_test(test_name)
    
    ! Placeholder test - will be implemented with Bowyer-Watson algorithm
    write(*,*) '  TODO: Implement after Bowyer-Watson algorithm'
    
    call end_test()
end subroutine

subroutine test_quality_constraints()
    character(len=*), parameter :: test_name = 'Quality Constraints (placeholder)'
    
    call start_test(test_name)
    
    ! Placeholder test - will be implemented with quality improvement
    write(*,*) '  TODO: Implement after quality improvement module'
    
    call end_test()
end subroutine

subroutine test_edge_cases()
    character(len=*), parameter :: test_name = 'Edge Cases'
    type(point_t) :: pa, pb, pc
    
    call start_test(test_name)
    
    ! Test collinear points
    pa = point_t(0.0_dp, 0.0_dp, 1, .true.)
    pb = point_t(1.0_dp, 0.0_dp, 2, .true.)
    pc = point_t(2.0_dp, 0.0_dp, 3, .true.)
    
    call assert_equal(orientation(pa, pb, pc), ORIENTATION_COLLINEAR, 'Collinear detection')
    
    ! Test very small triangle
    pa = point_t(0.0_dp, 0.0_dp, 1, .true.)
    pb = point_t(1e-15_dp, 0.0_dp, 2, .true.)
    pc = point_t(0.0_dp, 1e-15_dp, 3, .true.)
    
    call assert_true(triangle_area(pa, pb, pc) < 1e-20_dp, 'Very small triangle area')
    
    call end_test()
end subroutine


! Helper subroutines and functions

subroutine start_test(test_name)
    character(len=*), intent(in) :: test_name
    test_count = test_count + 1
    write(*,'(A,I0,A,A)') 'Test ', test_count, ': ', test_name
end subroutine

subroutine end_test()
    passed_tests = passed_tests + 1
    write(*,*) '  PASSED'
end subroutine

subroutine assert_equal(actual, expected, description)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: description
    if (actual /= expected) then
        write(*,'(A,A,A,I0,A,I0)') '  FAILED: ', description, ' - got ', actual, ', expected ', expected
        stop 1
    end if
end subroutine

subroutine assert_true(condition, description)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: description
    if (.not. condition) then
        write(*,'(A,A)') '  FAILED: ', description
        stop 1
    end if
end subroutine

subroutine assert_false(condition, description)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: description
    if (condition) then
        write(*,'(A,A)') '  FAILED: ', description
        stop 1
    end if
end subroutine

subroutine assert_approx_equal(actual, expected, tolerance, description)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: description
    if (abs(actual - expected) > tolerance) then
        write(*,'(A,A,A,ES15.8,A,ES15.8)') '  FAILED: ', description, ' - got ', actual, ', expected ', expected
        stop 1
    end if
end subroutine

end program test_triangulation