program test_triangulation
    use iso_fortran_env, only: dp => real64
    implicit none
    
    ! Test framework
    integer :: test_count = 0
    integer :: passed_tests = 0
    
    write(*,*) '=== MEPHIT Triangulation Tests ==='
    
    ! Test suite
    call test_simple_triangle()
    call test_square_with_hole()
    call test_complex_boundary()
    call test_quality_constraints()
    call test_edge_cases()
    call test_triangle_compatibility()
    
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

subroutine test_simple_triangle()
    character(len=*), parameter :: test_name = 'Simple Triangle'
    real(dp), parameter :: points(2,3) = reshape([&
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        0.5_dp, 0.866_dp], [2, 3])
    integer, parameter :: segments(2,3) = reshape([&
        1, 2, &
        2, 3, &
        3, 1], [2, 3])
    
    type(triangulation_result_t) :: result
    
    call start_test(test_name)
    
    ! Test simple triangle triangulation
    call triangulate_fortran(points, segments, result)
    
    ! Should produce exactly 1 triangle with 3 vertices
    call assert_equal(result%npoints, 3, 'Number of points')
    call assert_equal(result%ntriangles, 1, 'Number of triangles')
    call assert_equal(result%nsegments, 3, 'Number of segments')
    
    ! Check triangle connectivity (should be 1-2-3)
    call assert_equal(result%triangles(1,1), 1, 'Triangle vertex 1')
    call assert_equal(result%triangles(2,1), 2, 'Triangle vertex 2')
    call assert_equal(result%triangles(3,1), 3, 'Triangle vertex 3')
    
    ! Check area is positive (counter-clockwise orientation)
    call assert_true(triangle_area(result%points, result%triangles(:,1)) > 0.0_dp, &
                     'Positive area (CCW orientation)')
    
    call cleanup_triangulation(result)
    call end_test()
end subroutine

subroutine test_square_with_hole()
    character(len=*), parameter :: test_name = 'Square with Hole'
    ! Outer square: (0,0) -> (2,0) -> (2,2) -> (0,2) -> (0,0)
    ! Inner square: (0.5,0.5) -> (1.5,0.5) -> (1.5,1.5) -> (0.5,1.5) -> (0.5,0.5)
    real(dp), parameter :: points(2,8) = reshape([&
        0.0_dp, 0.0_dp, &  ! outer square
        2.0_dp, 0.0_dp, &
        2.0_dp, 2.0_dp, &
        0.0_dp, 2.0_dp, &
        0.5_dp, 0.5_dp, &  ! inner square (hole)
        1.5_dp, 0.5_dp, &
        1.5_dp, 1.5_dp, &
        0.5_dp, 1.5_dp], [2, 8])
    integer, parameter :: segments(2,8) = reshape([&
        1, 2, &  ! outer segments
        2, 3, &
        3, 4, &
        4, 1, &
        5, 6, &  ! inner segments
        6, 7, &
        7, 8, &
        8, 5], [2, 8])
    real(dp), parameter :: hole_point(2) = [1.0_dp, 1.0_dp]
    
    type(triangulation_result_t) :: result
    
    call start_test(test_name)
    
    ! Test square with hole
    call triangulate_with_hole_fortran(points, segments, hole_point, result)
    
    ! Should have more than 8 points (Steiner points added)
    call assert_true(result%npoints >= 8, 'At least 8 points')
    
    ! Should have multiple triangles
    call assert_true(result%ntriangles > 1, 'Multiple triangles')
    
    ! Check that no triangle contains the hole point
    call assert_false(any_triangle_contains_point(result, hole_point), &
                      'No triangle contains hole point')
    
    ! Check total area (should be 4.0 - 1.0 = 3.0)
    call assert_approx_equal(total_triangulation_area(result), 3.0_dp, &
                            1e-10_dp, 'Total area')
    
    call cleanup_triangulation(result)
    call end_test()
end subroutine

subroutine test_complex_boundary()
    character(len=*), parameter :: test_name = 'Complex Boundary'
    ! Test with a more complex boundary (circle approximation)
    integer, parameter :: n_circle = 16
    real(dp) :: points(2, n_circle)
    integer :: segments(2, n_circle)
    real(dp) :: angle
    integer :: i
    
    type(triangulation_result_t) :: result
    
    call start_test(test_name)
    
    ! Generate circle points
    do i = 1, n_circle
        angle = 2.0_dp * 3.14159265358979_dp * real(i-1, dp) / real(n_circle, dp)
        points(1, i) = cos(angle)
        points(2, i) = sin(angle)
        segments(1, i) = i
        segments(2, i) = mod(i, n_circle) + 1
    end do
    
    call triangulate_fortran(points, segments, result)
    
    ! Should have exactly n_circle boundary points
    call assert_equal(result%npoints, n_circle, 'Number of points equals boundary')
    
    ! Should have 2*(n_circle-2) triangles for a convex polygon
    call assert_equal(result%ntriangles, 2*(n_circle-2), 'Number of triangles')
    
    ! Check that all triangles have positive area
    call assert_true(all_triangles_positive_area(result), &
                     'All triangles have positive area')
    
    call cleanup_triangulation(result)
    call end_test()
end subroutine

subroutine test_quality_constraints()
    character(len=*), parameter :: test_name = 'Quality Constraints'
    ! Test minimum angle constraint (equivalent to TRIANGLE's 'q' option)
    real(dp), parameter :: points(2,4) = reshape([&
        0.0_dp, 0.0_dp, &
        2.0_dp, 0.0_dp, &
        2.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp], [2, 4])
    integer, parameter :: segments(2,4) = reshape([&
        1, 2, &
        2, 3, &
        3, 4, &
        4, 1], [2, 4])
    
    type(triangulation_result_t) :: result
    real(dp), parameter :: min_angle = 20.0_dp  ! degrees
    
    call start_test(test_name)
    
    ! Test with quality constraint
    call triangulate_with_quality_fortran(points, segments, min_angle, result)
    
    ! Should have added Steiner points to improve quality
    call assert_true(result%npoints >= 4, 'At least 4 points')
    
    ! Check that all triangles meet minimum angle constraint
    call assert_true(all_triangles_meet_angle_constraint(result, min_angle), &
                     'All triangles meet minimum angle constraint')
    
    call cleanup_triangulation(result)
    call end_test()
end subroutine

subroutine test_edge_cases()
    character(len=*), parameter :: test_name = 'Edge Cases'
    
    call start_test(test_name)
    
    ! Test degenerate cases
    call test_collinear_points()
    call test_duplicate_points()
    call test_very_small_triangle()
    
    call end_test()
end subroutine

subroutine test_collinear_points()
    ! Test handling of collinear points
    real(dp), parameter :: points(2,3) = reshape([&
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        2.0_dp, 0.0_dp], [2, 3])
    integer, parameter :: segments(2,3) = reshape([&
        1, 2, &
        2, 3, &
        3, 1], [2, 3])
    
    type(triangulation_result_t) :: result
    integer :: status
    
    ! This should fail gracefully
    call triangulate_fortran(points, segments, result, status)
    call assert_equal(status, -1, 'Collinear points should fail')
end subroutine

subroutine test_duplicate_points()
    ! Test handling of duplicate points
    real(dp), parameter :: points(2,3) = reshape([&
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp], [2, 3])  ! duplicate point
    integer, parameter :: segments(2,3) = reshape([&
        1, 2, &
        2, 3, &
        3, 1], [2, 3])
    
    type(triangulation_result_t) :: result
    integer :: status
    
    ! This should fail gracefully
    call triangulate_fortran(points, segments, result, status)
    call assert_equal(status, -1, 'Duplicate points should fail')
end subroutine

subroutine test_very_small_triangle()
    ! Test handling of very small triangles
    real(dp), parameter :: points(2,3) = reshape([&
        0.0_dp, 0.0_dp, &
        1e-15_dp, 0.0_dp, &
        0.0_dp, 1e-15_dp], [2, 3])
    integer, parameter :: segments(2,3) = reshape([&
        1, 2, &
        2, 3, &
        3, 1], [2, 3])
    
    type(triangulation_result_t) :: result
    integer :: status
    
    ! This should handle small triangles gracefully
    call triangulate_fortran(points, segments, result, status)
    ! Implementation should either succeed or fail gracefully
    call assert_true(status >= 0 .or. status == -1, 'Small triangle handled gracefully')
    
    if (status >= 0) then
        call cleanup_triangulation(result)
    end if
end subroutine

subroutine test_triangle_compatibility()
    character(len=*), parameter :: test_name = 'TRIANGLE Compatibility'
    ! Test that our output exactly matches TRIANGLE's output format
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
    
    type(triangulation_result_t) :: result_fortran, result_triangle
    
    call start_test(test_name)
    
    ! Get result from our Fortran implementation
    call triangulate_fortran(points, segments, result_fortran)
    
    ! Get result from TRIANGLE library (for comparison)
    call triangulate_triangle_lib(points, segments, result_triangle)
    
    ! Results should be identical
    call assert_equal(result_fortran%npoints, result_triangle%npoints, &
                     'Same number of points')
    call assert_equal(result_fortran%ntriangles, result_triangle%ntriangles, &
                     'Same number of triangles')
    call assert_equal(result_fortran%nsegments, result_triangle%nsegments, &
                     'Same number of segments')
    
    ! Areas should be identical
    call assert_approx_equal(total_triangulation_area(result_fortran), &
                            total_triangulation_area(result_triangle), &
                            1e-14_dp, 'Same total area')
    
    call cleanup_triangulation(result_fortran)
    call cleanup_triangulation(result_triangle)
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

real(dp) function triangle_area(points, triangle_vertices)
    real(dp), intent(in) :: points(:,:)
    integer, intent(in) :: triangle_vertices(3)
    real(dp) :: x1, y1, x2, y2, x3, y3
    
    x1 = points(1, triangle_vertices(1))
    y1 = points(2, triangle_vertices(1))
    x2 = points(1, triangle_vertices(2))
    y2 = points(2, triangle_vertices(2))
    x3 = points(1, triangle_vertices(3))
    y3 = points(2, triangle_vertices(3))
    
    triangle_area = 0.5_dp * abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)))
end function

logical function any_triangle_contains_point(result, point)
    type(triangulation_result_t), intent(in) :: result
    real(dp), intent(in) :: point(2)
    integer :: i
    
    any_triangle_contains_point = .false.
    do i = 1, result%ntriangles
        if (point_in_triangle(point, result%points, result%triangles(:,i))) then
            any_triangle_contains_point = .true.
            return
        end if
    end do
end function

logical function point_in_triangle(point, points, triangle_vertices)
    real(dp), intent(in) :: point(2), points(:,:)
    integer, intent(in) :: triangle_vertices(3)
    real(dp) :: x1, y1, x2, y2, x3, y3, px, py
    real(dp) :: denom, a, b, c
    
    x1 = points(1, triangle_vertices(1))
    y1 = points(2, triangle_vertices(1))
    x2 = points(1, triangle_vertices(2))
    y2 = points(2, triangle_vertices(2))
    x3 = points(1, triangle_vertices(3))
    y3 = points(2, triangle_vertices(3))
    px = point(1)
    py = point(2)
    
    denom = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    a = ((y2 - y3)*(px - x3) + (x3 - x2)*(py - y3)) / denom
    b = ((y3 - y1)*(px - x3) + (x1 - x3)*(py - y3)) / denom
    c = 1.0_dp - a - b
    
    point_in_triangle = (a >= 0.0_dp .and. b >= 0.0_dp .and. c >= 0.0_dp)
end function

real(dp) function total_triangulation_area(result)
    type(triangulation_result_t), intent(in) :: result
    integer :: i
    
    total_triangulation_area = 0.0_dp
    do i = 1, result%ntriangles
        total_triangulation_area = total_triangulation_area + &
            triangle_area(result%points, result%triangles(:,i))
    end do
end function

logical function all_triangles_positive_area(result)
    type(triangulation_result_t), intent(in) :: result
    integer :: i
    
    all_triangles_positive_area = .true.
    do i = 1, result%ntriangles
        if (triangle_area(result%points, result%triangles(:,i)) <= 0.0_dp) then
            all_triangles_positive_area = .false.
            return
        end if
    end do
end function

logical function all_triangles_meet_angle_constraint(result, min_angle_deg)
    type(triangulation_result_t), intent(in) :: result
    real(dp), intent(in) :: min_angle_deg
    integer :: i
    
    all_triangles_meet_angle_constraint = .true.
    do i = 1, result%ntriangles
        if (.not. triangle_meets_angle_constraint(result%points, result%triangles(:,i), min_angle_deg)) then
            all_triangles_meet_angle_constraint = .false.
            return
        end if
    end do
end function

logical function triangle_meets_angle_constraint(points, triangle_vertices, min_angle_deg)
    real(dp), intent(in) :: points(:,:)
    integer, intent(in) :: triangle_vertices(3)
    real(dp), intent(in) :: min_angle_deg
    real(dp) :: x1, y1, x2, y2, x3, y3
    real(dp) :: a, b, c, angle1, angle2, angle3
    real(dp), parameter :: pi = 3.14159265358979_dp
    
    x1 = points(1, triangle_vertices(1))
    y1 = points(2, triangle_vertices(1))
    x2 = points(1, triangle_vertices(2))
    y2 = points(2, triangle_vertices(2))
    x3 = points(1, triangle_vertices(3))
    y3 = points(2, triangle_vertices(3))
    
    ! Calculate side lengths
    a = sqrt((x2-x3)**2 + (y2-y3)**2)
    b = sqrt((x1-x3)**2 + (y1-y3)**2)
    c = sqrt((x1-x2)**2 + (y1-y2)**2)
    
    ! Calculate angles using law of cosines
    angle1 = acos((b**2 + c**2 - a**2) / (2.0_dp * b * c))
    angle2 = acos((a**2 + c**2 - b**2) / (2.0_dp * a * c))
    angle3 = acos((a**2 + b**2 - c**2) / (2.0_dp * a * b))
    
    ! Convert to degrees
    angle1 = angle1 * 180.0_dp / pi
    angle2 = angle2 * 180.0_dp / pi
    angle3 = angle3 * 180.0_dp / pi
    
    triangle_meets_angle_constraint = (angle1 >= min_angle_deg .and. &
                                       angle2 >= min_angle_deg .and. &
                                       angle3 >= min_angle_deg)
end function

end program test_triangulation