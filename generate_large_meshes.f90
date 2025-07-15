program generate_large_meshes
    use iso_fortran_env, only: dp => real64
    use delaunay_types
    use triangulation_fortran
    implicit none
    
    write(*,*) '=== Generating Large Mesh Files for Plotting ==='
    
    call generate_random_10_points()
    call generate_random_25_points()
    call generate_random_50_points()
    call generate_random_100_points()
    
    write(*,*) 'All large mesh files generated successfully!'
    
contains

subroutine generate_random_10_points()
    real(dp), parameter :: points(2,10) = reshape([&
        0.1_dp, 0.3_dp, 0.7_dp, 0.9_dp, 0.2_dp, 0.8_dp, 0.4_dp, 0.6_dp, 0.5_dp, 0.35_dp, &
        0.2_dp, 0.8_dp, 0.1_dp, 0.5_dp, 0.9_dp, 0.3_dp, 0.7_dp, 0.4_dp, 0.6_dp, 0.75_dp], [2, 10])
    
    type(triangulation_result_t) :: result
    
    write(*,*) 'Generating random 10 points mesh...'
    
    call triangulate_fortran(points, reshape([0], [2, 0]), result)
    call write_mesh_file(result, 'large_random_10.msh')
    
    write(*,'(A,I0,A,I0,A)') '  → large_random_10.msh: ', result%npoints, ' points, ', result%ntriangles, ' triangles'
    
    call cleanup_triangulation(result)
end subroutine

subroutine generate_random_25_points()
    real(dp) :: points(2,25)
    integer :: i
    type(triangulation_result_t) :: result
    
    ! Generate pseudo-random points
    do i = 1, 25
        points(1, i) = mod(i * 17 + 7, 100) / 100.0_dp
        points(2, i) = mod(i * 23 + 13, 100) / 100.0_dp
    end do
    
    write(*,*) 'Generating random 25 points mesh...'
    
    call triangulate_fortran(points, reshape([0], [2, 0]), result)
    call write_mesh_file(result, 'large_random_25.msh')
    
    write(*,'(A,I0,A,I0,A)') '  → large_random_25.msh: ', result%npoints, ' points, ', result%ntriangles, ' triangles'
    
    call cleanup_triangulation(result)
end subroutine

subroutine generate_random_50_points()
    real(dp) :: points(2,50)
    integer :: i
    type(triangulation_result_t) :: result
    
    ! Generate pseudo-random points
    do i = 1, 50
        points(1, i) = mod(i * 31 + 11, 200) / 200.0_dp
        points(2, i) = mod(i * 37 + 19, 200) / 200.0_dp
    end do
    
    write(*,*) 'Generating random 50 points mesh...'
    
    call triangulate_fortran(points, reshape([0], [2, 0]), result)
    call write_mesh_file(result, 'large_random_50.msh')
    
    write(*,'(A,I0,A,I0,A)') '  → large_random_50.msh: ', result%npoints, ' points, ', result%ntriangles, ' triangles'
    
    call cleanup_triangulation(result)
end subroutine

subroutine generate_random_100_points()
    real(dp) :: points(2,100)
    integer :: i
    type(triangulation_result_t) :: result
    
    ! Generate pseudo-random points
    do i = 1, 100
        points(1, i) = mod(i * 41 + 17, 300) / 300.0_dp
        points(2, i) = mod(i * 47 + 23, 300) / 300.0_dp
    end do
    
    write(*,*) 'Generating random 100 points mesh...'
    
    call triangulate_fortran(points, reshape([0], [2, 0]), result)
    call write_mesh_file(result, 'large_random_100.msh')
    
    write(*,'(A,I0,A,I0,A)') '  → large_random_100.msh: ', result%npoints, ' points, ', result%ntriangles, ' triangles'
    
    call cleanup_triangulation(result)
end subroutine

subroutine write_mesh_file(result, filename)
    type(triangulation_result_t), intent(in) :: result
    character(len=*), intent(in) :: filename
    
    integer :: unit, i
    
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
end subroutine

end program generate_large_meshes