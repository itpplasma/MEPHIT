program test_fortran_interface_simple
    use iso_fortran_env, only: dp => real64
    use mephit_triangulation_fortran
    implicit none
    
    write(*,*) '=== Testing Direct Fortran Triangulation Interface ==='
    
    call test_direct_fortran_interface()
    call test_complex_geometry()
    
    write(*,*) 'Direct Fortran interface tests completed successfully!'
    
contains

subroutine test_direct_fortran_interface()
    logical :: file_exists
    integer, parameter :: npt_inner = 0, npt_outer = 4
    real(dp), parameter :: node_R(4) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: node_Z(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: R_O = 0.5_dp, Z_O = 0.5_dp
    character(len=30) :: filename = 'test_direct_fortran.msh'
    
    write(*,*) 'Test 1: Direct Fortran Interface (Square)'
    
    call FEM_triangulate_external_fortran(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, filename)
    
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
        call verify_mesh_file(filename)
        write(*,*) '  ✓ Direct Fortran interface test passed'
    else
        write(*,*) '  ✗ Direct Fortran interface test failed'
    end if
    
end subroutine test_direct_fortran_interface

subroutine test_complex_geometry()
    logical :: file_exists
    integer, parameter :: npt_inner = 1, npt_outer = 6
    real(dp), parameter :: node_R(7) = [0.0_dp, 1.0_dp, 0.5_dp, -0.5_dp, -1.0_dp, -0.5_dp, 0.5_dp]
    real(dp), parameter :: node_Z(7) = [0.0_dp, 0.0_dp, 0.866_dp, 0.866_dp, 0.0_dp, -0.866_dp, -0.866_dp]
    real(dp), parameter :: R_O = 0.0_dp, Z_O = 0.0_dp
    character(len=30) :: filename = 'test_complex_geometry.msh'
    
    write(*,*) 'Test 2: Complex Geometry (Hexagon with inner point)'
    
    call FEM_triangulate_external_fortran(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, filename)
    
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
        call verify_mesh_file(filename)
        write(*,*) '  ✓ Complex geometry test passed'
    else
        write(*,*) '  ✗ Complex geometry test failed'
    end if
    
end subroutine test_complex_geometry

subroutine verify_mesh_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, npoints, ntriangles, dummy
    real(dp) :: x, y
    integer :: v1, v2, v3
    integer :: i
    
    open(newunit=unit, file=filename, status='old', action='read')
    
    read(unit, *) npoints, ntriangles, dummy
    
    write(*,'(A,I0,A,I0,A)') '    Mesh: ', npoints, ' points, ', ntriangles, ' triangles'
    
    do i = 1, npoints
        read(unit, *) x, y, dummy
    end do
    
    do i = 1, ntriangles
        read(unit, *) v1, v2, v3, dummy
        if (v1 < 1 .or. v1 > npoints .or. v2 < 1 .or. v2 > npoints .or. v3 < 1 .or. v3 > npoints) then
            write(*,'(A,I0,A,3I0)') '    Warning: Invalid triangle ', i, ' vertices: ', v1, v2, v3
        end if
    end do
    
    close(unit)
    
end subroutine verify_mesh_file

end program test_fortran_interface_simple