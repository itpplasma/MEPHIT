module mephit_triangulation_wrapper
    !> Configurable triangulation wrapper for MEPHIT
    !> Allows switching between TRIANGLE library and native Fortran implementation
    
    use iso_fortran_env, only: dp => real64
    use iso_c_binding, only: c_char, c_int, c_double, c_null_char
    implicit none
    
    private
    public :: FEM_triangulate_external_wrapper
    public :: set_triangulation_method
    public :: TRIANGULATION_TRIANGLE, TRIANGULATION_FORTRAN
    
    ! Triangulation method constants
    integer, parameter :: TRIANGULATION_TRIANGLE = 1
    integer, parameter :: TRIANGULATION_FORTRAN = 2
    
    ! Current triangulation method (default to Fortran)
    integer :: current_method = TRIANGULATION_FORTRAN
    
    ! C interface to original TRIANGLE library
    interface
        subroutine FEM_triangulate_external_c(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname) &
            bind(C, name='FEM_triangulate_external')
            use iso_c_binding, only: c_char, c_int, c_double
            integer(c_int), intent(in), value :: npt_inner, npt_outer
            real(c_double), intent(in), dimension(*) :: node_R, node_Z
            real(c_double), intent(in), value :: R_O, Z_O
            character(c_char), intent(in) :: fname(*)
        end subroutine FEM_triangulate_external_c
    end interface
    
contains

subroutine set_triangulation_method(method)
    !> Set the triangulation method
    integer, intent(in) :: method
    
    if (method == TRIANGULATION_TRIANGLE .or. method == TRIANGULATION_FORTRAN) then
        current_method = method
        select case (method)
        case (TRIANGULATION_TRIANGLE)
            write(*,*) 'Triangulation method set to: TRIANGLE library'
        case (TRIANGULATION_FORTRAN)
            write(*,*) 'Triangulation method set to: Native Fortran'
        end select
    else
        write(*,*) 'Warning: Invalid triangulation method, keeping current method'
    end if
end subroutine set_triangulation_method

subroutine FEM_triangulate_external_wrapper(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    !> Wrapper that calls either TRIANGLE or Fortran implementation
    integer, intent(in) :: npt_inner, npt_outer
    real(dp), intent(in) :: node_R(:), node_Z(:)
    real(dp), intent(in) :: R_O, Z_O
    character(len=*), intent(in) :: fname
    
    select case (current_method)
    case (TRIANGULATION_TRIANGLE)
        call use_triangle_library(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    case (TRIANGULATION_FORTRAN)
        call use_fortran_implementation(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    case default
        write(*,*) 'Error: Unknown triangulation method'
        stop 1
    end select
end subroutine FEM_triangulate_external_wrapper

subroutine use_triangle_library(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    !> Call original TRIANGLE library through C interface
    integer, intent(in) :: npt_inner, npt_outer
    real(dp), intent(in) :: node_R(:), node_Z(:)
    real(dp), intent(in) :: R_O, Z_O
    character(len=*), intent(in) :: fname
    
    character(len=len(fname)+1) :: c_fname
    
    ! Convert Fortran string to C string
    c_fname = fname // c_null_char
    
    ! Call C interface
    call FEM_triangulate_external_c(int(npt_inner, c_int), int(npt_outer, c_int), &
        real(node_R, c_double), real(node_Z, c_double), &
        real(R_O, c_double), real(Z_O, c_double), c_fname)
    
end subroutine use_triangle_library

subroutine use_fortran_implementation(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    !> Call native Fortran implementation
    use mephit_triangulation_fortran, only: FEM_triangulate_external_fortran
    
    integer, intent(in) :: npt_inner, npt_outer
    real(dp), intent(in) :: node_R(:), node_Z(:)
    real(dp), intent(in) :: R_O, Z_O
    character(len=*), intent(in) :: fname
    
    call FEM_triangulate_external_fortran(npt_inner, npt_outer, node_R, node_Z, R_O, Z_O, fname)
    
end subroutine use_fortran_implementation

end module mephit_triangulation_wrapper