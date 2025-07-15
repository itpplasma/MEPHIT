! Patch for mephit_mesh.F90 to use configurable triangulation
! Add this to the module imports section:

use mephit_triangulation_wrapper, only: FEM_triangulate_external_wrapper, &
    set_triangulation_method, TRIANGULATION_TRIANGLE, TRIANGULATION_FORTRAN

! Replace the existing FEM_triangulate_external interface (around line 272) with:

! Configuration subroutine to set triangulation method
subroutine set_mesh_triangulation_method(method)
    integer, intent(in) :: method
    call set_triangulation_method(method)
end subroutine set_mesh_triangulation_method

! Replace the call in write_FreeFem_mesh (around line 3080) from:
!   call FEM_triangulate_external(npt_inner, npt_outer, bdry_R, bdry_Z, R_mid, Z_mid, &
!     decorate_filename('outer.msh', '', basename_suffix) // c_null_char)

! To:
!   call FEM_triangulate_external_wrapper(npt_inner, npt_outer, bdry_R, bdry_Z, R_mid, Z_mid, &
!     decorate_filename('outer.msh', '', basename_suffix))

! Also add the configuration subroutine to the public interface:
! public :: set_mesh_triangulation_method