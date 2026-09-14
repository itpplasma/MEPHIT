program test_mesh

  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr, c_null_char

  implicit none (external)

  interface
    function MFEM_init(tor_mode, mesh_file, edgemap_file) result(maxwell_solver) &
      bind(C, name = 'MFEM_init')
      use iso_c_binding, only: c_char, c_int, c_ptr
      integer(c_int), intent(in), value :: tor_mode
      character(c_char), intent(in) :: mesh_file(*)
      character(c_char), intent(in) :: edgemap_file(*)
      type(c_ptr) :: maxwell_solver
    end function MFEM_init

    function test_map_edges(maxwell_solver, test_edgemap_file) result(status) &
      bind(C, name = 'test_map_edges')
      import, only: c_char, c_int, c_ptr
      type(c_ptr), intent(in), value :: maxwell_solver
      character(c_char), intent(in) :: test_edgemap_file(*)
      integer(c_int) :: status
    end function test_map_edges

    subroutine MFEM_deinit(maxwell_solver) bind(C, name = 'MFEM_deinit')
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(in), value :: maxwell_solver
    end subroutine MFEM_deinit
  end interface

  type(c_ptr) :: maxwell_solver
  character(len = *), parameter :: mesh_file = 'core_plasma.mesh', &
    edgemap_file = 'edgemap.dat', test_edgemap_file = 'test_edgemap.dat'
  integer, parameter :: n_tor = 2
  integer status

  maxwell_solver = MFEM_init(n_tor, mesh_file // c_null_char, edgemap_file // c_null_char)
  status = test_map_edges(maxwell_solver, test_edgemap_file // c_null_char)
  call MFEM_deinit(maxwell_solver)
  if (status /= 0) error stop status

end program test_mesh
