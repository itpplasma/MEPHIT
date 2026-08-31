execute_process(
  COMMAND
    "${C_COMPILER}"
    -std=c11
    -D_GNU_SOURCE
    -DREAL=double
    -DUSE_FORTFEM
    "-I${MEPHIT_INCLUDE}"
    "-I${FORTFEM_INCLUDE}"
    -c
    "${MEPHIT_FEM_SOURCE}"
    -o
    "${TEST_OBJECT}"
  RESULT_VARIABLE compile_status
  OUTPUT_VARIABLE compile_output
  ERROR_VARIABLE compile_error)

if(NOT compile_status EQUAL 0)
  message(FATAL_ERROR
    "FortFEM FEM dispatch must compile without Triangle headers:\n"
    "${compile_output}${compile_error}")
endif()
