set(backend_was_launched FALSE)
foreach(attempt RANGE 1 5)
  execute_process(
    COMMAND
      "${MEPHIT_RUN}" 0
      "${TEST_TEMP}/fortfem-intentionally-missing.in"
      "probe"
      "${TEST_TEMP}"
      "/usr/bin/printf"
    OUTPUT_VARIABLE run_output
    ERROR_VARIABLE run_error
    TIMEOUT 10)

  string(FIND "${run_output}${run_error}" "MEPHIT_0x" backend_marker)
  if(NOT backend_marker EQUAL -1)
    set(backend_was_launched TRUE)
    break()
  endif()
endforeach()

if(backend_was_launched)
  message(FATAL_ERROR "FortFEM run launched the external FEM backend")
endif()

execute_process(
  COMMAND
    "${MEPHIT_RUN}" 0
    "${TEST_TEMP}/fortfem-intentionally-missing.in"
    "probe"
    "${TEST_TEMP}"
  OUTPUT_VARIABLE four_argument_output
  ERROR_VARIABLE four_argument_error
  TIMEOUT 10)

string(FIND
  "${four_argument_output}${four_argument_error}"
  "FreeFem script file"
  obsolete_argument_error)
if(NOT obsolete_argument_error EQUAL -1)
  message(FATAL_ERROR "FortFEM run still requires a FreeFem script argument")
endif()
