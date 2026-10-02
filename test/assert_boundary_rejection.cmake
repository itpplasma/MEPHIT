execute_process(COMMAND "${PROGRAM}" escaped
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors
  TIMEOUT 60)
if("${status}" STREQUAL "0")
  message(FATAL_ERROR "Contour crossing the known rectangle was accepted")
endif()
if(NOT "${output}${errors}" MATCHES
    "Closed-contour ODE evaluation leaves EQDSK rectangle")
  message(FATAL_ERROR "Expected domain rejection, got: ${status}\n${output}${errors}")
endif()
message(STATUS "Known circular excursion rejected by surface-domain check")
