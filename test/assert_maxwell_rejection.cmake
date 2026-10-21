execute_process(COMMAND "${PROGRAM}" "${PROBE}"
                RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0)
  message(FATAL_ERROR "${PROBE}: invalid exterior was accepted")
endif()
string(CONCAT combined "${output}" "${error}")
if(NOT combined MATCHES "${EXPECTED}")
  message(FATAL_ERROR "${PROBE}: rejection lacked required reason: ${combined}")
endif()
