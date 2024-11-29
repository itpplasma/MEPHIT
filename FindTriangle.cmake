find_path(TRIANGLE_INCLUDE_DIR
  NAMES triangle.h
  PATHS ENV TRIANGLE_DIR
)

if (TRIANGLE_INCLUDE_DIR)
  message(STATUS "Found triangle.h in ${TRIANGLE_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "triangle.h not found! Set the TRIANGLE_DIR environment variable.")
endif()

set(TRIANGLE_INCLUDE_DIR ${TRIANGLE_INCLUDE_DIR} CACHE INTERNAL "Path to triangle include")
