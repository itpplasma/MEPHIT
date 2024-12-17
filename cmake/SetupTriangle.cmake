if (DEFINED ENV{TRIANGLE_DIR} OR DEFINED TRIANGLE_DIR)
  add_custom_target(TRIANGLE ALL)
else()
  message(STATUS "Downloading and building external libraries")
endif()

set(TRIANGLE_INCLUDE_DIR ${TRIANGLE_DIR} CACHE STRING "Path to triangle include")
set(TRIANGLE_LIB_DIR ${TRIANGLE_DIR} CACHE STRING "Path to triangle lib")

set(triangle_lib ${TRIANGLE_LIB_DIR}/libtriangle${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE STRING "TRIANGLE lib")
