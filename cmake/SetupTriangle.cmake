if(NOT "${TRIANGLE_DIR}" STREQUAL "")
  add_custom_target(TRIANGLE ALL)
elseif(DEFINED ENV{TRIANGLE_DIR})
  add_custom_target(TRIANGLE ALL)
  set(TRIANGLE_DIR $ENV{TRIANGLE_DIR} CACHE STRING "Triangle directory")
else()# Define Triangle repository
  set(TRIANGLE_REPO "https://salsa.debian.org/science-team/triangle.git")
  set(TRIANGLE_PATCHES_DIR "debian/patches")
  set(TRIANGLE_PATCH_SERIES "${TRIANGLE_PATCHES_DIR}/series")

  # External Project to fetch, patch, and build Triangle
  ExternalProject_Add(TRIANGLE
      PREFIX ${CMAKE_BINARY_DIR}/triangle
      GIT_REPOSITORY ${TRIANGLE_REPO}
      GIT_TAG master
      CONFIGURE_COMMAND "" # Skip configure step
      BUILD_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/setup_triangle.sh <SOURCE_DIR>
      INSTALL_COMMAND ""
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
      BUILD_BYPRODUCTS <SOURCE_DIR>/libtriangle.so
  )
endif()

set(TRIANGLE_INCLUDE_DIR ${TRIANGLE_DIR} CACHE STRING "Path to triangle include")
set(TRIANGLE_LIB_DIR ${TRIANGLE_DIR} CACHE STRING "Path to triangle lib")

set(triangle_lib ${TRIANGLE_LIB_DIR}/libtriangle${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE STRING "TRIANGLE lib")
