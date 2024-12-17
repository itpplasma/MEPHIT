if(NOT "${MFEM_DIR}" STREQUAL "")
  add_custom_target(MFEM ALL)
elseif(DEFINED ENV{MFEM_DIR})
  add_custom_target(MFEM ALL)
  set(MFEM_DIR $ENV{MFEM_DIR})
else()
  message(STATUS "Downloading and building MFEM")

  # Define version and URL for MFEM
  set(MFEM_VERSION 4.7)
  set(MFEM_URL "https://github.com/mfem/mfem/archive/refs/tags/v${MFEM_VERSION}.tar.gz")

  cmake_host_system_information(RESULT NUM_CORES QUERY NUMBER_OF_PHYSICAL_CORES)

  # Specify external project
  include(ExternalProject)
  ExternalProject_Add(
    MFEM
    PREFIX ${CMAKE_BINARY_DIR}/mfem
    URL ${MFEM_URL}
    CMAKE_ARGS
      -DMFEM_USE_SUITESPARSE=1
      -DCMAKE_CXX_FLAGS=-fPIC
      -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/mfem/install
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_BYPRODUCTS
      ${CMAKE_BINARY_DIR}/mfem/src/MFEM-build/libmfem${CMAKE_STATIC_LIBRARY_SUFFIX}
  )
  set(MFEM_DIR ${CMAKE_BINARY_DIR}/mfem/src/MFEM-build)
endif()

set(MFEM_INCLUDE_DIRS ${MFEM_DIR} CACHE STRING "MFEM include")
set(MFEM_LIBRARIES ${MFEM_DIR}/libmfem${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE STRING "MFEM lib")
