find_library(gsl_lib gsl ${GSL_LIB})

if (NOT DEFINED ENV{FGSL_DIR} AND NOT DEFINED FGSL_DIR)
  ExternalProject_Add(
    FGSL
    URL https://github.com/reinh-bader/fgsl/archive/refs/tags/v1.6.0.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/fgsl
    CONFIGURE_COMMAND mkdir m4 && autoreconf -i && ./configure
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/fgsl/src/FGSL/.libs/libfgsl${CMAKE_STATIC_LIBRARY_SUFFIX}
  )
  set(FGSL_DIR ${CMAKE_BINARY_DIR}/fgsl/src/FGSL CACHE STRING "FGSL directory")
else ()
  add_custom_target(FGSL ALL)
  set(FGSL_DIR $ENV{FGSL_DIR} CACHE STRING "FGSL directory")
endif ()

set(FGSL_LIB ${FGSL_DIR}/.libs CACHE STRING "FGSL library path")
set(FGSL_INC ${FGSL_DIR} CACHE STRING "FGSL include")
set(fgsl_lib ${FGSL_DIR}/.libs/libfgsl${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE STRING "FGSL lib")
