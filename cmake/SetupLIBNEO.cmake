if(DEFINED ENV{LIBNEO_DIR})
  set(LIBNEO_DIR $ENV{LIBNEO_DIR} CACHE STRING "libneo directory")
  add_custom_target(LIBNEO ALL)
else()
  ExternalProject_Add(
    LIBNEO
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG main
    PREFIX ${CMAKE_BINARY_DIR}/libneo
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/libneo${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/libmagfie${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/src/hdf5_tools/libhdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX}
  )
  set(LIBNEO_DIR ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build CACHE STRING "libneo directory")
endif()

set(magfie_include_dir ${LIBNEO_DIR}/include)
set(magfie_lib_dir ${LIBNEO_DIR})
set(HDF5_tools_lib_dir ${LIBNEO_DIR}/src/hdf5_tools)

set(neo_lib ${LIBNEO_DIR}/libneo${CMAKE_SHARED_LIBRARY_SUFFIX})
set(magfie_lib ${magfie_lib_dir}/libmagfie${CMAKE_SHARED_LIBRARY_SUFFIX})
set(HDF5_tools_lib ${HDF5_tools_lib_dir}/libhdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX})
