### Basic settings
cmake_minimum_required(VERSION 3.13.4)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

### Define the project
project(MEPHIT LANGUAGES C CXX Fortran)

### Define options
option(WITH_MFEM "Compile with MFEM" FALSE)

### Specify paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIBRARY_OUTPUT_PATH})

### Set compiler flags
if ("$ENV{FC}" MATCHES "ifort")
  # ifort configuration
  set(CMAKE_Fortran_FLAGS "-fpp -g -assume norealloc_lhs -traceback -check all -fpe-all=0")
  set(WARNING_FLAGS "-warn all")
  set(LEGACY_FLAGS "-check noarg_temp_created")
else ()
  # gfortran configuration
  set(CMAKE_Fortran_FLAGS "-O2 -march=native -cpp -g -fno-realloc-lhs -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow")
  set(WARNING_FLAGS "-Wall -Wextra -pedantic -fmax-errors=1 -Wno-unused-label")
  set(LEGACY_FLAGS "-fcheck=no-array-temps")
  set(CMAKE_C_FLAGS "-O2 -g -march=native -Wconversion -Wfloat-equal -Wshadow -Wredundant-decls -Winline ${WARNING_FLAGS}")
endif ()

# Prevent linker error with hdf5tools on Ubuntu
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,--copy-dt-needed-entries")
endif()

### Define external paths
if ("${CMAKE_SYSTEM}" MATCHES "Darwin")
  set(tmpdir $TMPDIR)
else ()
  if ("$ENV{XDG_RUNTIME_DIR}" STREQUAL "")
    set(tmpdir /tmp/runtime-$USER)
  else ()
    set(tmpdir $XDG_RUNTIME_DIR)
  endif ()
endif ()

### Find dependencies in standard directories
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(SuiteSparse)
if(SUITESPARSE_FOUND)
  set(LEGACY_SUITESPARSE ON)
  message(STATUS "SuiteSparse found.")
else()
  message(STATUS "SuiteSparse not found. Trying new cmake search.")
  find_package(UMFPACK REQUIRED)
endif()

cmake_policy(SET CMP0074 NEW)
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
find_package(Boost 1.74 REQUIRED)
find_package(GSL REQUIRED)

### NetCDF
find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
				OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

message(STATUS "NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "NetCDF Fortran libs: " ${NETCDF_FLIBS})

# Replace space by semicolon in the Fortran libs
string(REPLACE " " ";" NETCDF_FLIBS ${NETCDF_FLIBS})

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

include(ExternalProject)
include(SetupCODE)
include(SetupTriangle)
if(WITH_MFEM)
  include(SetupMFEM)
  add_definitions(-DUSE_MFEM)
endif()
include(Util)
find_or_fetch(libneo)

### Project source files
# SuiteSparse
set(SUITESPARSE_SRC_FILES
  src/external/umf4_f77wrapper.c
  src/external/umf4_f77zwrapper.c
)
if(LEGACY_SUITESPARSE)
  set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
else()
  set(SUITESPARSE_FLAGS "")
endif()
set_source_files_properties(src/external/umf4_f77wrapper.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(src/external/umf4_f77zwrapper.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")

set(COMMON_FORTRAN_SRC_FILES
  src/external/netlib_mod.f90
  src/external/d1mach.f
  src/fftw3.f90
  src/field_line_integration_for_SYNCH.f90
  src/magdata_in_symfluxcoord.f90
  src/odeint_allroutines.f
  src/points_2d.f90
  src/sparse_mod.f90
)
set_source_files_properties(${COMMON_FORTRAN_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${LEGACY_FLAGS}")

set(MEPHIT_FORTRAN_SRC_FILES
  src/mephit_conf.f90
  src/mephit_util.f90
  src/mephit_mesh.F90
  src/mephit_pert.f90
  src/mephit_flr2.f90
  src/mephit_iter.F90
)
set_source_files_properties(src/mephit_test.f90 ${MEPHIT_FORTRAN_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${WARNING_FLAGS}")

set(MEPHIT_C_SRC_FILES
  src/mephit_util.c
  src/mephit_fem.c
  src/hyper1F1.c
)
set(MEPHIT_CPP_SRC_FILES
  src/mephit_fem.cpp
)
if(WITH_MFEM)
  set(MEPHIT_CPP_SRC_FILES ${MEPHIT_CPP_SRC_FILES}
    src/magnetic_differential_equation.cpp
    src/miscellaneous.cpp
  )
endif()
set_source_files_properties(src/mephit_run.c ${MEPHIT_C_SRC_FILES} ${MEPHIT_CPP_SRC_FILES}
	PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -DREAL=double -I${TRIANGLE_INCLUDE_DIR} -I${SUITESPARSE_INCLUDE_DIRS} -L${TRIANGLE_LIB_DIR}")


### Define library
add_library(mephit SHARED
  ${SUITESPARSE_SRC_FILES}
  ${COMMON_FORTRAN_SRC_FILES}
  ${MEPHIT_FORTRAN_SRC_FILES}
  ${MEPHIT_C_SRC_FILES}
  ${MEPHIT_CPP_SRC_FILES}
)
add_dependencies(mephit TRIANGLE)
target_include_directories(mephit PUBLIC
  ${magfie_include_dir}
  ${NETCDF_INCLUDES}
  ${NETCDF_INCLUDES}
  ${TRIANGLE_INCLUDE_DIR}
)
set(MEPHIT_LIBS
  hdf5::hdf5
  hdf5::hdf5_hl
  hdf5::hdf5_fortran
  hdf5::hdf5_hl_fortran
  GSL::gsl
  BLAS::BLAS
  LAPACK::LAPACK
  ${FFTW_LIBRARIES}
  ${triangle_lib}
  LIBNEO::magfie
  LIBNEO::neo
  LIBNEO::hdf5_tools
)
target_link_libraries(mephit ${MEPHIT_LIBS})
if(LEGACY_SUITESPARSE)
  target_link_libraries(mephit ${SUITESPARSE_LIBRARIES})
else()
  target_link_libraries(mephit SuiteSparse::umfpack)
endif()
if(WITH_MFEM)
  add_dependencies(mephit MFEM)
  target_include_directories(mephit PUBLIC ${MFEM_INCLUDE_DIRS})
  target_link_libraries(mephit ${MFEM_LIBRARIES})
endif()

### Define executables
add_executable(mephit_run.x src/mephit_run.c)
add_executable(mephit_test.x src/mephit_test.f90)
target_link_libraries(mephit_run.x mephit)
target_link_libraries(mephit_test.x mephit)

### Add helper script to binaries
configure_file(scripts/mephit.bash ${PROJECT_BINARY_DIR}/scripts/mephit.bash @ONLY)
configure_file(scripts/ff-mephit.bash ${PROJECT_BINARY_DIR}/scripts/ff-mephit.bash COPYONLY)
configure_file(scripts/maxwell_daemon.edp ${PROJECT_BINARY_DIR}/scripts/maxwell_daemon.edp COPYONLY)
configure_file(scripts/extmesh.edp ${PROJECT_BINARY_DIR}/scripts/extmesh.edp COPYONLY)
configure_file(scripts/mephit.gdb ${PROJECT_BINARY_DIR}/scripts/mephit.gdb COPYONLY)
configure_file(data/convexwall_asdex.dat ${PROJECT_BINARY_DIR}/data/convexwall_asdex.dat COPYONLY)
configure_file(data/convexwall_kilca.dat ${PROJECT_BINARY_DIR}/data/convexwall_kilca.dat COPYONLY)
configure_file(data/convexwall_mastu.dat ${PROJECT_BINARY_DIR}/data/convexwall_mastu.dat COPYONLY)
configure_file(data/field_divB0.inp ${PROJECT_BINARY_DIR}/data/field_divB0.inp COPYONLY)
configure_file(data/preload_for_SYNCH.inp ${PROJECT_BINARY_DIR}/data/preload_for_SYNCH.inp COPYONLY)
