### Basic settings
cmake_minimum_required(VERSION 3.13.4)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

### Define the project
project(MEPHIT LANGUAGES Fortran C CXX)

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
  set(BLA_VENDOR OpenBLAS)
endif ()

# Prevent linker error with hdf5tools on Ubuntu
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--copy-dt-needed-entries")

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

if ("$ENV{LIBNEO_DIR}" STREQUAL "")
  set(LIBNEO_DIR "${CMAKE_SOURCE_DIR}/../libneo/build")
else ()
  set(LIBNEO_DIR "$ENV{LIBNEO_DIR}")
endif ()
set(magfie_include_dir ${LIBNEO_DIR}/include)
set(magfie_lib_dir ${LIBNEO_DIR})
set(HDF5_tools_include_dir ${LIBNEO_DIR})
set(HDF5_tools_lib_dir ${LIBNEO_DIR}/src/hdf5_tools)

if ("$ENV{MFEM_DIR}" STREQUAL "")
  set(MFEM_DIR "${CMAKE_SOURCE_DIR}/../mfem/build")
else ()
  set(MFEM_DIR "$ENV{MFEM_DIR}")
endif ()

if ("$ENV{FGSL_DIR}" STREQUAL "")
  set(FGSL_DIR ${CMAKE_SOURCE_DIR}/../fgsl)
else ()
  set(FGSL_DIR $ENV{FGSL_DIR})
endif ()

### Find dependencies in standard directories
find_package(LAPACK)
find_package(SuiteSparse REQUIRED)
find_package(SuperLU)
cmake_policy(SET CMP0074 NEW)
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)
find_package(MFEM REQUIRED)
find_package(Boost 1.74 REQUIRED)
find_package(Triangle REQUIRED)

message(STATUS "TRIANGLE_INCLUDE_DIR: ${TRIANGLE_INCLUDE_DIR}")

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_DIR}/.libs CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_DIR} CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)

# hdf5_tools
find_library(HDF5_tools_lib hdf5_tools ${HDF5_tools_lib_dir} NO_DEFAULT_PATH)

# magfie
find_library(magfie_lib magfie ${magfie_lib_dir} NO_DEFAULT_PATH)

### Project source files
# SuperLU
set(SUPERLU_SRC_FILES
  src/external/c_fortran_dgssv.c
  src/external/c_fortran_zgssv.c
)
set(SUPERLU_FLAGS "${CFLAGS_DEBUG} ${CFLAGS} -I${SUPERLU_INCLUDES}")
set_source_files_properties(${SUPERLU_SRC_FILES} PROPERTIES COMPILE_FLAGS ${SUPERLU_FLAGS})

# SuiteSparse
set(SUITESPARSE_SRC_FILES
  src/external/umf4_f77wrapper.c
  src/external/umf4_f77zwrapper.c
)
set(SUITESPARSE_FLAGS -I${SUITESPARSE_INCLUDE_DIRS})
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
  src/getIfunc.f90
  src/progonka.f90
  src/response_current.f90
  src/W2_arr.f90
)
set_source_files_properties(${COMMON_FORTRAN_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${LEGACY_FLAGS}")

set(MEPHIT_FORTRAN_SRC_FILES
  src/mephit_conf.f90
  src/mephit_util.f90
  src/mephit_mesh.f90
  src/mephit_pert.f90
  src/mephit_iter.f90
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
  src/magnetic_differential_equation.cpp
  src/miscellaneous.cpp
)
set_source_files_properties(src/mephit_run.c ${MEPHIT_C_SRC_FILES} ${MEPHIT_CPP_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -I${TRIANGLE_INCLUDE_DIR}")

### Define library
add_library(mephit SHARED
  ${SUITESPARSE_SRC_FILES}
  # ${SUPERLU_SRC_FILES}
  ${COMMON_FORTRAN_SRC_FILES}
  ${MEPHIT_FORTRAN_SRC_FILES}
  ${MEPHIT_C_SRC_FILES}
  ${MEPHIT_CPP_SRC_FILES}
)
target_include_directories(mephit PUBLIC
  ${HDF5_INCLUDE_DIRS}
  ${HDF5_tools_include_dir}
  ${magfie_include_dir}
  ${FGSL_INC}
  ${NETCDF_INCLUDES}
  ${MFEM_INCLUDE_DIRS}
  ${TRIANGLE_INCLUDE_DIR}
)
set(MEPHIT_LIBS
  # ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${magfie_lib}
  ${HDF5_tools_lib}
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${NETCDF_LIBRARIES}
  ${gsl_lib}
  ${fgsl_lib}
  ${MFEM_LIBRARIES}
  blas
  lapack
  fftw3
  triangle
)
target_link_libraries(mephit ${MEPHIT_LIBS})

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
