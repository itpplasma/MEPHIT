### Basic settings
cmake_minimum_required(VERSION 3.13.4)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

### Define the project
project(MEPHIT LANGUAGES Fortran C CXX)

### Specify paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

### Set compiler flags
if ("$ENV{FC}" MATCHES "ifort")
  # ifort configuration
  set(CMAKE_Fortran_FLAGS "-fpp -g -assume norealloc_lhs -traceback -check all -fpe-all=0")
  set(WARNING_FLAGS "-warn all")
  set(LEGACY_FLAGS "-check noarg_temp_created")
else ()
  # gfortran configuration
  set(CMAKE_Fortran_FLAGS "-O2 -march=native -cpp -g -fno-realloc-lhs -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow")
  set(WARNING_FLAGS "-Wall -Wextra -pedantic -fmax-errors=1")
  set(LEGACY_FLAGS "-fcheck=no-array-temps")
  set(CMAKE_C_FLAGS "-O2 -g -march=native -Wconversion -Wfloat-equal -Wshadow -Wredundant-decls -Winline ${WARNING_FLAGS}")
endif ()

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
  set(LIBNEO_DIR "${CMAKE_SOURCE_DIR}/../libneo")
else ()
  set(LIBNEO_DIR "$ENV{LIBNEO_DIR}")
endif ()
set(HDF5_tools_include_dir ${LIBNEO_DIR}/build)
set(HDF5_tools_lib_dir ${LIBNEO_DIR}/build/src/hdf5_tools)

if ("$ENV{FGSL_PATH}" STREQUAL "")
  if ("${CMAKE_SYSTEM}" MATCHES "Darwin")
    set(FGSL_PATH /usr/local)
  else ()
    set(FGSL_PATH /usr)
  endif ()
else ()
  set(FGSL_PATH $ENV{FGSL_PATH})
endif ()

### Find dependencies in standard directories
find_package(SuiteSparse REQUIRED)
find_package(SuperLU)
cmake_policy(SET CMP0074 NEW)
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_PATH}/lib/ CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_PATH}/include/fgsl/ CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)

# hdf5_tools
find_library(HDF5_tools_lib hdf5_tools ${HDF5_tools_lib_dir} NO_DEFAULT_PATH)

### External source files
include(${CMAKE_SOURCE_DIR}/CMakeSources.txt)
include_directories(${HDF5_INCLUDE_DIRS})

### Project source files
set(COMMON_FORTRAN_SRC_FILES
  src/bdivfree.f90
  src/binsrc.f90
  src/field_divB0.f90
  src/field_line_integration_for_SYNCH.f90
  src/from_nrtype.f90
  src/leastsq.f90
  src/lubksb.f
  src/ludcmp.f
  src/magdata_in_symfluxcoord.f90
  src/odeint_allroutines.f
  src/plag_coeff.f90
  src/points_2d.f90
  src/preload_for_SYNCH.f90
  src/spl_three_to_five_mod.f90
  src/spline5_RZ.f90
  src/external/netlib_mod.f90
  src/external/d1mach.f
  src/fftw3.f90
  src/sparse_mod.f90
)
set(MEPHIT_FORTRAN_SRC_FILES
  src/mephit_iter.f90
  src/mephit_conf.f90
  src/mephit_util.f90
  src/mephit_mesh.f90
  src/mephit_pert.f90
)

set(MEPHIT_C_SRC_FILES
  src/mephit_util.c
  src/mephit_fem.c
)

### Set compiler flags
set_source_files_properties(${COMMON_FORTRAN_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${LEGACY_FLAGS}")
set_source_files_properties(${MEPHIT_FORTRAN_SRC_FILES}
  src/mephit_mesher.f90
  src/mephit_test.f90
  src/vacfield.f90
  src/standardise_equilibrium.f90
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${WARNING_FLAGS}")
set_source_files_properties(src/mephit_run.c ${MEPHIT_C_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS}")

### Define library
add_library(mephit SHARED
  ${SUITESPARSE_SRC_FILES}
  # ${SUPERLU_SRC_FILES}
  ${COMMON_FORTRAN_SRC_FILES}
  ${MEPHIT_FORTRAN_SRC_FILES}
  ${MEPHIT_C_SRC_FILES}
)
target_include_directories(mephit PUBLIC ${HDF5_tools_include_dir} ${FGSL_INC} ${NETCDF_INCLUDES})
set(MEPHIT_LIBS
  # ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${HDF5_tools_lib}
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${NETCDF_LIBRARIES}
  ${gsl_lib}
  ${fgsl_lib}
  blas
  lapack
  fftw3
)
target_link_libraries(mephit ${MEPHIT_LIBS})

### Define executables
add_executable(mephit_run.x src/mephit_run.c)
add_executable(mephit_test.x src/mephit_test.f90)
add_executable(standardise_equilibrium.x src/standardise_equilibrium.f90)
add_executable(vacfield.x src/vacfield.f90)
target_link_libraries(mephit_run.x mephit)
target_link_libraries(mephit_test.x mephit)
target_link_libraries(standardise_equilibrium.x mephit)
target_link_libraries(vacfield.x mephit)

### Add helper script to binaries
configure_file(scripts/mephit.bash ${EXECUTABLE_OUTPUT_PATH}/mephit @ONLY)
