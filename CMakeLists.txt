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

if ("$ENV{KILCA_DIR}" STREQUAL "")
  set(KILCA_DIR "${CMAKE_SOURCE_DIR}/../KiLCA/build")
else ()
  set(KILCA_DIR "$ENV{KILCA_DIR}")
endif ()

if ("$ENV{LIBNEO_DIR}" STREQUAL "")
  set(LIBNEO_DIR "${CMAKE_SOURCE_DIR}/../libneo/build")
else ()
  set(LIBNEO_DIR "$ENV{LIBNEO_DIR}")
endif ()
set(HDF5_tools_include_dir ${LIBNEO_DIR})
set(HDF5_tools_lib_dir ${LIBNEO_DIR}/src/hdf5_tools)

if ("$ENV{MFEM_DIR}" STREQUAL "")
  set(MFEM_DIR "${CMAKE_SOURCE_DIR}/../mfem/build")
else ()
  set(MFEM_DIR "$ENV{MFEM_DIR}")
endif ()

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
find_package(LAPACK)
find_package(SuiteSparse REQUIRED)
find_package(SuperLU)
cmake_policy(SET CMP0074 NEW)
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)
find_package(MFEM REQUIRED)
find_package(Boost 1.74 REQUIRED)

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_PATH}/lib/ CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_PATH}/include/fgsl/ CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)

# KiLCA
find_library(kilca_lib libkilca.a ${KILCA_DIR} NO_DEFAULT_PATH)

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
  src/getIfunc.f90
  src/progonka.f90
  src/response_current.f90
  src/W2_arr.f90
)
set(MEPHIT_FORTRAN_SRC_FILES
  src/mephit_conf.f90
  src/mephit_util.f90
  src/mephit_mesh.f90
  src/mephit_pert.f90
  src/mephit_iter.f90
)

set(MEPHIT_C_SRC_FILES
  src/mephit_util.c
  src/mephit_fem.c
)

set(MEPHIT_CPP_SRC_FILES
  src/mephit_fem.cpp
  src/magnetic_differential_equation.cpp
  src/miscellaneous.cpp
)

### Set compiler flags
set_source_files_properties(${COMMON_FORTRAN_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${LEGACY_FLAGS}")
set_source_files_properties(${MEPHIT_FORTRAN_SRC_FILES}
  src/mephit_test.f90
  src/vacfield.f90
  src/standardise_equilibrium.f90
  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${WARNING_FLAGS} ${LEGACY_FLAGS}")
set_source_files_properties(src/mephit_run.c ${MEPHIT_C_SRC_FILES} ${MEPHIT_CPP_SRC_FILES}
  PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS}")

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
  ${HDF5_tools_include_dir}
  ${FGSL_INC}
  ${NETCDF_INCLUDES}
  ${MFEM_INCLUDE_DIRS}
)
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
  ${MFEM_LIBRARIES}
  ${kilca_lib}
  blas
  lapack
  fftw3
  triangle
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
configure_file(scripts/mephit.bash ${PROJECT_BINARY_DIR}/scripts/mephit.bash @ONLY)
configure_file(scripts/ff-mephit.bash ${PROJECT_BINARY_DIR}/scripts/ff-mephit.bash COPYONLY)
configure_file(scripts/maxwell_daemon.edp ${PROJECT_BINARY_DIR}/scripts/maxwell_daemon.edp COPYONLY)
configure_file(scripts/extmesh.edp ${PROJECT_BINARY_DIR}/scripts/extmesh.edp COPYONLY)
configure_file(scripts/mephit.gdb ${PROJECT_BINARY_DIR}/scripts/mephit.gdb COPYONLY)
configure_file(scripts/magdifplot.tex ${PROJECT_BINARY_DIR}/scripts/magdifplot.tex COPYONLY)
configure_file(scripts/mephit_plot.py ${PROJECT_BINARY_DIR}/scripts/mephit_plot.py COPYONLY)
configure_file(scripts/test_plot.py ${PROJECT_BINARY_DIR}/scripts/test_plot.py COPYONLY)
configure_file(data/convexwall_asdex.dat ${PROJECT_BINARY_DIR}/data/convexwall_asdex.dat COPYONLY)
configure_file(data/convexwall_kilca.dat ${PROJECT_BINARY_DIR}/data/convexwall_kilca.dat COPYONLY)
configure_file(data/field_divB0.inp ${PROJECT_BINARY_DIR}/data/field_divB0.inp COPYONLY)
configure_file(data/preload_for_SYNCH.inp ${PROJECT_BINARY_DIR}/data/preload_for_SYNCH.inp COPYONLY)
