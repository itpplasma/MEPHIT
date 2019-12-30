### Basic settings
cmake_minimum_required (VERSION 2.8.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

### Define the project
project(NEO-EQ)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in.sample)


### Specify paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

### SuiteSparse and SuperLU
find_package(SuiteSparse REQUIRED)
find_package(SuperLU REQUIRED)

### External source files
include(${CMAKE_SOURCE_DIR}/CMakeSources.txt)

### Set compiler flags
set(CMAKE_Fortran_FLAGS "-O2 -cpp -g -fbounds-check -Wall")  # -ffpe-trap=invalid,zero,overflow
# set(CMAKE_Fortran_FLAGS "-O2 -fpp -g -check all -warn all")  # ifort configuration

### Define libraries
add_library(magdif SHARED
  PRELOAD/SRC/from_nrtype.f90
  PRELOAD/SRC/spline5_RZ.f90
  PRELOAD/SRC/field_divB0.f90
  PRELOAD/SRC/bdivfree.f90
  PRELOAD/SRC/orbit_mod.f90
  PRELOAD/SRC/mesh_mod.f90
  PRELOAD/SRC/binsrc.f90
  PRELOAD/SRC/plag_coeff.f90
  RUN/SRC/arnoldi.f90
  MHD/SRC/sparse_mod.f90
  MHD/SRC/magdif.f90
  MHD/SRC/magdif_config.f90
  MHD/SRC/magdif_util.f90
  ${SUITESPARSE_SRC_FILES}
  ${SUPERLU_SRC_FILES}
)
add_library(geomint_mesh SHARED
  PRELOAD/SRC/field_divB0.f90
  PRELOAD/SRC/binsrc.f90
  PRELOAD/SRC/field_line_integration_for_SYNCH.f90
  PRELOAD/SRC/magdata_in_symfluxcoord.f90
  PRELOAD/SRC/odeint_allroutines.f
  PRELOAD/SRC/plag_coeff.f90
  PRELOAD/SRC/points_2d.f90
  PRELOAD/SRC/preload_for_SYNCH.f90
  PRELOAD/SRC/spl_three_to_five_mod.f90
)

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_PATH}/lib/ CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_PATH}/include/fgsl/ CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
set(FGSL_FLAGS "${CMAKE_Fortran_FLAGS} -I${FGSL_INC}")

### Define executables
add_executable(magdif_test.x
  MHD/SRC/magdif_test.f90
)
add_executable(minimal_example.x
  PRELOAD/SRC/minimal_example.f90
)
add_executable(geomint_mesh.x
  PRELOAD/SRC/geomint_mesh.f90
)

set_target_properties(minimal_example.x PROPERTIES COMPILE_FLAGS ${FGSL_FLAGS})

### Linking
target_link_libraries(magdif_test.x
  magdif
  ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  blas
  lapack
)
target_link_libraries(minimal_example.x
  magdif
  ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  blas
  lapack
  ${gsl_lib} ${fgsl_lib}
)
target_link_libraries(geomint_mesh.x
  magdif
  geomint_mesh
  ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  blas
  lapack
)
