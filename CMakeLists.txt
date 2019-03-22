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
set(CMAKE_Fortran_FLAGS "-O2 -g -fbounds-check -Wall")

### Define libraries
add_library(magdif SHARED
  PRELOAD/SRC/from_nrtype.f90
  PRELOAD/SRC/spline5_RZ.f90
  PRELOAD/SRC/field_divB0.f90
  PRELOAD/SRC/bdivfree.f90
  PRELOAD/SRC/orbit_mod.f90
  PRELOAD/SRC/mesh_mod.f90
  MHD/SRC/sparse_mod.f90
  MHD/SRC/magdif.f90
  MHD/SRC/magdif_config.f90
  ${SUITESPARSE_SRC_FILES}
  ${SUPERLU_SRC_FILES}
)

### Define executables
add_executable(magdif_test.x
  MHD/SRC/magdif_test.f90
)

### Linking
target_link_libraries(magdif_test.x
  magdif
  ${SUPERLU_LIBRARIES}
	${SUITESPARSE_LIBRARIES}
  blas
  lapack
)
