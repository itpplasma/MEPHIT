### Basic settings
cmake_minimum_required (VERSION 2.8.6)

### Define the project
project(NEO-EQ)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

### External source files
include(${CMAKE_SOURCE_DIR}/CMakeSources.txt)

### Specify paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

### Set compiler flags
set(CMAKE_Fortran_FLAGS "-O2 -g")

### Define libraries
add_library(magdif
  PRELOAD/SRC/from_nrtype.f90
  PRELOAD/SRC/spline5_RZ.f90
  PRELOAD/SRC/field_divB0.f90
  PRELOAD/SRC/bdivfree.f90
  PRELOAD/SRC/orbit_mod.f90
  PRELOAD/SRC/mesh_mod.f90
  MHD/SRC/sparse_mod.f90
  MHD/SRC/magdif.f90
)

### Define executables
add_executable(magdif_test
  MHD/SRC/magdif_test.f90
)

### Linking
target_link_libraries(magdif_test
  ${SuperLU_lib}
  ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
  ${SuiteSparseConfig_lib}
  magdif
  blas
  lapack
)

