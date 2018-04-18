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

### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(SuperLU_lib superlu ${SUPERLU_DIR}lib NO_DEFAULT_PATH)
find_library(Umfpack_lib umfpack ${SUITESPARSE_DIR}/UMFPACK/Lib NO_DEFAULT_PATH)
find_library(Amd_lib amd ${SUITESPARSE_DIR}/AMD/Lib NO_DEFAULT_PATH)
find_library(Cholmod_lib cholmod ${SUITESPARSE_DIR}/CHOLMOD/Lib NO_DEFAULT_PATH)
find_library(Colamd_lib colamd ${SUITESPARSE_DIR}/COLAMD/Lib NO_DEFAULT_PATH)
find_library(Camd_lib camd ${SUITESPARSE_DIR}/CAMD/Lib NO_DEFAULT_PATH)
#find_library(Metis_lib metis ${SUITESPARSE_DIR}/metis-5.1.0/Lib ${METIS_PATH} NO_DEFAULT_PATH)
find_library(Ccolamd_lib ccolamd ${SUITESPARSE_DIR}/CCOLAMD/Lib NO_DEFAULT_PATH)
# required with newer SuiteSparse versions:
find_library(SuiteSparseConfig_lib suitesparseconfig ${SUITESPARSE_DIR}/SuiteSparse_config NO_DEFAULT_PATH)

### Set compiler flags
set(CMAKE_Fortran_FLAGS "-O2 -g -fbounds-check")

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
  ${SuperLU_lib}
  ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Ccolamd_lib}
  ${SuiteSparseConfig_lib} # required with newer SuiteSparse versions
  #${Metis_lib}
  blas
  lapack
)

