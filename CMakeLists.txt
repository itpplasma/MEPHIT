### Basic settings
cmake_minimum_required (VERSION 2.8.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

### Define the project
project(NEO-EQ)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)


### Specify paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

### SuiteSparse and SuperLU
find_package(SuiteSparse REQUIRED)
find_package(SuperLU REQUIRED)

### External source files
include(${CMAKE_SOURCE_DIR}/CMakeSources.txt)

### Set compiler flags
# gfortran configuration
set(CMAKE_Fortran_FLAGS "-O2 -cpp -g -fcheck=all -ffpe-trap=invalid,zero,overflow")
set(WARNING_FLAGS "-Wall -Wextra")
# ifort configuration
# set(CMAKE_Fortran_FLAGS "-O2 -fpp -g -check all")
# set(WARNING_FLAGS "-warn all")

### Define libraries
add_library(magdif SHARED
  src/bdivfree.f90
  src/binsrc.f90
  src/field_divB0.f90
  src/field_line_integration_for_SYNCH.f90
  src/from_nrtype.f90
  src/leastsq.f90
  src/lubksb.f
  src/ludcmp.f
  src/magdata_in_symfluxcoord.f90
  src/mesh_mod.f90
  src/odeint_allroutines.f
  src/orbit_mod.f90
  src/plag_coeff.f90
  src/points_2d.f90
  src/preload_for_SYNCH.f90
  src/spl_three_to_five_mod.f90
  src/spline5_RZ.f90
  src/external/netlib_mod.f90
  src/external/d1mach.f
  src/arnoldi.f90
  src/sparse_mod.f90
  src/magdif.f90
  src/magdif_config.f90
  src/magdif_util.f90
  src/magdif_mesh_mod.f90
  ${SUITESPARSE_SRC_FILES}
  ${SUPERLU_SRC_FILES}
)

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_PATH}/lib/ CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_PATH}/include/fgsl/ CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
set(FGSL_FLAGS "-I${FGSL_INC}")

### Add helper script to binaries
configure_file(scripts/magdif.bash bin/magdif COPYONLY)

### Define executables
add_executable(magdif_test.x src/magdif_test.f90)
add_executable(standardise_equilibrium.x src/standardise_equilibrium.f90)
add_executable(magdif_mesher.x src/magdif_mesher.f90)
add_executable(vacfield.x src/vacfield.f90)

### Add compiler flags
set_target_properties(magdif PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${FGSL_FLAGS}")
set_source_files_properties(src/magdif.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/magdif_config.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/magdif_util.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/magdif_mesh_mod.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/magdif_mesher.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/magdif_test.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/vacfield.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")
set_source_files_properties(src/standardise_equilibrium.f90 PROPERTIES COMPILE_FLAGS "${WARNING_FLAGS}")

### Linking
set(MAGDIF_LIBS ${SUPERLU_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${gsl_lib} ${fgsl_lib} blas lapack)
target_link_libraries(magdif ${MAGDIF_LIBS})
target_link_libraries(magdif_test.x magdif)
target_link_libraries(standardise_equilibrium.x magdif)
target_link_libraries(magdif_mesher.x magdif)
target_link_libraries(vacfield.x magdif)
