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
  src/bdivfree.f90
  src/binsrc.f90
  src/brent.f90
  src/field_divB0.f90
  src/field_line_integration_for_SYNCH.f90
  src/from_nrtype.f90
  src/in_out_cw.f90
  src/inside_of_wall.f90
  src/leastsq.f90
  src/lubksb.f
  src/ludcmp.f
  src/magdata_in_symfluxcoord.f90
  src/mesh_mod.f90
  src/mods4ax.f90
  src/odeint_allroutines.f
  src/orbit_mod.f90
  src/plag_coeff.f90
  src/points_2d.f90
  src/polylag.f
  src/preload_for_SYNCH.f90
  src/rtbis.f90
  src/spl_three_to_five_mod.f90
  src/spline5_RZ.f90
  src/twostreights.f90
  src/arnoldi.f90
  src/check_in_tri.f90
  src/sub_orbit_fast.f90
  src/sparse_mod.f90
  src/magdif.f90
  src/magdif_config.f90
  src/magdif_util.f90
  ${SUITESPARSE_SRC_FILES}
  ${SUPERLU_SRC_FILES}
)

# GSL
find_library(gsl_lib gsl ${GSL_LIB})

# FGSL
set(FGSL_LIB ${FGSL_PATH}/lib/ CACHE STRING "FGSL lib")
set(FGSL_INC ${FGSL_PATH}/include/fgsl/ CACHE STRING "FGSL include")
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
set(FGSL_FLAGS "${CMAKE_Fortran_FLAGS} -I${FGSL_INC}")

### Define executables
add_executable(magdif_test.x src/magdif_test.f90)
add_executable(minimal_example.x src/minimal_example.f90)
add_executable(geomint_mesh.x src/geomint_mesh.f90)
add_executable(axis.x src/axis.f90)
add_executable(tri_mesh.x src/tri_mesh.f90)
add_executable(readcarre_m.x src/readcarre_m.f90)
add_executable(writeout.x src/writeout.f90)
add_executable(vacfield.x src/vacfield.f90)
add_executable(result_spectrum.x src/result_spectrum.f90)

set_target_properties(minimal_example.x PROPERTIES COMPILE_FLAGS ${FGSL_FLAGS})

### Linking
set(MAGDIF_LIBS magdif ${SUPERLU_LIBRARIES} ${SUITESPARSE_LIBRARIES} blas lapack)
target_link_libraries(magdif_test.x ${MAGDIF_LIBS})
target_link_libraries(minimal_example.x ${MAGDIF_LIBS} ${gsl_lib} ${fgsl_lib})
target_link_libraries(geomint_mesh.x ${MAGDIF_LIBS})
target_link_libraries(axis.x ${MAGDIF_LIBS})
target_link_libraries(tri_mesh.x ${MAGDIF_LIBS})
target_link_libraries(readcarre_m.x ${MAGDIF_LIBS})
target_link_libraries(writeout.x ${MAGDIF_LIBS})
target_link_libraries(vacfield.x ${MAGDIF_LIBS})
target_link_libraries(result_spectrum.x ${MAGDIF_LIBS})
