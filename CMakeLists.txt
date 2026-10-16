### Basic settings
cmake_minimum_required(VERSION 3.22)
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
  set(CMAKE_Fortran_FLAGS "-O2 -march=native -g -fno-realloc-lhs -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow")
  set(WARNING_FLAGS "-Wall -Wextra -pedantic -fmax-errors=1 -Wno-unused-label")
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

### Find dependencies in standard directories
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include(SetupFortSparse)
include(SetupFreeFemNative)

cmake_policy(SET CMP0074 NEW)
find_package(Boost 1.74 REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)

include(ExternalProject)
include(SetupCODE)
include(SetupTriangle)
if(WITH_MFEM)
  include(SetupMFEM)
endif()
include(Util)
find_or_fetch(libneo)

### Project source files
set(COMMON_FORTRAN_SRC_FILES
  src/external/netlib_mod.f90
  src/external/d1mach.f
  src/fftw3.f90
  src/field_line_integration_for_SYNCH.f90
  src/magdata_in_symfluxcoord.f90
  src/points_2d.f90
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
set_source_files_properties(src/mephit_test.f90 src/mephit_test_refine.f90 ${MEPHIT_FORTRAN_SRC_FILES}
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
    src/mephit_native_mfem.cpp
  )
endif()
set_source_files_properties(src/mephit_run.c ${MEPHIT_C_SRC_FILES} ${MEPHIT_CPP_SRC_FILES}
	PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -DREAL=double -I${TRIANGLE_INCLUDE_DIR} -L${TRIANGLE_LIB_DIR}")


### Define library
add_library(mephit_sparse STATIC src/sparse_mod.f90)
set_target_properties(mephit_sparse PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mephit_sparse PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(mephit_sparse PUBLIC fortsparse)

add_library(mephit SHARED
  ${COMMON_FORTRAN_SRC_FILES}
  ${MEPHIT_FORTRAN_SRC_FILES}
  ${MEPHIT_C_SRC_FILES}
  ${MEPHIT_CPP_SRC_FILES}
)
add_dependencies(mephit TRIANGLE)
target_include_directories(mephit PUBLIC
  ${magfie_include_dir}
  ${TRIANGLE_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIRS}
)
set(MEPHIT_LIBS
  mephit_sparse
  GSL::gsl
  BLAS::BLAS
  LAPACK::LAPACK
  ${FFTW_LIBRARIES}
  ${triangle_lib}
  LIBNEO::magfie
  LIBNEO::odeint
  LIBNEO::hdf5_tools
)
target_link_libraries(mephit ${MEPHIT_LIBS})
if(WITH_MFEM)
  target_compile_features(mephit PRIVATE cxx_std_17)
  target_compile_definitions(mephit PRIVATE USE_MFEM)
  add_dependencies(mephit MFEM)
  target_include_directories(mephit PUBLIC ${MFEM_INCLUDE_DIRS})
  target_link_libraries(mephit ${MFEM_LIBRARIES})
endif()

### Define executables
add_executable(mephit_run.x src/mephit_run.c)
add_executable(mephit_test.x src/mephit_test.f90)
add_executable(mephit_post.x src/mephit_post.f90)
add_executable(mephit_test_refine.x src/mephit_test_refine.f90)
target_link_libraries(mephit_run.x mephit)
target_link_libraries(mephit_test.x mephit)
target_link_libraries(mephit_post.x mephit)
target_link_libraries(mephit_test_refine.x mephit)
include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_executable(test_native_sparse test/test_native_sparse.f90)
  target_link_libraries(test_native_sparse PRIVATE mephit_sparse)
  add_test(NAME native_sparse_adapter COMMAND test_native_sparse)
  add_test(NAME native_sparse_matched_probe COMMAND test_native_sparse matched)
  set_tests_properties(native_sparse_matched_probe PROPERTIES TIMEOUT 3
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  set_tests_properties(native_sparse_adapter PROPERTIES TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  add_test(NAME native_sparse_failures COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/test/test_native_sparse.py $<TARGET_FILE:test_native_sparse>)
  set_tests_properties(native_sparse_failures PROPERTIES TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  add_test(NAME native_sparse_linkage COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/test/test_native_sparse.py --linkage $<TARGET_FILE:mephit>)
  set_tests_properties(native_sparse_linkage PROPERTIES TIMEOUT 1)
  add_test(NAME native_freefem_adapter COMMAND ${FREEFEM_EXECUTABLE} -nw -nowait
    ${CMAKE_SOURCE_DIR}/test/test_native_freefem.edp)
  set_tests_properties(native_freefem_adapter PROPERTIES TIMEOUT 2
    ENVIRONMENT "FF_LOADPATH=${PROJECT_BINARY_DIR}/scripts;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  add_test(NAME native_freefem_singular COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/test/test_native_sparse.py --freefem-singular
    ${FREEFEM_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/test_native_freefem_singular.edp)
  set_tests_properties(native_freefem_singular PROPERTIES TIMEOUT 2
    ENVIRONMENT "FF_LOADPATH=${PROJECT_BINARY_DIR}/scripts;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  add_test(NAME native_freefem_linkage COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/test/test_native_sparse.py --linkage $<TARGET_FILE:fortsparse_freefem>)
  set_tests_properties(native_freefem_linkage PROPERTIES TIMEOUT 1)
  if(WITH_MFEM)
    add_executable(test_native_mfem test/test_native_mfem.cpp
      src/mephit_native_mfem.cpp src/miscellaneous.cpp)
    target_compile_features(test_native_mfem PRIVATE cxx_std_17)
    add_dependencies(test_native_mfem MFEM)
    target_include_directories(test_native_mfem PRIVATE src ${MFEM_INCLUDE_DIRS})
    target_link_libraries(test_native_mfem PRIVATE fortsparse ${MFEM_LIBRARIES})
    add_test(NAME native_mfem_adapter COMMAND test_native_mfem)
    set_tests_properties(native_mfem_adapter PROPERTIES TIMEOUT 1
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endif()
  add_executable(test_process_supervisor
    src/mephit_run.c src/mephit_util.c test/test_process_stub.c)
  target_include_directories(test_process_supervisor PRIVATE src)
  target_link_libraries(test_process_supervisor GSL::gsl)
  add_executable(test_pipe_client test/test_pipe_client.c)
  target_include_directories(test_pipe_client PRIVATE src)
  target_link_libraries(test_pipe_client mephit)
  add_test(NAME freefem_failure COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/test/test_freefem_failure.py
    $<TARGET_FILE:test_process_supervisor> $<TARGET_FILE:test_pipe_client>)
  set_tests_properties(freefem_failure PROPERTIES TIMEOUT 30)
  add_test(NAME refine_unit_partition_edge COMMAND mephit_test_refine.x)
  add_executable(test_mesh_quadrature_cache test/test_mesh_quadrature_cache.f90)
  target_link_libraries(test_mesh_quadrature_cache mephit)
  add_test(NAME mesh_quadrature_cache COMMAND test_mesh_quadrature_cache)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/mesh-quadrature-cache)
  set_tests_properties(mesh_quadrature_cache PROPERTIES
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/mesh-quadrature-cache)
  add_executable(test_helical_current test/test_helical_current.f90)
  target_link_libraries(test_helical_current mephit)
  add_test(NAME helical_current COMMAND test_helical_current)
  add_executable(test_layer_width test/test_layer_width.f90)
  target_link_libraries(test_layer_width mephit)
  add_test(NAME layer_width COMMAND test_layer_width)
  add_executable(test_resonant_current test/test_resonant_current.f90)
  target_link_libraries(test_resonant_current mephit)
  add_test(NAME resonant_current COMMAND test_resonant_current)
  # sheet (shielding) current response against the analytic cylinder; needs FreeFem++,
  # python3 with h5py and scipy, and runs the meshing, one iteration and the sheet phase
  find_program(FREEFEM_EXECUTABLE NAMES FreeFem++ FreeFem++-mpi)
  if(FREEFEM_EXECUTABLE)
    add_test(NAME sheet_cylinder COMMAND ${CMAKE_SOURCE_DIR}/test/sheet_cylinder/run_test.sh
      ${PROJECT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${PROJECT_BINARY_DIR}/sheet-cylinder)
    set_tests_properties(sheet_cylinder PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1800
      LABELS freefem)
    # boundary sheet of the moved equilibrium current against the analytic cylinder
    add_test(NAME edge_cylinder COMMAND ${CMAKE_SOURCE_DIR}/test/edge_cylinder/run_test.sh
      ${PROJECT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${PROJECT_BINARY_DIR}/edge-cylinder)
    set_tests_properties(edge_cylinder PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1800
      LABELS freefem)
  endif()
endif()

### Add helper script to binaries
configure_file(scripts/mephit.bash ${PROJECT_BINARY_DIR}/scripts/mephit.bash @ONLY)
configure_file(scripts/ff-mephit.bash ${PROJECT_BINARY_DIR}/scripts/ff-mephit.bash @ONLY)
configure_file(scripts/maxwell_daemon.edp ${PROJECT_BINARY_DIR}/scripts/maxwell_daemon.edp COPYONLY)
configure_file(scripts/extmesh.edp ${PROJECT_BINARY_DIR}/scripts/extmesh.edp COPYONLY)
configure_file(scripts/mephit.gdb ${PROJECT_BINARY_DIR}/scripts/mephit.gdb COPYONLY)
configure_file(data/convexwall_asdex.dat ${PROJECT_BINARY_DIR}/data/convexwall_asdex.dat COPYONLY)
configure_file(data/convexwall_iter.dat ${PROJECT_BINARY_DIR}/data/convexwall_iter.dat COPYONLY)
configure_file(data/convexwall_kilca.dat ${PROJECT_BINARY_DIR}/data/convexwall_kilca.dat COPYONLY)
configure_file(data/convexwall_mastu.dat ${PROJECT_BINARY_DIR}/data/convexwall_mastu.dat COPYONLY)
configure_file(data/field_divB0.inp ${PROJECT_BINARY_DIR}/data/field_divB0.inp COPYONLY)
configure_file(data/preload_for_SYNCH.inp ${PROJECT_BINARY_DIR}/data/preload_for_SYNCH.inp COPYONLY)
