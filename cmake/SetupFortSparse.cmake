# The application selects project-owned native factors only. Dense BLAS/LAPACK
# remain allowed; sparse-library factorization and helper fallbacks do not.
set(FORTSPARSE_REF "d5ab98a8b5d12550eb9aacf795eba97d782f1735" CACHE STRING
    "Pinned FortSparse revision")
set(FORTSPARSE_PATH "" CACHE PATH "Optional local FortSparse source override")

function(mephit_add_fortsparse)
  foreach(_backend SUPERLU UMFPACK_IPC MUMPS)
    if(FORTSPARSE_ENABLE_${_backend})
      message(FATAL_ERROR "MEPHIT native-only build requires FORTSPARSE_ENABLE_${_backend}=OFF")
    endif()
    set(FORTSPARSE_ENABLE_${_backend} OFF CACHE BOOL "Disabled by native-only MEPHIT" FORCE)
  endforeach()
  set(FORTSPARSE_USE_SYSTEM_DEPS ON CACHE BOOL "Use installed dense dependencies" FORCE)
  set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "Dependency tests are separate" FORCE)
  set(FORTSPARSE_BUILD_BENCHMARKS OFF CACHE BOOL "No dependency benchmarks" FORCE)
  set(FORTSPARSE_INSTALL_HELPER OFF CACHE BOOL "No sparse helper" FORCE)
  set(FORTSPARSE_BUILD_FREEFEM_PLUGIN ON CACHE BOOL "Native Maxwell solver module" FORCE)
  # Two-dimensional Maxwell/MFEM matrices need fill-reducing ordering. AMD
  # supplies only a BSD ordering permutation; all numerical factors stay native.
  set(FORTSPARSE_ENABLE_AMD ON CACHE BOOL "AMD ordering for native factors" FORCE)
  set(FORTSPARSE_ENABLE_METIS OFF CACHE BOOL "Native-only dependency baseline" FORCE)
  # MEPHIT's legacy global flag is not valid for a dependency that uses standard
  # allocatable assignment. Keep the existing application compiler policy local.
  string(REPLACE "-fno-realloc-lhs" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  string(REPLACE "-assume norealloc_lhs" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  if(FORTSPARSE_PATH)
    set(_fortsparse_source "${FORTSPARSE_PATH}")
  else()
    include(FetchContent)
    FetchContent_Declare(fortsparse
      GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
      GIT_TAG ${FORTSPARSE_REF})
    FetchContent_GetProperties(fortsparse)
    if(NOT fortsparse_POPULATED)
      FetchContent_Populate(fortsparse)
    endif()
    set(_fortsparse_source "${fortsparse_SOURCE_DIR}")
  endif()
  add_subdirectory("${_fortsparse_source}" "${CMAKE_BINARY_DIR}/fortsparse" EXCLUDE_FROM_ALL)
endfunction()
mephit_add_fortsparse()
