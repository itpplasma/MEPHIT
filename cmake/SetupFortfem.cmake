# Configure the native FortFEM backend from a reproducible source revision.
# A local source tree is useful for development; disabling the fetch keeps the
# installed-package workflow available for downstream builds.
include(FetchContent)

set(FORTFEM_REF "605dc7f056aa9b15cf8e0673eaa13b13ce12a273" CACHE STRING
    "fortfem commit to fetch")
set(FORTFEM_PATH "" CACHE PATH
    "local fortfem source directory (empty = fetch or find package)")
option(FORTFEM_FETCH "Fetch the pinned FortFEM source when no local path is set" ON)

if(FORTFEM_PATH)
    add_subdirectory(${FORTFEM_PATH} ${CMAKE_CURRENT_BINARY_DIR}/fortfem EXCLUDE_FROM_ALL)
elseif(FORTFEM_FETCH)
    if(CMAKE_VERSION VERSION_LESS 3.22)
        message(FATAL_ERROR
            "WITH_FORTFEM=ON with source fetching requires CMake 3.22 or newer")
    endif()
    FetchContent_Declare(fortfem
        GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
        GIT_TAG ${FORTFEM_REF}
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortfem)
else()
    find_package(fortfem CONFIG REQUIRED)
endif()

if(NOT TARGET fortfem::capi)
    message(FATAL_ERROR "FortFEM was configured without the fortfem::capi target")
endif()

# MEPHIT uses -fno-realloc-lhs for its legacy Fortran sources.  FortFEM and
# its in-tree FortSparse dependency use standard allocatable assignment, so
# restore the standard GNU behavior for those targets only.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(TARGET fortfem_capi)
        target_compile_options(fortfem_capi PRIVATE
            "$<$<COMPILE_LANGUAGE:Fortran>:-frealloc-lhs>")
    endif()
    if(TARGET fortsparse)
        target_compile_options(fortsparse PRIVATE
            "$<$<COMPILE_LANGUAGE:Fortran>:-frealloc-lhs>")
    endif()
endif()
