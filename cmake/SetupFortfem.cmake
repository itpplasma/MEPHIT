# Fetch fortfem for its Triangle-compatible mesher (target
# fortfem_triangle_compat, MIT-licensed). This is the current fortfem main
# commit as of 2026-08-17, pinned so the MEPHIT branch remains reproducible.
include(FetchContent)

set(FORTFEM_REF "605dc7f056aa9b15cf8e0673eaa13b13ce12a273" CACHE STRING
    "fortfem commit to fetch")
set(FORTFEM_PATH "" CACHE PATH
    "local fortfem source directory (empty = fetch from git)")
set(FORTFEM_TRIANGLE_COMPAT_ONLY ON CACHE BOOL "" FORCE)

if(FORTFEM_PATH)
    add_subdirectory(${FORTFEM_PATH} ${CMAKE_CURRENT_BINARY_DIR}/fortfem EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(fortfem
        GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
        GIT_TAG ${FORTFEM_REF}
    )
    FetchContent_MakeAvailable(fortfem)
endif()
