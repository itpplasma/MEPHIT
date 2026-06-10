# Fetch fortfem for its Triangle-compatible mesher (target
# fortfem_triangle_compat, MIT-licensed). Pinned to the triangle-compat
# branch until lazy-fortran/fortfem PR #56 is merged and tagged; retarget
# GIT_TAG to that release afterwards.
include(FetchContent)

set(FORTFEM_TRIANGLE_COMPAT_ONLY ON CACHE BOOL "" FORCE)

FetchContent_Declare(fortfem
    GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
    GIT_TAG triangle-compat
)
FetchContent_MakeAvailable(fortfem)
