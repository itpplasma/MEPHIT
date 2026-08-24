# Dependency overrides for a hermetic build. With nothing set, libneo, Triangle,
# and (optionally) MFEM are fetched and built from source.
#
#   -DLIBNEO_REF=<ref>   libneo branch, tag, or SHA to fetch
#   -DLIBNEO_PATH=<dir>  local libneo source directory (empty = fetch from git)
#   -DTRIANGLE_DIR=<dir>  prebuilt Triangle directory (see SetupTriangle)
#   -DMFEM_DIR=<dir>      prebuilt MFEM directory (see SetupMFEM)
set(LIBNEO_REF "ff51ec3b7dcd410a0b7f88defbecb492590c0268" CACHE STRING
    "libneo branch, tag, or SHA to fetch")
set(LIBNEO_PATH "" CACHE PATH "local libneo source directory (empty = fetch from git)")
