# Dependency overrides for a hermetic build. With nothing set, libneo, fortfem,
# and (optionally) MFEM are fetched and built from source.
#
#   -DLIBNEO_REF=<ref>   libneo branch, tag, or SHA to fetch (empty = auto)
#   -DLIBNEO_PATH=<dir>  local libneo source directory (empty = fetch from git)
#   -DMFEM_DIR=<dir>      prebuilt MFEM directory (see SetupMFEM)
set(LIBNEO_REF "" CACHE STRING "libneo branch, tag, or SHA to fetch (empty = auto)")
set(LIBNEO_PATH "" CACHE PATH "local libneo source directory (empty = fetch from git)")
