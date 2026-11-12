# Dependency overrides for a hermetic build. With nothing set, libneo, Triangle,
# and (optionally) MFEM are fetched and built from source.
#
#   -DITP_LIBNEO_REF=<ref>                 another libneo branch, tag, or SHA
#   -DFETCHCONTENT_SOURCE_DIR_LIBNEO=<dir> local libneo source directory
#   -DTRIANGLE_DIR=<dir>  prebuilt Triangle directory (see SetupTriangle)
#   -DMFEM_DIR=<dir>      prebuilt MFEM directory (see SetupMFEM)
set(MEPHIT_LIBNEO_REF d30174fb9a9c6520301463ef2406bd2cbfbcd55b)
