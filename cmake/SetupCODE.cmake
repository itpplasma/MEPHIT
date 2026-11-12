# Dependency overrides for a hermetic build. With nothing set, libneo, Triangle,
# and (optionally) MFEM are fetched and built from source.
#
#   -DITP_LIBNEO_REF=<ref>                 another libneo branch, tag, or SHA
#   -DFETCHCONTENT_SOURCE_DIR_LIBNEO=<dir> local libneo source directory
#   -DTRIANGLE_DIR=<dir>  prebuilt Triangle directory (see SetupTriangle)
#   -DMFEM_DIR=<dir>      prebuilt MFEM directory (see SetupMFEM)
set(MEPHIT_LIBNEO_REF aac3f9ef324d5ebab5a4ead26c504cd6548c1bee)
# The fortio pinned by that libneo stops on any HDF5 dataset above 2 GiB because
# it never completes a short POSIX write; libneo adopted this fortio later.
set(FORTIO_REF "314872622fb28c259fb3c11fd200f67e7dc96227" CACHE STRING
    "fortio SHA that libneo builds against")
