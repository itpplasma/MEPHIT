# The optional SDK module belongs to FortSparse; this application owns only its
# script location and launcher. No private solver bridge or external fallback.
if(NOT TARGET fortsparse_freefem)
  message(FATAL_ERROR "Pinned FortSparse must provide its native FreeFem module")
endif()
set_target_properties(fortsparse_freefem PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/scripts)
find_program(FREEFEM_EXECUTABLE NAMES FreeFem++ REQUIRED)

# SDK/runtime layouts vary: Debian places BinaryIO in lib/freefem++ and idp
# beneath the SDK include directory, while source installs use lib/ff++/VERSION.
# Resolve these resources for the configured prefix, preserving that choice in
# both the launcher and tests. Never substitute another solver/runtime.
set(FREEFEM_BINARYIO_DIR "" CACHE PATH "Optional matching FreeFem BinaryIO directory")
set(FREEFEM_IDP_DIR "" CACHE PATH "Optional matching FreeFem idp directory")
get_filename_component(_ff_bin_dir "${FREEFEM_EXECUTABLE}" DIRECTORY)
get_filename_component(_ff_prefix "${_ff_bin_dir}" DIRECTORY)
get_target_property(_ff_sdk_dirs fortsparse_freefem INCLUDE_DIRECTORIES)
set(_ff_binaryio_hints "${_ff_prefix}/lib/freefem++")
set(_ff_idp_hints)
foreach(_ff_include IN LISTS _ff_sdk_dirs)
  list(PREPEND _ff_binaryio_hints "${_ff_include}/../lib")
  list(APPEND _ff_idp_hints "${_ff_include}/idp" "${_ff_include}/../idp")
endforeach()
if(FREEFEM_BINARYIO_DIR)
  set(_ff_binaryio_hints "${FREEFEM_BINARYIO_DIR}")
endif()
if(FREEFEM_IDP_DIR)
  set(_ff_idp_hints "${FREEFEM_IDP_DIR}")
endif()
find_file(_mephit_binaryio NAMES BinaryIO.so HINTS ${_ff_binaryio_hints}
  REQUIRED NO_CACHE NO_DEFAULT_PATH)
find_file(_mephit_getargv NAMES getARGV.idp HINTS ${_ff_idp_hints}
  REQUIRED NO_CACHE NO_DEFAULT_PATH)
get_filename_component(MEPHIT_FREEFEM_PLUGIN_DIR "${_mephit_binaryio}" DIRECTORY)
get_filename_component(MEPHIT_FREEFEM_IDP_DIR "${_mephit_getargv}" DIRECTORY)
message(STATUS "FreeFem BinaryIO directory: ${MEPHIT_FREEFEM_PLUGIN_DIR}")
message(STATUS "FreeFem idp directory: ${MEPHIT_FREEFEM_IDP_DIR}")
