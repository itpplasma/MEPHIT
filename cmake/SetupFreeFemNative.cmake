# The optional SDK module belongs to FortSparse; this application owns only its
# script location and launcher. No private solver bridge or external fallback.
if(NOT TARGET fortsparse_freefem)
  message(FATAL_ERROR "Pinned FortSparse must provide its native FreeFem module")
endif()
set_target_properties(fortsparse_freefem PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/scripts)
find_program(FREEFEM_EXECUTABLE NAMES FreeFem++ REQUIRED)
