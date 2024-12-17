if(DEFINED ENV{MFEM_DIR} OR DEFINED MFEM_DIR)
  add_custom_target(MFEM ALL)
else()
  message(STATUS "Downloading and building MFEM")
endif()

set(MFEM_INCLUDE_DIRS ${MFEM_DIR} CACHE STRING "MFEM include")
set(MFEM_LIBRARIES ${MFEM_DIR}/libmfem${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE STRING "MFEM lib")
