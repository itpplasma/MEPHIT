include(FetchContent)

set(DELAUNAY32_REF "ae2d0cf9167deefd0d180e39a30be3cc8f0b0698" CACHE STRING
    "Delaunay32 commit to fetch")
set(DELAUNAY32_PATH "" CACHE PATH
    "Local Delaunay32 source directory (empty = fetch from git)")

set(DELAUNAY32_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(DELAUNAY32_BUILD_EXTRAS OFF CACHE BOOL "" FORCE)
set(DELAUNAY32_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(DELAUNAY32_BUILD_TESTS OFF CACHE BOOL "" FORCE)

if(DELAUNAY32_PATH)
  set(delaunay32_SOURCE_DIR "${DELAUNAY32_PATH}")
else()
  FetchContent_Declare(
    delaunay32
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/morishuz/delaunay32.git
    GIT_TAG ${DELAUNAY32_REF}
  )
  FetchContent_Populate(delaunay32)
endif()

add_subdirectory(
  ${delaunay32_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/delaunay32
  EXCLUDE_FROM_ALL
)
set_target_properties(delaunay32 PROPERTIES POSITION_INDEPENDENT_CODE ON)
