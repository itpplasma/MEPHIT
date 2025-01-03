
# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2013 Google Inc. All rights reserved.
# http://code.google.com/p/ceres-solver/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: alexs.mac@gmail.com (Alex Stewart)
#
# FindSuiteSparse.cmake - Find SuiteSparse libraries & dependencies.
#
# This module defines the following variables:
#
# SUITESPARSE_FOUND: TRUE iff SuiteSparse and all dependencies have been found.
# SUITESPARSE_INCLUDE_DIRS: Include directories for all SuiteSparse components.
# SUITESPARSE_LIBRARIES: Libraries for all SuiteSparse component libraries and
#                        dependencies.
# SUITESPARSE_VERSION: Extracted from UFconfig.h (<= v3) or
#                      SuiteSparse_config.h (>= v4).
# SUITESPARSE_MAIN_VERSION: Equal to 4 if SUITESPARSE_VERSION = 4.2.1
# SUITESPARSE_SUB_VERSION: Equal to 2 if SUITESPARSE_VERSION = 4.2.1
# SUITESPARSE_SUBSUB_VERSION: Equal to 1 if SUITESPARSE_VERSION = 4.2.1
#
# SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION: TRUE iff running
#     on Ubuntu, SUITESPARSE_VERSION is 3.4.0 and found SuiteSparse is a system
#     install, in which case found version of SuiteSparse cannot be used to link
#     a shared library due to a bug (static linking is unaffected).
#
# The following variables control the behaviour of this module:
#
# SUITESPARSE_INCLUDE_DIR_HINTS: List of additional directories in which to
#                                search for SuiteSparse includes,
#                                e.g: /timbuktu/include.
# SUITESPARSE_LIBRARY_DIR_HINTS: List of additional directories in which to
#                                search for SuiteSparse libraries,
#                                e.g: /timbuktu/lib.
#
# The following variables define the presence / includes & libraries for the
# SuiteSparse components searched for, the SUITESPARSE_XX variables are the
# union of the variables for all components.
#
# == Symmetric Approximate Minimum Degree (AMD)
# AMD_FOUND
# AMD_INCLUDE_DIR
# AMD_LIBRARY
#
# == Constrained Approximate Minimum Degree (CAMD)
# CAMD_FOUND
# CAMD_INCLUDE_DIR
# CAMD_LIBRARY
#
# == Column Approximate Minimum Degree (COLAMD)
# COLAMD_FOUND
# COLAMD_INCLUDE_DIR
# COLAMD_LIBRARY
#
# Constrained Column Approximate Minimum Degree (CCOLAMD)
# CCOLAMD_FOUND
# CCOLAMD_INCLUDE_DIR
# CCOLAMD_LIBRARY
#
# == Sparse Supernodal Cholesky Factorization and Update/Downdate (CHOLMOD)
# CHOLMOD_FOUND
# CHOLMOD_INCLUDE_DIR
# CHOLMOD_LIBRARY
#
# == Multifrontal LU Factorization (UMFPACK)
# UMFPACK_FOUND
# UMFPACK_INCLUDE_DIR
# UMFPACK_LIBRARY
#
# == Multifrontal Sparse QR (SuiteSparseQR)
# SUITESPARSEQR_FOUND
# SUITESPARSEQR_INCLUDE_DIR
# SUITESPARSEQR_LIBRARY
#
# == Common configuration for all but CSparse (SuiteSparse version >= 4).
# SUITESPARSE_CONFIG_FOUND
# SUITESPARSE_CONFIG_INCLUDE_DIR
# SUITESPARSE_CONFIG_LIBRARY
#
# == Common configuration for all but CSparse (SuiteSparse version < 4).
# UFCONFIG_FOUND
# UFCONFIG_INCLUDE_DIR
#
# Optional SuiteSparse Dependencies:
#
# == Serial Graph Partitioning and Fill-reducing Matrix Ordering (METIS)
# METIS_FOUND
# METIS_LIBRARY
#
# == Intel Thread Building Blocks (TBB)
# TBB_FOUND
# TBB_LIBRARIES
# Called if we failed to find SuiteSparse or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
MACRO(SUITESPARSE_REPORT_NOT_FOUND REASON_MSG)
  UNSET(SUITESPARSE_FOUND)
  UNSET(SUITESPARSE_INCLUDE_DIRS)
  UNSET(SUITESPARSE_LIBRARIES)
  UNSET(SUITESPARSE_VERSION)
  UNSET(SUITESPARSE_MAIN_VERSION)
  UNSET(SUITESPARSE_SUB_VERSION)
  UNSET(SUITESPARSE_SUBSUB_VERSION)
  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  IF (SuiteSparse_FIND_QUIETLY)
    MESSAGE(STATUS "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  ELSEIF (SuiteSparse_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  ELSE()
    # Neither QUIETLY nor REQUIRED, use WARNING which emits a message
    # but continues configuration and allows generation.
    MESSAGE(WARNING "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  ENDIF (SuiteSparse_FIND_QUIETLY)
ENDMACRO(SUITESPARSE_REPORT_NOT_FOUND)
# Specify search directories for include files and libraries (this is the union
# of the search directories for all OSs).  Search user-specified hint
# directories first if supplied, and search user-installed locations first
# so that we prefer user installs to system installs where both exist.
LIST(APPEND SUITESPARSE_CHECK_INCLUDE_DIRS
  ${SUITESPARSE_INCLUDE_DIR_HINTS}
  /opt/local/include
  /opt/local/include/ufsparse # Mac OS X
  /usr/local/homebrew/include # Mac OS X
  /usr/local/include
  /usr/local/include/suitesparse
  /usr/include/suitesparse # Ubuntu
  /usr/include)
LIST(APPEND SUITESPARSE_CHECK_LIBRARY_DIRS
  ${SUITESPARSE_LIBRARY_DIR_HINTS}
  /opt/local/lib
  /opt/local/lib/ufsparse # Mac OS X
  /usr/local/homebrew/lib # Mac OS X
  /usr/local/lib
  /usr/local/lib/suitesparse
  /usr/lib/suitesparse # Ubuntu
  /usr/lib)
# BLAS.
FIND_PACKAGE(BLAS QUIET)
IF (NOT BLAS_FOUND)
  SUITESPARSE_REPORT_NOT_FOUND(
    "Did not find BLAS library (required for SuiteSparse).")
ENDIF (NOT BLAS_FOUND)
# LAPACK.
FIND_PACKAGE(LAPACK QUIET)
IF (NOT LAPACK_FOUND)
  SUITESPARSE_REPORT_NOT_FOUND(
    "Did not find LAPACK library (required for SuiteSparse).")
ENDIF (NOT LAPACK_FOUND)
# AMD.
SET(AMD_FOUND TRUE)
FIND_LIBRARY(AMD_LIBRARY NAMES amd
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${AMD_LIBRARY})
  MESSAGE(STATUS "Found AMD library: ${AMD_LIBRARY}")
ELSE (EXISTS ${AMD_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find AMD library.")
  SET(AMD_FOUND FALSE)
ENDIF (EXISTS ${AMD_LIBRARY})
MARK_AS_ADVANCED(AMD_LIBRARY)
FIND_PATH(AMD_INCLUDE_DIR NAMES amd.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${AMD_INCLUDE_DIR})
  MESSAGE(STATUS "Found AMD header in: ${AMD_INCLUDE_DIR}")
ELSE (EXISTS ${AMD_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find AMD header.")
  SET(AMD_FOUND FALSE)
ENDIF (EXISTS ${AMD_INCLUDE_DIR})
MARK_AS_ADVANCED(AMD_INCLUDE_DIR)
# CAMD.
SET(CAMD_FOUND TRUE)
FIND_LIBRARY(CAMD_LIBRARY NAMES camd
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${CAMD_LIBRARY})
  MESSAGE(STATUS "Found CAMD library: ${CAMD_LIBRARY}")
ELSE (EXISTS ${CAMD_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CAMD library.")
  SET(CAMD_FOUND FALSE)
ENDIF (EXISTS ${CAMD_LIBRARY})
MARK_AS_ADVANCED(CAMD_LIBRARY)
FIND_PATH(CAMD_INCLUDE_DIR NAMES camd.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${CAMD_INCLUDE_DIR})
  MESSAGE(STATUS "Found CAMD header in: ${CAMD_INCLUDE_DIR}")
ELSE (EXISTS ${CAMD_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CAMD header.")
  SET(CAMD_FOUND FALSE)
ENDIF (EXISTS ${CAMD_INCLUDE_DIR})
MARK_AS_ADVANCED(CAMD_INCLUDE_DIR)
# COLAMD.
SET(COLAMD_FOUND TRUE)
FIND_LIBRARY(COLAMD_LIBRARY NAMES colamd
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${COLAMD_LIBRARY})
  MESSAGE(STATUS "Found COLAMD library: ${COLAMD_LIBRARY}")
ELSE (EXISTS ${COLAMD_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find COLAMD library.")
  SET(COLAMD_FOUND FALSE)
ENDIF (EXISTS ${COLAMD_LIBRARY})
MARK_AS_ADVANCED(COLAMD_LIBRARY)
FIND_PATH(COLAMD_INCLUDE_DIR NAMES colamd.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${COLAMD_INCLUDE_DIR})
  MESSAGE(STATUS "Found COLAMD header in: ${COLAMD_INCLUDE_DIR}")
ELSE (EXISTS ${COLAMD_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find COLAMD header.")
  SET(COLAMD_FOUND FALSE)
ENDIF (EXISTS ${COLAMD_INCLUDE_DIR})
MARK_AS_ADVANCED(COLAMD_INCLUDE_DIR)
# CCOLAMD.
SET(CCOLAMD_FOUND TRUE)
FIND_LIBRARY(CCOLAMD_LIBRARY NAMES ccolamd
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${CCOLAMD_LIBRARY})
  MESSAGE(STATUS "Found CCOLAMD library: ${CCOLAMD_LIBRARY}")
ELSE (EXISTS ${CCOLAMD_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CCOLAMD library.")
  SET(CCOLAMD_FOUND FALSE)
ENDIF (EXISTS ${CCOLAMD_LIBRARY})
MARK_AS_ADVANCED(CCOLAMD_LIBRARY)
FIND_PATH(CCOLAMD_INCLUDE_DIR NAMES ccolamd.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${CCOLAMD_INCLUDE_DIR})
  MESSAGE(STATUS "Found CCOLAMD header in: ${CCOLAMD_INCLUDE_DIR}")
ELSE (EXISTS ${CCOLAMD_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CCOLAMD header.")
  SET(CCOLAMD_FOUND FALSE)
ENDIF (EXISTS ${CCOLAMD_INCLUDE_DIR})
MARK_AS_ADVANCED(CCOLAMD_INCLUDE_DIR)
# CHOLMOD.
SET(CHOLMOD_FOUND TRUE)
FIND_LIBRARY(CHOLMOD_LIBRARY NAMES cholmod
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${CHOLMOD_LIBRARY})
  MESSAGE(STATUS "Found CHOLMOD library: ${CHOLMOD_LIBRARY}")
ELSE (EXISTS ${CHOLMOD_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CHOLMOD library.")
  SET(CHOLMOD_FOUND FALSE)
ENDIF (EXISTS ${CHOLMOD_LIBRARY})
MARK_AS_ADVANCED(CHOLMOD_LIBRARY)
FIND_PATH(CHOLMOD_INCLUDE_DIR NAMES cholmod.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${CHOLMOD_INCLUDE_DIR})
  MESSAGE(STATUS "Found CHOLMOD header in: ${CHOLMOD_INCLUDE_DIR}")
ELSE (EXISTS ${CHOLMOD_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find CHOLMOD header.")
  SET(CHOLMOD_FOUND FALSE)
ENDIF (EXISTS ${CHOLMOD_INCLUDE_DIR})
MARK_AS_ADVANCED(CHOLMOD_INCLUDE_DIR)
# KLU.
SET(KLU_FOUND TRUE)
FIND_LIBRARY(KLU_LIBRARY NAMES klu
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${KLU_LIBRARY})
  MESSAGE(STATUS "Found KLU library: ${KLU_LIBRARY}")
ELSE (EXISTS ${KLU_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find KLU library.")
  SET(KLU_FOUND FALSE)
ENDIF (EXISTS ${KLU_LIBRARY})
MARK_AS_ADVANCED(KLU_LIBRARY)
FIND_PATH(KLU_INCLUDE_DIR NAMES klu.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${KLU_INCLUDE_DIR})
  MESSAGE(STATUS "Found KLU header in: ${KLU_INCLUDE_DIR}")
ELSE (EXISTS ${KLU_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find KLU header.")
  SET(KLU_FOUND FALSE)
ENDIF (EXISTS ${KLU_INCLUDE_DIR})
MARK_AS_ADVANCED(KLU_INCLUDE_DIR)
# UMFPACK.
SET(UMFPACK_FOUND TRUE)
FIND_LIBRARY(UMFPACK_LIBRARY NAMES umfpack
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${UMFPACK_LIBRARY})
  MESSAGE(STATUS "Found UMFPACK library: ${UMFPACK_LIBRARY}")
ELSE (EXISTS ${UMFPACK_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find UMFPACK library.")
  SET(UMFPACK_FOUND FALSE)
ENDIF (EXISTS ${UMFPACK_LIBRARY})
MARK_AS_ADVANCED(UMFPACK_LIBRARY)
FIND_PATH(UMFPACK_INCLUDE_DIR NAMES umfpack.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${UMFPACK_INCLUDE_DIR})
  MESSAGE(STATUS "Found UMFPACK header in: ${UMFPACK_INCLUDE_DIR}")
ELSE (EXISTS ${UMFPACK_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find UMFPACK header.")
  SET(UMFPACK_FOUND FALSE)
ENDIF (EXISTS ${UMFPACK_INCLUDE_DIR})
MARK_AS_ADVANCED(UMFPACK_INCLUDE_DIR)
# SuiteSparseQR.
SET(SUITESPARSEQR_FOUND TRUE)
FIND_LIBRARY(SUITESPARSEQR_LIBRARY NAMES spqr
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${SUITESPARSEQR_LIBRARY})
  MESSAGE(STATUS "Found SuiteSparseQR library: ${SUITESPARSEQR_LIBRARY}")
ELSE (EXISTS ${SUITESPARSEQR_LIBRARY})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find SUITESPARSEQR library.")
  SET(SUITESPARSEQR_FOUND FALSE)
ENDIF (EXISTS ${SUITESPARSEQR_LIBRARY})
MARK_AS_ADVANCED(SUITESPARSEQR_LIBRARY)
FIND_PATH(SUITESPARSEQR_INCLUDE_DIR NAMES SuiteSparseQR.hpp
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
  MESSAGE(STATUS "Found SuiteSparseQR header in: ${SUITESPARSEQR_INCLUDE_DIR}")
ELSE (EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
  SUITESPARSE_REPORT_NOT_FOUND("Did not find SUITESPARSEQR header.")
  SET(SUITESPARSEQR_FOUND FALSE)
ENDIF (EXISTS ${SUITESPARSEQR_INCLUDE_DIR})
MARK_AS_ADVANCED(SUITESPARSEQR_INCLUDE_DIR)
IF (SUITESPARSEQR_FOUND)
  # SuiteSparseQR may be compiled with Intel Threading Building Blocks,
  # we assume that if TBB is installed, SuiteSparseQR was compiled with
  # support for it, this will do no harm if it wasn't.
  SET(TBB_FOUND TRUE)
  FIND_LIBRARY(TBB_LIBRARIES NAMES tbb
    PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
  IF (EXISTS ${TBB_LIBRARIES})
    MESSAGE(STATUS "Found Intel Thread Building Blocks (TBB) library: "
      "${TBB_LIBRARIES}, assuming SuiteSparseQR was compiled with TBB.")
  ELSE (EXISTS ${TBB_LIBRARIES})
    MESSAGE(STATUS "Did not find Intel TBB library, assuming SuiteSparseQR was "
      "not compiled with TBB.")
    SET(TBB_FOUND FALSE)
  ENDIF (EXISTS ${TBB_LIBRARIES})
  MARK_AS_ADVANCED(TBB_LIBRARIES)
  IF (TBB_FOUND)
    FIND_LIBRARY(TBB_MALLOC_LIB NAMES tbbmalloc
      PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
    IF (EXISTS ${TBB_MALLOC_LIB})
      MESSAGE(STATUS "Found Intel Thread Building Blocks (TBB) Malloc library: "
        "${TBB_MALLOC_LIB}")
      # Append TBB malloc library to TBB libraries list whilst retaining
      # any CMake generated help string (cache variable).
      LIST(APPEND TBB_LIBRARIES ${TBB_MALLOC_LIB})
      GET_PROPERTY(HELP_STRING CACHE TBB_LIBRARIES PROPERTY HELPSTRING)
      SET(TBB_LIBRARIES "${TBB_LIBRARIES}" CACHE STRING ${HELP_STRING})
      # Add the TBB libraries to the SuiteSparseQR libraries (the only
      # libraries to optionally depend on TBB).
      LIST(APPEND SUITESPARSEQR_LIBRARY ${TBB_LIBRARIES})
    ELSE (EXISTS ${TBB_MALLOC_LIB})
      # If we cannot find all required TBB components do not include it as
      # a dependency.
      MESSAGE(STATUS "Did not find Intel Thread Building Blocks (TBB) Malloc "
        "Library, discarding TBB as a dependency.")
      SET(TBB_FOUND FALSE)
    ENDIF (EXISTS ${TBB_MALLOC_LIB})
    MARK_AS_ADVANCED(TBB_MALLOC_LIB)
  ENDIF (TBB_FOUND)
ENDIF(SUITESPARSEQR_FOUND)
# UFconfig / SuiteSparse_config.
#
# If SuiteSparse version is >= 4 then SuiteSparse_config is required.
# For SuiteSparse 3, UFconfig.h is required.
SET(SUITESPARSE_CONFIG_FOUND TRUE)
SET(UFCONFIG_FOUND TRUE)
FIND_LIBRARY(SUITESPARSE_CONFIG_LIBRARY NAMES suitesparseconfig
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
  MESSAGE(STATUS "Found SuiteSparse_config library: "
    "${SUITESPARSE_CONFIG_LIBRARY}")
ELSE (EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
  MESSAGE(STATUS "Did not find SuiteSparse_config library (should be present "
    "in SuiteSparse >= v4 installs).")
ENDIF (EXISTS ${SUITESPARSE_CONFIG_LIBRARY})
MARK_AS_ADVANCED(SUITESPARSE_CONFIG_LIBRARY)
FIND_PATH(SUITESPARSE_CONFIG_INCLUDE_DIR NAMES SuiteSparse_config.h
  PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
IF (EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
  MESSAGE(STATUS "Found SuiteSparse_config header in: "
    "${SUITESPARSE_CONFIG_INCLUDE_DIR}")
  SET(UFCONFIG_FOUND FALSE)
ELSE (EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
  MESSAGE(STATUS "Did not find SuiteSparse_config header (should be present "
    "in SuiteSparse >= v4 installs).")
ENDIF (EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
MARK_AS_ADVANCED(SUITESPARSE_CONFIG_INCLUDE_DIR)
IF (EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND
    EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
  # SuiteSparse_config (SuiteSparse version >= 4) requires librt library for
  # timing by default when compiled on Linux or Unix, but not on OSX (which
  # does not have librt).
  IF (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
    FIND_LIBRARY(LIBRT_LIBRARY NAMES rt
      PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
    IF (LIBRT_LIBRARY)
      MESSAGE(STATUS "Adding librt: ${LIBRT_LIBRARY} to "
        "SuiteSparse_config libraries (required on Linux & Unix [not OSX] if "
        "SuiteSparse is compiled with timing).")
    ELSE (LIBRT_LIBRARY)
      MESSAGE(STATUS "Could not find librt, but found SuiteSparse_config, "
        "assuming that SuiteSparse was compiled without timing.")
    ENDIF (LIBRT_LIBRARY)
    MARK_AS_ADVANCED(LIBRT_LIBRARY)
    LIST(APPEND SUITESPARSE_CONFIG_LIBRARY ${LIBRT_LIBRARY})
  ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux" OR UNIX AND NOT APPLE)
ELSE (EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND
      EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
  SET(SUITESPARSE_CONFIG_FOUND FALSE)
  FIND_PATH(UFCONFIG_INCLUDE_DIR NAMES UFconfig.h
    PATHS ${SUITESPARSE_CHECK_INCLUDE_DIRS})
  IF (EXISTS ${UFCONFIG_INCLUDE_DIR})
    MESSAGE(STATUS "Found UFconfig header in: ${UFCONFIG_INCLUDE_DIR}")
  ELSE (EXISTS ${UFCONFIG_INCLUDE_DIR})
    MESSAGE(STATUS "Did not find UFconfig header (should be present "
      "in SuiteSparse < v4 installs)")
    SET(UFCONFIG_FOUND FALSE)
  ENDIF (EXISTS ${UFCONFIG_INCLUDE_DIR})
  MARK_AS_ADVANCED(UFCONFIG_INCLUDE_DIR)
ENDIF (EXISTS ${SUITESPARSE_CONFIG_LIBRARY} AND
       EXISTS ${SUITESPARSE_CONFIG_INCLUDE_DIR})
# Extract the SuiteSparse version from the appropriate header (UFconfig.h for
# <= v3, SuiteSparse_config.h for >= v4).
IF (UFCONFIG_FOUND)
  # SuiteSparse version <= 3.
  SET(SUITESPARSE_VERSION_FILE ${UFCONFIG_INCLUDE_DIR}/UFconfig.h)
  IF (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
    SUITESPARSE_REPORT_NOT_FOUND(
      "Could not find file: ${SUITESPARSE_VERSION_FILE} containing version "
      "information for <= v3 SuiteSparse installs, but UFconfig was found "
      "(only present in <= v3 installs).")
  ELSE (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
    FILE(READ ${SUITESPARSE_VERSION_FILE} UFCONFIG_CONTENTS)
    STRING(REGEX MATCH "#define SUITESPARSE_MAIN_VERSION [0-9]+"
      SUITESPARSE_MAIN_VERSION "${UFCONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_MAIN_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_MAIN_VERSION "${SUITESPARSE_MAIN_VERSION}")
    STRING(REGEX MATCH "#define SUITESPARSE_SUB_VERSION [0-9]+"
      SUITESPARSE_SUB_VERSION "${UFCONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_SUB_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_SUB_VERSION "${SUITESPARSE_SUB_VERSION}")
    STRING(REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION [0-9]+"
      SUITESPARSE_SUBSUB_VERSION "${UFCONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_SUBSUB_VERSION "${SUITESPARSE_SUBSUB_VERSION}")
    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 4.;2.;1 nonsense.
    SET(SUITESPARSE_VERSION
      "${SUITESPARSE_MAIN_VERSION}.${SUITESPARSE_SUB_VERSION}.${SUITESPARSE_SUBSUB_VERSION}")
  ENDIF (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
ENDIF (UFCONFIG_FOUND)
IF (SUITESPARSE_CONFIG_FOUND)
  # SuiteSparse version >= 4.
  SET(SUITESPARSE_VERSION_FILE
    ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h)
  IF (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
    SUITESPARSE_REPORT_NOT_FOUND(
      "Could not find file: ${SUITESPARSE_VERSION_FILE} containing version "
      "information for >= v4 SuiteSparse installs, but SuiteSparse_config was "
      "found (only present in >= v4 installs).")
  ELSE (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
    FILE(READ ${SUITESPARSE_VERSION_FILE} SUITESPARSE_CONFIG_CONTENTS)
    STRING(REGEX MATCH "#define SUITESPARSE_MAIN_VERSION [0-9]+"
      SUITESPARSE_MAIN_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_MAIN_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_MAIN_VERSION "${SUITESPARSE_MAIN_VERSION}")
    STRING(REGEX MATCH "#define SUITESPARSE_SUB_VERSION [0-9]+"
      SUITESPARSE_SUB_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_SUB_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_SUB_VERSION "${SUITESPARSE_SUB_VERSION}")
    STRING(REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION [0-9]+"
      SUITESPARSE_SUBSUB_VERSION "${SUITESPARSE_CONFIG_CONTENTS}")
    STRING(REGEX REPLACE "#define SUITESPARSE_SUBSUB_VERSION ([0-9]+)" "\\1"
      SUITESPARSE_SUBSUB_VERSION "${SUITESPARSE_SUBSUB_VERSION}")
    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 4.;2.;1 nonsense.
    SET(SUITESPARSE_VERSION
      "${SUITESPARSE_MAIN_VERSION}.${SUITESPARSE_SUB_VERSION}.${SUITESPARSE_SUBSUB_VERSION}")
  ENDIF (NOT EXISTS ${SUITESPARSE_VERSION_FILE})
ENDIF (SUITESPARSE_CONFIG_FOUND)
# METIS (Optional dependency).
FIND_LIBRARY(METIS_LIBRARY NAMES libmetis.so.5
  PATHS ${SUITESPARSE_CHECK_LIBRARY_DIRS})
IF (EXISTS ${METIS_LIBRARY})
  MESSAGE(STATUS "Found METIS library: ${METIS_LIBRARY}.")
ELSE (EXISTS ${METIS_LIBRARY})
  MESSAGE(STATUS "Did not find METIS library (optional SuiteSparse dependency)")
ENDIF (EXISTS ${METIS_LIBRARY})
MARK_AS_ADVANCED(METIS_LIBRARY)
# Only mark SuiteSparse as found if all required dependencies have been found.
SET(SUITESPARSE_FOUND FALSE)
IF (AMD_FOUND AND
    CAMD_FOUND AND
    COLAMD_FOUND AND
    CCOLAMD_FOUND AND
    CHOLMOD_FOUND AND
    UMFPACK_FOUND AND
    SUITESPARSEQR_FOUND AND
    (SUITESPARSE_CONFIG_FOUND OR UFCONFIG_FOUND) AND
    BLAS_FOUND AND
    LAPACK_FOUND)
  SET(SUITESPARSE_FOUND TRUE)
  LIST(APPEND SUITESPARSE_INCLUDE_DIRS
    ${AMD_INCLUDE_DIR}
    ${CAMD_INCLUDE_DIR}
    ${COLAMD_INCLUDE_DIR}
    ${CCOLAMD_INCLUDE_DIR}
    ${CHOLMOD_INCLUDE_DIR}
    ${UMFPACK_INCLUDE_DIR}
    ${SUITESPARSEQR_INCLUDE_DIR})
  # Handle config separately, as otherwise at least one of them will be set
  # to NOTFOUND which would cause any check on SUITESPARSE_INCLUDE_DIRS to fail.
  IF (SUITESPARSE_CONFIG_FOUND)
    LIST(APPEND SUITESPARSE_INCLUDE_DIRS
      ${SUITESPARSE_CONFIG_INCLUDE_DIR})
  ENDIF (SUITESPARSE_CONFIG_FOUND)
  IF (UFCONFIG_FOUND)
    LIST(APPEND SUITESPARSE_INCLUDE_DIRS
      ${UFCONFIG_INCLUDE_DIR})
  ENDIF (UFCONFIG_FOUND)
  # As SuiteSparse includes are often all in the same directory, remove any
  # repetitions.
  LIST(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)
  # Important: The ordering of these libraries is *NOT* arbitrary, as these
  # could potentially be static libraries their link ordering is important.
  LIST(APPEND SUITESPARSE_LIBRARIES
    ${SUITESPARSEQR_LIBRARY}
    ${UMFPACK_LIBRARY}
    ${CHOLMOD_LIBRARY}
    ${CCOLAMD_LIBRARY}
    ${CAMD_LIBRARY}
    ${COLAMD_LIBRARY}
    ${AMD_LIBRARY})
  IF (SUITESPARSE_CONFIG_FOUND)
    LIST(APPEND SUITESPARSE_LIBRARIES
      ${SUITESPARSE_CONFIG_LIBRARY})
  ENDIF (SUITESPARSE_CONFIG_FOUND)
  IF (METIS_FOUND)
    LIST(APPEND SUITESPARSE_LIBRARIES
      ${METIS_LIBRARY})
  ENDIF (METIS_FOUND)
  IF (KLU_FOUND)
    LIST(APPEND SUITESPARSE_LIBRARIES
      ${KLU_LIBRARY})
  ENDIF (KLU_FOUND)
ELSE()
  SUITESPARSE_REPORT_NOT_FOUND(
    "Failed to find some/all required components of SuiteSparse.")
ENDIF()
# Determine if we are running on Ubuntu with the package install of SuiteSparse
# which is broken and does not support linking a shared library.
SET(SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION FALSE)
IF (CMAKE_SYSTEM_NAME MATCHES "Linux" AND
    SUITESPARSE_VERSION VERSION_EQUAL 3.4.0)
  FIND_PROGRAM(LSB_RELEASE_EXECUTABLE lsb_release)
  IF (LSB_RELEASE_EXECUTABLE)
    # Any even moderately recent Ubuntu release (likely to be affected by
    # this bug) should have lsb_release, if it isn't present we are likely
    # on a different Linux distribution (should be fine).
    EXECUTE_PROCESS(COMMAND ${LSB_RELEASE_EXECUTABLE} -si
      OUTPUT_VARIABLE LSB_DISTRIBUTOR_ID
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    IF (LSB_DISTRIBUTOR_ID MATCHES "Ubuntu" AND
        SUITESPARSE_LIBRARIES MATCHES "/usr/lib/libamd")
      # We are on Ubuntu, and the SuiteSparse version matches the broken
      # system install version and is a system install.
      SET(SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION TRUE)
      MESSAGE(STATUS "Found system install of SuiteSparse "
        "${SUITESPARSE_VERSION} running on Ubuntu, which has a known bug "
        "preventing linking of shared libraries (static linking unaffected).")
    ENDIF (LSB_DISTRIBUTOR_ID MATCHES "Ubuntu" AND
      SUITESPARSE_LIBRARIES MATCHES "/usr/lib/libamd")
  ENDIF (LSB_RELEASE_EXECUTABLE)
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux" AND
  SUITESPARSE_VERSION VERSION_EQUAL 3.4.0)
# Handle REQUIRED and QUIET arguments to FIND_PACKAGE
INCLUDE(FindPackageHandleStandardArgs)
# A change to CMake after release 2.8.10.2 means that
# FindPackageHandleStandardArgs() unsets <LibraryName>_FOUND without checking
# if it is one of the variables passed whose existence & validity is verified
# by FindPackageHandleStandardArgs() in conjunction with handling the REQUIRED
# and QUIET optional arguments, as such we use an intermediary variable.
SET(SUITESPARSE_FOUND_COPY ${SUITESPARSE_FOUND})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SuiteSparse
  REQUIRED_VARS SUITESPARSE_INCLUDE_DIRS SUITESPARSE_LIBRARIES SUITESPARSE_FOUND_COPY
  VERSION_VAR SUITESPARSE_VERSION)
