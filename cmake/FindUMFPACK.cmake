# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[.rst:
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
FindUMFPACK
-----------

Find UMFPACK libraries and headers

Usage :

find_package(UMFPACK REQUIRED)
target_link_libraries(yourlib PRIVATE UMFPACK::UMFPACK)

Set UMFPACK_ROOT=<where UMFPACK is installed>
if it's not in a "classic" place or if you want a specific version

header : umfpack.h
lib : <prefi>umfpack.<suffix>

Note : umfpack lib usually requires linking to a blas library.
It is up to the user of this module to find a BLAS and link to it.

Warning : this routine tries first to find umfpack in SuiteSparse.

#]=======================================================================]


# Try suitesparse
if(NOT SuiteSparse_UMFPACK_FOUND)
  find_package(SuiteSparse COMPONENTS UMFPACK)
endif()

if(SuiteSparse_UMFPACK_FOUND)
  set(UMFPACK_FOUND TRUE CACHE INTERNAL "")
  return()
endif()

include(FindPackageHandleStandardArgs)

if(NOT UMFPACK_ROOT)
  set(UMFPACK_ROOT $ENV{UMFPACK_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME UMFPACK MODULE umfpack)

if(NOT UMFPACK_INCLUDE_DIR)
  find_path(UMFPACK_INCLUDE_DIR NAMES umfpack.h
    PATH_SUFFIXES include ufsparse
    ${_UMFPACK_INC_SEARCH_OPTS}
    )
endif()

if(NOT UMFPACK_LIBRARIES)
  find_library(UMFPACK_LIBRARIES NAMES umfpack
    ${_UMFPACK_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)

endif()

# -- Library setup --
find_package_handle_standard_args(UMFPACK REQUIRED_VARS UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIR)

if(UMFPACK_FOUND)
  if(NOT TARGET UMFPACK::UMFPACK)
    add_library(UMFPACK::UMFPACK IMPORTED INTERFACE)
    set_property(TARGET UMFPACK::UMFPACK PROPERTY INTERFACE_LINK_LIBRARIES ${UMFPACK_LIBRARIES})
    if(UMFPACK_INCLUDE_DIR)
      set_target_properties(UMFPACK::UMFPACK PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}")
    endif()
  endif()
endif()
