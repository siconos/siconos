#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2019 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --
#[=======================================================================[.rst:
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

if(UMFPACK_ROOT)
  set(_UMFPACK_SEARCH_OPTS
    HINTS ${UMFPACK_ROOT}
    NO_DEFAULT_PATH)
else()
  # Try pkgconfig
  find_package(PkgConfig QUIET)
  pkg_check_modules(PKGC_UMFPACK umfpack QUIET)
  if(PKGC_UMFPACK_FOUND)
    set(UMFPACK_LIBRARIES "${PKGC_UMFPACK_LINK_LIBRARIES}")
  endif()
  set(_UMFPACK_SEARCH_OPTS
    HINTS ${PKGC_UMFPACK_INCLUDE_DIRS} ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
endif()

if(NOT UMFPACK_INCLUDE_DIR)
  find_path(UMFPACK_INCLUDE_DIR NAMES umfpack.h
    PATH_SUFFIXES include ufsparse
    ${_UMFPACK_SEARCH_OPTS}
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

