#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2018 INRIA.
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
# Find MUMPS libraries and headers
#
# Usage :
# 
# find_package(MUMPS REQUIRED)
# target_link_libraries(yourlib PRIVATE MUMPS::MUMPS)
#
# It will handles both linking and include directories for your target.
#
# This module sets the following variables:
#
# MUMPS_FOUND - set to true if a library implementing the MUMPS interface
#    is found
# MUMPS_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use MUMPS
# MUMPS_INCLUDE_DIR - location of mumps headers found by cmake

# Set MUMPS_DIR=<where mumps is installed> if it's not in a "classic" place or if you want a specific version

include(FindPackageHandleStandardArgs)

# On debian and Suse, both mpi and mpi-free version of MUMPS can be installed in parallel
# The latter have a "_seq" suffix
if(WITH_MPI)
  set(__MUMPS_NAMES dmumps_ptscotch dmumps_scotch dmumps)
else()
  set(__MUMPS_NAMES dmumps_seq dmumps)
  set(__SUFFIX seq)
  # on debian systems one may have mumps+[pt]scotch packages also
endif()


# ---- Search for mumps header(s) ---

# -- First try : use MUMPS_DIR provided explicitely at cmake call
if(MUMPS_DIR)
  find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
    PATHS ${MUMPS_DIR}
    PATH_SUFFIXES include
    PATH_SUFFIXES MUMPS mumps mumps${__SUFFIX}
    )
endif()

# -- Second try : pkg-config
# 
find_package(PkgConfig)
pkg_check_modules(PKGC_MUMPS mumps QUIET)
if(PKGC_MUMPS_FOUND)
  set(MUMPS_LIBRARIES "${PKGC_MUMPS_LIBRARIES}")
endif()

find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
  PATHS ${PKGC_MUMPS_INCLUDE_DIRS})

# -- Last try : default behavior ...
find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h)

# --- Search for mumps libraries ---

if(NOT MUMPS_LIBRARY)
  if(MUMPS_DIR)
    find_library(MUMPS_LIBRARY NAMES ${__MUMPS_NAMES}
    PATHS ${MUMPS_DIR} ENV LIBRARY_PATH ENV LD_LIBRARY_PATH
    PATH_SUFFIXES lib lib64 
    )
    
  endif()
  find_library(MUMPS_LIBRARY NAMES ${__MUMPS_NAMES})
  
endif()

# -- Search for extra libraries, required to link with mumps --
# Try to be smart and detect whether the MUMPS lib file has a "_seq" suffix. If yes, we add it to all the other library
if(MUMPS_LIBRARY MATCHES "dmumps_seq")
  set(__MUMPS_COMMON_NAMES mumps_common_seq mumps_common)
else()
  set(__MUMPS_COMMON_NAMES mumps_common_ptscotch mumps_common_scotch mumps_common mumps_common_seq)
endif()

if(MUMPS_DIR)
  find_library(MUMPS_COMMON_LIBRARY NAMES ${__MUMPS_COMMON_NAMES}
    PATHS ${MUMPS_DIR} ENV LIBRARY_PATH ENV LD_LIBRARY_PATH
    PATH_SUFFIXES lib lib64 
    )
  
endif()
find_library(MUMPS_COMMON_LIBRARY NAMES ${__MUMPS_COMMON_NAMES})

set(MUMPS_LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})

# Extras
set(extras_libs ptscotch scotch metis pord${__SUFFIX})
foreach(extra IN LISTS extras_libs)
  if(MUMPS_DIR)
    find_library(${extra}_LIBRARY NAMES ${extra}
      PATHS ${MUMPS_DIR} ENV LIBRARY_PATH ENV LD_LIBRARY_PATH
      PATH_SUFFIXES lib lib64 
      )
    
  endif()
  find_library(${extra}_LIBRARY NAMES ${extra})
  if(${extra}_LIBRARY)
    list(APPEND ${MUMPS_LIBRARIES} ${extra}_LIBRARY)
  endif()
    
endforeach()

# -- Library setup --
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)

if(MUMPS_FOUND)
  
  if(NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS IMPORTED INTERFACE)
    set_property(TARGET MUMPS::MUMPS PROPERTY INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARIES}")
    if(MUMPS_INCLUDE_DIR)
      set_target_properties(MUMPS::MUMPS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
    endif()
    if(EXISTS "${MUMPS_LIBRARIES}")
      set_target_properties(MUMPS::MUMPS PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${MUMPS_LIBRARIES}")
    endif()
  endif()
endif()
mark_as_advanced(MUMPS_LIBRARIES MUMPS_INCLUDE_DIR)

