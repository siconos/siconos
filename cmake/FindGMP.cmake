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
# Find GMP libraries and headers
#
# Usage :
# 
# find_package(GMP REQUIRED)
# target_link_libraries(yourlib PRIVATE GMP::GMP)
#
# It will handles both linking and include directories for your target.
#
# This module sets the following variables:
#
# GMP_FOUND - set to true if a library implementing the GMP interface
#    is found
#  GMP_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use GMP
#  GMP_INCLUDE_DIR - location of blas headers found by cmake

# Set GMP_DIR=where gmp is installed if it's not in a "classic" place or if you want a specific version
#

include(FindPackageHandleStandardArgs)

# -- First try : use GMP_DIR provided explicitely at cmake call
if(GMP_DIR)
  find_path(GMP_INCLUDE_DIR NAMES gmp.h
    PATHS ${GMP_DIR}
    PATH_SUFFIXES include
    )
endif()

# --- Second try : pkg-config
# 
find_package(PkgConfig)
pkg_check_modules(PKGC_GMP gmp QUIET)
if(PKGC_GMP_FOUND)
  set(GMP_LIBRARIES "${PKGC_GMP_LINK_LIBRARIES}")
endif()

find_path(GMP_INCLUDE_DIR NAMES gmp.h
  PATHS ${PKGC_GMP_INCLUDE_DIRS})

# --- Last try : default behavior ...
find_path(GMP_INCLUDE_DIR NAMES gmp.h)


if(NOT GMP_LIBRARIES)
  if(GMP_DIR)
    find_library(GMP_LIBRARIES NAMES gmp
    PATHS ${GMP_DIR} ENV LIBRARY_PATH ENV LD_LIBRARY_PATH
    PATH_SUFFIXES lib lib64 
    )
    
  endif()
  find_library(GMP_LIBRARIES NAMES gmp)

endif()

# -- Library setup --
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARIES GMP_INCLUDE_DIR)

if(GMP_FOUND)
  
  if(NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_property(TARGET GMP::GMP PROPERTY INTERFACE_LINK_LIBRARIES ${GMP_LIBRARIES})
    if(GMP_INCLUDE_DIR)
      set_target_properties(GMP::GMP PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
    endif()
    if(EXISTS "${GMP_LIBRARIES}")
      set_target_properties(GMP::GMP PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${GMP_LIBRARIES}")
    endif()
  endif()
endif()
mark_as_advanced(GMP_LIBRARIES GMP_INCLUDE_DIR)

