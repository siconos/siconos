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
FindGMP
-----------

Find gmp libraries and headers

Usage :
 
find_package(GMP REQUIRED)
target_link_libraries(yourlib PRIVATE GMP::GMP)

Set GMP_ROOT=<where gmp is installed>
if it's not in a "classic" place or if you want a specific version

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT GMP_ROOT)
    set(GMP_ROOT $ENV{GMP_ROOT})
endif()

if(GMP_ROOT)
  set(_GMP_SEARCH_OPTS
    "HINTS ${GMP_ROOT} NO_DEFAULT_PATH")
else()
  # Try pkgconfig
  find_package(PkgConfig QUIET)
  pkg_check_modules(PKGC_GMP gmp QUIET)
  if(PKGC_GMP_FOUND)
    set(GMP_LIBRARIES "${PKGC_GMP_LIBRARIES}")
  endif()
  set(_GMP_SEARCH_OPTS
    "HINTS ${PKGC_GMP_INCLUDE_DIRS} ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH")
endif()

find_path(GMP_INCLUDE_DIR NAMES gmp.h
  PATH_SUFFIXES include
  ${_GMP_SEARCH_OPTS}
  )

if(NOT GMP_LIBRARIES)
  find_library(GMP_LIBRARIES NAMES gmp
    ${_GMP_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  
endif()

# -- Library setup --
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARIES GMP_INCLUDE_DIR)

if(GMP_FOUND)
  
  if(NOT TARGET GMP::GMP)
    add_library(GMP::GMP IMPORTED INTERFACE)
    set_property(TARGET GMP::GMP PROPERTY INTERFACE_LINK_LIBRARIES ${GMP_LIBRARIES})
    if(GMP_INCLUDE_DIR)
      set_target_properties(GMP::GMP PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
    endif()
    # if(EXISTS "${GMP_LIBRARIES}")
    #   set_target_properties(GMP::GMP PROPERTIES
    #     IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    #     IMPORTED_LOCATION "${GMP_LIBRARIES}")
    # endif()
  endif()
endif()

