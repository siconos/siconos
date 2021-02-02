#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2020 INRIA.
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
FindLPSOLVE
-----------

Find lpsolve libraries and headers

Usage :
 
find_package(LPSOLVE REQUIRED)
target_link_libraries(yourlib PRIVATE LPSOLVE::LPSOLVE)

Set LPSOLVE_ROOT=<where lpsolve is installed>
if it's not in a "classic" place or if you want a specific version

Header : lp_lib.h
Lib : <prefix>lpsolve55.<suffix>

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT LPSOLVE_ROOT)
  set(LPSOLVE_ROOT $ENV{LPSOLVE_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME LPSOLVE MODULE lpsolve)

find_path(LPSOLVE_INCLUDE_DIR NAMES lp_lib.h
  PATH_SUFFIXES lpsolve include
  ${_LPSOLVE_INC_SEARCH_OPTS}
  )

if(NOT LPSOLVE_LIBRARIES)
  # Fix debian nonsense:
  # see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=503314
  if(EXISTS "/etc/debian_version")
    find_library(LPSOLVE_LIBRARIES NAMES liblpsolve55.so
      ${_LPSOLVE_SEARCH_OPTS}
      PATH_SUFFIXES lib/lp_solve lib)
  else()
    find_library(LPSOLVE_LIBRARIES NAMES lpsolve55
      ${_LPSOLVE_SEARCH_OPTS}
      PATH_SUFFIXES lib lib64)
  endif()
endif()



# -- Library setup --
find_package_handle_standard_args(LPSOLVE
  REQUIRED_VARS LPSOLVE_LIBRARIES LPSOLVE_INCLUDE_DIR)

if(LPSOLVE_FOUND)
  
  if(NOT TARGET LPSOLVE::LPSOLVE)
    add_library(LPSOLVE::LPSOLVE IMPORTED INTERFACE)
    set_property(TARGET LPSOLVE::LPSOLVE PROPERTY INTERFACE_LINK_LIBRARIES ${LPSOLVE_LIBRARIES})
    if(LPSOLVE_INCLUDE_DIR)
      set_target_properties(LPSOLVE::LPSOLVE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LPSOLVE_INCLUDE_DIR}")
    endif()
  endif()
endif()


