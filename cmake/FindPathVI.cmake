#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
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
FindPathVI
-----------

Find Path library

Usage :
 
find_package(PathVI REQUIRED)
target_link_libraries(yourlib PRIVATE PathVI::PathVI)

Set PathVI_ROOT=<where cppunit is installed>
if it's not in a "classic" place or if you want a specific version

Header : vi_desc.h
Lib : <prefix>pathvi<PathVI_VERSION>.<suffix>

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT PathVI_ROOT)
    set(PathVI_ROOT $ENV{PathVI_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME PathVI MODULE pathvi)

if(NOT PathVI_LIBRARIES)
  find_library(PathVI_LIBRARIES NAMES pathvi${PathVI_VERSION}
    ${_PathVI_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
endif()

# -- Library setup --
find_package_handle_standard_args(PathVI
  REQUIRED_VARS PathVI_LIBRARIES)

if(PathVI_FOUND)
  
  if(NOT TARGET PathVI::PathVI)
    add_library(PathVI::PathVI IMPORTED INTERFACE)
    set_property(TARGET PathVI::PathVI PROPERTY INTERFACE_LINK_LIBRARIES ${PathVI_LIBRARIES})
    if(PathVI_INCLUDE_DIR)  # FP: should we set this with externals/PATH_SDK ??
      set_target_properties(PathVI::PathVI PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PathVI_INCLUDE_DIR}")
    endif()
  endif()
endif()

