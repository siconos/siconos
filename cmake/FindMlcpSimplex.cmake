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
FindMlcpSimplex
-----------

Find the mlcp_simplex library and header.

Usage :
 
find_package(MlcpSimplex REQUIRED)
target_link_libraries(yourlib PRIVATE MlcpSimplex::MlcpSimplex)

Set MlcpSimplex_ROOT=<where mlcp simplex is installed>
if it's not in a "classic" place or if you want a specific version

Header : external_mlcp_simplex.h
Lib : <prefix>MlcpSimplex.<suffix>

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT MlcpSimplex_ROOT)
  set(MlcpSimplex_ROOT $ENV{MlcpSimplex_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME MlcpSimplex MODULE MlcpSimplex)

find_path(MlcpSimplex_INCLUDE_DIR NAMES external_mlcp_simplex.h
  PATH_SUFFIXES include
  ${_MlcpSimplex_INC_SEARCH_OPTS}
  )

if(NOT MlcpSimplex_LIBRARIES)
  find_library(MlcpSimplex_LIBRARIES NAMES MlcpSimplex
    ${_MlcpSimplex_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
endif()

# -- Library setup --
find_package_handle_standard_args(MlcpSimplex
  REQUIRED_VARS MlcpSimplex_LIBRARIES MlcpSimplex_INCLUDE_DIR)

if(MlcpSimplex_FOUND)
  if(NOT TARGET MlcpSimplex::MlcpSimplex)
    add_library(MlcpSimplex::MlcpSimplex IMPORTED INTERFACE)
    set_property(TARGET MlcpSimplex::MlcpSimplex PROPERTY INTERFACE_LINK_LIBRARIES ${MlcpSimplex_LIBRARIES})
    if(MlcpSimplex_INCLUDE_DIR)
      set_target_properties(MlcpSimplex::MlcpSimplex PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MlcpSimplex_INCLUDE_DIR}")
    endif()
  endif()
endif()
