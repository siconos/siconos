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
FindCPLEX
-----------

Find CPLEX libraries and headers

See https://www.ibm.com/analytics/cplex-optimizer

Usage :
 
find_package(CPLEX REQUIRED)
target_link_libraries(yourlib PRIVATE CPLEX::CPLEX)

Set CPLEX_ROOT=<where CPLEX is installed>
if it's not in a "classic" place or if you want a specific version

lib : <prefix>cplex.<suffix>

# XXX we do not look for include now ... --xhub

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT CPLEX_ROOT)
    set(CPLEX_ROOT $ENV{CPLEX_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME CPLEX MODULE cplex)

find_path(CPLEX_INCLUDE_DIR NAMES CPLEX.h
  PATH_SUFFIXES include
  ${_CPLEX_INC_SEARCH_OPTS}
  )

if(NOT CPLEX_LIBRARIES)
  find_library(CPLEX_LIBRARIES NAMES CPLEX
    ${_CPLEX_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  
endif()

# -- Library setup --
find_package_handle_standard_args(CPLEX
  REQUIRED_VARS CPLEX_LIBRARIES CPLEX_INCLUDE_DIR)

if(CPLEX_FOUND)
  
  if(NOT TARGET CPLEX::CPLEX)
    add_library(CPLEX::CPLEX IMPORTED INTERFACE)
    set_property(TARGET CPLEX::CPLEX PROPERTY INTERFACE_LINK_LIBRARIES ${CPLEX_LIBRARIES})
    if(CPLEX_INCLUDE_DIR)
      set_target_properties(CPLEX::CPLEX PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIR}")
    endif()
  endif()
endif()

