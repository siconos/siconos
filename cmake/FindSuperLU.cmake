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
FindSuperLU
-----------

Find SuperLU libraries and headers

Usage :
 
find_package(SuperLU REQUIRED COMPONENTS <name>)
target_link_libraries(yourlib PRIVATE SuperLU::SuperLU)

with name among :
- STD (standard version of superlu),
- MT (openmp/threads interface to superlu),
- DIST (manycore heterogeous node architecture interface).

If COMPONENTS is not set, STD version is searched.
Warning : do not use several components in your build.

Check details here : https://portal.nersc.gov/project/sparse/superlu/#superlu_mt

Set SuperLU_ROOT=<where SuperLU is installed>
if it's not in a "classic" place or if you want a specific version

header : slu_ddefs.h
libs : <prefix>superlu.<suffix>

#]=======================================================================]

include(FindPackageHandleStandardArgs)

list(LENGTH SuperLU_FIND_COMPONENTS nb_comp)
if(nb_comp GREATER 1)
  message(FATAL_ERROR "Try to find several SuperLU conflicting components : ${SuperLU_FIND_COMPONENTS}.")
endif()
if(NOT SuperLU_FIND_COMPONENTS)
  set(SuperLU_FIND_COMPONENTS STD)
endif()

if(${SuperLU_FIND_COMPONENTS} STREQUAL STD)
  set(SuperLU_HEADER slu_ddefs.h)
  set(SuperLU_LIBNAME superlu)
elseif(${SuperLU_FIND_COMPONENTS} STREQUAL MT)
  set(SuperLU_HEADER slu_mt_ddefs.h)
  set(SuperLU_LIBNAME superlu_mt)
elseif(${SuperLU_FIND_COMPONENTS} STREQUAL DIST)
  set(SuperLU_HEADER superlu_ddefs.h)
  set(SuperLU_LIBNAME superlu_dist)
else()
  message(FATAL_ERROR "Unknown SuperLU component ${SuperLU_FIND_COMPONENTS}. Search failed.")
endif()
    
if(NOT SuperLU_ROOT)
  set(SuperLU_ROOT $ENV{SuperLU_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME SuperLU MODULE ${SuperLU_LIBNAME})


find_path(SuperLU_INCLUDE_DIR NAMES ${SuperLU_HEADER}
  PATH_SUFFIXES include superlu
  ${_SuperLU_INC_SEARCH_OPTS}
  )

if(NOT SuperLU_LIBRARIES)
  find_library(SuperLU_LIBRARIES NAMES ${SuperLU_LIBNAME}
    ${_SuperLU_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)
  
endif()

# -- Library setup --
find_package_handle_standard_args(SuperLU
  REQUIRED_VARS SuperLU_LIBRARIES SuperLU_INCLUDE_DIR)

if(SuperLU_FOUND)
  if(NOT TARGET SuperLU::SuperLU)
    add_library(SuperLU::SuperLU IMPORTED INTERFACE)
    set_property(TARGET SuperLU::SuperLU PROPERTY INTERFACE_LINK_LIBRARIES ${SuperLU_LIBRARIES})
    if(SuperLU_INCLUDE_DIR)
      set_target_properties(SuperLU::SuperLU PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${SuperLU_INCLUDE_DIR}")
    endif()
  endif()
endif()

