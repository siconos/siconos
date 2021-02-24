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
FindCPPUNIT
-----------

Find cppunit libraries and headers

Usage :
 
find_package(CPPUNIT REQUIRED)
target_link_libraries(yourlib PRIVATE CPPUNIT::CPPUNIT)

Set CPPUNIT_ROOT=<where cppunit is installed>
if it's not in a "classic" place or if you want a specific version

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT CPPUNIT_ROOT)
    set(CPPUNIT_ROOT $ENV{CPPUNIT_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME CPPUNIT MODULE cppunit)

find_path(CPPUNIT_INCLUDE_DIR NAMES TestCase.h
  PATH_SUFFIXES include cppunit
  ${_CPPUNIT_INC_SEARCH_OPTS}
  )

if(NOT CPPUNIT_LIBRARIES)
  if(WIN32)
    find_library(CPPUNIT_LIBRARIES NAMES cppunit_dll
      ${_CPPUNIT_SEARCH_OPTS} 
      PATH_SUFFIXES lib lib64)
  else()
    find_library(CPPUNIT_LIBRARIES NAMES cppunit
      ${_CPPUNIT_SEARCH_OPTS}
      PATH_SUFFIXES lib lib64)
  endif()
 
endif()

# -- Library setup --
find_package_handle_standard_args(CPPUNIT
  REQUIRED_VARS CPPUNIT_LIBRARIES CPPUNIT_INCLUDE_DIR)

if(CPPUNIT_FOUND)
  
  if(NOT TARGET CPPUNIT::CPPUNIT)
    add_library(CPPUNIT::CPPUNIT IMPORTED INTERFACE)
    set_property(TARGET CPPUNIT::CPPUNIT PROPERTY INTERFACE_LINK_LIBRARIES ${CPPUNIT_LIBRARIES})
    if(CPPUNIT_INCLUDE_DIR)
      set_target_properties(CPPUNIT::CPPUNIT PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CPPUNIT_INCLUDE_DIR}")
    endif()
  endif()
endif()
