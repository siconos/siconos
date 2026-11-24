# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[.rst:
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
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
