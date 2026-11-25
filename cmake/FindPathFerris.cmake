# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[.rst:
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
FindPathFerris
-----------

Find Path library

Usage :

find_package(PathFerris REQUIRED)
target_link_libraries(yourlib PRIVATE PathFerris::PathFerris)

Set PathFerris_ROOT=<where pathferris is installed>
if it's not in a "classic" place or if you want a specific version

Lib : <prefix>path<PathFerris_VERSION>.<suffix>

#]=======================================================================]

include(FindPackageHandleStandardArgs)

if(NOT PathFerris_ROOT)
  set(PathFerris_ROOT $ENV{PathFerris_ROOT})
endif()

if(NOT PathFerris_VERSION)
  set(PathFerris_VERSION 47)
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME PathFerris MODULE path${PathFerris_VERSION})

if(NOT PathFerris_LIBRARIES)
  find_library(PathFerris_LIBRARIES NAMES path${PathFerris_VERSION}
    ${_PathFerris_SEARCH_OPTS}
    #HINTS ENV LD_LIBRARY_PATH)
    PATH_SUFFIXES lib lib64
    )
 endif()

# -- Library setup --
find_package_handle_standard_args(PathFerris
  REQUIRED_VARS PathFerris_LIBRARIES)

if(PathFerris_FOUND)

  if(NOT TARGET PathFerris::PathFerris)
    add_library(PathFerris::PathFerris IMPORTED INTERFACE)
    set_property(TARGET PathFerris::PathFerris PROPERTY INTERFACE_LINK_LIBRARIES ${PathFerris_LIBRARIES})
    if(PathFerris_INCLUDE_DIR)  # FP: should we set this with externals/PATH_SDK ??
      set_target_properties(PathFerris::PathFerris PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PathFerris_INCLUDE_DIR}")
    endif()
  endif()
endif()
