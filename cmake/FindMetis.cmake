# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[.rst:
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
FindMetis
-----------

Find Metis libraries and headers

Usage :

find_package(Metis REQUIRED)
target_link_libraries(yourlib PRIVATE Metis::Metis)

Set Metis_ROOT=<where Metis is installed>
if it's not in a "classic" place or if you want a specific version

header : metis.h
lib : <prefi>metis.<suffix>

#]=======================================================================]


include(FindPackageHandleStandardArgs)

if(NOT Metis_ROOT)
  set(Metis_ROOT $ENV{Metis_ROOT})
endif()

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME Metis MODULE metis)

if(NOT Metis_INCLUDE_DIR)
  find_path(Metis_INCLUDE_DIR NAMES metis.h
    PATH_SUFFIXES include
    ${_UMFPACK_INC_SEARCH_OPTS}
    )
endif()

if(NOT Metis_LIBRARIES)
  find_library(Metis_LIBRARIES NAMES metis
    ${_Metis_SEARCH_OPTS}
    PATH_SUFFIXES lib lib64)

endif()

# -- Library setup --
find_package_handle_standard_args(Metis REQUIRED_VARS Metis_LIBRARIES Metis_INCLUDE_DIR)

if(Metis_FOUND)
  if(NOT TARGET Metis::Metis)
    add_library(Metis::Metis IMPORTED INTERFACE)
    set_property(TARGET Metis::Metis PROPERTY INTERFACE_LINK_LIBRARIES ${Metis_LIBRARIES})
    if(Metis_INCLUDE_DIR)
      set_target_properties(Metis::Metis PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Metis_INCLUDE_DIR}")
    endif()
  endif()
endif()
