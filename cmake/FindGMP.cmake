# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
#[=======================================================================[.rst:
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
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

# Try to help find_package process (pkg-config ...)
set_find_package_hints(NAME GMP MODULE gmp)


find_path(GMP_INCLUDE_DIR NAMES gmp.h
  PATH_SUFFIXES include
  ${_GMP_INC_SEARCH_OPTS}
  )

if(NOT GMP_LIBRARIES)
if(WIN32)
    set(GMPLIB mpir)
else()
    set(GMPLIB gmp)
endif()
  find_library(GMP_LIBRARIES NAMES ${GMPLIB}
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
