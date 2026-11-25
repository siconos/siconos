# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
#/

# Use this file to define functions or macros that
# are not available in older version of cmake, but
# used in Siconos cmake files.
#
# e.g. : target_link_options exists from cmake 3.13
# and could be replace by a call to set_target_properties for older version.
#

# if(${CMAKE_VERSION} VERSION_LESS "3.14")
#   # https://cmake.org/cmake/help/latest/prop_gbl/CMAKE_ROLE.html
#   set_property(GLOBAL PROPERTY CMAKE_ROLE PROJECT)

# endif()


# if(${CMAKE_VERSION} VERSION_LESS "3.13")

#   function(target_link_options target prop value)
#     set_target_properties(${target} PROPERTIES LINK_OPTIONS  ${value})
#   endfunction()

# endif()
