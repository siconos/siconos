# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2020 INRIA.
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
#/

# Use this file to define functions or macros that
# are not available in older version of cmake, but
# used in Siconos cmake files.
#
# e.g. : target_link_options exists from cmake 3.13
# and could be replace by a call to set_target_properties for older version.
# 

if(${CMAKE_VERSION} VERSION_LESS "3.14")
  # https://cmake.org/cmake/help/latest/prop_gbl/CMAKE_ROLE.html
  set_property(GLOBAL PROPERTY CMAKE_ROLE PROJECT)
  
endif()


if(${CMAKE_VERSION} VERSION_LESS "3.13")

  function(target_link_options target prop value)
    set_target_properties(${target} PROPERTIES LINK_OPTIONS  ${value})
  endfunction()
  
endif()
