# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2024 INRIA.
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
Check if nlohmann is required by Siconos.

See https://github.com/nlohmann/json

A library to use json from c++.

If so:
- look for the proper nlohmann version 
- configure siconos to link with nlohmann

The following options control if and how JSON is used

* WITH_JSON: activate json (based on nlohmann) and look for it. Stop if not found.
* JSON_ROOT=<some_path>: activate json/nlohmann and look for it in <some_path>. Stop if not found
* JSON_INSTALL: activate json/nlohmann, download and install the proper version at the same place as Siconos (CMAKE_INSTALL_PREFIX). This
location might be used later as input to nlohmann_ROOT.

#]=======================================================================]

# function(set_bullet_target)
#   if(BULLET_DEFINITIONS)
#     string(REPLACE "-D" "" BULLET_DEFINITIONS ${BULLET_DEFINITIONS})
#     set(BULLET_DEFINITIONS ${BULLET_DEFINITIONS} PARENT_SCOPE)
#   endif()
#   create_target(NAME BULLET::BULLET
#     LIBRARIES "${BULLET_LIBRARIES}"
#     LIBRARY_DIRS "${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS}"
#     INCLUDE_DIRS "${BULLET_INCLUDE_DIRS}"
#     COMPILE_DEFINITIONS "${BULLET_DEFINITIONS}")
#   # Add bullet headers and libs to the build.
#   # first draft ... turn this to private later
#   target_link_libraries(${COMPONENT} PUBLIC $<BUILD_INTERFACE:BULLET::BULLET>)
#   set(SICONOS_HAS_BULLET TRUE CACHE INTERNAL "True if Bullet API has been found and is activated.")
# endfunction()

# Three ways:
  # - JSON_INSTALL=ON : use fetchcontent to download and install json/nlohmann as a siconos part;
  # - WITH_JSON is ON, nothing more: look for json/nlohmann, check version and link with siconos components.
  # - User asks explicitely for a specific (already installed) version of json/nlohmann
  #   by providing JSON_ROOT on cmake command line.
  #   => find it and check the version

  
# Full config :
# - Download, build and install Bullet
# - Create a target JSON::NLOHMANN
# - Link mechanics with this target
if(WITH_JSON_INSTALL)
  include(FetchContent)
  message(STATUS "nlohmann for json will be downloaded from github repository and installed as a siconos component")

  # cmake_policy(SET CMP0077 NEW) # option() honors normal variables
  # cmake_policy(PUSH)
  # cmake_policy(SET CMP0072 NEW) # FindOpenGL prefers GLVND
  # cmake_policy(PUSH)
  # if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20")
  #   cmake_policy(SET CMP0115 OLD) # 
  #   cmake_policy(PUSH)
  # endif()
  set(FETCHCONTENT_QUIET OFF) # verbose mode for fetchcontent. Comment/uncomment according to your needs.

  
  set(JSON_Install ON) # To install nlohmman targets
  FetchContent_Declare(json
    GIT_REPOSITORY    https://github.com/nlohmann/json.git
    GIT_TAG          v3.11.2 # Tag 3.11.2  corresponds to version 3.2, so it seems ...
    GIT_SHALLOW TRUE
    UPDATE_DISCONNECTED TRUE # Do not update git repo at each run
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    )
  FetchContent_MakeAvailable(json)
    
  message(STATUS "Built, installed and used nlohmann  version ${nlohmann_json_VERSION_STRING} in ${nlohmann_json_ROOT_DIR}.")
  set(nlohmann_json_VERSION 3.2 CACHE INTERNAL "Json version") 
  set(SICONOS_HAS_JSON TRUE CACHE INTERNAL "Json activated and found")
  set(nlohmann_json_DIR ${CMAKE_INSTALL_PREFIX} CACHE INTERNAL "") # for siconos-config generation, to help finding nlhomman at runtime.
  
elseif(WITH_JSON OR JSON_ROOT)
  # Up to cmake 3.22 there is no way to get Bullet version using find_package standard.
  # It's then mandatory to use the 'config' versionb of find_package.
  # Anyway, find_package is not able to check the version since Bullet does not provide a bullet-config-version or BulettConfigVersion file.
  find_package(nlohmann_json 3.2 REQUIRED)
  set(SICONOS_HAS_JSON TRUE CACHE INTERNAL "Json activated and found") 
  set(nlohmann_json_VERSION 3.2 CACHE INTERNAL "Json version") 

endif()

